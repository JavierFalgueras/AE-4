#!/usr/bin/env python3
# encoding: utf-8
# ae4.py
# Carlos Villagrasa, Javier Falgueras
# juanfc 2019-02-16
# ae4 2020-02-11

__version__ = 0.001 # 2020-02-12

import os
import sys
import builtins
from pathlib import Path
import json
import argparse
import textwrap
from pprint import pprint

from numpy import *
from numpy.random import seed, randint, shuffle, sample, permutation
# format for printing numpy arrays
set_printoptions(formatter={'int': '{: 7d}'.format})

import xlsxwriter
# sys.path.append('/usr/local/lib/python3.7/site-packages')
# import xlrd
from datetime import datetime



# ######### #
# CONSTANTS #
# ######### #
INDIVIDUAL = 0   # they consume df
ACTOR      = 1   # they consume df+if of recipient
RECIPIENT  = 2   # they consume df (it doesn't consume directly)
RECIPROCAL = 3   # they consume (intra 2x df+if) (inter df1+if1+df2+if2)
NUMBER_OF_FORMS = 4
# Types of distribution
NEIGHBOURS_DISTRIBUTION = 0 # randomly distribute among neighbours around
RANDOM_GLOBAL_AVG       = 1 # random change around global average

# extension for (re)new conf output files
CONT_FILE_NAME_SUFIX    = "_cont"

# ############################################################################# #
#                                MAIN SUBPROGRAMS                               #
# ############################################################################# #

def newWorld():
    """Build our world.
     First index:  cell index
     Second index: species index
    So we have:
        gWorld[iCell, iSpecies]
    """
    return (zeros((gNumberOfCells,  gNumberOfSpecies), dtype=int),
            zeros((gNumberOfSpecies, NUMBER_OF_FORMS), dtype=int),
            zeros((gNumberOfSpecies, NUMBER_OF_FORMS), dtype=int)
            )

def randomDist(nitems):
    """Returns a list of nitems numbers randomly distributed in gNumberOfCells"""

    r = sort(randint(0, nitems+1, gNumberOfCells-1))
    return concatenate((r,[nitems])) - concatenate(([0],r))

def doInitialDistribution():
    """NOT anymore: In fact this is unnecessary as first step in the generation process
    will do it (again) distribute all equally
    doDistribute() left for the end of generations loop"""

    # If a world file is there, take it!, but only when the number of cells fits

    if gNumberOfCells == checkAndCountPrevWorldNumberOfCells(gWorldCompFileName):
        # adding the possibility of changing the number of cells
        # First: compute the previous number of them
        # OJO: traspuesto el fichero World: ahora cada fila es una celda, y cada columna especie
        #
        gWorld[:,:] = [list(map(int, line.split())) for line in open(gWorldCompFileName)]
    else:
        # Simplest, averaged distribution
        # Should we call doDistribute() at the end of this?
        for iSpecies in range(gNumberOfSpecies):
            n = gConf["species"][iSpecies]["NumberOfItems"]
            meanPerCell = n // gNumberOfCells
            remainder   = n %  gNumberOfCells
            gWorld[ :, iSpecies] = [meanPerCell for _ in range(gNumberOfCells)]
            # the excess of items are equally given out to the first cells
            for iCell in range(remainder):
                gWorld[iCell, iSpecies] += 1
        doDistribute()

def doDistribute():
    """Distribute each each species in cells
    There are two ways of distribution:
        - NEIGHBOURS_DISTRIBUTION
        - RANDOM_GLOBAL_AVG
    Each species (iSpecies) can have either of them
    """

    global gStdDevDistr, gWorld
    noNegger = vectorize(noNeg)
    gWorld = noNegger(gWorld)
    for iSpecies in range(gNumberOfSpecies):
        distType, distVal = getDist(iSpecies)

        #             NEIGHBOURS_DISTRIBUTION
        if distType == NEIGHBOURS_DISTRIBUTION:
            tempDist = zeros((gNumberOfCells), dtype=int)
            distLen = int(gNumberOfCells*distVal/100/2)
            for iCell in range(gNumberOfCells):
                nitems = gWorld[iCell, iSpecies]
                # TODO: if the number of items nitems
                # is large, we could consider directly writing the
                # average on each cell
                for _ in range(nitems):
                    dist = randint(-distLen, distLen+1)
                    p = iCell+dist
                    if p >= gNumberOfCells:
                        p %= gNumberOfCells
                    tempDist[p] += 1
            gWorld[ :, iSpecies] = tempDist


        else:       # RANDOM_GLOBAL_AVG: x_a \sigma + xm (1-\sigma)
            sigma = distVal/100
            nitems = gWorld[ :, iSpecies].sum()
            average = nitems / gNumberOfCells
            wildDist = randomDist(nitems)
            # by Javi :)
            wildDist = around(wildDist * sigma + (1-sigma) * average).astype(int)
            # compensate rounding simply adding the difference
            # to the first cell
            dif = nitems - wildDist.sum()
            i = 0
            while dif != 0:
                if wildDist[i] + dif >= 0:
                    wildDist[i] += dif
                    dif = 0
                else:
                    dif += wildDist[i]
                    wildDist[i] = 0
                i += 1

            gWorld[ :, iSpecies] = wildDist
    gStdDevDistr = around(std(gWorld[:,:], axis=0), decimals=1)

def doGrouping():
    """Form the groups following the group partner"""
    for iCell in range(gNumberOfCells):
        n = 0
        # printv("Cell grouping %d" % n)
        n += 1
        stillPossibleGroupInCell = False
        permutedListOrigGroups = permutation(gWithPartnerList)
        # printv("Grouping permutation list:", permutedListOrigGroups)
        for iSpecies in permutedListOrigGroups:
            phenFlex = gConf["species"][iSpecies]["PhenotypicFlexibility"]
            i_partnerList = iGetPartnerList(iSpecies)      # list of i partners for the groups
            i_groupedList = iGetGroupStartingInList(iSpecies)  # i already formed group that iSpecies starts
            # number of iSpecies items in that iCell
            ni = gWorld[iCell, iSpecies]

            # item to group with
            ni_partnerList = array([ gWorld[iCell, i_partner] for i_partner in i_partnerList ])

            # Phenotypic complexity tells the % (0-1) from the total
            # existing particular species, that should be grouped, so we
            # compare the current amount of items (ungrouped)
            # and check if there is room for more grouping
            ni_toGroup = int(ni * phenFlex)
            ni_total_partners = ni_partnerList.sum()
            if ni_toGroup < ni_total_partners:
                ni_toGroupList = trunc(ni_partnerList/ni_total_partners * ni_toGroup).astype(int)
            else:
                ni_toGroupList = ni_partnerList

            ni_feasible = ni_toGroupList.sum()

            #printv("ni_toGroup: %d with %d, ni_feasible: %d" % (ni_toGroup, i_grouped, ni_feasible))
            if ni_feasible > 0:
                for i in range(len(i_partnerList)):
                    n = ni_toGroupList[i]
                    if i_partnerList[i]  == iSpecies:
                        n //= 2
                    gWorld[ iCell,iSpecies] -= n
                    gWorld[iCell, i_partnerList[i]] -= n

                    gWorld[iCell, i_groupedList[i]] += n

def doAssociationAndQueue():
    """ performs the association and immediately moves the new associated
    to the queue ready to shuffle and then eat
    Participants in the association leave the cell both. Then individuals """

    # ######### #
    # CONSTANTS #
    # ######### #
    # act in queue
    AS_INDIVIDUAL = -1
    AS_RECIPIENT  = -2

    class CQueue():
        """Specific queue for each pair orig->associate to eat"""
        def __init__(self, size):
            self.size = size
            self.q = zeros((size, 2), dtype=int)
            self._top = 0

        def push(self, item):
            l = len(item)
            if l > 0:
                self.q[self._top:self._top+l] = item
                self._top += l

        def top(self):
            if self._top == 0:
                raise IndexError
            else:
                return self.q[self._top-1]

        def pop(self):
            if self._top == 0:
                raise IndexError
            else:
                self._top -= 1
                return self.q[self._top]

        def shuffle(self):
            shuffle(self.q)

        def __len__(self):
            return self._top

        def __getitem__(self, position):
            if position <= self._top:
                return self.q[position]
            else:
                raise IndexError

        def __str__(self):
            return str(self.q)

        def __iter__(self):
            self.curi = 0
            return self

        def __next__(self):
            if self.curi < self.size:
                self.curi += 1
                return self.q[self.curi-1]
            else:
                raise StopIteration


    if gArgs["varia"]:
        newDirFit  = zeros((gNumberOfSpecies), dtype=[("sum", "int"), ("N", "int")])



    for iCell in range(gNumberOfCells):
        # ENQUEUEING
        # queue the index of item and its receptor
        # the second number next to the index is:
        #     -2 for recipient AS_RECIPIENT
        #     -1 for individual AS_INDIVIDUAL
        # if second number >= 0
        #   If second number == first => receptor == orig => reciprocal intraspecific
        #   If second number != first => receptor != orig => reciprocal interspecific
        queue = CQueue(gWorld[iCell,:].sum())   # size of the queue, larger than necessary
        printv("number items in     Cell: ", gWorld[iCell,:].sum())


        listOfAssocOrigs = permutation(gListOfAssociationOrigs)
        for iOrig in listOfAssocOrigs:
            ni_Orig = gWorld[iCell, iOrig]
            iAssocList = iAssociatedTargetGetList(iOrig)
            # build list of populations
            ni_assocList = array([ gWorld[iCell, i] for i in iAssocList ])
            ni_assocListSum = ni_assocList.sum()

            if ni_Orig >= ni_assocListSum:
               ni_toAssocList = ni_assocList
            else:
               ni_toAssocList = trunc(ni_assocList * (ni_Orig/ni_assocListSum) ).astype(int)

            # Kind of association
            # IRAR: Individual, Recipient, Actor, Reciprocal
            for iAssoc, ni_Assoc in zip(iAssocList, ni_toAssocList):
                if iOrig == iAssoc: # is Reciprocal intraspecific  A><A
                    ni_Single = ni_Assoc % 2
                    n = ni_Assoc - ni_Single

                    queue.push(n * [[iOrig, iOrig]])              # RECIPROCAL INTRA
                    gWorld[iCell, iOrig] = ni_Single  # the only one left in cell, if
                    gStatsAnt[iOrig, RECIPROCAL] += n

                else:

                    # either A><B or A>B>C, A>B
                    # ni_feasible = min(ni_Orig, gWorld[iCell, iAssoc])
                    # OJJOJOSJOOOO OJO

                    queue.push(ni_Assoc * [[iOrig, iAssoc]])
                    gStatsAnt[iOrig, ACTOR] += ni_Assoc

                    # Here we tie the associated with ourself, if there was that association
                    if iOrig in iAssociatedTargetGetList(iAssoc): # is iOrig among the assoc of assoc
                        queue.push(ni_Assoc * [[iAssoc, iOrig]])  # RECIPROCAL INTER
                        gStatsAnt[iOrig, RECIPROCAL] += ni_Assoc
                    else:
                        queue.push(ni_Assoc * [[iAssoc, AS_RECIPIENT]])     # RECIPIENT
                        gStatsAnt[iOrig, RECIPIENT] += ni_Assoc

                    gWorld[iCell, iOrig] -= ni_Assoc
                    gWorld[iCell, iAssoc] -= ni_Assoc

                    # if iOrig in iAssociatedTargetGetList(iAssoc): # is iOrig among the assoc of assoc
                    #     #  A><B reciprocal interspecific
                    #     gWorld[iOrig, iCell, RECIPROCAL] += ni_feasible
                    #     gWorld[iAssoc, iCell, RECIPROCAL] += ni_feasible
                    # else:  # A>B>C or A>B  simple association
                    #     gWorld[iOrig, iCell,     ACTOR] += ni_feasible
                    #     gWorld[iAssoc, iCell, RECIPIENT] += ni_feasible

        # ADD INDIVIDUALS TO THE QUEUE
        for iSpecies in range(gNumberOfSpecies):
            n = gWorld[iCell, iSpecies]
            if n > 0:
                queue.push( n * [[iSpecies, AS_INDIVIDUAL]] )
                gWorld[iCell, iSpecies] = 0
                gStatsAnt[iSpecies, INDIVIDUAL] += n

        # CONSUME!
        printv("queue size before eating: ", len(queue))
        queue.shuffle()
        rsrc = gConf["NumberOfRsrcsInEachCell"]
        for couple in iter(queue):
            iOrig = couple[0]
            iDest = couple[1]

            dirFit   = gConf["species"][iOrig]["DirectOffspring"]

            toEat = 0
            if iDest < 0:  # iOrig is INDIVIDUAL or RECIPIENT
                toEat = dirFit
                ##TODO? toEat = dirFit + abs(indirFit)
                #       taking advantage of the indirFit for himself
                if rsrc >= toEat:
                    gWorld[iCell, iOrig] += toEat
                    if iDest == AS_INDIVIDUAL:
                        gStatsPost[iOrig, INDIVIDUAL] += 1
                    else:
                        gStatsPost[iOrig, RECIPIENT] += 1
            else:
                # giving indirect offspring
                indirFit = gConf["species"][iOrig]["IndirectOffspring"]
                if gArgs["varia"]:
                    ##TODO only for those with list of assoc and after
                    #      seeing the form
                    fitVarLimit   = gConf["species"][iOrig]["FitnessVariationLimit"]
                    dirFit, indirFit = fitnessVariations(dirFit, indirFit, fitVarLimit)
                    newDirFit[iOrig]["sum"] += dirFit * dirFit
                    newDirFit[iOrig]["N"] += dirFit

                # The others (ACTOR, and RECIPROCAL) have to give IndirectOffspring
                toEat = dirFit + abs(indirFit)
                if rsrc >= toEat:
                    gWorld[iCell, iOrig] += dirFit
                    gWorld[iCell, iDest] += indirFit
                    # is iOrig among the assoc of assoc
                    if iOrig == iDest or iOrig in iAssociatedTargetGetList(iDest):
                        gStatsPost[iOrig, RECIPROCAL] += 1
                    else:
                        gStatsPost[iOrig, ACTOR] += 1

            rsrc -= toEat

    if gArgs["varia"]:
        for iSpecies in range(gNumberOfSpecies):
            # CHANGE THE GLOBAL DirectOffspring gConf parameter
            if newDirFit[iSpecies]["N"]:
                incFit = gConf["species"][iSpecies]["DirectOffspring"] + gConf["species"][iSpecies]["IndirectOffspring"]
                dirFit = int(round(newDirFit[iSpecies]["sum"] / newDirFit[iSpecies]["N"]))
                gConf["species"][iSpecies]["DirectOffspring"] = dirFit
                gConf["species"][iSpecies]["IndirectOffspring"] = incFit - dirFit

def doUngroup():
    for iCell in range(gNumberOfCells):
        for iOrig in gWithPartnerList:
            iGroupList = iGetGroupStartingInList(iOrig)
            for iGroup in iGroupList:
                phenFlex = gConf["species"][iGroup]["PhenotypicFlexibility"]

                # number of iGroup items in that iCell
                ni = gWorld[iCell,iGroup] # form index is 0 always

                ni_unGroup = int(round(ni * (1.0 - phenFlex)))

                if ni_unGroup > 0:
                    iPartner = iGetPartnerFromOrigAndGroup(iOrig, iGroup)
                    gWorld[iCell,   iOrig] += ni_unGroup
                    gWorld[iCell,iPartner] += ni_unGroup
                    gWorld[iCell,  iGroup] -= ni_unGroup


# ############################################################################# #
#                                     TOOLS                                     #
# ############################################################################# #

def iFrom_id(id):
    """Returns index of species from its id, or -1"""
    i = gNumberOfSpecies - 1
    while i >= 0 and gConf["species"][i]["id"] != id:
        i -= 1
    return i

def iListFrom_idList(idList):
    """Returns index of species from its id, or -1"""
    return [iFrom_id(toFindId) for toFindId in idList]

def noNeg(n):
    if n < 0: n = 0
    return n

# ------------ Groups

def iGetPartnerList(iSpecies):
    "returns a list of indexes of the partners for grouping"
    toFindIdList = gConf["species"][iSpecies]["GroupPartner"]
    return [iFrom_id(toFindId) for toFindId in toFindIdList]

def iGetGroupStartingInList(iSpecies):
    idOrig = gConf["species"][iSpecies]["id"]
    toFindIList = [idOrig + '|' + idPartner for idPartner in gConf["species"][iSpecies]["GroupPartner"]]
    return [iFrom_id(toFindId) for toFindId in toFindIList]

def iGetPartnerFromOrigAndGroup(iOrig, iGroup):
    idOrig = gConf["species"][iOrig]["id"]
    idGroup = gConf["species"][iGroup]["id"]
    return iFrom_id(idGroup[len(idOrig)+1:])

def getListOfOrigGroups():
    theList = [i for i in range(gNumberOfSpecies) if len(gConf["species"][i]["GroupPartner"]) ]

    # Puts complex groups first
    theList.sort(key = lambda i:
            gConf["species"][i]["id"].count('|') +
            str(gConf["species"][i]["GroupPartner"]).count('|'), reverse=True)
    return theList

# ------------ Association

def collectAssociationActors():
    """Return a list of indexes of origs that are actually associated with
    some other Now when no associates, is [] """

    return [i for i in range(gNumberOfSpecies) if len(gConf["species"][i]["AssociatedSpecies"]) ]

def iAssociatedTargetGetList(iSpecies):
    '''list of is of its associate '''
    return iListFrom_idList(gConf["species"][iSpecies]["AssociatedSpecies"])

# ------------ Others

def getDist(iSpecies):
    """from 45n or 50n returns the integer and the type
    n: NEIGHBOURS_DISTRIBUTION
    r: RANDOM_GLOBAL_AVG
    """
    if "Distribution" in gConf["species"][iSpecies] or "Distribution" in gConf:
        if "Distribution" in gConf["species"][iSpecies]:
            dist = gConf["species"][iSpecies]["Distribution"]
        else:
            dist = gConf["Distribution"]

        if dist.endswith("n"): # average neighbours distribution
            distType = NEIGHBOURS_DISTRIBUTION
        else:
            distType = RANDOM_GLOBAL_AVG

        distVal = int(dist[:-1])
        if distVal > 100 or distVal < 0:
            print("Error in Distribution val (should be 0<=x<=100. Is:", distVal)
            sys.exit(1)

    if "DistType" in gConf["species"][iSpecies] or "DistType" in gConf:
        if "DistType" in gConf["species"][iSpecies]:
            dist = gConf["species"][iSpecies]["DistType"]
            distVal = int(gConf["species"][iSpecies]["DistVal"])
        else:
            dist = gConf["DistType"]
            distVal = int(gConf["DistVal"])

        if dist == "n": # average neighbours distribution
            distType = NEIGHBOURS_DISTRIBUTION
        else:
            distType = RANDOM_GLOBAL_AVG

        oldStyle = str(distVal)+dist
        if "DistType" in gConf["species"][iSpecies]:
            gConf["species"][iSpecies]["Distribution"] = oldStyle
        else:
            gConf["Distribution"] = oldStyle

        if distVal > 100 or distVal < 0:
            print("Error in Distribution val (should be 0<=x<=100. Is:", distVal)
            sys.exit(1)


    return distType, distVal

def fitnessVariations(direct, indirect, fitVarLimit):
    """Returns a pair of new rand int direc-indirect fitnesses inside the
    fitVarLimit"""

    tot = direct + abs(indirect)

    # the pairs ordered in the smoothest possible way
    directRange = list(range(tot, -1, -1)) + list(range(0, tot))
    indireRange = list(range(0, tot + 1)) + list(range(-tot, 0, 1))
    pairsOfFitness = list(zip(directRange, indireRange))

    currentI = pairsOfFitness.index((direct, indirect))
    lpairsOfFitness = len(pairsOfFitness)

    # 4 cases: fitVarLimit > Maxfitness, upper-overflow, lower-underflow, inside
    if 2*fitVarLimit >= lpairsOfFitness:
        bottom = 0
        topp = lpairsOfFitness
    elif currentI + fitVarLimit >= lpairsOfFitness:
        extraOver = currentI + fitVarLimit - (lpairsOfFitness - 1)
        bottom = currentI - fitVarLimit - extraOver
        topp = lpairsOfFitness
    elif currentI - fitVarLimit < 0:
        extraBelow = fitVarLimit - currentI
        bottom = 0
        topp = currentI + fitVarLimit + extraBelow
    else:
        bottom = currentI - fitVarLimit
        topp = currentI + fitVarLimit + 1

    newI = randint(bottom, topp)
    return pairsOfFitness[newI]

def ranking(l, value):
    rank = unique(l.flatten())
    return 1+where(rank==value)[0][0]

def checkAndCountPrevWorldNumberOfCells(fName):
    if not os.path.isfile(fName):
        return 0
    with open(fName) as f:
        return builtins.sum(1 for _ in f)

def initExcel():
    global gExcelCellHeader, gExcelCellID, gExcelWorkbook, gExcelWorksheet, gGlobalExcel

    # # Check if there already is a collective Excel
    # if gGlobalExcel:
    #     gExcelWorkbook = xlsxwriter.Workbook(gGlobalExcel)
    #     if os.path.isfile(gGlobalExcel):
    #         prexcel = xlrd.open_workbook(gGlobalExcel)
    #         sheets = prexcel.sheets()
    #         # run through the sheets and store sheets in workbook
    #         # this still doesn't write to the file yet
    #         for sheet in sheets: # write data from old file
    #             newSheet = gExcelWorkbook.add_worksheet(sheet.name)
    #             for row in range(sheet.nrows):
    #                 for col in range(sheet.ncols):
    #                     newSheet.write(row, col, sheet.cell(row, col).value)

    #     gExcelWorksheet = gExcelWorkbook.add_worksheet(gOutFNameBase)
    # else:
    #     excelOut = os.path.join(gOutDir, gOutFNameBase + ".xlsx")
    #     gExcelWorkbook = xlsxwriter.Workbook(excelOut)
    #     gExcelWorksheet = gExcelWorkbook.add_worksheet(gOutFNameBase)
    excelOut = os.path.join(gOutDir, gOutFNameBase + ".xlsx")
    gExcelWorkbook = xlsxwriter.Workbook(excelOut)
    gExcelWorksheet = gExcelWorkbook.add_worksheet()

    gExcelCellHeader = gExcelWorkbook.add_format({
                                           'align':'center', 'bg_color': '#CCCCFF', 'bold': True})
    gExcelCellID = gExcelWorkbook.add_format({'align':'center', 'bg_color': '#33FF99', 'bold': True})

    gExcelWorkbook.set_properties({
        'title':    'Evolutionary Automata',
        'subject':  gInitConfFile + '_' + gThedatetime,
        'author':   'Javier Falgueras, Juan Falgueras, Santiago Elena',
        'manager':  'Javier Falgueras, Juan Falgueras, Santiago Elena',
        'company':  'Univ. Valencia + Univ. Málaga',
        'category': 'Research Thesis',
        'keywords': 'Evolution Automata',
        'comments': " ".join(sys.argv)})

def saveExcel(numGen):
    """Saving in Excel"""
    # id, NumberOfItems, DirectOffspring, IndirectOffspring,
    # AssociatedSpecies, FitnessVariationLimit, INDIVIDUAL, ACTOR, RECIPIENT,
    # RECIPROCAL
    if 'gExcelCellHeader' not in globals():
        initExcel()
    txtOutName = os.path.join(gOutDir, gOutFNameBase + ".txt")
    txtOut = open(txtOutName, "a")

    globalsHeader =  ["NCel", "RsCel", "Dst"]
    iSpecHeader =  ["ID", "dst", "D", "I", ">", "Fv", "Gr", "Ph", "Drσ", "IND", "2", "ACT", "2", "RNT", "2", "RCL", "2"]
    globalsHeaderLen = len(globalsHeader)
    iSpecHeaderLen = len(iSpecHeader)

    if numGen == 1:
        gExcelWorksheet.write_row(0,0, globalsHeader, gExcelCellHeader)
    gExcelWorksheet.write_row(numGen,0, [gConf["NumberOfCells"], gConf["NumberOfRsrcsInEachCell"], gConf["Distribution"]])
    print(gConf["NumberOfCells"], "\t", gConf["NumberOfRsrcsInEachCell"],"\t", gConf["Distribution"], "\t", file=txtOut, end="", sep="")


    if numGen == 1:
        for iSpecies in range(gNumberOfSpecies):
            gExcelWorksheet.write_row(0,globalsHeaderLen+iSpecies*iSpecHeaderLen, iSpecHeader, gExcelCellHeader)


    for iSpecies in range(gNumberOfSpecies):
        gExcelWorksheet.write(numGen, globalsHeaderLen+iSpecies*iSpecHeaderLen,
                 gConf["species"][iSpecies]["id"], gExcelCellID)
        print(gConf["species"][iSpecies]["id"] + "\t", file=txtOut, end="", sep="")

        if "Distribution" in gConf["species"][iSpecies]:
            speciesDist = gConf["species"][iSpecies]["Distribution"]
        else:
            speciesDist = ""

        toWrite =[
             speciesDist,
             gConf["species"][iSpecies]["DirectOffspring"],
             gConf["species"][iSpecies]["IndirectOffspring"],
             gConf["species"][iSpecies]["AssociatedSpecies"],
             gConf["species"][iSpecies]["FitnessVariationLimit"],
             ",".join(gConf["species"][iSpecies]["GroupPartner"]),
             gConf["species"][iSpecies]["PhenotypicFlexibility"],
             gStdDevDistr[iSpecies],
             gStatsAnt[ iSpecies,INDIVIDUAL],
             gStatsPost[iSpecies,INDIVIDUAL],
             gStatsAnt[ iSpecies,ACTOR],
             gStatsPost[iSpecies,ACTOR],
             gStatsAnt[ iSpecies,RECIPIENT],
             gStatsPost[iSpecies,RECIPIENT],
             gStatsAnt[ iSpecies,RECIPROCAL],
             gStatsPost[iSpecies,RECIPROCAL]
            ]
        gExcelWorksheet.write_row(numGen, globalsHeaderLen+iSpecies*iSpecHeaderLen+1, toWrite)
        print("\t".join(map(str,toWrite)), file=txtOut, end="", sep="")
        if iSpecies < gNumberOfSpecies-1:
            print("\t", file=txtOut, end="")

        if numGen == 1:
            ori = iSpecies*iSpecHeaderLen + 1 + globalsHeaderLen
            end = ori + iSpecHeaderLen - 2
            gExcelWorksheet.set_column(ori, end, None, None, {'level': 1, 'hidden': True})

    print(file=txtOut)

def saveConf(genNumber):
    """Save conf in a new _cont.json file
    if not there, or rewrite previous _cont.json file if was the initial conf loaded
    Added: save corresponding final gWorld state"""

    # gThedatetime = datetime.now().strftime("%Y%m%d-%H%M%S.%f")


    # save the conf, but ala!, the number of items is totally different
    for iSpecies in range(gNumberOfSpecies):
        gConf["species"][iSpecies]["NumberOfItems"] = int(gWorld[:, iSpecies].sum())
    with open(gNewConfCompFileName, 'w') as outfile:
        json.dump(gConf, outfile, sort_keys = True, indent = 4,
                   ensure_ascii = False)
    # save the matrix of cells state for possible re-reading
    with open(gNewWorldCompFileName, 'w') as outfile:
        for iCell in range(gNumberOfCells):
            print("\t".join(map(str, gWorld[iCell, :])), file=outfile)

def checkConf(conf):
    """Verifies conf returning inconsistencies or an empty string"""
    # TO DO
    # verify gConf so
    #   - groups agree
    #
    problems = ""

    items = {"NumberOfCells", "NumberOfRsrcsInEachCell",
        "MultilevelDeath1Percent", "LambdaForEgoism", "species"}

    itemsspecies = {"id", "NumberOfItems", "DirectOffspring",
        "GroupPartner", "PhenotypicFlexibility", "AssociatedSpecies",
        "IndirectOffspring", "FitnessVariationLimit"}

    if not items.issubset(conf.keys()):
        problems += "Some essential item(s) is not defined\n\t Check the next items are all there:\n\t%s\n" % str(list(items))

    # Check values for Distribution, if there
    if "Distribution" in conf and conf["Distribution"][-1] not in "rn":
        problems += "Distribution global value must end either in r or n\n"

    # Check values for DistType, if there
    if "DistType" in conf and conf["DistType"] not in "rn":
        problems += "DistType global value must be either r or n. You gave %s\n" % conf["DistType"]


    if "DistVal" in conf and (100 < conf["DistVal"] or 0 > conf["DistVal"]):
        problems += "DistVal global value given: %d. It must be between 0..100\n" % conf["DistVal"]

    if ("DistType" in conf) and ("DistVal" not in conf) or ("DistType" not in conf) and ("DistVal" in conf):
        problems += "Global. Both DistType/DisVal values must be provided when one of them is provided\n"

    previousIds = set([])
    longListSpecies = len(conf["species"])
    if type(conf["species"]) != list or longListSpecies == 0:
        problems += "Species must be a list with species inside\n"
    for i in range(longListSpecies):
        ##TODO check the IndirectOffspring and variability are 0 when
        #      the list of associated is empty
        confi = conf["species"][i]
        if "Distribution" in confi and confi["Distribution"][-1] not in "rn":
            problems += "Distribution for species %d must end either in r or n\n" % i

        if "DistType" in confi and confi["DistType"] not in "rn":
            problems += "DistType for species %d must end either in r or n\n" % i

        if "DistVal" in confi and (100 < confi["DistVal"] or \
                                                0   > confi["DistVal"]):
            problems += "DistVal for species %d must between 0..100\n" % i

        if ("DistType" in confi) and ("DistVal" not in confi) or ("DistType" not in confi) and ("DistVal" in confi):
            problems += "Species %d. Both DistType/DisVal values must be provided when one of them is provided\n" % i

        theId = confi["id"]
        # partnerId = confi["GroupPartner"]
        if theId in previousIds:
            problems += "Id '%s' REPEATED\n" % theId
        else:
            previousIds.add(theId)
        # if partnerId != "":
        #     if iGetPartnerList(i) == -1:
        #         problems += "Species %d, id: '%s' has not the partner species '%s' in the conf\n" % \
        #             (i, theId, partnerId)
        #     if iGetGroupStartingInList(i) == -1:
        #         problems += "Group id: '%s' is not configured\n" % \
        #             (theId + '|' + partnerId)
        for aId in theId.split('|'):
            if iFrom_id(aId) == -1:
                problems += "Component id: '%s' from group '%s' is not configured\n" % \
                    (aId, theId)


        if not itemsspecies.issubset(confi.keys()):
            problems += "Some essential item(s) of species %s is not defined\n\t Check the next items for all species are all there:\n\t%s\n" % (theId, str(list(itemsspecies)))

    if problems:
        print("\nERROR(s)\nIn the configuration file '%s' or in the command line arguments added there are inconsistencies:\n\n%s\n\n" % (gInitConfFile, problems))
        sys.exit(1)

def completeDataFileNames(fromfilename, inDirectory):
    conf = os.path.join(inDirectory, fromfilename + ".json")
    world = os.path.join(inDirectory, fromfilename + ".world.txt")
    return conf, world

def printv(*args, **kwargs):
    if gArgs["verbose"]:
        if type(args[0]) == str:
            print(*args)
        else:
            pprint(args)

# ######################## #
# ARGS AND CONF PROCESSING #
# ######################## #

def defineAndGetCommandLineArgs():
    """Parses and returns command line arguments"""

    theArgParser = argparse.ArgumentParser(description="* Evolutionary Automata *",
                    formatter_class=argparse.RawTextHelpFormatter)


    # initial configuration to load from a file
    theArgParser.add_argument(
        "initFile",
        nargs="?",
        default="defaultInit",
        help=textwrap.dedent("""\
        When provided it sets the file with the initial configuration.
        All the configuration files are inside the './data' directory.
        Do not add the .extension to the file name.
        If not provided, it defaults to 'defaultInit'
        HINT:
            You can use
                egoism/test1
            as init file and then, _cont files will be generated inside that data/egoism/
            folder""")
    )

    # Number of generations to run
    theArgParser.add_argument(
        "--numGen", type=int,
        default=10,
        metavar="int",
        help="Sets the number of generations to run, default: 10")

    # We want phenotypic variability
    theArgParser.add_argument(
        "--varia", help=textwrap.dedent("""\
        It takes the FitnessVariationLimit for each species
        and changes in each generation the Direct/Indirect fitnesses
        inside that limits: currentValue+- FitnessVariationLimit
        following a list of pairs.
        False if not provided"""),
        action='store_true')

    # Setting random seed to 0 to repeat pseudorandom values
    theArgParser.add_argument(
        "--setRandomSeed", type=int,
        default=-1,
        metavar="int",
        help=textwrap.dedent("""\
        Sets the random seed to a fixed initial value so
          to repeat same random sequences. If not, each
          running will start with different seeds"""),
        )

    # We want phenotypic variability
    theArgParser.add_argument(
        "--verbose", help=textwrap.dedent("""\
        Gives as much detailed information as it can"""),
        default=False,
        action='store_true')

    # # We want egoism multilevel selection
    # theArgParser.add_argument(
    #     "--egoism", help=textwrap.dedent("""\
    #     Considers egoism of each item in multilevel selection"""),
    #     action='store_true')

    # We want to save final status
    theArgParser.add_argument(
        "--saveExcel", help=textwrap.dedent("""\
        Save stats in 'Excel' file and in a txt file.
        See Excel iSpecHeader for meaning of txt columns
        It is set anyway to True
            if any --outFName or outDir is set"""),
        action='store_true')

    # Directory for the results.
    # If no /start relative to ./results
    theArgParser.add_argument(
        "--outDir", type=str,
        metavar="'str'",
        help=textwrap.dedent("""\
        Specifies another than the default output directory
        where to save the .txt and .xlsx files with global outputs.
        If the path of the output dir starts in /
            is considered an ABSOLUTE path.
        else
            is considered a path relative to ./results/
        It sets --saveExcel to True
        Example:
           --outDir=assocTests
           --outDir=/absoluteDir
           --outDir=more/than/one/level""")
    )

    # Filename for the results.
    # If none provided, the default original
    # conf+date is taken
    theArgParser.add_argument(
        "--outFName", type=str,
        metavar="'str'",
        help=textwrap.dedent("""\
        Specifies another than the default output filename
        (initFilename+date)
        where to save the .txt and .xlsx files with global outputs.
        It sets --saveExcel to True
        It appends a new sheet after to the Excel file. The name
        of the Excel file is the name of the output directory
        Example:
            --outFNname=assocTest20""")
    )

    theArgParser.add_argument(
        "--NumberOfCells", type=int,
        default=argparse.SUPPRESS, metavar="int",
        help="Sets the total number of cells in our world")

    theArgParser.add_argument(
        "--NumberOfRsrcsInEachCell", type=float,
        default=argparse.SUPPRESS, metavar="float",
        help="Sets the items of resources in each cell for each generation")
    # theArgParser.add_argument(
    #     "--MultilevelDeath1Percent", type=int,
    #     default=argparse.SUPPRESS, metavar="int",
    #     help="Sets MultilevelDeath1Percent to. Range 0.0 to 1.0")
    # theArgParser.add_argument(
    #     "--LambdaForEgoism", type=float,
    #     default=argparse.SUPPRESS, metavar="float",
    #     help="Sets the coef. for egoism application")


    theArgParser.add_argument(
        "--Distribution", type=str,
        metavar="'str'",
        default=argparse.SUPPRESS,
        help=textwrap.dedent("""\
        It allows specify the desired type of distribution between:
          a) RANDOM_GLOBAL_AVG           from 0r to 100r (with suffix r)
                           to distribute all the elements among every cell
                           0r    means in each cell the same value averaged
                           100r: means a total random number in each cell
                           Items are taken from cells forgetting their
                           previous position and then distributed randomly
          b) NEIGHBOURS_DISTRIBUTION     from 0n to 100n (with suffix n)
                           to distribute each cell items among the
                           n%%/2 cells around

          0r   means null random, it is equivalent to 100n
                distributed in case (b)
          100r means _totally_ random numbers en each cell

        When distributed through neighbours, the number of cells given
          is taken from the left, from the right and the
          [i] cell to distribute randomly among them
             [i-n][i-(n-1)][…][i][i+1][…][i+n]
          The list is considered circular
          0 means no cells around are considered,
            so cells are isolated
          3 for example, means 3x2+1 cells
            are used for distribution: the central one plus the 3 at
            the left plus the 3 at the right""")
    )


    theArgParser.add_argument(
        "--DistType", type=str,
        metavar="'str'",
        default=argparse.SUPPRESS,
        help=textwrap.dedent("""\
        It allows specify the desired type of distribution between:
          a) RANDOM_GLOBAL_AVG        r
                           to distribute all the elements among every cell
                           0 r    means in each cell the same value averaged
                           100 r: means a total random number in each cell
                           Items are taken from cells forgetting their
                           previous position and then distributed randomly
          b) NEIGHBOURS_DISTRIBUTION  n
                           to distribute each cell items among the
                           n%%/2 cells around

        --DistType=r
        --DistType=n
        If used, this parameter cancels --Distribution parameter""")
    )

    theArgParser.add_argument(
        "--DistVal", type=int,
        default=argparse.SUPPRESS, metavar="int",
        help=textwrap.dedent("""\
        It specifies the amount (integer) between 0-100
        for the active kind of distribution (r or n)
        If used, this parameter cancels --Distribution parameter""")
    )


    theArgParser.add_argument(
        "--species", type=str,
        default=argparse.SUPPRESS, metavar="'str'",
        help=textwrap.dedent("""\
    Sets one or more species parameters
    Surround everything with " "
    No spaces

    Put first the 0-index of the species:
        --species="0;DirectOffspring=8;1;IndirectOffspring=1"
    Changes the DirectOffspring of the first species and the
    IndirectOffspring of the second.

    A series of indexes separated by spaces is valid:
        --species="0;3;1;DirectOffspring=8"
    changes DirectOffspring of species 0, 3 and 1

    A Python standard range (no 'steps' here) is valid:
        --species="0;8:10;DirectOffspring=8"
    as with Python 8:10 means 8,9 (10 is not included)

    When an index is negative, is suppose from the end:
        --species="0;3;1;8:-2;DirectOffspring=8"
    where 8:-2  is  8,9,10,11,12,13  if there were 15 species.

    An empty index in ranges means a extreme:
        --species=":;DirectOffspring=8"
    means change all species
    """))


    return vars(theArgParser.parse_args())

def readInitConfFile(fileName):
    """Load config from init file"""
    with open(fileName) as f:
        conf = json.load(f)

    return conf

def replaceArgsInConf(conf, args):
    """Replace the conf args with the args provided in the command line"""
    def parseSpeciesArgs(s):
        def isInt(anyNumberOrString):
            try:
                int(anyNumberOrString)
                return True
            except ValueError:
                return False

        lenSpecies = len(conf["species"])
        l = s.split(';')
        indxsRead = True
        for arg in l:
            if isInt(arg) or ':' in arg:   # int (+-) or range
                if indxsRead:
                    ns = []
                    indxsRead = False
                if isInt(arg):
                    ns.append(int(arg))
                else:
                    if arg == ':':
                        first = 0
                        last = lenSpecies
                    elif arg.startswith(':'):
                        first = 0
                        last = int(arg[1:])
                    elif arg.endswith(':'):
                        first = int(arg[1:])
                        last = lenSpecies
                    else:
                        first, last = map(int,arg.split(':'))
                    if first < 0: first += lenSpecies
                    if last  < 0: last  += lenSpecies
                    ns += list(range(first, last))
                    printv(ns)
            else:
                indxsRead = True
                k, val=arg.split('=')
                for n in ns:
                    if k in conf["species"][n]:
                        # if it was an existing key, adapt
                        t = type(conf["species"][n][k])
                        conf["species"][n][k] = t(val)
                    else:
                        # if not, create
                        conf["species"][n][k] = val

    # substitute init file conf specific parameters
    # provided through the args in the command line
    for param in args:
        if param in conf:
            if param == "species":
                parseSpeciesArgs(args["species"])
            else:
                if conf[param] != None:
                    t = type(conf[param])
                    conf[param] = t(args[param])
        else:
            conf[param] = args[param]
    return conf


# ############################################################################## #
#                                    MAIN                                        #
# ############################################################################## #

# ####### #
# GLOBALS #
# ####### #




# INPUT
# GET CONF FROM .json INIT FILE AND THEN FROM ARGS
# init conf files (json and world)
gArgs = defineAndGetCommandLineArgs()
gInitConfFile = gArgs["initFile"]  # base file name
gInDir = "data"
# A -> A.json, A.world.txt
gInitConfCompName, gWorldCompFileName = completeDataFileNames(gInitConfFile, gInDir)
# 1. A.json -> A_cont.json A_cont.json
if gInitConfFile.endswith(CONT_FILE_NAME_SUFIX):
    newInitFileNameBase = gInitConfFile
else:
    newInitFileNameBase = gInitConfFile + CONT_FILE_NAME_SUFIX

# 2. …and newWorld and newConf cont out
gNewConfCompFileName, gNewWorldCompFileName = completeDataFileNames(newInitFileNameBase, gInDir)

gConf = replaceArgsInConf(readInitConfFile(gInitConfCompName), gArgs)

gNumberOfSpecies = len(gConf["species"])

checkConf(gConf)


# SET UP CONVENIENT GLOBALS
#
gNumberOfCells   = gConf["NumberOfCells"]

gEgoism          = []
gToSaveExcel     = gArgs["saveExcel"]
gGlobalExcel = None
gListOfAssociationOrigs = collectAssociationActors() # list of species starting association


gWithPartnerList = getListOfOrigGroups()
printv("gWithPartnerList:", gWithPartnerList)


# OUTPUT DIRS AND FILENAMES
#
gOutDir = "results"
if gArgs["outDir"]:
    gToSaveExcel = True
    if gArgs["outDir"].startswith("/"):
        gOutDir = gArgs["outDir"]
    else:
        gOutDir = os.path.join(gOutDir, gArgs["outDir"])
    gGlobalExcel = os.path.join(gOutDir, "global.xlsx")

Path(gOutDir).mkdir(parents=True, exist_ok=True)

gThedatetime = datetime.now().strftime("%Y%m%d-%H%M%S")
gOutFNameBase = gInitConfFile + '_' + gThedatetime
if gArgs["outFName"]:
    gToSaveExcel = True
    gGlobalExcel = os.path.join(gOutDir, gOutDir + ".xlsx")
    gOutFNameBase = gArgs["outFName"]

printv(gConf)
printv(gArgs)

# exit(0)

# ########## #
# LETS DO IT #
# ########## #

# To avoid total random, we admit setting a specific seed
if gArgs["setRandomSeed"] != -1:
    seed(int(gArgs["setRandomSeed"]))

gWorld, gStatsAnt, gStatsPost = newWorld()

# if "egoism" in gArgs:
#     calcEgoism()

doInitialDistribution()

# print("     %s" % " ".join([gConf["species"][i]["id"] for i in range(gNumberOfSpecies)]))

print(" "*5, end='')
for i in range(gNumberOfSpecies):
    print("%8s" % (gConf["species"][i]["id"],), end='')
print()

for genNumber in range(1, gArgs["numGen"]+1):
    gStatsAnt.fill(0)
    gStatsPost.fill(0)
    doGrouping()
    print("%3d: %s Tot: %6d" % (genNumber, gWorld[:,:].sum(axis=0),  gWorld[:,:].sum(axis=0).sum()))
    doAssociationAndQueue()
    doUngroup()
    doDistribute()

    if gToSaveExcel:
        saveExcel(genNumber)

    saveConf(genNumber)

if gToSaveExcel:
    gExcelWorkbook.close()

