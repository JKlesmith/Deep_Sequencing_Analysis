#!/usr/bin/python

#Copyright (c) 2016, Justin R. Klesmith
#All rights reserved.
#QuickStats: Get the statistics from a enrich run

from __future__ import division
from subprocess import check_output
from math import log
import StringIO
import argparse
import time
import os

__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2016, Justin R. Klesmith"
__credits__ = ["Justin R. Klesmith", "Timothy A. Whitehead"]
__license__ = "BSD-3"
__version__ = "1.4X, Build: 201507X"
__maintainer__ = "Justin R. Klesmith"
__email__ = "klesmit3@msu.edu"

#Build Notes:
#1.3 - 20150616 - Fixed counting bug at the end of the tile in CodonSubs such that it's just less than and not equal to

#Get commandline arguments
parser = argparse.ArgumentParser(description='Quick Enrich Stats - Note: you must pre-normalize the data using QuickNormalize.py.')
parser.add_argument('-f', dest='file', action='store', help='File of your already normalized dataset')
parser.add_argument('-p', dest='path', action='store', help='What is the path to the enrich tile directory? ie: ./tile/')
args = parser.parse_args()

#Verify inputs
if args.file == None:
    print "No normalized file given"
    quit()

if args.path == None:
    print "No enrich path given"
    quit()

#Global vars
AA_Table = '*ACDEFGHIKLMNPQRSTVWY'
Mutations = {}
NumResi = None #Tile length
NormData = ""
StartResidue = None

def Build_Matrix():
    #Populate Mutation Dictionary with None Data
    for j in xrange(0+StartResidue,NumResi+StartResidue):
        for i in enumerate(AA_Table):
            try:
                #Mutations[ResID][MutID[1]][0 = NormLog2, 1 = Unselected, 2 = Selected]
                Mutations[j][i[1]] = [None, None, None]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None, None]

    return Mutations

def ImportNormData():
    global NumResi
    global NormData
    global StartResidue
    
    lines = 0
    normdata = ""

    #Import the previously normalized data
    with open(args.file) as infile:
        copy = False
        for line in infile:
            if line.strip() == "Location,Mutation,Normalized_ER,Unselected_Reads,Selected_Reads,RawLog2":
                copy = True
            elif line.strip() == "Normalized Heatmap":
                copy = False
            elif line.startswith("Tile Length: "):
                print "Tile length: "+line.strip()[13:]
                NumResi = int(line.strip()[13:])
            elif line.startswith("Start residue (-s): "):
                split = line.split(" ")
                StartResidue = int(split[3]) #Set the start residue
            elif copy:
                NormData = NormData + line
                lines = lines + 1
    
    #NumResi = int(lines / 21) #Set the tile length
    return normdata

def PopulateMutArrays():
    #Loop through the output
    for line in StringIO.StringIO(NormData):
        split = line.split(",")
        
        location = int(split[0])
        identity = str(split[1])
        
        Mutations[location][identity][0] = split[2]
        Mutations[location][identity][1] = split[3]
        Mutations[location][identity][2] = split[4].rstrip('\n')

    return Mutations

def DNAReads():
    reads = {} #Initialize the variable for the number of reads 0=unsel, 1=sel
    
    SC = 0
    UC = 0
    selectedcounts = ""
    unselectedcounts = ""
    
    if os.path.isfile(args.path+'data/output/counts_sel_example_F_N_include_filtered_B_DNA_qc'):
        selectedcounts = check_output(["awk", 'FNR>1{ print $9 }', args.path+'data/output/counts_sel_example_F_N_include_filtered_B_DNA_qc'])
    elif os.path.isfile(args.path+'data/output/counts_sel_example_F_N_include_filtered_R1_DNA_qc'):
        selectedcounts = check_output(["awk", 'FNR>1{ print $9 }', args.path+'data/output/counts_sel_example_F_N_include_filtered_R1_DNA_qc'])
    else:
        print "Can't find selected DNA counts"
        quit()
    
    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc'):
        unselectedcounts = check_output(["awk", 'FNR>1{ print $9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc'):
        unselectedcounts = check_output(["awk", 'FNR>1{ print $9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc'])
    else:
        print "Can't find unselected DNA counts"
        quit()
    
    #Loop through the output
    for line in StringIO.StringIO(selectedcounts):
        split = line.split(" ")
        SC = SC + int(split[0].rstrip('\n'))

    for line in StringIO.StringIO(unselectedcounts):
        split = line.split(" ")
        UC = UC + int(split[0].rstrip('\n'))
    
    reads[0] = str(UC) #Set the unselected reads
    reads[1] = str(SC) #Set the selected reads

    return reads
    
def MutationCounts():
    muts = {}
    
    NM00 = 0
    NM10 = 0
    NM15 = 0
    NM30 = 0
    NM50 = 0
    NM100 = 0
    
    FiveThreshold = 0
    Retained = 0
    
    for j in xrange(0+StartResidue,NumResi+StartResidue):
        for i in enumerate(AA_Table):
            if Mutations[j][i[1]][0] != "NS":
                if float(Mutations[j][i[1]][0]) > 0.00:
                    NM00 += 1
                
                if float(Mutations[j][i[1]][0]) > 0.10:
                    NM10 += 1
                
                if float(Mutations[j][i[1]][0]) > 0.15:
                    NM15 += 1            
                
                if float(Mutations[j][i[1]][0]) > 0.30:
                    NM30 += 1                
                
                if float(Mutations[j][i[1]][0]) > 0.50:
                    NM50 += 1
                
                if float(Mutations[j][i[1]][0]) > 1.00:
                    NM100 += 1
            
            
            if Mutations[j][i[1]][1] != "None":
                if int(Mutations[j][i[1]][1]) >= 5:
                    FiveThreshold += 1
                    
            if Mutations[j][i[1]][2] != "None":
                Retained += 1
    
    muts[0] = NM00
    muts[1] = NM10
    muts[2] = NM15
    muts[3] = NM30
    muts[4] = NM50
    muts[5] = NM100

    muts[6] = FiveThreshold
    muts[7] = Retained

    return muts

def Nonsynonymous():
    reads = {}
    
    Total = 0
    Single = 0
    WT = 0
    
    ALL = ""
    M1 = ""

    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_PRO_qc'):
        ALL = check_output(["awk", 'FNR>1{ print $1,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_PRO_qc'])    
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_PRO_qc'):
        ALL = check_output(["awk", 'FNR>1{ print $1,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_PRO_qc'])    
    else:
        print "Unsel protein counts not found"
        quit()
    
    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_PRO_qc.m1'):
        M1 = check_output(["awk", 'FNR>1{ print $9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_PRO_qc.m1'])    
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_PRO_qc.m1'):
        M1 = check_output(["awk", 'FNR>1{ print $9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_PRO_qc.m1'])
    else:
        print "Unsel protein counts.m1 not found"
        quit()

    #Loop through the output
    for line in StringIO.StringIO(ALL):
        split = line.split(" ")
        
        if split[0] == "NA-NA":
            WT = int(split[1])

        Total = Total + int(split[1].rstrip('\n'))

    for line in StringIO.StringIO(M1):
        split = line.split(" ")
        Single = Single + int(split[0].rstrip('\n'))

    reads[0] = WT #Wild-type
    reads[1] = Single #.m1
    reads[2] = (Total - Single - WT) #all - .m1 - WT
    return reads

def CodonSubs():
    codons = {}

    One = 0
    Two = 0
    Three = 0
    
    #Get the start of translation
    TranslateStart = 0
    TranslateEnd = 0
    with open(args.path+'input/example_local_config') as infile:
        for line in infile:
            if line.startswith("<translate_start>"):
                TSLen = len(line)
                TranslateStart = int(line[17:TSLen-20])
                TranslateEnd = TranslateStart+(3*NumResi)

    ALL = ""
    M1 = ""
    M2 = ""
    
    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc'):
        ALL = check_output(["awk", 'FNR>1{ print $4,$5 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc'):
        ALL = check_output(["awk", 'FNR>1{ print $4,$5 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc'])
    else:
        print "Counts unsel DNA not found."
        quit()
    
    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m1'):
        M1 = check_output(["awk", 'FNR>1{ print $5 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m1'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m1'):
        M1 = check_output(["awk", 'FNR>1{ print $5 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m1'])
    else:
        print "Counts unsel DNA.m1 not found."
        quit()

    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m2'):
        M2 = check_output(["awk", 'FNR>1{ print $5 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m2'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m2'):
        M2 = check_output(["awk", 'FNR>1{ print $5 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m2'])
    else:
        print "Counts unsel DNA.m2 not found."
        quit()
    
    #Check for single base mutations
    for line in StringIO.StringIO(M1):
        split = line.split(" ")
        if int(split[0]) >= TranslateStart and int(split[0]) < TranslateEnd: #Check to see that the base is in our tile
            One = One + 1

    #Check for double base mutations
    for line in StringIO.StringIO(M2):
        split2 = line.split(" ")
        location = split2[0].split(",") #Get the individual mutation locations
        
        if int(location[0]) >= TranslateStart and int(location[0]) < TranslateEnd: #Check to see that the base is in our tile
            if int(location[1]) >= TranslateStart and int(location[1]) < TranslateEnd: #Check to see that the base is in our tile
                codon1 = int((int(location[0]) - int(TranslateStart))/3)
                codon2 = int((int(location[1]) - int(TranslateStart))/3)
                if codon1 == codon2:
                    Two = Two + 1
    
    #Check for triple base mutations
    for line in StringIO.StringIO(ALL):
        split3 = line.split(" ")
        
        if split3[0] == "3": #Test to see that there are three mutations
            location = split3[1].split(",") #Get the individual mutation locations
            
            if int(location[0]) >= TranslateStart and int(location[0]) < TranslateEnd: #Check to see that the base is in our tile
                if int(location[1]) >= TranslateStart and int(location[1]) < TranslateEnd: #Check to see that the base is in our tile
                    if int(location[2]) >= TranslateStart and int(location[2]) < TranslateEnd: #Check to see that the base is in our tile
                        codon1 = int((int(location[0]) - int(TranslateStart))/3)
                        codon2 = int((int(location[1]) - int(TranslateStart))/3)
                        codon3 = int((int(location[2]) - int(TranslateStart))/3)
                        if codon1 == codon2 and codon2 == codon3:
                            Three = Three + 1

    codons[0] = One #1-base sub
    codons[1] = Two #2-base sub
    codons[2] = Three #3-base sub

    return codons
    
def RunStats():
    print "Stat run parameters:"
    print time.strftime("%H:%M:%S")
    print time.strftime("%m/%d/%Y")
    print "Nomalized file: "+args.file
    print "Data path: "+args.path

    reads = DNAReads()
    print "Unselected DNA sequences (reads) from Enrich: "+reads[0]
    print "Selected DNA sequences (reads) from Enrich: "+reads[1]
    
    mutations = MutationCounts()
    print "Number of mutations above 0.00: "+str(mutations[0])
    print "Number of mutations above 0.10: "+str(mutations[1])
    print "Number of mutations above 0.15: "+str(mutations[2])
    print "Number of mutations above 0.30: "+str(mutations[3])
    print "Number of mutations above 0.50: "+str(mutations[4])
    print "Number of mutations above 1.00: "+str(mutations[5])

    print "Number unselected mutants above threshold of 5: "+str(mutations[6])
    print "Number of mutations retained in the selected population (not given a 1 if significant in unsel): "+str(mutations[7])
    
    codons = CodonSubs()
    print "Percent of possible codon subsititions observed in the unselected population:"
    print "1-base substitution (#codons*9): {0:.1f}".format((codons[0]/(9*NumResi)*100))+"% "+str(codons[0])+"/"+str(9*NumResi)
    print "2-base substitutions (#codons*27): {0:.1f}".format((codons[1]/(27*NumResi)*100))+"% "+str(codons[1])+"/"+str(27*NumResi)
    print "3-base substitutions (#codons*27): {0:.1f}".format((codons[2]/(27*NumResi)*100))+"% "+str(codons[2])+"/"+str(27*NumResi)
    print "Total base substitutions: "+str(codons[0]+codons[1]+codons[2])+"/"+str(63*NumResi)
    
    nonsynonymous = Nonsynonymous()
    print "Percent of unselected reads with: "
    print "No nonsynonymous mutations: {0:.1f}".format((nonsynonymous[0]/int(reads[0]))*100)+"% "+str(nonsynonymous[0])+"/"+reads[0]
    print "One nonsynonymous mutation: {0:.1f}".format((nonsynonymous[1]/int(reads[0]))*100)+"% "+str(nonsynonymous[1])+"/"+reads[0]
    print "Multiple nonsynonymous mutations: {0:.1f}".format((nonsynonymous[2]/int(reads[0]))*100)+"% "+str(nonsynonymous[2])+"/"+reads[0]
    
    print "Coverage of possible single nonsynonymous amino acid mutations: {0:.1f}".format((mutations[6]/(NumResi*20))*100)+"% "+str(mutations[6])+"/"+str(NumResi*20)
    return
    
def main():
    #Write out preample
    print "QuickStats"
    print "Author: "+__author__
    print "Contact: "+__email__
    print __copyright__
    print "License: "+__license__
    print "Credits: "+__credits__[0]+", "+__credits__[1]
    print ""
    print "Please cite:"
    print "Github [user: JKlesmith] (www.github.com)"
    print "Klesmith JR, Bacik J-P, Michalczyk R, Whitehead TA. 2015. Comprehensive sequence-flux mapping of metabolic pathways in living cells."
    print "Kowalsky CA, Klesmith JR, Stapleton JA, Kelly V, Reichkitzer N, Whitehead TA. 2015. High-Resolution Sequence-Function Mapping of Full-Length Proteins. PLoS ONE 10(3):e0118193. doi:10.1371/journal.pone.0118193."
    print ""
    
    #Print out run stats
    ImportNormData()
    Build_Matrix()
    PopulateMutArrays()
    RunStats()

if __name__ == '__main__':
    main()