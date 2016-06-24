#!/usr/bin/python

#Copyright (c) 2016, Justin R. Klesmith
#All rights reserved.
#UnselectedDNAStats - Get the codon stats of the unselected libaries

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
__version__ = "1.0, Build: 20160624"
__maintainer__ = "Justin R. Klesmith"
__email__ = "klesmit3@msu.edu"

#Get commandline arguments
parser = argparse.ArgumentParser(description='Unselected DNA stats from Enrich')
parser.add_argument('-p', dest='path', action='store', help='What is the path to the enrich tile directory? ie: ./tile/')
parser.add_argument('-l', dest='tilelength', action='store', help='Tile Length in codons')
parser.add_argument('-t', dest='threshold', action='store', help='Significant threshold')
args = parser.parse_args()

if args.path == None:
    print "No enrich path given"
    quit()

#Global vars
NumResi = int(args.tilelength) #Tile length
Threshold = int(args.threshold)

def CodonCounts():
    codons = {}

    One = 0
    Two = 0
    Three = 0
    
    CountsSingle = 0
    CountsDouble = 0
    CountsTriple = 0
    
    SigSingle = 0
    SigDouble = 0
    SigTriple = 0
    
    TotalUnselCounts = 0
    
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
        ALL = check_output(["awk", 'FNR>1{ print $4,$5,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc'):
        ALL = check_output(["awk", 'FNR>1{ print $4,$5,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc'])
    else:
        print "Counts unsel DNA not found."
        quit()
    
    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m1'):
        M1 = check_output(["awk", 'FNR>1{ print $5,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m1'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m1'):
        M1 = check_output(["awk", 'FNR>1{ print $5,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m1'])
    else:
        print "Counts unsel DNA.m1 not found."
        quit()

    if os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m2'):
        M2 = check_output(["awk", 'FNR>1{ print $5,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m2'])
    elif os.path.isfile(args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m2'):
        M2 = check_output(["awk", 'FNR>1{ print $5,$9 }', args.path+'data/output/counts_unsel_example_F_N_include_filtered_R1_DNA_qc.m2'])
    else:
        print "Counts unsel DNA.m2 not found."
        quit()
    
    #Check for single base mutations
    for line in StringIO.StringIO(M1):
        split = line.split(" ")
        if int(split[0]) >= TranslateStart and int(split[0]) < TranslateEnd: #Check to see that the base is in our tile
            count = int(split[1].rstrip('\n'))
            One = One + 1
            CountsSingle = CountsSingle + count
            if count >= Threshold:
                SigSingle = SigSingle + 1

    #Check for double base mutations
    for line in StringIO.StringIO(M2):
        split2 = line.split(" ")
        location = split2[0].split(",") #Get the individual mutation locations
        
        if int(location[0]) >= TranslateStart and int(location[0]) < TranslateEnd: #Check to see that the base is in our tile
            if int(location[1]) >= TranslateStart and int(location[1]) < TranslateEnd: #Check to see that the base is in our tile
                codon1 = int((int(location[0]) - int(TranslateStart))/3)
                codon2 = int((int(location[1]) - int(TranslateStart))/3)
                if codon1 == codon2:
                    count = int(split2[1].rstrip('\n'))
                    Two = Two + 1
                    CountsDouble = CountsDouble + count
                    if count >= Threshold:
                        SigDouble = SigDouble + 1
    
    #Check for triple base mutations
    for line in StringIO.StringIO(ALL):
        split3 = line.split(" ")
        
        TotalUnselCounts = TotalUnselCounts + int(split3[2].rstrip('\n'))
        
        if split3[0] == "3": #Test to see that there are three mutations
            location = split3[1].split(",") #Get the individual mutation locations
            
            if int(location[0]) >= TranslateStart and int(location[0]) < TranslateEnd: #Check to see that the base is in our tile
                if int(location[1]) >= TranslateStart and int(location[1]) < TranslateEnd: #Check to see that the base is in our tile
                    if int(location[2]) >= TranslateStart and int(location[2]) < TranslateEnd: #Check to see that the base is in our tile
                        codon1 = int((int(location[0]) - int(TranslateStart))/3)
                        codon2 = int((int(location[1]) - int(TranslateStart))/3)
                        codon3 = int((int(location[2]) - int(TranslateStart))/3)
                        if codon1 == codon2 and codon2 == codon3:
                            count = int(split3[2].rstrip('\n'))
                            Three = Three + 1
                            CountsTriple = CountsTriple + count
                            if count >= Threshold:
                                SigTriple = SigTriple + 1

    codons[0] = One #1-base sub
    codons[1] = Two #2-base sub
    codons[2] = Three #3-base sub
    
    codons[3] = CountsSingle #Single Codons
    codons[4] = CountsDouble #Double Codons
    codons[5] = CountsTriple #Triple Codons
    
    codons[6] = SigSingle
    codons[7] = SigDouble
    codons[8] = SigTriple
    
    codons[9] = TotalUnselCounts

    return codons
    
def RunStats():

    codons = CodonCounts()
    print "Unselected DNA sequences (reads) from Enrich: "+str(codons[9])
    
    print "Percent of possible codon subsititions observed in the unselected population:"
    print "1-base substitution (#codons*9): {0:.1f}".format((codons[0]/(9*NumResi)*100))+"% "+str(codons[0])+"/"+str(9*NumResi)
    print "2-base substitutions (#codons*27): {0:.1f}".format((codons[1]/(27*NumResi)*100))+"% "+str(codons[1])+"/"+str(27*NumResi)
    print "3-base substitutions (#codons*27): {0:.1f}".format((codons[2]/(27*NumResi)*100))+"% "+str(codons[2])+"/"+str(27*NumResi)
    print "Total base substitutions: {0:.1f}".format((codons[0]+codons[1]+codons[2])/(63*NumResi)*100)+"% "+str(codons[0]+codons[1]+codons[2])+"/"+str(63*NumResi)
    
    print "Percent of reads of mutant codons with:"
    
    print "1-base: {0:.1f}".format((codons[3])/(codons[3]+codons[4]+codons[5])*100)+"% "
    print "2-base: {0:.1f}".format((codons[4])/(codons[3]+codons[4]+codons[5])*100)+"% "
    print "3-base: {0:.1f}".format((codons[5])/(codons[3]+codons[4]+codons[5])*100)+"% " 
    
    print "Significant coverage of possible programmed mutant codons: {0:.1f}".format((codons[6]+codons[7]+codons[8])/(63*NumResi)*100)+"% " 

    return
    
def main():
    #Write out preample
    print "Stats for unselected DNA libaries from Enrich"
    print __author__
    print "Contact: "+__email__
    print __copyright__
    print "License: "+__license__
    print "Credits: "+__credits__[0]+", "+__credits__[1]
    
    print "Path: "+args.path
    print "Tile Length: "+args.tilelength
    print "Significant count threshold: "+args.threshold
    
    #Print out run stats
    RunStats()

if __name__ == '__main__':
    main()