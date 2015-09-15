#!/usr/bin/python

#Copyright (c) 2015, Justin R. Klesmith
#All rights reserved.
#QuickNormalize: Normalize tiles from a growth or FACS selection

from __future__ import division
from subprocess import check_output
from math import log, sqrt, pow, e
from scipy import special
import numpy as np
import StringIO
import argparse
import time
import os

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2015, Justin R. Klesmith"
__credits__ = ["Justin R. Klesmith", "Jim Stapleton", "Timothy A. Whitehead"]
__license__ = "BSD-3"
__version__ = "1.9, Build: 20150615"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["klesmit3@msu.edu", "justinklesmith@gmail.com", "justinklesmith@evodyn.com"]

#Get commandline arguments
parser = argparse.ArgumentParser(description='QuickNormalize for Growth or FACS')
parser.add_argument('-u', dest='unselected', action='store', required=True, help='Unselected File')
parser.add_argument('-s', dest='selected', action='store', required=True, help='Selected File')
parser.add_argument('-c', dest='selectedwtcounts', action='store', required=True, help='Selected Wild-Type Counts')
parser.add_argument('-d', dest='unselectedwtcounts', action='store', required=True, help='Unselected Wild-Type Counts')
parser.add_argument('-t', dest='sigthreshold', action='store', nargs='?', const=1, default=5, help='Unselected counts for significance. Default = 5')
parser.add_argument('-w', dest='wildtype', action='store', nargs='?', const=1, default='./WTSeq', help='File with the wild-type amino acid sequence. Default = ./WTSeq')
args = parser.parse_args()

#Global Variables
if os.path.isfile(args.wildtype):
    with open(args.wildtype, 'r') as infile: #Open the file with the wild-type protein sequence
        WTSeq = infile.readline() #Read the first line of the WT sequence file
else:
    print "Wild-type sequence file not found...exit"
    quit()

SignificantThreshold = int(args.sigthreshold) #Number of counts in the unselected library and selected library to be significant
AA_Table = '*FWYPMILVAGCSTNQDEHKR'

Mutations = {} #Mutations matrix 
#UnselectedWTCounts = 31017421
UnselectedWTCounts = args.unselectedwtcounts
TotalSelectedCounts = 0
TotalUnselectedCounts = 0
TileLen = len(WTSeq)

def Build_Matrix():
    #Populate mutation matrix with None data
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            try:
                #Mutations[ResID][MutID[1]][0 = RawLog2, 1 = Log2-AddedDepletes, 2 = Unselected, 3 = Selected, ]
                Mutations[j][i[1]] = [None, None, None, None]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None, None, None]

    return Mutations
  
def Get_Unsel_Counts():
    global TotalUnselectedCounts

	#Get the unselected counts for a variant
    with open(args.unselected) as infile:
        for line in infile:
            split = line.split(",")
            ID = split[0].split("-")
            location = int(ID[0])
            identity = str(ID[1])
            counts = int(split[1].rstrip('\n'))
            Mutations[location][identity][2] = counts #Set the unselected counts
            TotalUnselectedCounts = TotalUnselectedCounts + counts
	
    return Mutations
	
def Get_Sel_Counts():
    global TotalSelectedCounts

	#Get the unselected counts for a variant
    with open(args.selected) as infile:
        for line in infile:
            split = line.split(" ")
            ID = split[0].split("-")
            location = int(ID[0])
            identity = str(ID[1])
            counts = int(split[1].rstrip('\n'))
            Mutations[location][identity][3] = counts #Set the selected counts
            TotalSelectedCounts = TotalSelectedCounts + counts

    return Mutations

def Enrich():
    print "Unselected Wild-Type Counts: "+str(UnselectedWTCounts)
    print "Selected Wild-Type Counts: "+str(args.selectedwtcounts)
    print "Total Unselected Counts: "+str(UnselectedWTCounts+TotalUnselectedCounts)
    print "Total Selected Counts: "+str(int(args.selectedwtcounts)+TotalSelectedCounts)
    
    #Calculate WT enrichment
    SelectedWT = int(args.selectedwtcounts)/(int(args.selectedwtcounts)+TotalSelectedCounts)
    UnselectedWT = UnselectedWTCounts/(UnselectedWTCounts+TotalUnselectedCounts)
    print "Wild-type enrichment: "+str(log((SelectedWT/UnselectedWT),2))
    
    print "Enriching the data - Above the unselected threshold and the selected read exists"
    print "Location,Mutation,Unselected_Reads,Selected_Reads,Enrichment"
    
    #Calculate the enrichment if values both exist and unsel is above threshold
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            if Mutations[j][i[1]][2] >= SignificantThreshold and Mutations[j][i[1]][3] != None:
                Selected = Mutations[j][i[1]][3]/(int(args.selectedwtcounts)+TotalSelectedCounts)
                Unselected = Mutations[j][i[1]][2]/(UnselectedWTCounts+TotalUnselectedCounts)
                log2 = log((Selected/Unselected),2)
                Mutations[j][i[1]][0] = log2

            print str(j)+","+i[1]+","+str(Mutations[j][i[1]][2])+","+str(Mutations[j][i[1]][3])+","+str(Mutations[j][i[1]][0])    
    
    #print "Enriching the data - Above 2x the unselected threshold and 1 count in selected for depleted values"    
    #print "Location,Mutation,Unselected_Reads,Selected_Reads,Enrichment"
    
    #Calculate the enrichment if values both exist and unsel is above threshold
    #for j in xrange(0,TileLen):
        #for i in enumerate(AA_Table):
            #if Mutations[j][i[1]][2] >= SignificantThreshold*2:
                #if Mutations[j][i[1]][3] == None:
                    #Selected = 1/(int(args.selectedwtcounts)+TotalSelectedCounts)
                    #Unselected = Mutations[j][i[1]][2]/(UnselectedWTCounts+TotalUnselectedCounts)
                    #log2d = log((Selected/Unselected),2)
                    #Mutations[j][i[1]][1] = log2d

            #print str(j)+","+i[1]+","+str(Mutations[j][i[1]][2])+","+str(Mutations[j][i[1]][3])+","+str(Mutations[j][i[1]][1])
    
    return Mutations

def Make_CSV():
    print "Enrichment Heatmap"
    #This makes a CSV style report of rows of letters and columns of residues
    
    #Print off the Number
    Numbering = " "
    for q in xrange(1,TileLen):
        Numbering = Numbering+","+str(q)
    print Numbering
        
    #Print off the WT Residue
    WTResi = " "
    for w in xrange(0,TileLen):
        WTResi = WTResi+","+WTSeq[w]
    print WTResi
    
    #Print off the mutations
    Output = ""
    for i in enumerate(AA_Table):
        Output = Output+i[1]+","
        for j in xrange(0,TileLen-1):
            Output = Output+str(Mutations[j][i[1]][0])+","
        Output = Output+"\n"
    print Output
    
    #Print off the mutations
    #Output2 = ""
    #for i in enumerate(AA_Table):
        # = Output2+i[1]+","
        #for j in xrange(0,TileLen-1):
            #Output2 = Output2+str(Mutations[j][i[1]][1])+","
        #Output2 = Output2+"\n"
    #print Output2
    
    #Write the heatmap to a newfile
    outfile = open('heatmap_'+args.selected+'.csv', 'w')
    outfile.write(Numbering+'\n')
    outfile.write(WTResi+'\n')
    outfile.write(Output)
    #outfile.write(Output2)
    
    return

def main():
    #Write out preamble
    print "QuickEnrich"
    print "Author: "+__author__
    print "Contact: "+__email__[0]+", "+__email__[1]+", "+__email__[2]
    print __copyright__
    print "Version: "+__version__
    print "License: "+__license__
    print "Credits: "+__credits__[0]+", "+__credits__[1]+", "+__credits__[2]
    print ""
    print "Please cite:"
    print "Github [user: JKlesmith] (www.github.com)"
    print ""
    print "Enrich run parameters:"
    print time.strftime("%H:%M:%S")
    print time.strftime("%m/%d/%Y")
    print "Unselected counts to be significant: "+str(args.sigthreshold)
    print "Wild-type sequence file: "+args.wildtype

    #Build Matrix
    Build_Matrix()
    
    #Get the selected counts
    Get_Unsel_Counts()
    Get_Sel_Counts()

    #Normalize the Data
    Enrich()

    #Print out a csv
    Make_CSV()

if __name__ == '__main__':
    main()