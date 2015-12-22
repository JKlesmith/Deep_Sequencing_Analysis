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
__credits__ = ["Justin R. Klesmith", "Caitlin A. Kowalsky", "Timothy A. Whitehead"]
__license__ = "BSD-3"
__version__ = "2.0, Build: 20150819"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["klesmit3@msu.edu", "justinklesmith@gmail.com", "justinklesmith@evodyn.com"]

#Get commandline arguments
parser = argparse.ArgumentParser(description='QuickNormalize '+__version__+' for Growth or FACS')
parser.add_argument('-n', dest='normtype', action='store', required=True, help='Normalization Type? Enter: growth or FACS')
parser.add_argument('-s', dest='startresidue', action='store', required=True, help='What is the start residue? ie: 0, 40, 80')
parser.add_argument('-l', dest='length', action='store', required=True, help='Length of your tile? ie: 40, 80')
parser.add_argument('-g', dest='gp', action='store', help='How many doublings/generations? (GROWTH) ie: 12.5')
parser.add_argument('-d', dest='stddev', action='store', help='Standard Deviation? (FACS) ie: 0.6')
parser.add_argument('-c', dest='percentcollected', action='store', help='Percent Collected? (FACS) ie: 0.05')
parser.add_argument('-p', dest='path', action='store', required=True, help='What is the path to the enrich output directory? ie: ./tile/data/output/')
parser.add_argument('-t', dest='sigthreshold', action='store', nargs='?', const=1, default=5, help='Unselected counts for significance. Default = 5')
parser.add_argument('-w', dest='wildtype', action='store', nargs='?', const=1, default='./WTSeq', help='File with the wild-type amino acid sequence. Default = ./WTSeq')
parser.add_argument('-o', dest='heatmap', action='store', nargs='?', const=1, default='True', help='Output a csv heatmap? Default = True') 
parser.add_argument('-y', dest='ewtenrichment', action='store', help='Manual Ewt enrichment value')
parser.add_argument('-z', dest='eiscalar', action='store', help='Manual Ei enrichment scalar')
args = parser.parse_args()

#Verify inputs
if args.normtype != "growth" and args.normtype != "FACS":
    print "Missing normalization type. Flag: -n"
    quit()
    
if args.startresidue == None:
    print "Missing start residue. Flag: -s"
    quit()

if args.length == None:
    print "Missing tile length. Flag: -l"
    quit()

if args.gp == None and args.normtype == "growth":
    print "Missing doublings. Flag: -g"
    quit()

if args.stddev == None and args.normtype == "FACS":
    print "Missing SD. Flag: -d"
    quit()

if args.percentcollected == None and args.normtype == "FACS":
    print "Missing percent collected. Flag: -c"
    quit() 

if args.path == None:
    print "Missing Enrich output path. Flag: -p"
    quit()

if args.ewtenrichment and args.eiscalar != None:
    #This section is only true if we want to provide our own WT enrichment and a scalar to add to Ei
    OverrideEwtEi = True
    ManualEwt = float(args.ewtenrichment)
    EiScalar = float(args.eiscalar)
else:
    OverrideEwtEi = False

#Global Variables
if os.path.isfile(args.wildtype):
    with open(args.wildtype, 'r') as infile: #Open the file with the wild-type protein sequence
        WTSeq = infile.readline() #Read the first line of the WT sequence file
else:
    print "Wild-type sequence file not found...exit"
    quit()

StartResidue = int(args.startresidue) #Starting residue for your tile
TileLen = int(args.length) #Length of your tile
Path = args.path #What is the path to the output directory
SignificantThreshold = int(args.sigthreshold) #Number of counts in the unselected library and selected library to be significant

if args.normtype == "growth":
    DoublingsGp = float(args.gp) #Number of doublings

if args.normtype == "FACS":
    SD = float(args.stddev) #Standard Deviation
    PC = float(args.percentcollected) #Percent collected
    THEOENRICHMENT = -log(PC, 2) #Theoretical maximum enrichment

#AA_Table = '*ACDEFGHIKLMNPQRSTVWY'
AA_Table = '*FWYPMILVAGCSTNQDEHKR'

Mutations = {} #Mutations matrix
Ewt = None #Initialize the variable for the wildtype enrichment

######################################################################################
#
#Main Program Functions
#
######################################################################################

######################################################################################
#Build_Matrix
#This does the initial population of the Mutations matrix that holds counts and
#enrichment values
######################################################################################

def Build_Matrix():
    #Populate mutation matrix with None data
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            try:
                #Mutations[ResID][MutID[1]][0 = RawLog2, 1 = Fitness, 2 = Unselected, 3 = Selected, 4=Unused, 5=WT]
                Mutations[j][i[1]] = [None, None, None, None, None, False]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None, None, None, None, False]

    return Mutations

######################################################################################
#Get_WT_Ewt
#This gets the wild-type enrichment from the enrich NA-NA output
######################################################################################
 
def Get_WT_Ewt():
    global Ewt
    #Extract NA-NA WT Ewt log2
    
    awk = ""
    
    if os.path.isfile(Path+'ratios_sel_example_F_N_include_filtered_B_PRO_qc_unsel_example_F_N_include_filtered_B_PRO_qc'):
        awk = check_output(["awk", '{ print $5,$6,$8 }', Path+'ratios_sel_example_F_N_include_filtered_B_PRO_qc_unsel_example_F_N_include_filtered_B_PRO_qc'])
    elif os.path.isfile(Path+'ratios_sel_example_F_N_include_filtered_R1_PRO_qc_unsel_example_F_N_include_filtered_R1_PRO_qc'):
        awk = check_output(["awk", '{ print $5,$6,$8 }', Path+'ratios_sel_example_F_N_include_filtered_R1_PRO_qc_unsel_example_F_N_include_filtered_R1_PRO_qc'])
    else:
        print "Selected protein ratios file not found...exit"
        quit()
    
    #Loop through the output
    for line in StringIO.StringIO(awk):
        split = line.split(" ")
        location = str(split[0])
        identity = str(split[1])
        if location == "NA" and identity == "NA":
            Ewt = float(split[2].rstrip('\n'))
            print "Wild-type log2 (Ewt): "+str(Ewt)

    return Ewt

######################################################################################
#Get_Mut_Ei
#This gets the enrichment of each mutation from the enrich output
######################################################################################    
    
def Get_Mut_Ei():
    #Extract Mut Ei log2
    
    awk = ""
    if os.path.isfile(Path+'ratios_sel_example_F_N_include_filtered_B_PRO_qc_unsel_example_F_N_include_filtered_B_PRO_qc.m1'):
        awk = check_output(["awk", 'FNR>1{ print $5,$6,$8 }', Path+'ratios_sel_example_F_N_include_filtered_B_PRO_qc_unsel_example_F_N_include_filtered_B_PRO_qc.m1'])
    elif os.path.isfile(Path+'ratios_sel_example_F_N_include_filtered_R1_PRO_qc_unsel_example_F_N_include_filtered_R1_PRO_qc.m1'):
        awk = check_output(["awk", 'FNR>1{ print $5,$6,$8 }', Path+'ratios_sel_example_F_N_include_filtered_R1_PRO_qc_unsel_example_F_N_include_filtered_R1_PRO_qc.m1'])
    else:
        print "Selected protein ratios .m1 file not found...exit"
        quit()
    
    #Loop through the output
    for line in StringIO.StringIO(awk):
        split = line.split(" ")
        location = int(split[0])
        identity = str(split[1])
        Ei = float(split[2].rstrip('\n'))
        
        #Check to see if we're above the tile length and go to next
        if location >= TileLen:
            continue
        
        #For FACS set a upper limit on enrichment, don't do anything for growth
        if args.normtype == "FACS":
            #Check to see if the enrichment is greater or equal than the theoretical
            if OverrideEwtEi == False: #Apply no scalar to the Ei
                if Ei >= THEOENRICHMENT:
                    Mutations[location][identity][0] = (THEOENRICHMENT - 0.001)
                else:
                    Mutations[location][identity][0] = Ei
            elif OverrideEwtEi == True: #Apply a scalar to the Ei
                if Ei >= (THEOENRICHMENT + EiScalar):
                    Mutations[location][identity][0] = ((THEOENRICHMENT + EiScalar) - 0.001)
                else:
                    Mutations[location][identity][0] = (Ei + EiScalar)
            
        elif args.normtype == "growth":
            Mutations[location][identity][0] = Ei

    return Mutations

######################################################################################
#Get_Unsel_Counts
#This gets the unselected counts for each mutation
######################################################################################

def Get_Unsel_Counts():
	#Get the unselected counts for a variant
    
    awk = ""
    if os.path.isfile(Path+'counts_unsel_example_F_N_include_filtered_B_PRO_qc.m1'):
        awk = check_output(["awk", 'FNR>1{ print $5,$6,$9 }', Path+'counts_unsel_example_F_N_include_filtered_B_PRO_qc.m1'])
    elif os.path.isfile(Path+'counts_unsel_example_F_N_include_filtered_R1_PRO_qc.m1'):
        awk = check_output(["awk", 'FNR>1{ print $5,$6,$9 }', Path+'counts_unsel_example_F_N_include_filtered_R1_PRO_qc.m1'])
    else:
        print "Unselected protein counts .m1 file not found...exit"
        quit()
        
    #Loop through the output
    for line in StringIO.StringIO(awk):
        split = line.split(" ")
        location = int(split[0])
        identity = str(split[1])
        counts = int(split[2].rstrip('\n'))
        
        #Check to see if we're above the tile length and go to next
        if location >= TileLen:
            continue
        
        Mutations[location][identity][2] = counts #Set the unselected counts
	
    return Mutations

######################################################################################
#Get_Sel_Counts
#This gets the selected counts for each mutation
######################################################################################

def Get_Sel_Counts():
	#Get the selected counts

    awk = ""
    if os.path.isfile(Path+'counts_sel_example_F_N_include_filtered_B_PRO_qc.m1'):
        awk = check_output(["awk", 'FNR>1{ print $5,$6,$9 }', Path+'counts_sel_example_F_N_include_filtered_B_PRO_qc.m1'])
    elif os.path.isfile(Path+'counts_sel_example_F_N_include_filtered_R1_PRO_qc.m1'):
        awk = check_output(["awk", 'FNR>1{ print $5,$6,$9 }', Path+'counts_sel_example_F_N_include_filtered_R1_PRO_qc.m1'])
    else:
        print "Selected protein counts .m1 file not found...exit"
        quit()
    
    #Loop through the output
    for line in StringIO.StringIO(awk):
        split = line.split(" ")
        location = int(split[0])
        identity = str(split[1])
        counts = int(split[2].rstrip('\n'))
        
        #Check to see if we're above the tile length
        if location >= TileLen:
            continue
        
        Mutations[location][identity][3] = counts #Set the selected counts

    return Mutations

######################################################################################
#Normalize
#This normalizes the enrichments to the wild-type using the fitness metric equations
######################################################################################

def Normalize():
    #Check to see if the wild-type enrichment is set
    if Ewt == None:
        print "Error: Wild-Type enrichment is not set...quit"
        quit()

    print ""
    print "Normalizing the data"
    print "Location,Mutation,Normalized_ER,Unselected_Reads,Selected_Reads,RawLog2"
    
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
        	#Check for a case where a significant variant fell out of the population
            if Mutations[j][i[1]][0] == None and Mutations[j][i[1]][2] >= SignificantThreshold and Mutations[j][i[1]][3] == None:
                Mutations[j][i[1]][0] = log((1/Mutations[j][i[1]][2]), 2) #Calculate the raw log2 for this variant and report it as less than this value
        
            #Calculate the fitness
            if Mutations[j][i[1]][0] != None and Mutations[j][i[1]][2] >= SignificantThreshold: #Report the significant fitness
                Ei = float(Mutations[j][i[1]][0])
                
                if args.normtype == "growth":
                    Mutant = (Ei/DoublingsGp)+1
                    WT = (Ewt/DoublingsGp)+1
                    if (Mutant/WT) < 0:
                        NE = -10 #Assign an extremely negative fitness for members who are greather than -8 raw log2 enrichment
                    else:
                        NE = log(Mutant/WT, 2)
                elif args.normtype == "FACS":              
                    WT = special.erfinv(1-PC*pow(2,(Ewt+1)))
                    Mutant = special.erfinv(1-PC*pow(2,(Ei+1)))
                    NE = (log(e, 2)*sqrt(2)*SD*(WT-Mutant))
                else:
                    print "Error: growth or FACS not set?"
                    quit()
                
                Mutations[j][i[1]][1] = "{0:.4f}".format(NE)
            elif Mutations[j][i[1]][2] < SignificantThreshold: #Report the insignificant NEs
                if WTSeq[j+StartResidue] == i[1]: #Check to see if it's wildtype else it's Not Significant
                    Mutations[j][i[1]][1] = "0.000"
                    Mutations[j][i[1]][5] = True #Set the WT flag
            	else:
                    Mutations[j][i[1]][1] = "NS" 
            elif Mutations[j][i[1]][2] == None and Mutations[j][i[1]][3] >= SignificantThreshold: #Error: Mutation with selected counts and no unselected
                Mutations[j][i[1]][1] = "Error: Sel with Zero Unsel"
            else:
                print "Error: unknown normalization problem."
                
            #Print out column data
            print str(j+StartResidue)+","+i[1]+","+Mutations[j][i[1]][1]+","+str(Mutations[j][i[1]][2])+","+str(Mutations[j][i[1]][3])+","+str(Mutations[j][i[1]][0])
    
    return Mutations

######################################################################################
#Make_CSV
#This outputs the fitness metrics to a CSV file to be imported into excel
######################################################################################
    
def Make_CSV():
    print "Normalized Heatmap"
    #This makes a CSV style report of rows of letters and columns of residues
    
    #Print off the Number
    Numbering = " "
    for q in xrange(1,TileLen+1):
        Numbering = Numbering+","+str(StartResidue+q)
    print Numbering
        
    #Print off the WT Residue
    WTResi = " "
    for w in xrange(0,TileLen):
        WTResi = WTResi+","+WTSeq[StartResidue+w]
    print WTResi
    
    #Print off the mutations
    Output = ""
    for i in enumerate(AA_Table):
        Output = Output+i[1]+","
        for j in xrange(0,TileLen):
            Output = Output+str(Mutations[j][i[1]][1])+","
        Output = Output+"\n"
    print Output
    
    if args.heatmap == "True":
        #Write the heatmap to a newfile
        outfile = open('heatmap_startresi_'+str(StartResidue)+'.csv', 'w')
        outfile.write(Numbering+'\n')
        outfile.write(WTResi+'\n')
        outfile.write(Output)
    
    return

######################################################################################
#main
#This is the main function that calls the sub-functions, it also outputs the run
#information including the command line parameters
######################################################################################    
    
def main():
    global Ewt
    
    #Write out preamble
    print "QuickNormalize"
    print "Author: "+__author__
    print "Contact: "+__email__[0]+", "+__email__[1]+", "+__email__[2]
    print __copyright__
    print "Version: "+__version__
    print "License: "+__license__
    print "Credits: "+__credits__[0]+", "+__credits__[1]+", "+__credits__[2]
    print ""
    print "Please cite:"
    print "Github [user: JKlesmith] (www.github.com)"
    print "Kowalsky CA, Klesmith JR, Stapleton JA, Kelly V, Reichkitzer N, Whitehead TA. 2015. High-Resolution Sequence-Function Mapping of Full-Length Proteins. PLoS ONE 10(3):e0118193. doi:10.1371/journal.pone.0118193."
    print "Klesmith JR, Bacik J-P, Michalczyk R, Whitehead TA. 2015. Comprehensive Sequence-Flux Mapping of a Levoglucosan Utilization Pathway in E. coli."
    print ""
    print "Normalization run parameters:"
    print time.strftime("%H:%M:%S")
    print time.strftime("%m/%d/%Y")
    print "Start residue (-s): "+args.startresidue
    print "Normalization type (-n): "+args.normtype
    if args.normtype == "growth":
        print "GROWTH: Doublings (gp) (-g): "+args.gp

    if args.normtype == "FACS":
        print "FACS: SD (-d): "+args.stddev
        print "FACS: Percent Collected (-c): "+args.percentcollected
        print "FACS: Theoretical max enrichment based off of percent collected: "+str(THEOENRICHMENT)
    
    print "Tile Length (-l): "+args.length
    print "Enrich output directory (-p): "+args.path
    print "Unselected counts to be significant (-t): "+str(args.sigthreshold)
    print "Wild-type sequence file (-w): "+args.wildtype
    
    #Build Matrix
    Build_Matrix()
    
    #Get the selected counts
    Get_Unsel_Counts()
    Get_Sel_Counts()
    
    #Get the raw log2 data
    if OverrideEwtEi == True:
        #Set the manual Ewt enrichment
        Ewt = ManualEwt
        print "Manually set Ewt (-y): "+str(Ewt)
        print "Ei scalar transform (-z): "+str(EiScalar)
    else:
        Get_WT_Ewt()

    Get_Mut_Ei()
    
    #Normalize the Data
    Normalize()

    #Print out a csv
    Make_CSV()

if __name__ == '__main__':
    main()