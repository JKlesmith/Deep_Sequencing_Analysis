#!/usr/bin/python

#Copyright (c) 2016, Justin R. Klesmith
#All rights reserved.
#EnrichmentExtract : Process the raw enrichments from the Enrich project dirs

from __future__ import division
from subprocess import check_output
from math import log, sqrt, pow, e
from scipy import special
import numpy as np
import StringIO, argparse, time, os

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2016, Justin R. Klesmith"
__credits__ = ["Justin R. Klesmith", "Timothy A. Whitehead"]
__license__ = "BSD-3"
__version__ = "1.X, Build: 201607X"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["klesmit3@msu.edu", "justinklesmith@gmail.com", "justinklesmith@evodyn.com"]

#Get commandline arguments
parser = argparse.ArgumentParser(description='QuickNormalize '+__version__+' for Growth or FACS')
parser.add_argument('-s', dest='startresidue', action='store', required=True, help='What is the start residue? ie: 0, 40, 80')
parser.add_argument('-q', dest='heatfilename', action='store', help='File name for the heatmap')
parser.add_argument('-p', dest='path', action='store', required=True, help='What is the path to the enrich project directory? ie: ./tile/')
parser.add_argument('-t', dest='sigthreshold', action='store', nargs='?', const=1, default=5, help='Unselected counts for significance. Default = 5')
args = parser.parse_args()

#Verify inputs
if args.path == None:
    print "Enrich Project Path Missing"
    quit()
    
StartResidue = int(args.startresidue) #Starting residue for your tile   
OutputFilename = args.heatfilename #Name for the heatmap csv
Path = args.path+"/data/output/" #What is the path to the output directory
ConfigPath = args.path+"/input/example_local_config" #Path to the config file
SignificantThreshold = int(args.sigthreshold) #Number of counts in the unselected library and selected library to be significant

with open(ConfigPath) as infile:
    for line in infile:
        if line.startswith("<wtPRO>"):
            Len = len(line)
            WTSeq = line[7:Len-10]
            TileLen = len(WTSeq)

print WTSeq
print TileLen            
            
#AA_Table = '*ACDEFGHIKLMNPQRSTVWY'
AA_Table = '*FWYPMILVAGCSTNQDEHKR'

Mutations = {} #Mutations matrix
Ewt = None #Initialize the variable for the wildtype enrichment
UCwt = None #Unselected WT counts
SCwt = None #Selected WT counts

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
                #Mutations[ResID][MutID[1]][0 = RawLog2, 1 = Depleted Enrichments, 2 = Unselected, 3 = Selected, 4=Unused, 5=WT]
                Mutations[j][i[1]] = [None, None, None, None, None, False]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None, None, None, None, False]

    return Mutations

######################################################################################
#Get_WT_Ewt
#This gets the wild-type enrichment from the enrich NA-NA output
######################################################################################
 
def Get_WT():
    global Ewt
    global UCwt
    global SCwt
    #Extract NA-NA WT Ewt log2
    
    awk = ""
    awk2 = ""
    awk3 = ""
    
    if os.path.isfile(Path+'ratios_sel_example_F_N_include_filtered_B_PRO_qc_unsel_example_F_N_include_filtered_B_PRO_qc'):
        awk = check_output(["awk", '{ print $5,$6,$8 }', Path+'ratios_sel_example_F_N_include_filtered_B_PRO_qc_unsel_example_F_N_include_filtered_B_PRO_qc'])
    elif os.path.isfile(Path+'ratios_sel_example_F_N_include_filtered_R1_PRO_qc_unsel_example_F_N_include_filtered_R1_PRO_qc'):
        awk = check_output(["awk", '{ print $5,$6,$8 }', Path+'ratios_sel_example_F_N_include_filtered_R1_PRO_qc_unsel_example_F_N_include_filtered_R1_PRO_qc'])
    else:
        print "Selected protein ratios file not found...exit"
        quit()
        
    if os.path.isfile(Path+'counts_sel_example_F_N_include_filtered_B_PRO_qc'):
        awk2 = check_output(["awk", '{ print $5,$6,$9 }', Path+'counts_sel_example_F_N_include_filtered_B_PRO_qc'])
    elif os.path.isfile(Path+'counts_sel_example_F_N_include_filtered_R1_PRO_qc'):
        awk2 = check_output(["awk", '{ print $5,$6,$9 }', Path+'counts_sel_example_F_N_include_filtered_R1_PRO_qc'])
    else:
        print "Sel protein counts file not found...exit"
        quit()
        
    if os.path.isfile(Path+'counts_unsel_example_F_N_include_filtered_B_PRO_qc'):
        awk3 = check_output(["awk", '{ print $5,$6,$9 }', Path+'counts_unsel_example_F_N_include_filtered_B_PRO_qc'])
    elif os.path.isfile(Path+'counts_unsel_example_F_N_include_filtered_R1_PRO_qc'):
        awk3 = check_output(["awk", '{ print $5,$6,$9 }', Path+'counts_unsel_example_F_N_include_filtered_R1_PRO_qc'])
    else:
        print "Unsel protein counts file not found...exit"
        quit()        

    
    #Loop through the output
    for line in StringIO.StringIO(awk):
        split = line.split(" ")
        location = str(split[0])
        identity = str(split[1])
        if location == "NA" and identity == "NA":
            Ewt = float(split[2].rstrip('\n'))
            print "Wild-type log2 (Ewt): "+str(Ewt)
            
    #Loop through the output
    for line in StringIO.StringIO(awk2):
        split = line.split(" ")
        location = str(split[0])
        identity = str(split[1])
        if location == "NA" and identity == "NA":
            SCwt = int(split[2].rstrip('\n'))
            print "Selected wild-type counts: "+str(SCwt)
            
    #Loop through the output
    for line in StringIO.StringIO(awk3):
        split = line.split(" ")
        location = str(split[0])
        identity = str(split[1])
        if location == "NA" and identity == "NA":
            UCwt = int(split[2].rstrip('\n'))
            print "Unselected wild-type counts: "+str(UCwt)

    return

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

def Enrich():
    #Check to see if the wild-type enrichment is set
    if Ewt == None:
        print "Error: Wild-Type enrichment is not set...quit"
        quit()

    print ""
    print "Raw Enrichment Data"
    print "Location,Mutation,Enrichment,Depleted_Enrichment,Unselected_Reads,Selected_Reads"
    
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            Enrichment = None
            DepEnrichment = None
            UCounts = None
            SCounts = None
            DepEnrichment = "False"
        
        	#Check for a case where a significant variant fell out of the population
            if Mutations[j][i[1]][0] == None and Mutations[j][i[1]][2] >= SignificantThreshold and Mutations[j][i[1]][3] == None:
                Enrichment = log((1/Mutations[j][i[1]][2]), 2) #Calculate the raw log2 for this variant and report it as less than this value
                Mutations[j][i[1]][1] = Enrichment
                DepEnrichment = "True"
                UCounts = Mutations[j][i[1]][2]
                SCounts = Mutations[j][i[1]][3]
                
            #Report the regular enrichment
            elif Mutations[j][i[1]][0] != None and Mutations[j][i[1]][2] >= SignificantThreshold: #Report the significant enrichments
                Enrichment = Mutations[j][i[1]][0]
                Mutations[j][i[1]][1] = Enrichment
                UCounts = Mutations[j][i[1]][2]
                SCounts = Mutations[j][i[1]][3]
            
            elif Mutations[j][i[1]][2] < SignificantThreshold: #Report the insignificant NEs
                if WTSeq[j] == i[1]: #Check to see if it's wildtype else it's Not Significant
                    Enrichment = Ewt
                    Mutations[j][i[1]][1] = Ewt
                    UCounts = UCwt
                    SCounts = SCwt
                    Mutations[j][i[1]][5] = True #Set the WT flag 
            
            elif Mutations[j][i[1]][2] == None and Mutations[j][i[1]][3] >= SignificantThreshold: #Error: Mutation with selected counts and no unselected
                Enrichment = "Error: Sel with Zero Unsel"
            
            else:
                print "Error: unknown enrichment problem."
                
                
            #Print out column data
            print str(j+StartResidue)+","+i[1]+","+str(Enrichment)+","+DepEnrichment+","+str(UCounts)+","+str(SCounts)
    
    return Mutations

######################################################################################
#Make_CSV
#This outputs the fitness metrics to a CSV file to be imported into excel
######################################################################################
    
def Make_CSV():
    print "Enrichment Heatmap"
    #This makes a CSV style report of rows of letters and columns of residues
    
    #Print off the Number
    Numbering = " "
    for q in xrange(0,TileLen):
        Numbering = Numbering+","+str(StartResidue+q)
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
        for j in xrange(0,TileLen):
            Output = Output+str(Mutations[j][i[1]][1])+","
        Output = Output+"\n"
    print Output

    #Write the heatmap to a newfile
    outfile = open('enrichmentheatmap_'+OutputFilename+'_'+str(StartResidue)+'.csv', 'w')
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
    #Build Matrix
    Build_Matrix()
    
    #Get the selected counts
    Get_Unsel_Counts()
    Get_Sel_Counts()
    
    #Get the enrichments
    Get_WT()
    Get_Mut_Ei()
    
    #Process the data
    Enrich()
    Make_CSV()

if __name__ == '__main__':
    main()