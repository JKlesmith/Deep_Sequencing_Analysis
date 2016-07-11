#!/usr/bin/python

#Copyright (c) 2016, Justin R. Klesmith
#All rights reserved.
#Entropy

from __future__ import division
from subprocess import check_output
from math import log, pow
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
parser = argparse.ArgumentParser(description='Shannon Entropy '+__version__)
parser.add_argument('-n', dest='normtype', action='store', required=True, help='Normalization Type? Enter: NonFACS or FACS')
parser.add_argument('-s', dest='startresidue', action='store', required=True, help='What is the start residue? ie: 0, 40, 80')
parser.add_argument('-d', dest='stddev', action='store', help='Standard Deviation? (FACS) ie: 0.6')
parser.add_argument('-c', dest='percentcollected', action='store', help='Percent Collected? (FACS) ie: 0.05')
parser.add_argument('-t', dest='sigthreshold', action='store', nargs='?', const=1, default=5, help='Unselected counts for significance. Default = 5')
parser.add_argument('-p', dest='path', action='store', required=True, help='What is the path to the enrich  directory? ie: ./tile/')
parser.add_argument('-y', dest='ewtenrichment', action='store', help='Manual Ewt enrichment value')
parser.add_argument('-z', dest='eiscalar', action='store', help='Manual Ei enrichment scalar')
args = parser.parse_args()

#Verify inputs
if args.path == None:
    print "Enrich Project Path Missing"
    quit()
    
if args.stddev == None and args.normtype == "FACS":
    print "Missing SD. Flag: -d"
    quit()
else:
    SD = float(args.stddev) #Standard Deviation

if args.percentcollected == None and args.normtype == "FACS":
    print "Missing percent collected. Flag: -c"
    quit()
else:
    PC = float(args.percentcollected) #Percent collected
    THEOENRICHMENT = -log(PC, 2) #Theoretical maximum enrichment 

if args.ewtenrichment and args.eiscalar != None:
    #This section is only true if we want to provide our own WT enrichment and a scalar to add to Ei
    OverrideEwtEi = True
    ManualEwt = float(args.ewtenrichment)
    EiScalar = float(args.eiscalar)
else:
    OverrideEwtEi = False    
    
StartResidue = int(args.startresidue) #Starting residue for your tile   
Path = args.path+"/data/output/" #What is the path to the output directory
ConfigPath = args.path+"/input/example_local_config" #Path to the config file
SignificantThreshold = int(args.sigthreshold) #Number of counts in the unselected library and selected library to be significant

with open(ConfigPath) as infile:
    for line in infile:
        if line.startswith("<wtPRO>"):
            Len = len(line)
            WTSeq = line[7:Len-10]
            TileLen = len(WTSeq)
  
AA_Table = 'ACDEFGHIKLMNPQRSTVWY'
Mutations = {} #Mutations matrix
Ewt = None #Initialize the variable for the wildtype enrichment
UCwt = None #Unselected WT counts
SCwt = None #Selected WT counts

def Build_Matrix():
    #Populate mutation array with None data
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            try:
                #Mutations[ResID][MutID[1]][0 = RawLog2, 1 = Unused, 2 = UnselectedCounts, 3 = SelectedCounts, 4=p-value, 5=WT]
                Mutations[j][i[1]] = [None, None, None, None, None, False]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None, None, None, None, False]

    return Mutations

def Get_WT_Ewt():
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
        
        #Skip stop codons
        if identity == "*":
            continue

        #Check to see if we're above the tile length and go to next
        if location >= TileLen:
            continue

        Ei = float(split[2].rstrip('\n'))
        
        #Check to see if the enrichment is greater or equal than the theoretical
        if args.normtype == "FACS":
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
        else:
            Mutations[location][identity][0] = Ei

    return Mutations

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
        
        #Skip stop codons
        if identity == "*":
            continue

        #Check to see if we're above the tile length and go to next
        if location >= TileLen:
            continue

        counts = int(split[2].rstrip('\n'))
        Mutations[location][identity][2] = counts #Set the unselected counts
	
    return Mutations

def AssignWT():
    #Assign the WT residues

    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            if i[1] == WTSeq[j]:
                Mutations[j][i[1]][5] = True

    return Mutations

def NumMutinCol(ID):
    #Returns the number of sig mutants at a residue location
    
    #Initialize our variable
    NUMSIGMUTS = 1 #Start at one to account for WT
    
    #Loop through the mut types
    for i in enumerate(AA_Table):
        if Mutations[ID][i[1]][2] >= SignificantThreshold:
            NUMSIGMUTS = NUMSIGMUTS + 1
    
    return NUMSIGMUTS
    
def Shannon():

    #Check to see if the wild-type enrichment is set
    if Ewt == None:
        print "Error: Wild-Type enrichment is not set...quit"
        quit()

    print ""
    print "Shannon Entropy"
    print "Location,WT Residue,Shannon Entropy,Number of Mutations Counted+WT"
    
    #Check to see if the wild-type enrichment is set
    if Ewt == None:
        print "Error: Wild-Type enrichment is not set...quit"
        quit()

    #Check for a case where a significant variant fell out of the population
    for j in xrange(0,TileLen):
        for i in enumerate(AA_Table):
            if Mutations[j][i[1]][0] == None and Mutations[j][i[1]][2] >= SignificantThreshold and Mutations[j][i[1]][3] == None:
                Mutations[j][i[1]][0] = log((1/Mutations[j][i[1]][2]), 2) #Calculate the raw log2 for this variant and report it as less than this value

    #Calculate p-values
    for j in xrange(0,TileLen):

        #First calculate the column sum
        pcol = 0
        for i in enumerate(AA_Table):
            if Mutations[j][i[1]][5] == False: #Check to see if it's Wild-Type
                if Mutations[j][i[1]][2] >= SignificantThreshold: #Check to see if the count is above the counting threshold
                    pcol = pcol + pow(2, float(Mutations[j][i[1]][0]))
            else:
                pcol = pcol + pow(2, float(Ewt))
        
        #Then calculate the individual p-value and store the shannon entropy
        for i in enumerate(AA_Table):
            if Mutations[j][i[1]][5] == False: #Check to see if it's Wild-Type
                if Mutations[j][i[1]][2] >= SignificantThreshold: #Check to see if the count is above the counting threshold
                    Mutations[j][i[1]][4] = (pow(2, float(Mutations[j][i[1]][0]))/pcol)
            else:
                Mutations[j][i[1]][4] = (pow(2, float(Ewt))/pcol)
                
    #Calculate the residue shannon entropy
    for j in xrange(0,TileLen):
        SE = 0
        for i in enumerate(AA_Table):
            if Mutations[j][i[1]][2] >= SignificantThreshold: #Check to see if the count is above the counting threshold
                SE = SE + -1*Mutations[j][i[1]][4]*log(Mutations[j][i[1]][4])
        
        #Normalize our Shannon Entropy
        try:
            SE = (SE*log(20))/log(NumMutinCol(j))
        except ZeroDivisionError:
            print "Your tile length is possibly too long or there is no mutations besides WT at position "+str(j)
        
        #Output the entropy values
        print str(StartResidue+j)+","+WTSeq[j]+","+str(SE)+","+str(NumMutinCol(j))

    return
    
def main():
    global Ewt
    
    #Write out preamble
    print "Shannon Entropy"
    print "Author: "+__author__
    print "Contact: "+__email__[0]+", "+__email__[1]+", "+__email__[2]
    print __copyright__
    print "Version: "+__version__
    print "License: "+__license__
    print "Credits: "+__credits__[0]+", "+__credits__[1]
    print ""
    print "Cite:"
    print "Github [user: JKlesmith] (www.github.com)"
    print ""
    print "Run parameters:"
    print time.strftime("%H:%M:%S")
    print time.strftime("%m/%d/%Y")
    print "Start residue (-s): "+str(StartResidue)
    print "Tile length: "+str(TileLen)
    print "Significant count threshold (-t): "+str(SignificantThreshold)
    print "Enrich directory (-p): "+Path

    #Make the internal matrix
    Build_Matrix()
    
    #Assign the WT residues
    AssignWT()
    
    #Get the counts
    Get_Unsel_Counts()
    
    #Get the log2 data
    if OverrideEwtEi == True:
        #Set the manual Ewt enrichment
        Ewt = ManualEwt
        print "Manually set Ewt (-y): "+str(Ewt)
        print "Ei scalar transform (-z): "+str(EiScalar)
    else:
        Get_WT_Ewt()
    
    Get_Mut_Ei()
    
    #Print out a csv
    Shannon()

if __name__ == '__main__':
    main()