#!/usr/bin/python

#Copyright (c) 2016, Justin R. Klesmith
#All rights reserved.
#Get the codon enrichment

from __future__ import division
from subprocess import check_output
from math import log, sqrt, pow, e
from scipy import special
import StringIO
import argparse
import time

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2016, Justin R. Klesmith"
__credits__ = ["Justin R. Klesmith", "Timothy A. Whitehead"]
__license__ = "BSD-3"
__version__ = "1.2, Build: 20160129"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["klesmit3@msu.edu", "justinklesmith@gmail.com", "justinklesmith@evodyn.com"]

parser = argparse.ArgumentParser(description='CodonNormalize for Growth or FACS')
parser.add_argument('-n', dest='normtype', action='store', required=True, help='Normalization Type? Enter: growth or FACS')
parser.add_argument('-s', dest='startresidue', action='store', required=True, help='What is the start of translation? ie: 25')
parser.add_argument('-l', dest='length', action='store', required=True, help='Length of your tile? ie: 40, 80 amino acids')
parser.add_argument('-g', dest='gp', action='store', help='How many doublings/generations? (GROWTH) ie: 12.5')
parser.add_argument('-d', dest='stddev', action='store', help='Standard Deviation? (FACS) ie: 0.6')
parser.add_argument('-c', dest='percentcollected', action='store', help='Percent Collected? (FACS) ie: 0.05')
parser.add_argument('-p', dest='path', action='store', required=True, help='What is the path to the tile directory? ie: ./tile/')
parser.add_argument('-t', dest='sigthreshold', action='store', nargs='?', const=1, default=5, help='Unselected counts for significance. Default = 5')
parser.add_argument('-o', dest='heatmap', action='store', nargs='?', const=1, default='False', help='Output a csv heatmap and file? Default = False')
parser.add_argument('-w', dest='synonwt', action='store', nargs='?', const=1, default='False', help='Output the wild-type codons in the synonymous mutation analysis? Default = False')
parser.add_argument('-x', dest='synon', action='store', nargs='?', const=1, default='False', help='Output the synonymous mutation analysis and file? Default = False')
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
    print "Missing data path. Flag: -p"
    quit()
    
if args.ewtenrichment and args.eiscalar != None:
    #This section is only true if we want to provide our own WT enrichment and a scalar to add to Ei
    OverrideEwtEi = True
    ManualEwt = float(args.ewtenrichment)
    EiScalar = float(args.eiscalar)
else:
    OverrideEwtEi = False
    
#Global Variables
StartResidue = int(args.startresidue) #Starting residue for your tile
TileLen = int(args.length) #Length of your tile
TranslateEnd = StartResidue+(3*TileLen) #Get the end of the DNA
Path = args.path #What is the path to the output directory
SignificantThreshold = int(args.sigthreshold) #Number of counts in the unselected library and selected library to be significant

if args.normtype == "growth":
    DoublingsGp = float(args.gp) #Number of doublings

if args.normtype == "FACS":
    SD = float(args.stddev) #Standard Deviation
    PC = float(args.percentcollected) #Percent collected
    THEOENRICHMENT = -log(PC, 2) #Theoretical maximum enrichment

Mutations = {} #Mutations matrix
WTCodon_Table = {}
WTAmino_Table = {}
Ewt = None #Initialize the variable for the wildtype enrichment
CountsUwt = None #Number of unselected counts
CountsSwt = None #Number of selected counts
WTSeq = None #Variable for the wild-type DNA sequence

DNA_Table = {'TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA',
'TTG', 'TCG', 'TAG', 'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC',
'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT',
'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG',
'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA',
'GTG', 'GCG', 'GAG', 'GGG'}
Translation_Table = {'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'}
Translation_Table_AtoD = {
'*' : ['TAA', 'TAG', 'TGA'],
'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
'C' : ['TGT', 'TGC'],
'D' : ['GAT', 'GAC'],
'E' : ['GAA', 'GAG'],
'F' : ['TTT', 'TTC'],
'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
'H' : ['CAT', 'CAC'],
'I' : ['ATT', 'ATC', 'ATA'],
'K' : ['AAA', 'AAG'],
'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
'M' : ['ATG'],
'N' : ['AAT', 'AAC'],
'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
'Q' : ['CAA', 'CAG'],
'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
'W' : ['TGG'],
'Y' : ['TAT', 'TAC']}

def Build_Mut_Arrays():
    #Populate Mutation Dictionary with None Data
    for j in xrange(0,TileLen):
        for i in enumerate(DNA_Table):
            try:
                #Mutations[Codon#][DNA[1]][0 = RawLog2, 1 = NormLog2, 2 = Unselected, 3 = Selected, 4 = DNA Sequence]
                Mutations[j][i[1]] = [None, None, None, None, None]
            except KeyError:
                Mutations[j] = {}
                Mutations[j][i[1]] = [None, None, None, None, None]
        
        #For each codon j set the wild-type DNA sequence
        WTCodon_Table[j] = WTSeq[j*3:(j*3)+3] #Set the WT codons
        
        #For each codon j set the amino acid
        try:
            WTAmino_Table[j] = Translation_Table[WTCodon_Table[j]]
        except KeyError:
            print "Error: Your tile length is probably set too high or the last codon is not three bases."
            exit(1)
            
    return Mutations, WTCodon_Table, WTAmino_Table

def Set_Enrichment(codonnum, codon, Ei):
    global Mutations

    #For FACS set a upper limit on enrichment, don't do anything for growth
    if args.normtype == "FACS":
        #Check to see if the enrichment is greater or equal than the theoretical
        if OverrideEwtEi == False: #Apply no scalar to the Ei
            if Ei >= THEOENRICHMENT:
                Mutations[codonnum][codon][0] = (THEOENRICHMENT - 0.001)
            else:
                Mutations[codonnum][codon][0] = Ei
        elif OverrideEwtEi == True: #Apply a scalar to the Ei
            if Ei >= (THEOENRICHMENT + EiScalar):
                Mutations[codonnum][codon][0] = ((THEOENRICHMENT + EiScalar) - 0.001)
            else:
                Mutations[codonnum][codon][0] = (Ei + EiScalar)
        
    elif args.normtype == "growth":
        Mutations[codonnum][codon][0] = Ei

    return

def Read_Enrichments():
    global Ewt
    global WTSeq
  
    awkout = ""
    try:
        awkout = check_output(["awk", 'FNR>1{ print $2,$4,$5,$6,$8 }', args.path+'ratios_sel_example_F_N_include_filtered_B_DNA_qc_unsel_example_F_N_include_filtered_B_DNA_qc'])
    except:
        print "AWK problem: probably incorrect path"
        exit(1)
  
    #First find the wild-type and build the arrays before we do anything else
    for line in StringIO.StringIO(awkout):
        split = line.split(" ")
        
        DNASeq = split[0]
        mutationcount = int(split[1])
        location = str(split[2])
        identity = str(split[3])
        enrichment = float(split[4].rstrip('\n'))
        
        if mutationcount == 0: #Wild-Type
            WTSeq = DNASeq[StartResidue:TranslateEnd] #Set the wild-type sequence that is being translated
            Ewt = enrichment #Set the wild-type sequence enrichment
            print "Wild-type log2 (Ewt): "+str(Ewt) #Print out the wild-type enrichment
            Build_Mut_Arrays() #Build the arrays now
 
    #Now populate the single, double, and triple mutations
    for line in StringIO.StringIO(awkout):
        split = line.split(" ")
        
        DNASeq = split[0]
        mutationcount = int(split[1])
        locationarr = split[2].split(",") #Note this will be splittable by ","
        identity = str(split[3])
        enrichment = float(split[4].rstrip('\n'))

        if mutationcount == 1: #Single mutations
            location = int(locationarr[0])
            
            if location >= StartResidue and location < TranslateEnd: #Check to see that the base is in our tile
                codonnum = int((location - StartResidue)/3) #Specific Codon Number
                codon = DNASeq[StartResidue+(codonnum*3):StartResidue+(codonnum*3)+3] #The mutated codon
                Set_Enrichment(codonnum, codon, enrichment) #Set the enrichment  
                Mutations[codonnum][codon][4] = DNASeq[StartResidue:TranslateEnd] #Set the sequence

        if mutationcount == 2: #Double mutations
            location1 = int(locationarr[0])
            location2 = int(locationarr[1]) 
            
            if location1 >= StartResidue and location1 < TranslateEnd: #Check to see that the base is in our tile
                if location2 >= StartResidue and location2 < TranslateEnd: #Check to see that the base is in our tile
                    codonnum1 = int((location1 - StartResidue)/3)
                    codonnum2 = int((location2 - StartResidue)/3)
                    if codonnum1 == codonnum2:
                        codon = DNASeq[StartResidue+(codonnum1*3):StartResidue+(codonnum1*3)+3] #The mutated codon
                        Set_Enrichment(codonnum1, codon, enrichment) #Set the enrichment         
                        Mutations[codonnum1][codon][4] = DNASeq[StartResidue:TranslateEnd] #Set the sequence
                
        if mutationcount == 3: #There are three mutations
            location1 = int(locationarr[0])
            location2 = int(locationarr[1])
            location3 = int(locationarr[2])
            
            if location1 >= StartResidue and location1 < TranslateEnd: #Check to see that the base is in our tile
                if location2 >= StartResidue and location2 < TranslateEnd: #Check to see that the base is in our tile
                    if location3 >= StartResidue and location3 < TranslateEnd: #Check to see that the base is in our tile
                        codonnum1 = int((location1 - StartResidue)/3)
                        codonnum2 = int((location2 - StartResidue)/3)
                        codonnum3 = int((location3 - StartResidue)/3)
                        if codonnum1 == codonnum2 and codonnum2 == codonnum3:
                            codon = DNASeq[StartResidue+(codonnum1*3):StartResidue+(codonnum1*3)+3] #The mutated codon
                            Set_Enrichment(codonnum1, codon, enrichment) #Set the enrichment 
                            Mutations[codonnum1][codon][4] = DNASeq[StartResidue:TranslateEnd] #Set the sequence

    return

def Get_Counts(Filename, ArrayKey):
    global Mutations
	#Get the unselected counts for a variant
    
    WildTypeCounts = None
    awkout = ""
    try:
        awkout = check_output(["awk", 'FNR>1{ print $2,$4,$5,$6,$9 }', args.path+Filename])
    except:
        print "AWK problem: file path probably incorrect"
        exit(1)
    
    #Loop through the file
    for line in StringIO.StringIO(awkout):
        split = line.split(" ")

        DNASeq = split[0]
        mutationcount = int(split[1])
        locationarr = split[2].split(",") #Note this will be splittable by ","
        identity = str(split[3])
        counts = int(split[4].rstrip('\n'))
  
        if mutationcount == 0: #Find the wild-type
            WildTypeCounts = counts

        if mutationcount == 1: #Single mutations
            location = int(locationarr[0])
            
            if location >= StartResidue and location < TranslateEnd: #Check to see that the base is in our tile
                codonnum = int((location - StartResidue)/3) #Specific Codon Number
                codon = DNASeq[StartResidue+(codonnum*3):StartResidue+(codonnum*3)+3] #The mutated codon
                Mutations[codonnum][codon][ArrayKey] = counts #Set the unselected counts

        if mutationcount == 2: #Double mutations
            location1 = int(locationarr[0])
            location2 = int(locationarr[1])
            
            if location1 >= StartResidue and location1 < TranslateEnd: #Check to see that the base is in our tile
                if location2 >= StartResidue and location2 < TranslateEnd: #Check to see that the base is in our tile
                    codonnum1 = int((location1 - StartResidue)/3)
                    codonnum2 = int((location2 - StartResidue)/3)
                    if codonnum1 == codonnum2:
                        codon = DNASeq[StartResidue+(codonnum1*3):StartResidue+(codonnum1*3)+3] #The mutated codon
                        Mutations[codonnum1][codon][ArrayKey] = counts #Set the unselected counts

        if mutationcount == 3: #Test to see that there are three mutations
            location1 = int(locationarr[0])
            location2 = int(locationarr[1])
            location3 = int(locationarr[2])            
            
            if location1 >= StartResidue and location1 < TranslateEnd: #Check to see that the base is in our tile
                if location2 >= StartResidue and location2 < TranslateEnd: #Check to see that the base is in our tile
                    if location3 >= StartResidue and location3 < TranslateEnd: #Check to see that the base is in our tile
                        codonnum1 = int((location1 - StartResidue)/3)
                        codonnum2 = int((location2 - StartResidue)/3)
                        codonnum3 = int((location3 - StartResidue)/3)
                        if codonnum1 == codonnum2 and codonnum2 == codonnum3:
                            codon = DNASeq[StartResidue+(codonnum1*3):StartResidue+(codonnum1*3)+3] #The mutated codon
                            Mutations[codonnum1][codon][ArrayKey] = counts #Set the unselected counts

    return WildTypeCounts

def Normalize():
    print "Normalizing the data"
    print "Location, DNA Sequence, Fitness Metric, Unselected_Reads, Selected_Reads, Raw Enrichment"
    
    for j in xrange(0,TileLen):
        for i in enumerate(DNA_Table):
        	#Check for a case where a significant variant fell out of the population
            if Mutations[j][i[1]][0] == None and Mutations[j][i[1]][2] >= SignificantThreshold and Mutations[j][i[1]][3] == None:
                Mutations[j][i[1]][0] = log((1/Mutations[j][i[1]][2]), 2) #Calculate the raw log2 for this variant and report it as less than this value
        
            #Calculate the fitness
            if Mutations[j][i[1]][0] != None and Mutations[j][i[1]][2] >= SignificantThreshold: #Report the significant fitness
                Ei = float(Mutations[j][i[1]][0])
                
                if args.normtype == "growth":
                    Mutant = (Ei/DoublingsGp)+1
                    WT = (Ewt/DoublingsGp)+1
                    NE = log(Mutant/WT, 2)
                elif args.normtype == "FACS":
                    WT = special.erfinv(1-PC*pow(2,(Ewt+1)))
                    Mutant = special.erfinv(1-PC*pow(2,(Ei+1)))
                    NE = (log(e, 2)*sqrt(2)*SD*(WT-Mutant))
                else:
                    print "Error: growth or FACS not set?"
                    quit()
                
                Mutations[j][i[1]][1] = "{0:.2f}".format(NE)
            elif Mutations[j][i[1]][2] < SignificantThreshold: #Report the insignificant NEs
                if WTCodon_Table[j] == i[1]: #Check to see if it's wildtype else it's Not Significant
                    Mutations[j][i[1]][1] = "0.0"
                    Mutations[j][i[1]][2] = "WT"
                    Mutations[j][i[1]][3] = "WT"
                    Mutations[j][i[1]][4] = WTSeq #Set the sequence
            	else:
                    Mutations[j][i[1]][1] = "NS" 
            elif Mutations[j][i[1]][2] == None and Mutations[j][i[1]][3] >= SignificantThreshold: #Error: Mutation with selected counts and no unselected
                Mutations[j][i[1]][1] = "Error: Sel with Zero Unsel"
            else:
                print "Error: unknown normalization problem."
                
            #Print out column data
            print str(j+1)+","+str(Mutations[j][i[1]][4])+","+Mutations[j][i[1]][1]+","+str(Mutations[j][i[1]][2])+","+str(Mutations[j][i[1]][3])+","+str(Mutations[j][i[1]][0])
    
    return Mutations

def Make_CSV():
    print "Normalized Heatmap"
    #This makes a CSV style report of rows of letters and columns of residues
    
    #Print off the Number
    Numbering = " "
    for q in xrange(1,TileLen+1):
        Numbering = Numbering+","+str(q)
    print Numbering
        
    #Print off the WT
    WTResi = " "
    for w in xrange(0,TileLen):
        WTResi = WTResi+","+WTCodon_Table[w]
    print WTResi
    
    #Print off the mutations
    Output = ""
    for i in enumerate(DNA_Table):
        Output = Output+i[1]+","
        for j in xrange(0,TileLen):
            Output = Output+str(Mutations[j][i[1]][1])+","
        Output = Output+"\n"
    print Output

    #Write the heatmap to a newfile
    outfile = open('CodonHeatmap.csv', 'w')
    outfile.write(Numbering+'\n')
    outfile.write(WTResi+'\n')
    outfile.write(Output)
    outfile.close()
    
    return

def Make_SynonCSV():
    #Output the synonymous mutations  
    
    print "Synonymous Mutation Analysis"
    print "ResID, Wild-Type Identity, Wild-Type Codon, Codon, Fitness Metric, Enrichment, Unselected Counts, Selected Counts"
    #q is a int that will loop through all of the codons
    Output = "ResID, Wild-Type Identity, Wild-Type Codon, Codon, Fitness Metric, Enrichment\n"
    for q in xrange(0, TileLen):
        #syn is the codon(s) for the amino acid specified from WTAmino_Table
        for syn in Translation_Table_AtoD[WTAmino_Table[q]]:
            if syn == WTCodon_Table[q]:
                if args.synonwt == 'True' or args.synonwt == 'true': 
                    #Output wild-type data
                    print str(q+1)+","+WTAmino_Table[q]+","+WTCodon_Table[q]+","+syn+","+Mutations[q][syn][1]+","+str(Ewt)+","+str(Mutations[q][syn][2])+","+str(Mutations[q][syn][3])
                    Output = Output + str(q+1)+","+WTAmino_Table[q]+","+WTCodon_Table[q]+","+syn+","+Mutations[q][syn][1]+","+str(Ewt)+","+str(Mutations[q][syn][2])+","+str(Mutations[q][syn][3])+"\n"
                else:
                    continue
            else:
                print str(q+1)+","+WTAmino_Table[q]+","+WTCodon_Table[q]+","+syn+","+Mutations[q][syn][1]+","+str(Mutations[q][syn][0])+","+str(Mutations[q][syn][2])+","+str(Mutations[q][syn][3])
                Output = Output + str(q+1)+","+WTAmino_Table[q]+","+WTCodon_Table[q]+","+syn+","+Mutations[q][syn][1]+","+str(Mutations[q][syn][0])+","+str(Mutations[q][syn][2])+","+str(Mutations[q][syn][3])+"\n"
    
    #Write an output file
    outfile = open('SynonCSV.csv', 'w')
    outfile.write(Output)
    outfile.close()
    
    return

def main():
    global CountsUwt
    global CountsSwt
    
    #Write out preamble
    print "CodonNormalize for growth or FACS selections"
    print __author__
    print "Contact: "+__email__[0]+", "+__email__[1]+", "+__email__[2]
    print __copyright__
    print "License: "+__license__
    print "Credits: "+__credits__[0]+", "+__credits__[1]

    print "Normalization run parameters:"
    print time.strftime("%H:%M:%S")
    print time.strftime("%m/%d/%Y")
    print "Translate Start: "+args.startresidue
    print "Normalization type: "+args.normtype
    if args.normtype == "growth":
        print "Doublings: "+args.gp

    if args.normtype == "FACS":
        print "SD: "+args.stddev
        print "Percent Collected: "+args.percentcollected
    
    print "Tile Length: "+args.length
    print "Data path: "+args.path
    print "Unselected counts to be significant: "+str(args.sigthreshold)

    #Build Mut Arrays and Get Enrichments
    Read_Enrichments()

    #Get the unselected counts
    CountsUwt = Get_Counts('counts_unsel_example_F_N_include_filtered_B_DNA_qc', 2)
    print "Number of unselected wild-type counts: "+str(CountsUwt)
    
    #Get the selected counts
    CountsSwt = Get_Counts('counts_sel_example_F_N_include_filtered_B_DNA_qc', 3)
    print "Number of selected wild-type counts: "+str(CountsSwt)

    #Normalize the Data
    Normalize()

    #Output a CSV
    if args.heatmap == 'True' or args.heatmap == 'true':
        Make_CSV()
    
    #Run the synonymous mutation analysis
    if args.synon == 'True' or args.synon == 'true':
        Make_SynonCSV()

if __name__ == '__main__':
    main()