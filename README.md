# Deep Sequencing Data Analysis
Software used for data analysis

# Script: QuickNormalize.py
This script normalizes a growth selection or a FACS screen.

Command line:
python QuickNormalize.py –n [growth or FACS] –s [start residue] –l [tile length] –g [growth: gp] –d [FACS: std dev] –c [FACS: percent collected] –p [path to enrich output directory] –t [significant unselected counts threshold (default: 5)] –w [path to wild-type AA sequence (default: ./WTSeq)] –o [Output a separate heatmap csv (default: True)]

Example inputs for flags:
-n String: growth or FACS
-s Integer: 0 40 80 …
-l Integer: 40 80 …
-g Float: 10.0
-d Float: 0.6
-c Float: 0.05 (For 5%)
-p String: ./output/ (Tailing slash is needed)
-t Integer: 5 (Default: 5)
-w String: ./wt_seq.txt (Default: ./WTSeq)
-o String: True or False (Default: True)

Example command line for growth:
python QuickNormalize.py –n growth –s 0 –l 40 –g 10.0 –p ./1/data/output/ > Tile1

Example command line for FACS:
python QuickNormalize.py –n FACS –s 0 –l 76 –d 0.6 –c 0.05 –p ./tile/data/output/ -w ./wt_seq.txt > Tile1

Help command:
python QuickNormalize.py –h

Notes:
It is highly recommended to direct the output to a file using > [file name] such that there is a saved copy of the normalization output. The output includes column and csv heatmap data used by other scripts for further analyses (stats, replicate errors). The output heatmap csv will be named heatmap_startresi_#.csv with # being the number given for the start residue. The wild-type amino acid sequence file is a single line ASCII file with the wild-type amino acid sequence. This file must be stripped of special characters hidden with rich-text editors. GNU nano can be used to edit this file and strip special hidden characters.

# Script: FACSEntropy.py
This script calculates the Shannon Entropy for a FACS screen.

Command line:
python FACSEntropy.py –s [start residue] –l [tile length] –d [FACS: std dev] –c [FACS: percent collected] –p [path to enrich output directory] –t [significant unselected counts threshold (default: 5)] –w [path to wild-type AA sequence (default: ./WTSeq)]

Example inputs for flags:
-s Integer: 0 40 80 …
-l Integer: 40 80 …
-d Float: 0.6
-c Float: 0.05 (For 5%)
-p String: ./output/ (Tailing slash is needed)
-t Integer: 5 (Default: 5)
-w String: ./wt_seq.txt (Default: ./WTSeq)

Example command line for FACS:
python FACSEntropy.py –s 0 –l 76 –d 0.6 –c 0.05 –p ./tile/data/output/ -w ./wt_seq.txt > Tile1

Help command:
python FACSEntropy.py –h

Notes:
Normalization of the dataset is not needed to run this script. It is highly recommended to direct the output to a file using > [file name] such that there is a saved copy of the entropy output. The script outputs residue number and entropy values. The wild-type amino acid sequence file is a single line ASCII file with the wild-type amino acid sequence. This file must be stripped of special characters hidden with rich-text editors. GNU nano can be used to edit this file and strip special hidden characters.

# Script: QuickStats.py
This script calculates the reportable statistics for a deep sequencing run.

Command line:
python QuickStats.py –f [path to file with normalized output] –p [path to root enrich tile directory]

Example inputs for flags:
-f String: ./Tile1Normed (File from QuickNormalize output)
-p String: ./tile/ (Tailing slash is needed)

Example command line:
python QuickStats.py –f ./Tile1Normed –p ./tile/ > Tile1Stats

Help command:
python QuickStats.py –h

Enrich files used:
data/output/counts_sel_example_F_N_include_filtered_B_DNA_qc
data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc
data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m1
data/output/counts_unsel_example_F_N_include_filtered_B_DNA_qc.m2
data/output/counts_unsel_example_F_N_include_filtered_B_PRO_qc
data/output/counts_unsel_example_F_N_include_filtered_B_PRO_qc.m1
input/example_local_config

Notes:
Normalization of the dataset is required to run this script. Additionally, the script uses the <translate_start> tag from the example_local_config file. Therefore, the enrich patch is needs to be applied. Unlike the normalization and other scripts, this script needs the root directory of the tile (i.e. the directory with the data and input directories). It is highly recommended to direct the output to a file using > [file name] such that there is a saved copy of the stats output. The script outputs all reportable statistics. This file assumes a certain naming scheme for the enrich output (listed above).