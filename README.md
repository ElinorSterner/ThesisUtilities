# ThesisUtilities
Tools for trees, etc
## cluster.py
Clusters by a list of similarities in a list in the top of the script. Iterates through Subtrees_unaligned from clade grabbing. 
Call in the command line, not within sublime text because sublime does not work from within the conda environment  (it won't find cd-hit)
Elinor 2/11/24

## plot_cladesize.py
Input is CladeSizesPerTaxon.csv output and a text file create by:
1. copying the 3rd row into sublime text. 
2. Replaced ( and ) with a space
3. Replaced space with a \n to make one clade size in each row. 
4. Named this file clade_sizes.txt
This script filters out clades with less than 2 sequences and plot with histograms and boxplots.

