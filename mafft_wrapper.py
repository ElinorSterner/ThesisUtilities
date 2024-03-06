import os
import sys
from Bio import SeqIO
from pathlib import Path


Path(f'Subtrees_At_Aligned').mkdir(parents=True, exist_ok=True)#makes output folder

indir = "SubTrees_unaligned_At_only"

for file in os.listdir(indir):
	if file.endswith("fasta"):
		os.system('mafft ' + indir + '/' + file + ' > ' + "Subtrees_At_Aligned/"  + file)