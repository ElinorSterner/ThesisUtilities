# purpose: grab sequences for eggnog mapper for functional analysis for allogromia gene expression
# copies the matching files into a directory (to_cat)
# after running this go to that dir and run in command line: cat *fasta > cat_ogs.fasta


import os
from pathlib import Path

Path(f'to_cat').mkdir(parents=True, exist_ok=True)#makes output folder


ogs = [line.strip('\n') for line in open("uniq_ogs.txt").readlines()]
print(len(ogs))

to_cat = []
for fasta in os.listdir("Clusters_0.99"):
	if fasta.endswith("fasta"):
		if fasta[8:18] in ogs:
			to_cat.append(fasta)
			os.system(f'cp Clusters_0.99/{fasta} to_cat/')
		else:
			print("no")

#print(to_cat)

