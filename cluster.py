# python cluster.py
# grabs the sequences of interest (taxa_tag) from clade grabbing output, and then runs CD-HIT on those at
# a given percent ID parameter. (edit "params = ['x', 'x']" at the end of the script)
# call in the command line, not within sublime text because sublime does not work from within the conda environment.
# Elinor 2/11/24




import os
from Bio import SeqIO
from pathlib import Path
import subprocess

params = ['0.99', '1.00']#enter the clustering parameters here




def cluster(percent_id):
	
	Path(f'Clusters_{percent_id}').mkdir(parents=True, exist_ok=True)#makes output folder
	in_dir = "SubTrees_unaligned_At_only"


	for subclade in os.listdir(in_dir):
		if subclade.endswith('fasta'):
			outfile = f'{subclade}_{percent_id}_output.fasta'
			print(outfile)

			command  = f'cd-hit -i {in_dir}/{subclade} -o Clusters_{percent_id}/{outfile} -c 0.99'
			
			subprocess.run([command], shell=True)


def select_taxa():

	Path(f'SubTrees_unaligned_At_only').mkdir(parents=True, exist_ok=True)#makes output folder
	in_dir = "Subtrees_unaligned"
	taxa_tag = "At_"

	for subclade in os.listdir("Subtrees_unaligned"): 
		if subclade.endswith('fasta'):

			to_keep = {}
			records = list(SeqIO.parse(f'Subtrees_unaligned/{subclade}', 'fasta'))

			for rec in records:
				if rec.id.startswith(taxa_tag):
					to_keep.update({rec.id : rec.seq})


			with open(f'SubTrees_unaligned_At_only/{taxa_tag}Only_{subclade}', 'w') as o:
				for idd, seq in to_keep.items():
					o.write(f'>{idd}\n{seq}\n')

if __name__ == '__main__':
	select_taxa()#this function goes into the subtrees unaligned files and takes out only taxa starting with "taxon_tag"
	
	params = ['0.99', '1.00']#enter the clustering parameters here
	for percent_id in params: # go through clustering parameters
		cluster(percent_id)



























