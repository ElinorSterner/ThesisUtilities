#call in the command line, not within sublime text because sublime does not work from within the conda environment.
#Elinor 2/11/24




import os
from pathlib import Path
import subprocess

params = ['0.99', '1.00']#enter the clustering parameters here




def cluster(percent_id):
	
	print(percent_id)

	Path(f'Clusters_{percent_id}').mkdir(parents=True, exist_ok=True)#makes output folder
	in_dir = "test_cluster"


	for subclade in os.listdir("test_cluster"):
		if subclade.endswith('fasta'):
			outfile = f'{subclade}_{percent_id}_output.fasta'
			print(outfile)
			command  = f'cd-hit -i {in_dir}/{subclade} -o Clusters_{percent_id}/{outfile} -c 0.99'
			

			subprocess.run([command], shell=True)


for percent_id in params: # go through clustering parameters
	cluster(percent_id)