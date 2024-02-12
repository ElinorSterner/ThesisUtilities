
#Written by Elinor 1/29/24 to rename ready to go files to rename them to include lifestage and nuclear type. Input raw group of ready to go files and file with information about each cell
# Should be a CSV file and the indicators need to be single characters.

import os
from Bio import SeqIO
from pathlib import Path


input_fasta = "unaligned"
Path(f'Renamed_{input_fasta}').mkdir(parents=True, exist_ok=True)#makes output folder


with open('tip_labels.csv', 'r') as o:
	cell_data = [i.strip('\n') for i in o.readlines()]
	cell_dict = {i.split(', ')[0] : i.split(', ')[1] for i in cell_data}


for file in os.listdir(input_fasta):
	new_rec = {}
	
	if file.endswith('fasta'):
		path = input_fasta +'/'+ file
		records = list(SeqIO.parse(path, 'fasta'))
	
		for record in records:
			taxa = record.id[0:10] 
			
			if taxa in cell_dict:
				new_id = cell_dict[taxa] + "_" + record.id
				new_seq = record.seq
				new_rec.update({new_id : new_seq})

			else:
				new_id = record.id
				new_seq = record.seq
				new_rec.update({new_id : new_seq})

	outpath = f"Renamed_{input_fasta}/{file}"

	with open(outpath, 'w') as o:
		for idd, seq in new_rec.items():
			o.write(f'>{idd}\n{seq}\n')







