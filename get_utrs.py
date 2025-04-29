#Elinor, April 2025
#Motivation: Get UTRs matching sequences post-clade grabbing
#Intent: Get full length transcripts (including UTRs) that match the sequences in a fasta or directory of fastas from EukPhylo 
#Dependencies:
#Inputs: fasta or directory of fastas of references, and the UTRs from here: EukPhylo_RO/Phylotol6_w_hook6.6/Raw_PTL1_Output/LKH/Plate*/Intermediate/XlaneBleeding/
#Outputs: UTRs that are renamed in the format of R2Gs
#Example: 

import os
from Bio import SeqIO
from pathlib import Path
import sys





'''
1. rename files to the life stages 
2. rename utr files: Ten digit code, length, "T" for transcript. change abreviations 
	 At_mu_Jv33_XX_0_Ct20_L5475_Cv199_E_Cn16_OG6_101144

3. match to the reference
4. Add the OG to the UTR
5. Write them out nicely


'''


def script_help():

	print('\nGet full length transcripts (including UTRs) that match the sequences in a fasta or directory of fastas formatted like R2G \n\nInput:\ntsv file of tab separated ten digit code and taxon info (taxonomy, lifestage, etc).\nAND\n\ndirectory of assemblies renamed in this format: ten_digit_code_assembledTranscripts.fasta\n\nOutput is multiple R plots faceted by taxon and a csv file of data. \n\nIt plots GC by length, and distributions of coverage, length and GC.\n\n To run: \n\n\tpython assess_transcriptomes.py -i <pathway to directory of assemblies>\n\nRun -h or --help for this message\n\n')

def get_args():
	#this parses user arguments. Checks if the files are renamed already or not (--renamed or --raw), and gets the directory of those files.
	
	if('--help' in sys.argv or '-h' in sys.argv):#check for help function in command line
		script_help()
		exit()


	if ('--input'in sys.argv or '-i' in sys.argv):#check for renamed parameter
		renamed = True
		try:
			if('--input' in sys.argv):
				input_dir = sys.argv[sys.argv.index('--input') + 1]
			else:
				input_dir = sys.argv[sys.argv.index('-i') + 1]
		except IndexError:
			print('\nSomething went wrong went parsing the arguments. Did you input a directory of assemblies?\n')

	if ('--reference'in sys.argv or '-r' in sys.argv):#check for renamed parameter
		renamed = True
		try:
			if('--reference' in sys.argv):
				reference_dir = sys.argv[sys.argv.index('--reference') + 1]
			else:
				reference_dir = sys.argv[sys.argv.index('-r') + 1]
		except IndexError:
			print('\nSomething went wrong went parsing the arguments. Did you input a directory of reference sequences?\n')


	rename_sequences(input_dir, reference_dir)


def rename_sequences(input_dir, reference_dir):

	# renames stuff, saves intermediate fastas 
	print('\n\nRenaming transcripts, saving all to intermediate folder\n\n')

	#life stage 
	Path(f'intermediate_{input_dir}').mkdir(parents=True, exist_ok=True)#makes output folder
	out = f'intermediate_{input_dir}'

	with open('tip_labels.csv', 'r') as o:
		cell_data = [i.strip('\n') for i in o.readlines()]
		cell_dict = {i.split(', ')[0] : i.split(', ')[1] for i in cell_data}


	for file in os.listdir(input_dir):
		new_rec = {}
		
		if file.endswith('fasta'):
			path = input_dir +'/'+ file
			records = list(SeqIO.parse(path, 'fasta'))
			taxa = file[0:10]

			for record in records:
				
				new_id = cell_dict[taxa] + "_" + record.id
			
				new_meta = record.id.split('_')
				contig = "Ct" +new_meta[1]
				length = 'L' +new_meta[2].strip('Len')
				if "Cov" in record.id:
					cov = 'Cv'+new_meta[3].strip('Cov')
				else:
					cov = 'CvNA'
				new_id = f'{cell_dict[taxa]}_XX_T_{contig}_{length}_{cov}' # replace the _XX_0_ with _XX_T_
				new_seq = record.seq
				new_rec.update({new_id : new_seq})

		new_file = f"{'_'.join(new_id.split('_')[0:3])}.200bp.renamed.fasta"
		outpath = f"intermediate_{input_dir}/{new_file}"

		with open(outpath, 'w') as o:
			for idd, seq in new_rec.items():
				o.write(f'>{idd}\n{seq}\n')

	matching_reference(out, reference_dir)


# match the renamed full-length transcripts to the ready-to-go format sequences
def matching_reference(outpath, reference_dir):


	print('\n\nMatching transcripts to reference...\n\n')



	query_record_ids = {}
	query_recs = {}

	for file in os.listdir(reference_dir):
		#print(file)
		if file.endswith('fasta'):
			records = list(SeqIO.parse(f'{reference_dir}/{file}', 'fasta'))
			#print(records)

			for rec in records:
				query_record_ids.update({rec.id : rec.seq})

	
	count = 0 # counting number of matches

	renamed_transcripts = []
	for file in os.listdir(outpath):
		if file.endswith('fasta'):
			records = list(SeqIO.parse(f'{outpath}/{file}', 'fasta'))

			for rec in records:
				taxon = rec.id[:10]
				contig = rec.id.split('_')[5]
				
				name_to_match = f'{taxon}_XX_0_{contig}' #match the format of the beginning of eukphylo sequences

				# match and rename sequences
				for key in query_record_ids.keys():
					if name_to_match in key: 
						count += 1						
						length = len(rec.seq)
					
						#for renaming: change 0 to T and add real length of transcript (they already are)
						name_pieces_to_keep = key.split('_Cv')[1]
						seq_id = f'{taxon}_XX_T_{contig}_L{length}_Cv{name_pieces_to_keep}'
						renamed_transcripts.append(f'>{seq_id}\n{rec.seq}\n')

	print(f'\n\nThere are {count} sequences in your reference that have matching transcripts\n\n')


	write_outputs(renamed_transcripts)



def write_outputs(renamed_transcripts):

	print('\n\nWriting matched transcripts to matched_transcripts directory\n\n')


	Path(f'matched_transcripts').mkdir(parents=True, exist_ok=True)#makes output folder

	# all sequences
	with open('matched_transcripts/matched_transcripts.fasta', 'w') as o:
		for i in renamed_transcripts:
			o.write(i)

	#by taxon
	records = list(SeqIO.parse('matched_transcripts/matched_transcripts.fasta', 'fasta'))
	for rec in records:
		cell = rec.id[:10]
		path = f'matched_transcripts/{cell}_transcripts.fasta'
		if os.path.exists(path) == True:
			with open(path, 'a') as o:
				o.write(f'>{rec.id}\n{rec.seq}\n')
		else:
			with open(path, 'w') as o:
				o.write(f'>{rec.id}\n{rec.seq}\n')
	print('\n\nDone! look for output in the matched_transcripts directory')



if __name__ == '__main__':
	all_data = {}
	
	try:
		f = open("tip_labels.csv")

	except FileNotFoundError:
		print('\n\n Did you include a tsv file of your cells? It should be called new_names.tsv and formatted like this: 10_digit_code\tnew 10 digit code\n\n')	
		exit()

	else:
		get_args()


