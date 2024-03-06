# written by Elinor March 2024 to get the sequences in grabbed clades from assembled transcripts to be able to look at UTRs. Writes them out by OG and by taxa.
# input: 
#		folder of relevant assembledTranscripts (assigned to "assemblies_dir" variable)
#		tip_labels.csv, a csv of 10 digit codes comma new 10 digit code (if they are unchagned make both columns identical)
#		folder of sequences that made it into trees/postguidance files, (called Subtrees_unaligned if you ran clade grabbing). (assigned to "grabbed_clades" variable)

import os
from Bio import SeqIO
from pathlib import Path

assemblies_dir = "assembled_transcripts"
grabbed_clades = "Subtrees_unaligned"


Path(f'RenamedTranscripts').mkdir(parents=True, exist_ok=True)#makes output folder for renamed transcripts
Path(f'MatchedTranscripts').mkdir(parents=True, exist_ok=True)#makes output folder for transcripts that are in subtrees

def rename_assemblies():

	with open('tip_labels.csv', 'r') as o:
		cell_data = [i.strip('\n') for i in o.readlines()]
		cell_dict = {i.split(', ')[0] : i.split(', ')[1] for i in cell_data}


	for taxa in os.listdir(assemblies_dir):
		all_recs = {}
		records = list(SeqIO.parse(f'{assemblies_dir}/{taxa}', 'fasta'))
		
		for record in records:
			contig = record.id.split("_")[1]
			length = record.id.split("_")[3]
			covrge = round(float(record.id.split("_")[5]))
			tdcode = taxa[:10]

			if tdcode in cell_dict:
				new_tdcode = cell_dict[tdcode]

			new_id =  f"{new_tdcode}_XX_0_Ct{contig}_L{length}_Cv{covrge}"
			all_recs.update({new_id:record.seq})

		file_name = f"{new_tdcode}_assembledTranscripts.fasta"
		with open(f'RenamedTranscripts/{file_name}', 'w') as o:
			for idd, seq in all_recs.items():
				o.write(f">{idd}\n{seq}\n")


def parse_subtrees():

	all_seqs_dicts = {}

	for og in os.listdir(grabbed_clades):
		og_code = og.split('.95gap')[0]
		clade = og.split('_')[-1].strip(".fasta")
		og_clade = og_code + "_" + clade

		records = list(SeqIO.parse(f'{grabbed_clades}/{og}', 'fasta'))
		for record in records:
			contig = record.id.split("Ct")[1].split("_")[0]
			taxa = record.id[:10]

			all_seqs_dicts.update({f'{taxa}-{contig}' : og_clade})

	find_matches(all_seqs_dicts)# find matches for everything in dictionary, which is all the seuqnces in trees



def find_matches(all_seqs_dicts):

	path = "MatchedTranscripts" #output path

	
	for assembly in os.listdir("RenamedTranscripts"):
		print(f"Finding matches in {assembly}....\n")
		transcripts_to_keep = {} # initiate dictionary of matching transcripts
		with open(f'RenamedTranscripts/{assembly}', 'r') as assem:
			
			records = list(SeqIO.parse(f'RenamedTranscripts/{assembly}', 'fasta'))
			for record in records:
				contig = record.id.split("Ct")[1].split("_")[0]# get the contig number of each sequence
				taxa = assembly[:10] # get the taxon

				for info, og in all_seqs_dicts.items(): #iterate though the dictionary of sequences in trees (made in parse_subtrees())
					#if the contig is in the sequence
					if contig == info.split('-')[1]:
						# if the taxa in trees with that contig matches the contig code in the transcripts (because the contig codes are not unique)
						if info.split('-')[0] == taxa:
							new_rec_id = f'{record.id}_{info}_{og}' # make a new name for the record that includes the OG assigned in phylotol part 2
							transcripts_to_keep.update({new_rec_id : record.seq}) # add these to dictionary

		# write out transcripts to keep for each cell
		outpath = f'{path}/{taxa}_selectedTranscripts.fasta'
		with open(outpath, 'w') as o:
			for idd, seq in transcripts_to_keep.items():
				o.write(f'>{idd}\n{seq}\n')


# makes another folder of the filtered assembled Transcripts but written out by OG
def transcripts_by_og():

	print('Writing filtered transcripts out by OG\n')
	Path(f'MatchedTranscriptsByOG').mkdir(parents=True, exist_ok=True)#makes output folder for renamed transcripts
	path = "MatchedTranscripts"

	seq_ids = [] # all the record ids  in all the filtered transcripts files
	all_seqs = {}# all the sequences in all the filtered transcripts files
	for file in os.listdir(path):
		if file.endswith(".fasta"):

			for rec in SeqIO.parse(path + '/' + file, 'fasta'):
				seq_ids.append(rec.id)
				all_seqs.update({ rec.id : str(rec.seq)})


		#Write out all of the seqs of interest from each postguidance file, but now separating out by taxon
	for og in list(dict.fromkeys([rec[-12:-2] for rec in seq_ids])):
		with open('MatchedTranscriptsByOG/' + og + '_Transcripts.fasta', 'w') as o:
			for rec in seq_ids:
				if (rec[-12:-2] == og):
					o.write('>' + rec + '\n' + all_seqs[rec] + '\n\n')


if __name__ == '__main__':
	#rename_assemblies()

	# this function calls files that have the sequences in trees, and then grabs those by contig from the assembled transcripts
	#parse_subtrees()
	transcripts_by_og()
	print("done")


