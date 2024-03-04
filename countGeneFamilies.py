#Written Feb 2024
#input: Subtrees unaligned

import os
from Bio import SeqIO


	
indir = "Subtrees_unaligned"

def TallyOGs():

	all_clade_data = []
	
	for file in os.listdir(indir):
		#print(file)

		if file.endswith('.fasta'):
			#initiate counts for all the cell types for each cluster that we parse through
			
			jv_count, ad_count, sz_count = 0, 0, 0
			multi_count, uni_count, no_data_count = 0, 0, 0

			og = file.split(".95")[0]
			clade = file.split("fasta_")[1].split(".fasta")[0]

			tree_code = f"{og}_{clade}"

			records = list(SeqIO.parse(f'{indir}/{file}', 'fasta'))


			for record in records:

				cell = record.id.split('_XX')[0]

				#tally lifestages
				if 'jv' in cell:
					jv_count += 1
				if 'ad' in cell:
					ad_count += 1
				if "sz" in cell:
					sz_count += 1

				#tally nuclear types
				if 'Mu' in cell:
					multi_count += 1
				if 'Un' in cell:
					uni_count += 1
				if 'Nd' in cell:
					no_data_count += 1

		porp_of_jv_in_cluster = str(round(jv_count/13, 2))
		porp_of_ad_in_cluster = str(round(ad_count/12, 2))
		porp_of_sz_in_cluster = str(round(sz_count/9, 2))

		porp_of_mu_in_cluster = str(round(multi_count/8, 2))
		porp_of_un_in_cluster = str(round(uni_count/6, 2))
		porp_of_nd_in_cluster = str(round(no_data_count/20, 2))

		to_write_per_tree = f'{tree_code}, , {porp_of_jv_in_cluster},  {porp_of_ad_in_cluster},  {porp_of_sz_in_cluster}, ,{porp_of_mu_in_cluster}, {porp_of_un_in_cluster}, {porp_of_nd_in_cluster}, ,{jv_count},{ad_count},{sz_count}, ,{multi_count}, {uni_count}, {no_data_count}'
		all_clade_data.append(to_write_per_tree)


	write_out(all_clade_data)

def write_out(all_clade_data):

	with open(f"{indir}_meta_data.csv", 'w') as o:
		
		header = ("OG, , porp_of_jv_in_cluster,  porp_of_ad_in_cluster,  porp_of_sz_in_cluster, ,porp_of_mu_in_cluster, porp_of_un_in_cluster, porp_of_nd_in_cluster, ,jv_count,ad_count,sz_count, ,multi_count, uni_count, no_data_count\n")


		o.write(header)

		for tree in all_clade_data:
			o.write(f'{tree}\n')



				


TallyOGs()

