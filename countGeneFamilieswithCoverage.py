#Written Feb 2024
#input: Subtrees unaligned

import os
from Bio import SeqIO


	
indir = "Subtrees_unaligned_At_only"

def TallyOGs():

	all_clade_data = []
	all_coverages = []
	
	for file in os.listdir(indir):

		if file.endswith('.fasta'):
			#initiate counts for all the cell types for each cluster that we parse through
			jv_count, ad_count, sz_count = 0, 0, 0
			multi_count, uni_count, no_data_count = 0, 0, 0
			jv_coverage, ad_coverage, sz_coverage = ([] for i in range(3))
			multi_coverage, uni_coverage, nd_coverage = ([] for i in range(3))


			og = file.split(".95")[0]
			clade = file.split("fasta_")[1].split(".fasta")[0]

			tree_code = f"{og}_{clade}"

			records = list(SeqIO.parse(f'{indir}/{file}', 'fasta'))


			for record in records:

				cell = record.id.split('_XX')[0]
				coverage = record.id.split('Cv')[1].split('_')[0]

				#tally lifestages
				if 'jv' in cell:
					jv_count += 1
					jv_coverage.append(int(record.id.split('Cv')[1].split('_')[0]))
				if 'ad' in cell:
					ad_count += 1
					ad_coverage.append(int(record.id.split('Cv')[1].split('_')[0]))

				if "sz" in cell:
					sz_count += 1
					sz_coverage.append(int(record.id.split('Cv')[1].split('_')[0]))

				#tally nuclear types
				if 'Mu' in cell:
					multi_count += 1
					multi_coverage.append(int(record.id.split('Cv')[1].split('_')[0]))
				if 'Un' in cell:
					uni_count += 1
					uni_coverage.append(int(record.id.split('Cv')[1].split('_')[0]))
				if 'Nd' in cell:
					no_data_count += 1
					nd_coverage.append(int(record.id.split('Cv')[1].split('_')[0]))

		porp_of_jv_in_cluster = str(round(jv_count/13, 2))
		porp_of_ad_in_cluster = str(round(ad_count/12, 2))
		porp_of_sz_in_cluster = str(round(sz_count/9, 2))

		porp_of_mu_in_cluster = str(round(multi_count/8, 2))
		porp_of_un_in_cluster = str(round(uni_count/6, 2))
		porp_of_nd_in_cluster = str(round(no_data_count/20, 2))


		avg_cov_jv = round(sum(jv_coverage) / len(jv_coverage))
		avg_cov_ad = round(sum(ad_coverage) / len(ad_coverage))
		avg_cov_sz = round(sum(sz_coverage) / len(sz_coverage))

		avg_cov_multi = round(sum(multi_coverage) / len(multi_coverage))
		avg_cov_uni = round(sum(uni_coverage) / len(uni_coverage))
		avg_cov_nd = round(sum(nd_coverage) / len(nd_coverage))


		to_write_per_tree = f'{tree_code}, , {porp_of_jv_in_cluster},  {porp_of_ad_in_cluster},  {porp_of_sz_in_cluster}, ,{porp_of_mu_in_cluster}, {porp_of_un_in_cluster}, {porp_of_nd_in_cluster}, ,{jv_count},{ad_count},{sz_count}, ,{multi_count}, {uni_count}, {no_data_count}'
		all_clade_data.append(to_write_per_tree)


		coverage_per_tree = f'{tree_code}, , {avg_cov_jv}, {avg_cov_ad}, {avg_cov_sz}, ,{avg_cov_multi}, {avg_cov_uni}, {avg_cov_nd}, , {jv_coverage}, {ad_coverage}, {sz_coverage}, {multi_coverage}, {uni_coverage}, {nd_coverage}'
		all_coverages.append(coverage_per_tree)



	#write_out_props(all_clade_data)
	write_out_covs(all_coverages)

def write_out_props(all_clade_data):

	with open(f"{indir}_meta_data.csv", 'w') as o:
		
		header = ("OG, , porp_of_jv_in_cluster,  porp_of_ad_in_cluster,  porp_of_sz_in_cluster, ,porp_of_mu_in_cluster, porp_of_un_in_cluster, porp_of_nd_in_cluster, ,jv_count, ad_count, sz_count, ,multi_count, uni_count, no_data_count\n")


		o.write(header)

		for tree in all_clade_data:
			o.write(f'{tree}\n')


def write_out_covs(all_coverages):

	with open(f"{indir}_coverages.csv", 'w') as o:
		
		header = ("OG, , avg_cov_jv, avg_cov_ad, avg_cov_sz, ,avg_cov_multi, avg_cov_uni, avg_cov_nd, , jv_coverage, ad_coverage, sz_coverage, multi_coverage, uni_coverage, nd_coverage\n")
		o.write(header)

		for tree in all_coverages:
			o.write(f'{tree}\n')


				


TallyOGs()

