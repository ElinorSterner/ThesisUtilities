#Written Feb 2024
#input: Subtrees unaligned

import os
from Bio import SeqIO
import statistics


	
indir = "Subtrees_unaligned_At_only"
outlabel = "epi_clades"

def TallyOGs():

	all_clade_data = []
	all_coverages = []
	all_data = []
	
	for file in os.listdir(indir):

		if file.endswith('.fasta'):
			#initiate counts for all the cell types for each cluster that we parse through
			seen = set() # initiate empty set for the seen cells to be added to for each clade
			jv_count, ad_count, sz_count = 0, 0, 0
			multi_count, uni_count, no_data_count = 0, 0, 0
			jv_coverage, ad_coverage, sz_coverage = ([] for i in range(3))
			multi_coverage, uni_coverage, nd_coverage = ([] for i in range(3))

			# file name format for each postguidance file: OG6_100120_postGuidance_preTrimAl_unaligned_masked_95.fas_0.fasta
			og = file.split("_postGuidance")[0]
			clade = file.split("fas_")[1].split(".fasta")[0]

			tree_code = f"{og}_{clade}"

			records = list(SeqIO.parse(f'{indir}/{file}', 'fasta'))


			for record in records:

				cell = record.id.split('_XX')[0]
				coverage = record.id.split('Cv')[1].split('_')[0]

				if cell not in seen:
					seen.add(cell)

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

					if "sz" not in cell:
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

		# proportions of each lifestage counting ALL cells, round off at 2
		prop_of_jv_in_cluster = str(round(jv_count/18, 2))
		prop_of_ad_in_cluster = str(round(ad_count/18, 2))
		prop_of_sz_in_cluster = str(round(sz_count/10, 2))

		# proportions of the totals in adults and juveniles, round off at 2
		prop_of_mu_in_cluster = str(round(multi_count/10, 2))
		prop_of_un_in_cluster = str(round(uni_count/13, 2))
		prop_of_nd_in_cluster = str(round(no_data_count/13, 2))


		mdn_cov_jv = statistics.median(jv_coverage)
		mdn_cov_ad = statistics.median(ad_coverage)
		mdn_cov_sz = statistics.median(sz_coverage)

		mdn_cov_multi = statistics.median(multi_coverage)
		mdn_cov_uni = statistics.median(uni_coverage)
		mdn_cov_nd = statistics.median(nd_coverage)


		to_write_per_tree = f'{tree_code}, , {prop_of_jv_in_cluster},  {prop_of_ad_in_cluster},  {prop_of_sz_in_cluster}, ,{prop_of_mu_in_cluster}, {prop_of_un_in_cluster}, {prop_of_nd_in_cluster}, ,{jv_count},{ad_count},{sz_count}, ,{multi_count}, {uni_count}, {no_data_count}'
		all_clade_data.append(to_write_per_tree)


		coverage_per_tree = f'{tree_code}, , {mdn_cov_jv}, {mdn_cov_ad}, {mdn_cov_sz}, ,{mdn_cov_multi}, {mdn_cov_uni}, {mdn_cov_nd}'
		all_coverages.append(coverage_per_tree)


		prop_cov_per_tree = f'{tree_code}, , {prop_of_jv_in_cluster}, {prop_of_ad_in_cluster}, {prop_of_sz_in_cluster}, {prop_of_mu_in_cluster}, {prop_of_un_in_cluster}, {prop_of_nd_in_cluster}, ,{mdn_cov_jv}, {mdn_cov_ad}, {mdn_cov_sz}, {mdn_cov_multi}, {mdn_cov_uni}, {mdn_cov_nd}'
		all_data.append(prop_cov_per_tree)



	write_out_props(all_clade_data)
	write_out_covs(all_coverages)
	write_out_all(all_data)

def write_out_props(all_clade_data):

	with open(f"{outlabel}_prop.csv", 'w') as o:
		
		header = ("OG, , prop_of_jv_in_cluster,  prop_of_ad_in_cluster,  prop_of_sz_in_cluster, ,prop_of_mu_in_cluster, prop_of_un_in_cluster, prop_of_nd_in_cluster, ,jv_count, ad_count, sz_count, ,multi_count, uni_count, no_data_count\n")


		o.write(header)

		for tree in all_clade_data:
			o.write(f'{tree}\n')


def write_out_covs(all_coverages):

	with open(f"{outlabel}_coverages.csv", 'w') as o:
		
		header = ("OG, , mdn_cov_jv, mdn_cov_ad, mdn_cov_sz, ,mdn_cov_multi, mdn_cov_uni, mdn_cov_nd\n")
		o.write(header)

		for tree in all_coverages:
			o.write(f'{tree}\n')


def write_out_all(all_data):

	with open(f"{outlabel}_prop_coverages.csv", 'w') as o:
		
		header = ('OG, , prop_of_jv_in_cluster, prop_of_ad_in_cluster, prop_of_sz_in_cluster, prop_of_mu_in_cluster, prop_of_un_in_cluster, prop_of_nd_in_cluster, ,mdn_cov_jv, mdn_cov_ad, mdn_cov_sz, mdn_cov_multi, mdn_cov_uni, mdn_cov_nd\n')

		o.write(header)

		for tree in all_data:
			o.write(f'{tree}\n')


				


TallyOGs()

