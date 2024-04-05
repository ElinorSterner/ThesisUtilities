# Written by Elinor Spring 24 to tally cells by lifestage and nuclear type (wirtten into fake ten digit codes)
# takes in .clstr files, tallies the number of cells up!
# ASK me a q! : esterner27@gmail.com

import os, itertools, re


input_dir = 'Clusters_0.99'

#iterate through .clsr files

all_clades = []


for clstr_file in os.listdir(input_dir):
	if clstr_file.endswith("clstr"):

		og = clstr_file.split('.')[0].split("Only_")[1]
		clade_number  = clstr_file.split('.')[2].split("fasta_")[1]
		clade_data = f'{og}_{clade_number}'
		#print(clade_data)

		clusters_list = []
		with open(f'{input_dir}/{clstr_file}', 'r') as clusters:
			for line in clusters.readlines():

				seq = ''#initiate seq list
				
				if line[0] == ">":#identify cluster labels with >
					cluster_name = line.strip('\n')
				
				else:
					line = line.strip('\n')#remove new line from the end so it can append the og and clade number at the end
					seq = f'{line}+{clade_data}' # join with a unique connector (+)to make spliting them easy later

				clusters_list.append(f'{cluster_name}, {seq}')#save cluster names with list of sequences into a dictionary

			grouped_clust_data = itertools.groupby(clusters_list, lambda x : x[:10]) #use groupby() to reformat into items connected by shared cluster number (x[:10]). this is an encoded data type, parsed out in the next few lines
  
			for cluster_number, sequences in grouped_clust_data: 
				whole_cluster = {cluster_number : list(sequences)} # make into a dictionary with sequences with the same cluster number together
				all_clades.append(whole_cluster) # save all clades together, will be easy to parse out next


all_clade_data = []#initiate list for data to write out to csv

for clade in all_clades:# for each clstr file in the directory
	for cluster, items in clade.items():

		#print(cluster)


		cluster_id = cluster[1:10] # get the name of the cluster from the dictionary key 

		seen = set() # initiate empty set for the seen cells to be added to for each cluster
		#initiate counts for all the cell types for each cluster that we parse through
		jv_count, ad_count, sz_count = 0, 0, 0
		multi_count, uni_count, no_data_count = 0, 0, 0

		for item in items:#for sequence found in each cluster (this is parsing through a list)

			item_list = item.replace('... ', '_').split(',') # break apart the item (a string), which contains the whole line of the cluster file. it has "..." and spaces so its awkward to work with
			if len(item_list) > 2:#if there are sequences in the line (the first line of each list is len = 1 because its ONLY the cluster id)
				
				cell = item_list[2][2:12]# ten digit code
				if cell not in seen:
					seen.add(cell)# add cell to "seen" now that it has been counted

					#add cell to lifestage tallies if it has not been counted yet 
					if 'jv' in cell:
						jv_count += 1
					if 'ad' in cell:
						ad_count += 1
					if "sz" in cell:
						sz_count += 1

					# for the nucleus counts, if it is NOT a schizont, do not count it (since there are so few schizonts with nuclear data)
					if "sz" not in cell:

						#tally nuclear types
						if 'Mu' in cell:
							multi_count += 1
						if 'Un' in cell:
							uni_count += 1
						if 'Nd' in cell:
							no_data_count += 1

				#get the og and clade number (from clade grabbing)
				og = item_list[2].split('+')[-1]
			
			else:
				cell = " "
				og = " "

		# proportions of each lifestage counting ALL cells, round off at 2
		porp_of_jv_in_cluster = str(round(jv_count/18, 2))
		porp_of_ad_in_cluster = str(round(ad_count/18, 2))
		porp_of_sz_in_cluster = str(round(sz_count/10, 2))

		# proportions of the totals in adults and juveniles, round off at 2
		porp_of_mu_in_cluster = str(round(multi_count/10, 2))
		porp_of_un_in_cluster = str(round(uni_count/13, 2))
		porp_of_nd_in_cluster = str(round(no_data_count/13, 2))




	to_write_per_cluster = f'{og}, {cluster_id}, , {porp_of_jv_in_cluster},  {porp_of_ad_in_cluster},  {porp_of_sz_in_cluster}, ,{porp_of_mu_in_cluster}, {porp_of_un_in_cluster}, {porp_of_nd_in_cluster}, ,{jv_count},{ad_count},{sz_count}, ,{multi_count}, {uni_count}, {no_data_count}'
	all_clade_data.append(to_write_per_cluster)



with open(f"{input_dir}_meta_data.csv", 'w') as o:
	
	header = ("OG, cluster, , porp_of_jv_in_cluster,  porp_of_ad_in_cluster,  porp_of_sz_in_cluster, ,porp_of_mu_in_cluster, porp_of_un_in_cluster, porp_of_nd_in_cluster, ,jv_count,ad_count,sz_count, ,multi_count, uni_count, no_data_count\n")


	o.write(header)

	for cluster in all_clade_data:
		#print(cluster)
		o.write(f'{cluster}\n')


		

					




