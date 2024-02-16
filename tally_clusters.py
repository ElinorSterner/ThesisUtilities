import os, itertools, re


input_dir = 'Clusters_0.99'

#iterate through .clsr files

all_clades = []


for clstr_file in os.listdir(input_dir):
	if clstr_file.endswith("clstr"):

		og = clstr_file.split('.')[0].split("Only_")[1]
		#print(og)

		clusters_list = []
		with open(f'{input_dir}/{clstr_file}', 'r') as clusters:
			for line in clusters.readlines():

				seq = ''
				
				if line[0] == ">":
					cluster_name = line.strip('\n')
				
				else:
					seq = line

				clusters_list.append(f'{cluster_name}, {seq}')

			grouped_clust_data = itertools.groupby(clusters_list, lambda x : x[:10]) 
  
			for cluster_number, sequences in grouped_clust_data: 
				whole_cluster = {cluster_number : list(sequences)} 
				all_clades.append(whole_cluster)


all_clade_data = []
for clade in all_clades:
	for cluster, items in clade.items():

		cluster_id = cluster

		jv_count = 0
		ad_count = 0
		sz_count = 0

		for item in items:
			cluster_id = item[1:10]


			item_list = item.replace('... ', '_').split(',')

			if len(item_list) > 2:
				
				cell = item_list[2][2:12]

				if 'jv' in cell:
					jv_count += 1
				if 'ad' in cell:
					ad_count += 1
				if "sz" in cell:
					sz_count += 1

				og = item_list[2].split('_')[10]  
			
			else:
				cell = " "
				og = " "
		
		porp_of_jv_in_cluster = str(round(jv_count/13, 2))
		porp_of_ad_in_cluster = str(round(ad_count/12, 2))
		porp_of_sz_in_cluster = str(round(sz_count/9, 2))


	to_write_per_cluster = f'{cluster_id}, OG6_{og}, ,{jv_count}, {porp_of_jv_in_cluster}, {ad_count}, {porp_of_ad_in_cluster}, {sz_count}, {porp_of_ad_in_cluster}'
	all_clade_data.append(to_write_per_cluster)
		#print(cluster_id, og, jv_count, porp_of_jv_in_cluster, ad_count, porp_of_ad_in_cluster)


						
with open("clusters_meta_data.csv", 'w') as o:
	
	header = ("cluster, OG, , jv_count, porp_of_jv_in_cluster, ad_count, porp_of_ad_in_cluster, sz_count, porp_of_ad_in_cluster\n")
	o.write(header)

	for cluster in all_clade_data:
		print(cluster)
		o.write(f'{cluster}\n')


		

					




