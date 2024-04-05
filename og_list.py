import os
from pathlib import Path



with open("epi_ogs.txt", 'r') as ogs:#txt file of all ogs
	ogs_list = [i.strip('\n') for i in ogs.readlines()]

print(ogs_list)



n=20

final = [ogs_list[i * n:(i + 1) * n] for i in range((len(ogs_list) + n - 1) // n )]  

print(final)

batch = 0
for i in final:
	batch += 1


	Path(f'ogs_{batch}_output').mkdir(parents=True, exist_ok=True)#makes output folder

	with open(f'ogs_{batch}.txt', 'w') as o:
		for og in i:
			o.write(f'{og}\n')


