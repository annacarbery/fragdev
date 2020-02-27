import os
from parse_functions import parse_residues


DATA_DIRECTORY = os.path.abspath('data')


for dir in os.listdir(DATA_DIRECTORY):
    pdb_files = []
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        pdb = os.path.join(DATA_DIRECTORY, dir, file)
        pdb_files.append(pdb)
    residues = parse_residues(pdb_files)
    clusters = {}
    for i in residues:
        this_residue_cluster = {}
        clusters[i.object_list[0].object_desc] = []
        for j in range(len(i.object_list[0].value_array)):
            this_residue_cluster[j] = [j]
            for k in range(len(i.object_list[0].value_array[j])):
                if i.object_list[0].value_array[j][k] < 2.5:
                    if j in this_residue_cluster:
                        if k < j:
                            this_residue_cluster[j].append(k)
                        else:
                            this_residue_cluster[j].append(k+1)
        if sum([len(this_residue_cluster[i]) for i in this_residue_cluster]) / 15 == 15.0:
           # print(i.object_list[0].object_desc, 'all in same cluster')
            clusters[i.object_list[0].object_desc].append(this_residue_cluster[0])
        else:
            for n in range(len(this_residue_cluster)):
                for m in range(len(this_residue_cluster)):
                    if sorted(this_residue_cluster[n]) == sorted(this_residue_cluster[m]):
                        if sorted(this_residue_cluster[n]) not in clusters[i.object_list[0].object_desc]:
                            clusters[i.object_list[0].object_desc].append(sorted(this_residue_cluster[n]))
    for i in clusters:
        print(i, clusters[i])



