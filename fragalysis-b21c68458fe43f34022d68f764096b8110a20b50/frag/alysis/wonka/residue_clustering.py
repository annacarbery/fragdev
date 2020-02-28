import os
from cluster_functions import cluster_dp

DATA_DIRECTORY = os.path.abspath('data')

def get_res(pdb):
    res_num = []
    all_res = {}
    for line in pdb:
        if line.startswith('ATOM'):
            res_num.append(int(line[23:26]))
    for i in range(min(res_num), max(res_num)+1):
        x = 0
        y = 0
        z = 0
        n = 0
        for line in pdb:
            if line.startswith('ATOM'):
                if int(line[23:26]) == i:
                    residue = line[17:26]
                    x += float(line[30:38].strip())
                    y += float(line[38:46].strip())
                    z += float(line[46:54].strip())
                    n += 1
        try:
            all_res[residue] = [x/n, y/n, z/n]
        except ZeroDivisionError:
            pass
    return all_res




for dir in os.listdir(DATA_DIRECTORY):
    coms = {}
    idents = []
    lam = 2.5
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        pdb = open(os.path.join(DATA_DIRECTORY, dir, file), 'r').readlines()
        coms[file] = get_res(pdb)
        idents.append(file)
    for i in coms[idents[0]]:
        print(i)
        vectors = []
        for j in coms:
            vectors.append(coms[j][i])
        clusters = cluster_dp(vectors, lam, idents)
        print(len(clusters), clusters)

