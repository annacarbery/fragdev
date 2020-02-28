import json

res = json.load(open('water_clusters.json', 'r'))
print(len(res['0']['mol_ids']))