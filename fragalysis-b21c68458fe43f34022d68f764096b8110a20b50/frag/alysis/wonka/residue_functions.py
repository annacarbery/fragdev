from get_functions import get_res, get_all_coords
import sys
from rdkit import Chem
from models import Protein

def register_res(prot, targ, prot_c, tot_prots):
    """Function to register all the resdiues for a protein"""
    try:
        chain = prot.code.split("_")[1]
    except IndexError:
        print("NO CHAIN ID", prot.code)
        chain = None
    # Now get a dict of res
    out_d = get_res(str(prot.pdb_info).split("\n"), chain.upper())
    # Now go through it
    old = -1
    tot = len(out_d)
    for i, item in enumerate(out_d):
        if not out_d[item]:
            continue
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rProtein %d of %d Registering residues %d%% complete..." % (prot_c + 1, tot_prots, old))
            sys.stdout.flush()
        res = Residue.objects.get_or_create(target_id=targ, prot_id=prot, res_num=item[3:], res_name=item[:3])[0]
        # Now get the min and max
        x_coords, y_coords, z_coords = get_all_coords(out_d[item])
        res.x_max = max(x_coords)
        res.x_min = min(x_coords)
        res.y_max = max(y_coords)
        res.y_min = min(y_coords)
        res.z_max = max(z_coords)
        res.z_min = min(z_coords)
        # Now add the pdb information - so we can use this later
        res.pdb_info = Chem.MolToPDBBlock(out_d[item])
        res.save()


def find_targ_rmsd(targ):
    """Function to find all RMSDs against all other RMSDs
    Takes a target
    Returns three dictionaries"""
    # Get all the proteins
    prots = [x for x in Protein.objects.filter(pdb_info__isnull=False, target_id=targ) if str(x.pdb_info) != ""]
    # Register all the residues
    # Print progress here for the loading of the proteins
    tot_prots = len(prots)
    for prot_counter, prot in enumerate(prots):
        register_res(prot, targ, prot_counter, tot_prots)
    print("Registered residues")
    # A dict to hold all the RMSD calcs in the end
    tot_d = {}
    # A dict to keep track of an individual proteins residues
    prot_dict = {}
    for p in prots:
        prot_dict[p.pk] = {}
    clust_d = {}
    res_set = set()
    # Loop through and calculate the RMSDs
    # Print a sensible protein counter
    for prot_c, prot_1 in enumerate(prots):
        sys.stdout.write("\rProtein %d of %d" % (prot_c + 1, tot_prots))
        sys.stdout.flush()
        for prot_2 in prots:
            # Function to find the rmsds between every residue in the protein
            out_d = find_all_rmsd(prot_1, prot_2)
            # Now tidy this up
            for res in out_d:
                res_set.add(res)
                if res in tot_d:
                    tot_d[res].append(out_d[res])
                else:
                    tot_d[res] = [out_d[res]]
                # Log this into the dictionary for all proteins
                if res in prot_dict[prot_1.pk]:
                    prot_dict[prot_1.pk][res].append(out_d[res])
                else:
                    prot_dict[prot_1.pk][res] = [out_d[res]]
                if res in prot_dict[prot_2.pk]:
                    prot_dict[prot_2.pk][res].append(out_d[res])
                else:
                    prot_dict[prot_2.pk][res] = [out_d[res]]
                # Log this into the dictionary for clustering
                if res in clust_d:
                    if prot_1.pk in clust_d[res]:
                        clust_d[res][prot_1.pk][prot_2.pk] = out_d[res]
                    else:
                        clust_d[res][prot_1.pk] = {prot_2.pk: out_d[res]}
                else:
                    clust_d[res] = {prot_1.pk: {prot_2.pk: out_d[res]}}
    # A dict to find the RMSDs beteeen
    out_d = {}
    # Loop through the residues
    for res in list(res_set):
        out_d[res] = {}
        for prot in prots:
            if prot.pk not in clust_d[res]:
                # If it doesn't exist give it zero for each case
                out_d[res][prot.pk] = [0.0 for x in prots]
                continue
            else:
                out_d[res][prot.pk] = []
            for prot_comp in prots:
                if prot_comp.pk not in clust_d[res][prot.pk]:
                    out_d[res][prot.pk].append(0.0)
                else:
                    if clust_d[res][prot.pk][prot_comp.pk]:
                        out_d[res][prot.pk].append(clust_d[res][prot.pk][prot_comp.pk])
                    else:
                        out_d[res][prot.pk].append(0.0)
    print("Got residue RMSDs")
    return tot_d, prot_dict, out_d


def calc_com(atoms):
    """Function to calculate the unweighted centre of mass"""
    numatoms = len(atoms)
    x_coord = 0.0
    y_coord = 0.0
    z_coord = 0.0
    # Assume all heavy atoms have the same mass
    for coords in atoms:
        x_coord += float(coords[0])
        y_coord += float(coords[1])
        z_coord += float(coords[2])
    return x_coord / numatoms, y_coord / numatoms, z_coord / numatoms


def mark_res_shift(out_d, prot_d, clus_d, target, lam=2.5):
    """Function to make the residue shifts for a target"""
    # Anaylsing reisude shifts
    print("Analysing residue shifts")
    tot = len(out_d)
    old = -1
    for i, res in enumerate(out_d):
        old = i * 100 / tot
        sys.stdout.write("\rSummarising residues %d%% complete..." % old)
        sys.stdout.flush()
        rp = ResShift()
        rp.target_id = target
        rp.res_num = int(res[3:])
        rp.res_name = res[:3]
        rp.max_shift = max(out_d[res])
        if [x for x in out_d[res] if x != 0.0]:
            rp.min_shift = min([x for x in out_d[res] if x != 0.0])
        else:
            rp.min_shift = 0.0
        rp.avg_shift = float(sum([x for x in out_d[res] if x])) / float(len(out_d[res]))
        try:
            rp.validate_unique()
            rp.save()
        except ValidationError:
            # Otherwise updateeverything
            rp = ResShift.objects.get(target_id=target, res_num=int(res[3:]), res_name=res[:3])
            rp.max_shift = max(out_d[res])
            if [x for x in out_d[res] if x != 0.0]:
                rp.min_shift = min([x for x in out_d[res] if x != 0.0])
            else:
                rp.min_shift = 0.0
            rp.avg_shift = float(sum([x for x in out_d[res] if x])) / float(len(out_d[res]))
            rp.save()
    # Now lets relate these to the residue objects
    print("General residues summarised")
    tot = len(prot_d)
    old = -1
    for i, prot in enumerate(prot_d):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rSummarising individual residues %d%% complete..." % old)
            sys.stdout.flush()
        for res in prot_d[prot]:
            res_num = int(res[3:])
            res_name = res[:3]
            # Pull it down and storethis info
            my_res = Residue.objects.filter(prot_id=prot, res_num=res_num, res_name=res_name)
            if my_res:
                my_res = my_res[0]
            else:
                continue
            my_res.max_shift = max(prot_d[prot][res])
            if [x for x in prot_d[prot][res] if x != 0.0]:
                my_res.min_shift = min([x for x in prot_d[prot][res] if x != 0.0])
            else:
                my_res.min_shift = 0.0
            my_res.avg_shift = float(sum([x for x in prot_d[prot][res] if x])) / float(len(prot_d[prot][res]))
            my_res.save()
    print("Individual residues summarised")
    # Now save the residues into clusters
    tot = len(clus_d)
    old = -1
    for i, res in enumerate(clus_d):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rClustering residues %d%% complete..." % old)
            sys.stdout.flush()
        # Now cluster this residue
        dp = dpmeans([clus_d[res][prot] for prot in clus_d[res]], lam, 0, False)
        dp.run()
        # Now save these clusters
        for j, prot in enumerate(clus_d[res]):
            my_res = Residue.objects.filter(prot_id__pk=prot, res_num=int(res[3:]), res_name=res[:3])
            if my_res:
                my_res = my_res[0]
                # Update the cluster id
                my_res.clust_id = dp.dataClusterId[j]
                my_res.save()
    print("Residues clustered")