import scipy.cluster.hierarchy as hc
from scipy.sparse import coo_matrix
from os.path import join as pjoin
from tqdm import tqdm
from scipy.spatial.distance import squareform
import json
import pandas
import os
import fastcluster
with open(pjoin("data", "dbs", "sra_data.r1.json")) as handle:
    sra_md = json.load(handle)
sra_md = {v['SRA_ID'] :v for v in sra_md.values()}

bins_per_sra = {s : 0 for s in sra_md}
for l in tqdm(os.listdir("all_bins")):
    bins_per_sra[l.split(".")[0]] += 1

with open("data/dbs/pulled_samples.json") as handle:
    pulled = set(json.load(handle))

with open("data/dbs/sra2rrnafreq.tsv") as handle:
    sra2rrna = { l.split()[0] : float(l.strip().split()[1]) for l in handle }

too_much_rrna = {k for k,v in sra2rrna.items() if v > 5}

with open(pjoin("data/dbs/all_similarities.tsv")) as handle:
    simis = list(handle.readlines())
    simis = { ( l.split()[0], l.split()[1] ) : float(l.split()[2][:-1]) for l in tqdm(simis)}

simis = {k : v for k, v in tqdm(simis.items()) if k[0] not in too_much_rrna and k[1] not in too_much_rrna and k[1] in sra_md and k[0] in sra_md}
all_samples = {a for k in simis.keys() for a in k}
for k in all_samples:
    simis[k,k] = 1

    #fixing a few missing simetric simis

for k in tqdm(list(simis)):
    if simis.get((k[1], k[0])) != simis[k] :
        simis[(k[1], k[0])] = simis[k]


index2sample = {i : a for i,a in enumerate(all_samples)}
sample2index = {v : k for k,v in index2sample.items()}
data = []
i = []
j = []
for k,v in simis.items():
    data += [float(v)]
    i += [sample2index[k[0]]]
    j += [sample2index[k[1]]]
sparse_mat = coo_matrix((data,(i,j)))
Y = sparse_mat.toarray()
X = squareform(1-Y)
Y = pandas.DataFrame(Y)

linked = fastcluster.linkage(X, method='ward')
rootnode, nodelist = hc.to_tree(linked, rd = True)

sub_tree_leaflists = dict()
sub_tree_Gbases = dict()
sub_tree_simis = dict()

def sub_tree_leaves(node):
    if node.id not in sub_tree_leaflists:
        if node.is_leaf() :
            sub_tree_leaflists[node.id] = [index2sample[node.id]]
        else :
            sub_tree_leaflists[node.id] = sub_tree_leaves(node.get_left()) + sub_tree_leaves(node.get_right())
    return sub_tree_leaflists[node.id]

def get_sub_tree_Gbases(node):
    if node.id not in sub_tree_Gbases:
        samples = sub_tree_leaves(node)
        total_bases = sum([int(sra_md[s]['nb_bases']) for s in samples])
        sub_tree_Gbases[node.id] = total_bases/1000000000
    return sub_tree_Gbases[node.id]

def get_sub_tree_mean_simi(node):
    if node.id not in sub_tree_simis:
        if len(sub_tree_leaves(node)) > 1:
            index = [sample2index[s] for s in sub_tree_leaves(node)]
            submatrix = Y.loc[index,index]
            submatrix.diag = 0
            sub_tree_simis[node.id] = (submatrix.sum().sum()-len(submatrix))/(len(submatrix)*len(submatrix)-len(submatrix))
        else :
            sub_tree_simis[node.id] = None
    return sub_tree_simis[node.id]

def split_tree(node, cutoff = 500):
    if node.is_leaf():
        return [node]
    left = node.get_left()
    right = node.get_right()
    lsize = get_sub_tree_Gbases(left)
    rsize = get_sub_tree_Gbases(right)
    if lsize < cutoff :
        left_clusters = [left]
    else :
        left_clusters = split_tree(left, cutoff)
    if rsize < cutoff :
        right_clusters = [right]
    else :
        right_clusters = split_tree(right, cutoff)

    return left_clusters + right_clusters

clsts = split_tree(rootnode)
filt_clsts = [ c for c in clsts if get_sub_tree_Gbases(c) > 30 and len(sub_tree_leaves(c)) > 1 and get_sub_tree_mean_simi(c) > 0]
unclusts = all_samples.difference(set(sum([ sub_tree_leaves(c) for c in filt_clsts],[])))
cs_meta = {c.id : {
                    'mean_simi' : get_sub_tree_mean_simi(c),
                    'raw_bases' : get_sub_tree_Gbases(c) ,
                    'nb_samples' : len(sub_tree_leaves(c)),
                    'bins_in_singles' : sum([bins_per_sra[cc] for cc in sub_tree_leaves(c)]),
                    'nb_SRA_projects' : len({sra_md[cc]['study_id'] for cc in sub_tree_leaves(c) if 'study_id' in sra_md[cc]}),
                    'contains_pulled' : any([cc in pulled for cc in sub_tree_leaves(c)]),
                    'clustering_type' : "Hierarchical clustering pruned at 500Gb nodes",
                    'samples'   : sub_tree_leaves(c),
                    } for c in filt_clsts}
i = 1
j = 1
full_cs_meta = {}
for v in cs_meta.values():
    if v['contains_pulled']:
        full_cs_meta['zEBG_HC500_coass_{:03}'.format(i)] = v
        i += 1
    else :
        full_cs_meta['HC500_coass_{:03}'.format(j)] = v
        j += 1


with open("data/dbs/coassemblies.json", "w") as handle:
    json.dump(full_cs_meta, handle, sort_keys = True, indent = 2)
