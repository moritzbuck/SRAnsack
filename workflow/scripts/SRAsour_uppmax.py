import json
from os.path import join as pjoin
from subprocess import call
import os
from tqdm import tqdm
from sourmash import MinHash
from sourmash.signature import SourmashSignature
from sourmash.signature import load_signatures
from multiprocessing import Pool
import sys

run_id = int(os.environ['run_id'])

def load_sig(f):
    with open(f) as handle:
        ll = "".join( handle.readlines()).replace('\n','')
    return list(load_signatures(ll)) [0]

def load_minhash(f):
    with open(f) as handle:
        ll = handle.read()
    ll = [int(t) for t in ll.split('"mins":[')[1].split('],"md5sum"')[0].split(",")]
    return ll


threash = 0.001
nb_matches = 10000

all_sigs = ["sigs/" + s for s in os.listdir('sigs') if s.endswith(".sig")]

block_size = 37

sig_blocks = [all_sigs[i:(i+block_size)] for i in list(range(0,len(all_sigs), block_size))]


block_order = [(i,j) for i,b in enumerate(sig_blocks) for j,v in enumerate(sig_blocks)]

bsize = int(len(block_order)/100)
print("Doing block {} to {} out of {}".format(run_id*bsize, ((run_id+1)*bsize), len(block_order)))
block_order = block_order[run_id*bsize:((run_id+1)*bsize)]
blocks = {a for aa in block_order for a in aa}
sigs = {bb for i, b in enumerate(sig_blocks) if i in blocks for bb in b}

print("copying sigs to temp-folder")
for s in sigs:
    shutil.copy(s, os.environ['SNIC_TMP'] + "/")

sig_blocks = [[bb.replace("sigs", os.environ['SNIC_TMP']) for bb in b]for b in sig_blocks]

def run_bloc(i):
    if os.path.exists("blocks/block_{}.csv".format(i)):
        return False
    sub_sigs_1 = {v.split("/")[-1][:-4] : load_sig(v) for v in sig_blocks[block_order[i][0]]}
    sub_sigs_2 = {v.split("/")[-1][:-4] : load_sig(v) for v in sig_blocks[block_order[i][1]]}
    dists = {(k,l) : v.similarity(w, ignore_abundance=True) for k,v in sub_sigs_1.items() for l,w in sub_sigs_2.items()}
    dists = {k : v for k,v in dists.items() if v >0.05 and  k[0] != k[1]}
    print("Done bloc:", i)
    with open("blocks/block_{}_{}.csv".format(run_id,i), "w") as handle:
        handle.writelines(["{}\t{}\t{}\n".format(k[0],k[1],v) for k,v in dists.items()] )
    return True


pool = Pool(processes=24)

p = pool.map(run_bloc, list(range(len(block_order))))
