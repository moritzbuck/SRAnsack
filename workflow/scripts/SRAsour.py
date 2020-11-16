import json
from os.path import join as pjoin
from subprocess import call
import os
from tqdm import tqdm
from sourmash import MinHash
from sourmash.signature import SourmashSignature
from sourmash.signature import load_signatures
from multiprocessing import Pool
from numba import jit


def load_sig(f):
    with open(f) as handle:
        ll = "".join( handle.readlines()).replace('\n','')
    return list(load_signatures(ll)) [0]

@jit
def test(i,j):
    return i+j

def load_minhash(f):
    with open(f) as handle:
        ll = handle.read()
    ll = [int(t) for t in ll.split('"mins":[')[1].split('],"md5sum"')[0].split(",")]
    return ll


root = "/home/moritz/data/SRAprov/data"
threash = 0.001
nb_matches = 10000

all_sigs = ["sigs/" + s for s in os.listdir('sigs')]

block_size = 100

sig_blocks = [all_sigs[i:(i+block_size)] for i in list(range(0,len(all_sigs), block_size))]


block_order = [(i,j) for i,b in enumerate(sig_blocks) for j,v in enumerate(sig_blocks)]


def run_bloc(i):
    sub_sigs_1 = {v.split("/")[-1][:-4] : load_sig(v) for v in sig_blocks[block_order[i][0]]}
    sub_sigs_2 = {v.split("/")[-1][:-4] : load_sig(v) for v in sig_blocks[block_order[i][1]]}
    dists = {(k,l) : v.similarity(w, ignore_abundance=True) for k,v in sub_sigs_1.items() for l,w in sub_sigs_2.items()}
    dists = {k : v for k,v in dists.items() if v >0.05 and  k[0] != k[1]}
    print("Done bloc:", i)
    return dists

@jit
def run_bloc_fast(i):
    b1 = {v.split("/")[-1][:-4] : load_minhash(v) for v in sig_blocks[block_order[i][0]]}
    b2 = {v.split("/")[-1][:-4] : load_minhash(v) for v in sig_blocks[block_order[i][1]]}

#    sub_sigs = {v : load_minhash("sigs/" + v + ".sig") for v in b1.union(b2)}

    dists = {(k,l) : len(v.intersection(w))/len(v.union(w))  for k,v in b1.items() for l,w in b2.items()}
    dists = {k : v for k,v in dists.items() if v >0.05 and  k[0] != k[1]}
    print("Done bloc:", i)
    return dists


def timer(fct, i):
    start = timeit.default_timer()
    ret = fct(i)
    stop = timeit.default_timer()
    print('Time: ', stop - start)
    return ret


pool = Pool(processes=6)

p = pool.map(run_bloc, list(range(len(block_order))))


def test1(ll):
