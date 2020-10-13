import json
import pandas
from os.path import join as pjoin
from subprocess import call
import os
from tqdm import tqdm
import tempfile

root = "/home/moritz/data/SRAprov/data"
threash = 0.001
nb_matches = 10000

sigs = []
for v in tqdm(os.walk(root)):
    for vv in v[2]:
        if vv.endswith(".sig.gz"):
            sigs += [pjoin(root,v[0], vv)]

with open(pjoin(root, "dbs", "all_sigs.txt"), "w") as handle:
    handle.writelines([l + "\n" for l in sigs])

index_cmd = "sourmash index   -k31  {root}/dbs/SRAaquas.sbt.json  --from-file {root}/dbs/all_sigs.txt {root}/library/SRR10176937/SRR10176937.sig.gz".format(root=root)
cmd = "sourmash search -q --threshold {threash} -o {tempfile} --num-results {nb_res}  -k31     {sig} {root}/dbs/SRAaquas.sbt.json"

call(index_cmd, shell=True)


def compare2sigs(sig):
    temp_out = tempfile.NamedTemporaryFile()
    idd = sig.split("/")[-1][:-7]

    imp_cmd = cmd.format(tempfile = temp_out.name, sig = sig, root = root, threash = threash, nb_res= nb_matches)
    call(imp_cmd, shell=True)

    vals = dict()
    with open(temp_out.name) as handle:
        head = handle.readline()[:-1]
        for l in handle.readlines():
            stuff = l[:-1].split(",")
            vals[(idd, stuff[1])] = float(stuff[0])
    return vals



compares = dict()
for q in tqdm(sigs):
    comapres.update(compare2sigs(q))
