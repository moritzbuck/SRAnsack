include : "workflow/rules/SRAsingle.smk"

from os.path import join as pjoin
import json

root = "/home/moritz/data/SRAprov/"
data = pjoin(root, "data")

with open(pjoin(data, "dbs", "sra_data.json")) as handle:
    sra_md = json.load(handle)

black_list = {'SRR2657549', 'SRR12682498', 'SRR12681850', 'SRR9617869', 'SRR10168434'}

#sras_to_do = [v['SRA_ID'] for v in sra_md.values() if float(v['nb_bases']) > 10000000]

sras_to_do = [v['SRA_ID'] for v in sra_md.values() if float(v['nb_bases']) > 100000000 and v['SRA_ID'] not in black_list]
#sras_to_do = black_list
print("#We have ", len(sras_to_do))

rule all :
    input : [pjoin(data, "libraries", g , g + ".sig.gz") for g in sras_to_do]
