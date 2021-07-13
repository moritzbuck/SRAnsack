import sys, os
from tqdm import tqdm
from subprocess import call
from Bio import SeqIO
from os.path import join as pjoin
import json
import re
from datetime import datetime
from sys import stderr
import shutil
from math import ceil, floor
import pandas


def title2log(title, llen = 90) :
    text_insert = "{title} started at : {time}".format(title = title, time = datetime.now())
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(pjoin(log_file), "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    print(text, file = stderr, flush = True)

def into_line(text, llen = 90) :
    text_insert =  text
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(pjoin(log_file), "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    print(text, file = stderr, flush = True)

def df_info():
    into_line("available space on {temp_folder} : {space}".format(space = shutil.disk_usage(temp_folder).free/1000000000, **attribs))


script , coass_id, params, final_location, threads = sys.argv

with open(params) as handle:
    params = json.load(handle)

len_cutoff = int(params['coass_len_cutoff'])
max_redundance = float(params['coass_max_redundance'])
min_completeness = float(params['coass_min_completeness'])
min_size = float(params['coass_min_bin_size'])
retries = int(params['retries'])
buffer_folder =  params['buffer']
max_buffer_size = params['seqio_buffer_size']

with open(params['coass_file']) as handle:
    coass_dat = json.load(handle)[coass_id]

stats = {}
stats.update(coass_dat)
coass = coass_dat['samples']

scratch = os.environ.get('SNIC_TMP', os.environ.get('SCRATCH_DIR', params['scratch']))
temp_folder = pjoin(scratch, coass_id)
os.makedirs(temp_folder, exist_ok = True)
os.makedirs(final_location, exist_ok = True)

library_loc = pjoin(final_location, "data","coassemblies", coass_id)
log_file = pjoin(final_location, "data","coassemblies", coass_id, coass_id + ".log")
attribs = { 'temp_folder' : temp_folder,
            'coass_id'    : coass_id,
            'threads'     : threads,
            'len_cutoff' : len_cutoff,
            'min_size' : int(min_size),
            'log_file' : log_file
}

dl_lib = "parallel-fastq-dump --tmpdir {temp_folder}  --threads {threads} -s {sraid} --split-e --skip-technical --outdir {temp_folder}  >> {log_file}  2>&1"

single_lib_prep = """
fastp -h /dev/null -j /dev/null  {fastp_inlibs} {fastp_outlibs} -w {threads}  >> {log_file} 2>&1
rm {raw_libs}
pigz {qc_libs}

bbnorm.sh {norm_libs_in} {norm_libs_out} t={threads} pigz=t  2>> {log_file}
cat {temp_folder}/normed_{lib}{paired}.fastq.gz >> {temp_folder}/{coass_id}{paired}.fastq.gz
rm {temp_folder}/normed_{lib}{paired}.fastq.gz
if [ -f {temp_folder}/normed_{lib}_2.fastq.gz ]
then
    cat {temp_folder}/normed_{lib}_2.fastq.gz >> {temp_folder}/{coass_id}_2.fastq.gz
    touch {temp_folder}/{coass_id}.fastq.gz
    rm {temp_folder}/normed_{lib}_2.fastq.gz
fi
"""

title2log("Starting to get grab and process libs")
df_info()
for i,lib in []:#enumerate(coass):
    title2log("Processing {}/{}: {}".format((i+1),len(coass), lib))
    df_info()
    call(dl_lib.format(sraid = lib, **attribs), shell = True)

    read_libs = [pjoin(temp_folder,i) for i in os.listdir(temp_folder) if i.split("/")[-1].startswith(lib) and i.endswith(".fastq")]
    paired = True if len(read_libs)  > 1 else False
    paired = True
    inlibs = "{flag1}{fwd} " + ("{flag2}{rev}" if paired else "")

    single_lib_dat = {
    'fastp_inlibs' : inlibs.format(flag1 = "--in1 ", flag2 = "--in2 ", fwd = pjoin(temp_folder,lib + "_1.fastq"), rev = pjoin(temp_folder,lib + "_2.fastq")),
    'fastp_outlibs' : inlibs.format(flag1 = "--out1 ", flag2 = "--out2 ", fwd = pjoin(temp_folder,"QCed_" + lib + "_1.fastq"), rev = pjoin(temp_folder,"QCed_" + lib + "_2.fastq")),
    'buffer' :buffer_folder,
    'qc_libs' : inlibs.format(flag1 = "", flag2 = "", fwd = pjoin(temp_folder,"QCed_" + lib + "_1.fastq"), rev = pjoin(temp_folder,"QCed_" + lib + "_2.fastq")),
    'raw_libs' : inlibs.format(flag1 = "", flag2 = "", fwd = pjoin(temp_folder,lib + "_1.fastq"), rev = pjoin(temp_folder,lib + "_2.fastq")),
    'qc_libs_gz' : inlibs.format(flag1 = "", flag2 = "", fwd = pjoin(temp_folder,"QCed_" + lib + "_1.fastq.gz"), rev = pjoin(temp_folder,"QCed_" + lib + "_2.fastq.gz")),
    'norm_libs_in' : inlibs.format(flag1 = "in=", flag2 = "in2=", fwd = pjoin(temp_folder,"QCed_" + lib + "_1.fastq.gz"), rev = pjoin(temp_folder,"QCed_" + lib + "_2.fastq.gz")),
    'norm_libs_out' : inlibs.format(flag1 = "out=", flag2 = "out2=", fwd = pjoin(temp_folder,"normed_" + lib + "_1.fastq.gz"), rev = pjoin(temp_folder,"normed_" + lib + "_2.fastq.gz")),
    'paired' : "_1" if paired else "",
    'lib' : lib
    }
    call(single_lib_prep.format(**single_lib_dat, **attribs), shell = True)


into_line("All libraries pre-processed")
df_info()

#computing compression rate
gziped_qced_size = sum([os.path.getsize(pjoin(temp_folder,f)) for f in os.listdir(temp_folder) if f.startswith("QCed_")])
gziped_normed_size = sum([os.path.getsize(pjoin(temp_folder,f)) for f in os.listdir(temp_folder) if f.startswith(coass_id) and f.endswith(".fastq.gz")])

stats['QC_gziped_size'] = gziped_qced_size
stats['QC_sample_normed_size'] = gziped_normed_size


megahit_line = """
if [ -s {temp_folder}/{coass_id}_1.fastq.gz ]
then
    bbnorm.sh in={temp_folder}/{coass_id}_1.fastq.gz in2={temp_folder}/{coass_id}_2.fastq.gz out={temp_folder}/pass2_{coass_id}_1.fastq.gz  out2={temp_folder}/pass2_{coass_id}_2.fastq.gz t=24 pigz=t 2>> {log_file}
fi
if [ -s {temp_folder}/{coass_id}.fastq.gz ]
then
    bbnorm.sh in={temp_folder}/{coass_id}.fastq.gz out={temp_folder}/pass2_{coass_id}.fastq.gz t=24 pigz=t 2>> {log_file}
fi
rm {temp_folder}/{coass_id}_1.fastq.gz {temp_folder}/{coass_id}_2.fastq.gz {temp_folder}/{coass_id}.fastq.gz 2>> {log_file}

du `ls | grep pass2` > {temp_folder}/du_double_normed.txt

if [ ! -s {temp_folder}/pass2_{coass_id}_1.fastq.gz ]
then
    megahit -r {temp_folder}/pass2_{coass_id}.fastq.gz -t {threads} -o {temp_folder}/assembly --min-contig-len {len_cutoff} 2>> {log_file}
else
    if [ ! -s {temp_folder}/pass2_{coass_id}.fastq.gz ]
    then
        megahit -1 {temp_folder}/pass2_{coass_id}_1.fastq.gz -2 {temp_folder}/pass2_{coass_id}_2.fastq.gz  -t {threads} -o {temp_folder}/assembly --min-contig-len {len_cutoff} 2>> {log_file}
    else
        megahit -1 {temp_folder}/pass2_{coass_id}_1.fastq.gz -2 {temp_folder}/pass2_{coass_id}_2.fastq.gz -r {temp_folder}/pass2_{coass_id}.fastq.gz -t {threads} -o {temp_folder}/assembly --min-contig-len {len_cutoff} 2>> {log_file}
    fi
fi
rm {temp_folder}/pass2_{coass_id}_1.fastq.gz {temp_folder}/pass2_{coass_id}_2.fastq.gz {temp_folder}/pass2_{coass_id}.fastq.gz
"""



title2log("Assembling the s*** out of it")

call(megahit_line.format(**attribs), shell = True)
with open(pjoin(temp_folder, "du_double_normed.txt")) as handle:
    stats['QC_conormed_size'] = sum([int(l.split()[0]) for l in handle])


title2log("Making bowtie2 index")
df_info()

call("bowtie2-build --threads {threads} {temp_folder}/assembly/final.contigs.fa {temp_folder}/index  >> {log_file} 2>&1".format(**attribs), shell=True)

bowtie_lines = """
bowtie2 -p{threads} -x {temp_folder}/index {reads} -S {temp_folder}/mapping.sam 2>> {log_file}
samtools view -b -S -@{threads}  {temp_folder}/mapping.sam >  {temp_folder}/mapping.bam 2>> {log_file}
samtools sort -@ 24 -o {temp_folder}/mapping.sorted.bam {temp_folder}/mapping.bam 2>> {log_file}
jgi_summarize_bam_contig_depths --outputDepth {temp_folder}/mapping.tsv --referenceFasta {temp_folder}/assembly/final.contigs.fa  {temp_folder}/mapping.sorted.bam >> {log_file} 2>&1
rm {temp_folder}/QCed_{lib}_1.fastq.gz {temp_folder}/QCed_{lib}_2.fastq.gz {temp_folder}/QCed_{lib}.fastq.gz
"""

title2log("Starting mapping libs to ass")
contig_lens = {}
coverages = {}
var = {}

for i,lib in enumerate(coass):
    title2log("Mapping {}/{}: {}".format((i+1),len(coass), lib))
    df_info()

    read_libs = [pjoin(temp_folder,i) for i in os.listdir(temp_folder) if i.endswith(".fastq.gz") and i.startswith("QCed_" + lib)]
    paired = True if len(read_libs)  > 1 else False
    if paired:
        inlibs = "-1 {temp_folder}/QCed_{lib}_1.fastq.gz -2 {temp_folder}/QCed_{lib}_2.fastq.gz".format(temp_folder = temp_folder, lib = lib)
    else :
        inlibs = "-1 {temp_folder}/QCed_{lib}.fastq.gz".format(temp_folder = temp_folder, lib = lib)
    single_lib_dat = {
    'lib' : lib,
    'reads' : inlibs
    }
    call(bowtie_lines.format(**single_lib_dat, **attribs), shell = True)

    with open(pjoin(temp_folder, "mapping.tsv")) as handle:
        handle.readline()
        coverages[lib] = {}
        var[lib] = {}
        for l in handle:
            ll = l.strip().split()
            if ll[0] not in contig_lens:
                contig_lens[ll[0]] = float(ll[1])
            coverages[lib][ll[0]] = float(ll[3])
            var[lib][ll[0]] = float(ll[4])
    for f in ['mapping.sam', 'mapping.bam', 'mapping.sorted.bam', 'mapping.tsv']:
        os.remove(pjoin(temp_folder, f))


title2log("Done mapping libs to ass")
into_line("regenerating mapping table and cleanup")
df_info()

for f in os.listdir(temp_folder):
    if f.startswith("index"):
        os.remove(pjoin(temp_folder, f))


head = ['contigName', 'contigLen', 'totalAvgDepth'] + [c + k for c in coass for k in ['', '-var']]
var = {k + "-var" : v for k,v in var.items()}
totalAvgDepth = {vv : 0 for v in coverages.values() for vv in v}
for samp,v in coverages.items():
    for contig, cov in v.items():
        totalAvgDepth[contig] += cov
clean = { 'contigLen' : contig_lens,  'totalAvgDepth' : totalAvgDepth}
clean.update(coverages)
clean.update(var)

stats['full_assembly_size'] = sum([v for v in contig_lens.values()])

mapping_table = pandas.DataFrame.from_dict(clean)[head[1:]]
mapping_table.to_csv(pjoin(temp_folder, "temp_mapping.tsv"), index_label = 'contigName', sep = "\t")

title2log("Making them bins")

binning_line = """
metabat2 --saveCls --noBinOut -i {temp_folder}/assembly/final.contigs.fa -o {temp_folder}/{coass_id}_bins.tsv -a {temp_folder}/temp_mapping.tsv -s {min_size} -t {threads}   >> {log_file} 2>&1
"""
call(binning_line.format(**attribs), shell = True)

title2log("Processing the bins, and renaming contigs accordingly")
df_info()

with open(pjoin(temp_folder, coass_id + "_bins.tsv")) as handle:
    ctg2bin = {l.split()[0] : l.strip().split()[1] for l in handle}

zeros = lambda ll : len(str(ll))

bins = set(ctg2bin.values())
bins.remove('0')
zeros_bin = zeros(len(bins))
bins2bin_id = {b : coass_id + "-bin_" + str(i+1).zfill(zeros_bin) for i,b in enumerate(bins)}
bins2bin_id['0'] = coass_id + "-unbinned"

ctg2bin = {k : bins2bin_id[v] for k,v in ctg2bin.items()}

count_ctgs = {b : 0 for b in bins2bin_id.values()}
for v in ctg2bin.values():
    count_ctgs[v] += 1

enum_ctgs = {b : 1 for b in bins2bin_id.values()}
ctg2new_ctg = {}

os.makedirs(pjoin(temp_folder, "bins"), exist_ok = True)

bin_handles = {b : open(pjoin(temp_folder, "bins", b +".fa"), "w") for b in bins2bin_id.values()}
buffer = []

binned_size = 0
with open(pjoin(temp_folder, coass_id + ".fna"), "w") as handle:
    for i,s in tqdm(enumerate(SeqIO.parse(pjoin(temp_folder, "assembly", "final.contigs.fa"), "fasta"))):
        bbin = ctg2bin[s.id]
        if bbin != coass_id + "-unbinned" :
            binned_size += len(s)
        ctg2new_ctg[s.id] = bbin + "-ctg_" + str(enum_ctgs[bbin]).zfill(zeros(count_ctgs[v]))
        s.id = ctg2new_ctg[s.id]
        enum_ctgs[bbin] += 1
        s.description = ""
        buffer += [s]
        SeqIO.write(s, bin_handles[bbin], "fasta")
        if len(buffer) > max_buffer_size:
            SeqIO.write(buffer, handle, "fasta")
            buffer = []
    SeqIO.write(buffer, handle, "fasta")

stats['binned_assembly_size'] = binned_size

for h in bin_handles.values():
    h.close()


mapping_table.index = [ctg2new_ctg[i] for i in mapping_table.index]
mapping_table.to_csv(pjoin(temp_folder, coass_id + "_mapping.tsv"), index_label = 'contigName', sep = "\t")


title2log("Running checkm")
df_info()

checkm_line = """
checkm taxonomy_wf life Prokaryote -x fa -t {threads} {temp_folder}/bins/ {temp_folder}/checkm > {temp_folder}/checkm.txt  2>> {log_file}
"""

call(checkm_line.format(**attribs), shell = True)


with open(pjoin(temp_folder, "checkm.txt")) as handle:
    all_lines = [l.strip() for l in  handle.readlines() if " INFO:" not in l]


if "ERROR" not in all_lines[0]:
    all_lines = [re.sub(r"  +","\t", a).split("\t") for a in all_lines]
    header_lines = [i for i,l in enumerate(all_lines) if 'Bin Id' in l and 'Completeness' in l and 'Contamination' in l]
    header_lines = header_lines[0]
    header_line = all_lines[header_lines]
    lines = [l for i,l in enumerate(all_lines) if i != header_lines and len(l) == len(header_line)]
    lines = [{a : b if a in ['Marker lineage', 'Bin Id'] else float(b) for a,b in zip(header_line,l) }for l in lines]

    chekm_out = {l['Bin Id'] : {k: l[k] for k in ('Completeness', 'Contamination')} for l in lines if l['Completeness'] > min_completeness and l['Contamination'] < max_redundance}
else :
    chekm_out = {}


with open(pjoin(temp_folder, coass_id + ".bins.json"), "w") as handle:
    json.dump(chekm_out, handle, indent=4, sort_keys=True)

with open(pjoin(temp_folder, coass_id + ".stats.json"), "w") as handle:
    json.dump(stats, handle, indent=4, sort_keys=True)


title2log("Moving to final location and cleaning-up")
into_line("Finished with {} bins".format(len(chekm_out)))
df_info()

os.makedirs(library_loc, exist_ok=True)
os.makedirs(pjoin(library_loc, "bins"), exist_ok=True)

shutil.move(pjoin(temp_folder, coass_id + ".fna"), library_loc)
shutil.move(pjoin(temp_folder, coass_id + ".bins.json"), library_loc)
shutil.move(pjoin(temp_folder, coass_id + ".stats.json"), library_loc)
shutil.move(pjoin(temp_folder, coass_id + "_mapping.tsv"), library_loc)

for b in chekm_out:
    shutil.copy(pjoin(temp_folder, "bins", b + ".fa"), pjoin(library_loc, "bins" ))

for path, subdirs, files in os.walk(library_loc):
    for name in files:
        call("pigz " + pjoin(path, name), shell = True)


shutil.rmtree(temp_folder)
df_info()
