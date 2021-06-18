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

def title2log(title, llen = 90) :
    text_insert = "{title} started at : {time}".format(title = title, time = datetime.now())
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(pjoin(temp_folder, coass_id + ".log"), "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    print(text, file = stderr, flush = True)

def into_line(text, llen = 90) :
    text_insert =  text
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(pjoin(temp_folder, coass_id + ".log"), "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    print(text, file = stderr, flush = True)


script , coass_id, params, final_location, threads = sys.argv

with open(params) as handle:
    params = json.load(handle)

len_cutoff = int(params['coass_len_cutoff'])
max_redundance = float(params['coass_max_redundance'])
min_completeness = float(params['coass_min_completeness'])
min_size = float(params['coass_min_bin_size'])
retries = int(params['retries'])
buffer_folder =  params['buffer']
with open(params['coass_file']) as handle:
    coass_dat = json.load(handle)[coass_id]

coass = coass_dat['samples']

scratch = os.environ.get('SNIC_TMP', os.environ.get('SCRATCH_DIR', params['scratch']))
temp_folder = pjoin(scratch, coass_id)
os.makedirs(temp_folder, exist_ok = True)

attribs = { 'temp_folder' : temp_folder,
            'coass_id'    : coass_id,
            'threads'     : threads,
            'len_cutoff' : len_cutoff,
            '' :
}

dl_lib = "parallel-fastq-dump --tmpdir {temp_folder}  --threads {threads} -s {sraid} --split-e --skip-technical --outdir {temp_folder}  >> {temp_folder}/{coass_id}.log  2>&1"

single_lib_prep = """
fastp -h /dev/null -j /dev/null  {fastp_inlibs} {fastp_outlibs} -w {threads}  >> {temp_folder}/{coass_id}.log 2>&1
rm {raw_libs}
pigz {qc_libs}

bbnorm.sh {norm_libs_in} {norm_libs_out} t={threads} pigz=t  2>> {temp_folder}/{coass_id}.log
cat {temp_folder}/normed_{lib}{paired}.fastq.gz >> {temp_folder}/{coass_id}{paired}.fastq.gz
rm {temp_folder}/normed_{lib}{paired}.fastq.gz
if [ -f {temp_folder}/normed_{lib}_2.fastq.gz ]
then
    cat {temp_folder}/normed_{lib}_2.fastq.gz >> {temp_folder}/{coass_id}_2.fastq.gz
    touch {temp_folder}/normed_{lib}.fastq.gz
    rm {temp_folder}/normed_{lib}_2.fastq.gz
fi
"""

title2log("Starting to get grab and process libs")
for i,lib in enumerate(coass):
    title2log("Processing {}/{}: {}".format((i+1),len(coass), lib))
    call(dl_lib.format(sraid = lib, **attribs), shell = True)

    read_libs = [pjoin(temp_folder,i) for i in os.listdir(temp_folder) if i.split("/")[-1].startswith(lib) and i.endswith(".fastq")]
    paired = True if len(read_libs)  > 1 else False

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

title2log("All libraries pre-processed")

#computing compression rate
gziped_qced_size = sum([os.path.getsize(pjoin(temp_folder,f)) for f in os.listdir(temp_folder) if f.startswith("QCed_")])
gziped_normed_size = sum([os.path.getsize(pjoin(temp_folder,f)) for f in os.listdir(temp_folder) if f.startswith(coass_id) and f.endswith(".fastq.gz")])

megahit_line = """
megahit -1 {temp_folder}/{coass_id}_1.fastq.gz -2 {temp_folder}/{coass_id}_2.fastq.gz -r {temp_folder}/{coass_id}.fastq.gz -t {threads} -o {temp_folder}/assembly --min-contig-len {len_cutoff} 2>> {temp_folder}/{coass_id}.log
rm {temp_folder}/{coass_id}_1.fastq.gz {temp_folder}/{coass_id}_2.fastq.gz {temp_folder}/{coass_id}.fastq.gz
"""

title2log("Assembling the s*** out of it")
call(megahit_line.format(**attribs), shell = True)

title2log("Making bowtie2 index")
call("bowtie2-build --threads {threads} {temp_folder}/assembly/final.contigs.fa {temp_folder}/index  >> {temp_folder}/{coass_id}.log 2>&1".format(**attribs), shell=True)

bowtie_lines = """
bowtie2 -p{threads} -x {temp_folder}/index {reads} -S {temp_folder}/mapping.sam 2>> {temp_folder}/{coass_id}.log
samtools view -b -S -@{threads}  {temp_folder}/mapping.sam >  {temp_folder}/mapping.bam 2>> {temp_folder}/{coass_id}.log
samtools sort -@ 24 -o {temp_folder}/mapping.sorted.bam {temp_folder}/mapping.bam 2>> {temp_folder}/{coass_id}.log
jgi_summarize_bam_contig_depths --outputDepth {temp_folder}/mapping.tsv --referenceFasta {temp_folder}/assembly.fna  {temp_folder}/mapping.sorted.bam >> {temp_folder}/{coass_id}.log 2>&1
"""

title2log("Starting mapping libs to ass")
contig_lens = {}
coverages = {}
var = {}

for i,lib in enumerate(coass):
    title2log("Mapping {}/{}: {}".format((i+1),len(coass), lib))

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
                contig_lens[ll[0]] = int(ll[1])
            coverages[lib][ll[0]] = float(ll[3])
            var[lib][ll[0]] = float(ll[3])


title2log("Mapping")
call(bowtie_lines.format(threads = threads, temp=temp_folder, sraid = SRA_ID, reads = reads), shell = True)

binning_line = """
metabat2 -i {temp_folder}/assembly/final.contigs.fna -o {temp_folder}/bins/{coass_id} -a {temp_folder}/mapping.tsv -s {min_size} -t {threads}   >> {temp_folder}/{coass_id}.log 2>&1
checkm taxonomy_wf life Prokaryote -x fa -t {threads} {temp_folder}/bins/ {temp_folder}/checkm > {temp_folder}/checkm.txt  2>> {temp_folder}/{coass_id}.log
"""

title2log("Binning")
call(binning_line.format(**attribs), shell = True)

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


reads2 = "-in1=" + clean_libs[0] + ((" -in2=" + clean_libs[1]) if paired else "")
sub_libs = reads2.replace(".clean", ".subclean").replace("-in", "-out")

with open(pjoin(temp_folder, SRA_ID + ".bins.json"), "w") as handle:
    json.dump(chekm_out, handle, indent=4, sort_keys=True)

read_proc_line = """
reformat.sh {in_reads} {out_reads} samplereadstarget={subcount} sampleseed=42 t={threads}  2>> {temp}/{sraid}.log
sourmash compute --track-abundance --merge {sraid} -k31 --scaled 1000 {reads} -o {temp}/{sraid}.sig -p{threads}  2>> {temp}/{sraid}.log
"""

title2log("Read subsetting and sketching")

reads3 = " ".join(clean_libs)
call(read_proc_line.format(sraid = SRA_ID, temp=temp_folder, threads = threads, reads = reads3, in_reads = reads2, out_reads=sub_libs, subcount = rarefaction), shell=True)

title2log("Moving to final location and cleaning-up")
into_line("Finished with {} bins".format(len(chekm_out)))

library_loc = pjoin(final_location, "data","libraries", SRA_ID)
os.makedirs(library_loc, exist_ok=True)
os.makedirs(pjoin(library_loc, "bins"), exist_ok=True)

shutil.move(pjoin(temp_folder, "assembly.fna"), pjoin(library_loc, SRA_ID + ".fna"))
shutil.move(pjoin(temp_folder, SRA_ID + ".bins.json"), library_loc)
shutil.move(pjoin(temp_folder, SRA_ID + ".fastp.json"), library_loc)
for lib in clean_libs:
    shutil.move(lib.replace(".clean", ".subclean"), library_loc)
shutil.move(pjoin(temp_folder, SRA_ID + ".log"), library_loc)
shutil.move(pjoin(temp_folder, SRA_ID + ".sig"), library_loc)

for b in chekm_out:
    shutil.copy(pjoin(temp_folder, "bins", b + ".fa"), pjoin(library_loc, "bins" ))

for path, subdirs, files in os.walk(library_loc):
    for name in files:
        call("pigz " + pjoin(path, name), shell = True)

shutil.rmtree(temp_folder)
