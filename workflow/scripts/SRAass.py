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
    with open(pjoin(temp_folder, SRA_ID + ".log"), "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    print(text, file = stderr, flush = True)

def into_line(text, llen = 90) :
    text_insert =  text
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(pjoin(temp_folder, SRA_ID + ".log"), "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    print(text, file = stderr, flush = True)


script , SRA_ID, temp_folder, final_location, threads, len_cutoff, max_redundance, min_completeness , rarefaction = sys.argv


len_cutoff = int(len_cutoff)
max_redundance = float(max_redundance)
min_completeness = float(min_completeness)

sratools_line = "parallel-fastq-dump --threads {threads} -s {sraid} --split-e --skip-technical --outdir {temp}  >> {temp}/{sraid}.log  2>&1"

temp_folder = pjoin(temp_folder, SRA_ID)
os.makedirs(temp_folder, exist_ok=True)

title2log("Starting SRA " + SRA_ID)

title2log("Downloading reads from SRA")
call(sratools_line.format(sraid = SRA_ID, temp = temp_folder, threads = threads), shell = True)

read_libs = [pjoin(temp_folder,i) for i in os.listdir(temp_folder) if i.split("/")[-1].startswith(SRA_ID) and i.endswith(".fastq")]
paired = True if len(read_libs)  > 1 else False

if paired :
    read_libs = sorted([l for l in read_libs if "_1.fastq" in l or "_2.fastq" in l])

fastp_line = "fastp -h /dev/null -j {temp}/{sraid}.fastp.json  --in1 {lib1} {pot_lib2} --out1 {out1} {pot_out2} -w {threads}  >> {temp}/{sraid}.log 2>&1"
clean_libs = [l.replace(".fastq", ".clean.fastq") for l in read_libs]

title2log("QCing reads")
call(fastp_line.format(lib1 = read_libs[0], out1 = clean_libs[0],
     pot_lib2 = "--in2 " + read_libs[1] if paired else "",
     pot_out2 = "--out2 " + clean_libs[1] if paired else "",
     threads = threads, sraid = SRA_ID, temp=temp_folder
), shell = True)

reads = ("-1 " + clean_libs[0] + " -2 " + clean_libs[1]) if paired else "-r " + clean_libs[0]
megahit_line = "megahit {reads} -t {threads} -o {temp}/assembly 2>> {temp}/{sraid}.log"

title2log("Assembling")
call(megahit_line.format(reads = reads, threads = threads, sraid = SRA_ID, temp = temp_folder), shell = True)

seqs = [s for s in SeqIO.parse(pjoin(temp_folder, "assembly", "final.contigs.fa"), "fasta") if len(s) > len_cutoff]

zeros = len(str(len(seqs)))

for i,s in enumerate(seqs):
    s.id = SRA_ID + "_" + str(i).zfill(zeros)
    s.description = ""

SeqIO.write(seqs, pjoin(temp_folder, "assembly.fna"), "fasta")


bowtie_lines = """
bowtie2-build --threads {threads} {temp}/assembly.fna {temp}/index  >> {temp}/{sraid}.log 2>&1
bowtie2 -p24 -x {temp}/index {reads} -S {temp}/mapping.sam 2>> {temp}/{sraid}.log
samtools view -b -S -@24  {temp}/mapping.sam >  {temp}/mapping.bam 2>> {temp}/{sraid}.log
samtools sort -@ 24 -o {temp}/mapping.sorted.bam {temp}/mapping.bam 2>> {temp}/{sraid}.log
"""

title2log("Mapping")
call(bowtie_lines.format(threads = threads, temp=temp_folder, sraid = SRA_ID, reads = reads), shell = True)

binning_line = """
jgi_summarize_bam_contig_depths --outputDepth {temp}/mapping.tsv --referenceFasta {temp}/assembly.fna  {temp}/mapping.sorted.bam >> {temp}/{sraid}.log 2>&1
metabat2 -i {temp}/assembly.fna -o {temp}/bins/{sraid} -a {temp}/mapping.tsv -s 500000 -t {threads}   >> {temp}/{sraid}.log 2>&1
checkm taxonomy_wf life Prokaryote -x fa -t 24 {temp}/bins/ {temp}/checkm > {temp}/checkm.txt  2>> {temp}/{sraid}.log
"""

title2log("Binning")
call(binning_line.format(threads = threads, temp=temp_folder, sraid = SRA_ID), shell = True)

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

library_loc = pjoin(final_location, "data","library", SRA_ID)
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
