---
title: "SRAnsack : metagenome signatures and assembled genomes for the exploration of SRA's aquatic data"
author: "Moritz Buck^1,\\*^Malhiheh?, Stefan?"
date: "^1^Swedish University of Agricultural Sciences, ^\\*^corresponding author"
bibliography: biblio.bib
---


The SRA (short read archive) is full of underused metagenomic data that is hard to put in context due to lacking metadata. Even for well characterized metagenomes, only selected assemblies and metagenome assembled genomes (MAGs) are made public. This makes it extremely difficult to put your own metagenomes in the right context of genomic and metagenomic diversity.

To remedy this, we have developed a simple snakemake-based pipeline that processes all aquatic metagenome[^1] from SRA independently, equally and reasonably efficiently, with the possibility to expand to other biomes at a later point. The pipeline produces a set of data-collection: a collection of assemblies and MAGs, a collection of min-hash signatures and a collection of rarefied libraries. which are all made publicly available, and easily accessible for the public's use.

Some of this data is also further processed to some degree to make it more usable and allow for easier comparison with your own data. The MAG collection is clustered into metagenomic Operation Taxonomic Units (mOTUS), and taxonomically annotated. The signature collection is used to make an 'all vs all' comparison of all the input SRA-libraries to compute similarity values, which allows for fast identification of highly similar samples. Additionally the rarefied libraries are screened for rRNA genes as many amplicon data-sets are wrongly annotated as metagenome in SRA. We processed a total of 20.653 libraries, obtaining ~1300Gb of assemblies, containing 144.222 MAGs which cluster in ~27.000 metagenomic Operation Taxonomic Units (mOTUs). This is similar in approach to Joint Genome Institute's genomic catalog of Earth’s microbiomes[@nayfach_2021], on a somewhat different scope, and a less restrictive data-set. Which itself is predated by efforts such as [@parks_2017]. Other more "read"-focused approaches are also around, such as [@mitchell_2020].

This is but the first step of this growing data-set. The similarities in samples will be used in the future to compute a large collection of co-assemblies to expand the MAG collection to the more rare fraction of these samples and additionally the rarefied libraries will be mapped to representative MAGs of the MAG collection to identify global distributions. This short paper is here to make this first iteration of the data available as the authors think it is valuable for it to reach the community. Additionally some controversy is currently rife about the 'misuse' of some the data submitted to SRA. I hereby present the community with my project and some preliminary results, out of which any actor that is unhappy to have his public data involved is invited to contact the corresponding author, who will remove the data-set from any future peer-reviewed publication, and down-stream analysis, and remark in the supplemental data on the removal of the data.

![](/home/moritz/temp/tSNE.pdf)
*tSNE transformation of the similarity matrix of the library MinHashes*

## Methods
The libraries of interest are identified with the `Entrez`-module[@bio_entrez] of the `biopython` python-package[@cock_2009], and then downloaded with `parallel-fastq-dump`[@valieris_2021] a conveniently parallelized wrapper for `fastq-dump`[@sratools]. The reads were quality filtered and trimmed with `fastp`[@chen_2018](version 0.20.1) due to it's agnosticism to adapters, and speed.
`MinHash`-signatures were computed for all libraries with `sourmash`[@brown_2016](version 3.5.0 and with the `--track-abundance -k31 --scaled 1000` options) and all of the obtained signatures were compaired pairwise to obtain similarities. The libraries were rarefied to 100.000 read-pairs (or reads if library is unpaired) with the `reformat.sh` program from the `bbtools` software suite[@bushnell_2014](version 38.18), to provide a handy exploration dataset of reads.  The rRNA content of all libraries was obtained from the rarefied libraries with `SortMeRNA`[@kopylova_2012] (version 4.2.0, in future version it will be obtained from whole library, it was an after though in this case).
The full QCed libraries are assembled with megahit [@li_2015](version 1.2.9), only the contigs larger than 2.5kb were kept to reduce size and increase average quality of the assembly data. Each library is mapped to its assembly with bowtie2[@langmead_2012](version 2.4.1), and the mapping are post-processed with samtools[@li_2009](version 1.10).
The mapping combined with the assembly is then used to bin the contigs with metabat2[@kang_2019](version 2.15), the obtained bins are then quality checked with CheckM[@parks_2015](version 1.1.3), and filtered based on a 30\% completeness and 5\% redundancy threashold. These are pretty lax filters based on [@olm_2017] and personal experience, but are selected for capturing a wide rather then particularly accurate biodiversity. All bins were taxonomically annotated with `GTDBtk` [@chaumeil_2020](version 1.4.0, with release 95 of the database), and clustered into metagenomic Operational Taxonomic Units (mOTUs) with `mOTUlizer`[@buck_2021](version 0.2.2), to speed up the process, only bins classified within the same family were clustered together.
This whole workflow is intended to be wrapped in a single `snakemake`-pipeline[@koster_2012], as of now it is not entirely so, only the processing of the single libraries is, and the rest is a number of separate scripts still to be properly integrated, nevertheless, all is available at `github.com/moritzbuck/SRAnsack`.


## Usage cases

This data can be used for a number of purposes. You can look for mOTUs of your favorite microorganism and back track to samples rich in these, or use the read library against your favorite genome to find samples where they are highly abundant.
Conversely you can compute a MinHash-signature of your freshly sequenced metagenome and compare it to the available signatures to find related samples even if the contextual data is missing or misleading. An other similar use would be to find your favorite metagenome from the SRA and check the library similarity file for related samples. Finally, what we intend to do it, first we will use the similarities between libraries to partition the data-set to generate coassemblies, and bins from these which will probably multiply around ~5 fold the number of MAGs, and we will use a combination of metadata and similarity between samples to generate a number of reference sets for certain biomes, so as to give a better genomic context for any metagenomic study done in them.

## Available data

* MAG collection: ~100-zipped Gb
* MAG metadata:
* MAG ANIs: ~1 Gb
* Assembly collection: ~400-zipped Gb
* Signature collection: ~300-zipped Gb
* Rarefied read collection: ~400-zipped Gb
* Library metadata:
* Library similarities:

[^1]: We call aquatic metagenomes any metagenomic dataset with their ncbi taxonomy as any of : peat metagenome, seawater metagenome, marine metagenome, marine plankton metagenome, marine sediment metagenome, karst metagenome, lagoon metagenome, lake water metagenome, glacier lake metagenome, freshwater metagenome, freshwater sediment metagenome, aquatic metagenome, bog metagenome, drinking water metagenome, estuary metagenome, Winogradsky column metagenome,  alkali sediment metagenome,  sediment metagenome or water metagenome.


## Bibliography
