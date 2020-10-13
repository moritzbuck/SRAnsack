import json

with open("params.json") as handle:
    params = json.load(handle)


rule SRAass:
    output : assembly = "{root}/data/libraries/{sra_id}/{sra_id}.fna.gz",
             bins_md = "{root}/data/libraries/{sra_id}/{sra_id}.bins.json.gz",
             signatures = "{root}/data/libraries/{sra_id}/{sra_id}.sig.gz",
             qc = "{root}/data/libraries/{sra_id}/{sra_id}.fastp.json.gz",
             log = "{root}/data/libraries/{sra_id}/{sra_id}.log.gz"
    params : script = "workflow/scripts/SRAass.py", scratch = params['scratch'], len_cutoff= params['ass_len_cutoff'], max_redundance = params['max_redundance'], min_completeness = params['min_completeness'], rarefaction = params['rarefaction_read_collection']
    threads : 24
    conda : "../envs/SRAass.yaml"
    shell : """
        python {params.script} {wildcards.sra_id} {params.scratch} {wildcards.root} {threads} {params.len_cutoff} {params.max_redundance} {params.min_completeness} {params.rarefaction}
        """
