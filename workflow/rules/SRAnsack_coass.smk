import json

rule SRAnsack_coass:
    output : assembly = "{root}/data/coassemblies/{coass_id}/{coass_id}.fna.gz",
             bins_md = "{root}/data/coassemblies/{coass_id}/{coass_id}.bins.json.gz",
             log = "{root}/data/coassemblies/{coass_id}/{coass_id}.log.gz"
    params : script = "workflow/scripts/SRAansack_coass.py", params = "params.json"
    threads : 24
    conda : "../envs/SRAnsack_ass.yaml"
    shell : """
        python {params.script} {wildcards.coass_id} {params.params} {wildcards.root} {threads}
        """
