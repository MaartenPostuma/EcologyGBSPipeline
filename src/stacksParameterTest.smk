configfile: "configParameter.yaml"
LARGEM = list(range(1,12))
rule all:
    input:
        MparameterPNG=expand("{dir}/LargeMparameter.png",dir=config["outdir"]),
        testLargeM=expand("{dir}/LargeM{largeM}/test{largeM}",largeM=LARGEM,dir=config["outdir"])

rule subset_popmap:
    input:
        config["popmap"]
    output:
        expand("{dir}/popmapSub.tsv",dir=config["outdir"])
    params:
        config["nInds"]
    shell:
        "paste <(shuf -n {params} {input} | cut -f1) <(yes opt | head -n {params}) > {output}"

rule mkdirLargeM:
    input:
        expand("{dir}/popmapSub.tsv",dir=config["outdir"])
    output:
        expand("{dir}/LargeM{{largeM}}/test{{largeM}}",dir=config["outdir"])
    params:
        largeM="{largeM}"
    shell:
        "cat {input} | head -n {params.largeM} > {output}"

rule runStacksLargeM:
    input:
        inDir=expand("{dir}/LargeM{{largeM}}/test{{largeM}}",dir=config["outdir"])
    output:
        outLog=expand("{dir}/LargeM{{largeM}}/denovo_map.log",dir=config["outdir"])
    params:
        largeM="{largeM}",
        popmapSub=expand("{dir}/popmapSub.tsv",dir=config["outdir"]),
        outDir=expand("{dir}/LargeM{{largeM}}/",dir=config["outdir"]),
        indir=config["indir"]
    conda:
        "src/env/stacks.yaml"
    threads: 4
    shell:
        """
        denovo_map.pl --samples {params.indir} --popmap {params.popmapSub} -T {threads} \
        -o {params.outDir} -n {params.largeM} -M {params.largeM}  \
        -X 'populations: -R 80'
        """

rule extractInfoLargeM:
    input:
        outLog=expand("{dir}/LargeM{largeM}/denovo_map.log",largeM=LARGEM,dir=config["outdir"])
    output:
        MparameterTSV=expand("{dir}/LargeMparameter.tsv",dir=config["outdir"])
    params:
        dir=config["outdir"]
    shell:
        """
        paste <(cat {params.dir}/LargeM*/denovo_map.log | grep "Kept" | cut -f2,14 -d " ") \
        <(cat {params.dir}/LargeM*/denovo_map.log | grep "per-sample coverage" | cut -f2 -d "=" | cut -f1 -d "x") \
        <(cat {params.dir}/LargeM*/denovo_map.log | grep "\-M" | grep "denovo_map.pl" | cut -f 13 -d " ") >\
        {params.dir}/LargeMparameter.tsv
        """

rule makePlotLargeM:
    input:
        MparameterTSV=expand("{dir}/LargeMparameter.tsv",dir=config["outdir"])
    output:
        MparameterPNG=expand("{dir}/LargeMparameter.png",dir=config["outdir"])
    params:
        dir=config["outdir"]
    shell:
        "Rscript src/LargeMparameterTest.R {params.dir}"
