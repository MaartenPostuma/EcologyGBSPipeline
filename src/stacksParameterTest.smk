param_nInds=config["param_StacksTest"]["nInds"]

LARGEM = list(range(1,12))
#rule all:
#    input:
#        MparameterPNG=expand("{dir}/stacksTestparameter.png",dir=config["outputDir"]),
#        testLargeM=expand("{dir}/stacksTest{largeM}/test{largeM}",largeM=LARGEM,dir=config["outputDir"])

rule subset_popmap:
    input:
        popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"]),
    output:
        expand("{dir}/popmapSub.tsv",dir=config["outputDir"])
    params:
        param_nInds
    shell:
        "paste <(shuf -n {params} {input} | cut -f1) <(yes opt | head -n {params}) > {output}"

rule mkdirLargeM:
    input:
        expand("{dir}/popmapSub.tsv",dir=config["outputDir"])
    output:
        expand("{dir}/stacksTest{{largeM}}/test{{largeM}}",dir=config["outputDir"])
    params:
        largeM="{largeM}"
    shell:
        "cat {input} | head -n {params.largeM} > {output}"

rule runStacksLargeM:
    input:
        samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
        samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
        inputDir=expand("{dir}/stacksTest{{largeM}}/test{{largeM}}",dir=config["outputDir"])
    output:
        outLog=expand("{dir}/stacksTest{{largeM}}/denovo_map.log",dir=config["outputDir"])
    params:
        largeM="{largeM}",
        popmapSub=expand("{dir}/popmapSub.tsv",dir=config["outputDir"]),
        outputDir=expand("{dir}/stacksTest{{largeM}}/",dir=config["outputDir"]),
        inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
    conda:
        "env/stacks.yaml"
    threads: 4
    shell:
        """
        denovo_map.pl --samples {params.inputDir} --popmap {params.popmapSub} -T {threads} \
        -o {params.outputDir} -n {params.largeM} -M {params.largeM}  \
        -X 'populations: -R 80'
        """

rule extractInfoLargeM:
    input:
        outLog=expand("{dir}/stacksTest{largeM}/denovo_map.log",largeM=LARGEM,dir=config["outputDir"])
    output:
        MparameterTSV=expand("{dir}/stacksTestparameter.tsv",dir=config["outputDir"])
    params:
        dir=config["outputDir"]
    shell:
        """
        paste <(cat {params.dir}/stacksTest*/denovo_map.log | grep "Kept" | cut -f2,14 -d " ") \
        <(cat {params.dir}/stacksTest*/denovo_map.log | grep "per-sample coverage" | cut -f2 -d "=" | cut -f1 -d "x") \
        <(cat {params.dir}/stacksTest*/denovo_map.log | grep "\-M" | grep "denovo_map.pl" | cut -f 13 -d " ") >\
        {params.dir}/stacksTestparameter.tsv
        """

rule makePlotLargeM:
    input:
        MparameterTSV=expand("{dir}/stacksTestparameter.tsv",dir=config["outputDir"])
    output:
        MparameterPNG=expand("{dir}/stacksTestparameter.png",dir=config["outputDir"])
    params:
        dir=config["outputDir"]
    conda:
        "env/R.yaml"
    shell:
        "Rscript src/stacksTest/parameterTest.R {params.dir}"
