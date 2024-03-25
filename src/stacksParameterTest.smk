param_nInds=config["param_StacksTest"]["nInds"]
def getParam_nInds(param_nInds):
    if param_nInds == "default" or param_oligo == "":
        id = 20
    else:
        id = param_nInds
    return id


LARGEM = list(range(1,12))
#rule all:
#    input:
#        MparameterPNG=expand("{dir}/stacksTestparameter.png",dir=config["outputDir"]),
#        testLargeM=expand("{dir}/stacksTest{largeM}/test{largeM}",largeM=LARGEM,dir=config["outputDir"])

rule subset_popmap:
    input:
        popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"]),
    output:
        expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"])
    params:
        getParam_nInds(param_nInds)
    shell:
        "paste <(shuf -n {params} {input} | cut -f1) <(yes opt | head -n {params}) > {output}"

rule mkdirLargeM:
    input:
        expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"])
    output:
        expand("{dir}/stacksTest/M{{largeM}}/test{{largeM}}",dir=config["outputDir"])
    params:
        largeM="{largeM}"
    shell:
        "cat {input} | head -n {params.largeM} > {output}"

rule runStacksLargeM:
    input:
        samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
        samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
        popmapSub= expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"]),
        inputDir=expand("{dir}/stacksTest/M{{largeM}}/test{{largeM}}",dir=config["outputDir"])
    output:
        outLog=expand("{dir}/stacksTest/M{{largeM}}/denovo_map.log",dir=config["outputDir"])
    params:
        largeM="{largeM}",
        outputDir=expand("{dir}/stacksTest/M{{largeM}}/",dir=config["outputDir"]),
        inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
    conda:
        "env/stacks.yaml"
    threads: 4
    shell:
        """
        denovo_map.pl --samples {params.inputDir} --popmap {input.popmapSub} -T {threads} \
        -o {params.outputDir} -n {params.largeM} -M {params.largeM}  \
        -X 'populations: -R 80'
        """

rule extractInfoLargeM:
    input:
        outLog=expand("{dir}/stacksTest/M{largeM}/denovo_map.log",largeM=LARGEM,dir=config["outputDir"])
    output:
        MparameterTSV=expand("{dir}/stacksTest/parameter.tsv",dir=config["outputDir"])
    params:
        dir=config["outputDir"]
    shell:
        """
        paste <(cat {params.dir}/stacksTest/*/denovo_map.log | grep "Kept" | cut -f2,14 -d " ") \
        <(cat {params.dir}/stacksTest/*/denovo_map.log | grep "per-sample coverage" | cut -f2 -d "=" | cut -f1 -d "x") \
        <(cat {params.dir}/stacksTest/*/denovo_map.log | grep "\-M" | grep "denovo_map.pl" | cut -f 13 -d " ") >\
        {params.dir}/stacksTest/parameter.tsv
        """

rule makePlotLargeM:
    input:
        MparameterTSV=expand("{dir}/stacksTest/parameter.tsv",dir=config["outputDir"])
    output:
        MparameterPNG=expand("{dir}/stacksTest/parameter.png",dir=config["outputDir"])
    params:
        dir=expand("{dir}/stacksTest/",dir=config["outputDir"])
    conda:
        "env/R.yaml"
    shell:
        "Rscript src/stacksTest/parameterTest.R {params.dir}"
