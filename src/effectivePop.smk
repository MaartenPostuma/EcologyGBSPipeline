configfile: "src/effectivePop/configEP.yaml"
import pandas as pd
import os
import random
df = pd.read_csv("fraaiHertshooi/results/stacksFiles/popmapFiltered.tsv", sep='\t', dtype="object")


rule all:
    EffPop=expand("{path}/effectivePopsize/{pop}/{pop}_NE.txt",path=config["inputDir"],pop=POPS)


rule splitPopmap:
    input:
        popmap=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["inputDir"]),
    output:
        popmap=expand("{path}/effectivePopsize/popmaps/{{pop}}map.tsv",path=config["inputDir"])
    conda:
        "env/R.yaml"
    params:
        outputDir=config["inputDir"]
    shell:
        """src/effectivePop/effectivePopSplit.R {input.popmap} {params.outputDir}"""

rule perPopStacks:
    input:
        popmap=expand("{path}/effectivePopsize/{{pop}}/{{pop}}map.tsv",path=config["inputDir"])
    output:
        vcf=expand("{path}/effectivePopsize/{{pop}}/populations.snps.genepop",path=config["inputDir"]),
    params:
        inputDir=expand("{path}/",path=config["inputDir"])
        parDir=expand("{path}/effectivePopsize/{{pop}}/",path=config["inputDir"])
        maf=config["maf"]
        max_missing=config["max_missing"]
    threads:
            4
    conda:
            "env/stacks.yaml"
    resources:
            mem_mb=10000,
            runtime=30,
            cpus_per_task=4
    shell:
        """
        populations -M {input.popmap} -P {params.inputDir}/filters/max_missing~{params.max_missing}/maf{params.maf} --min-maf {params.maf} --vcf -O {params.parDir} --threads {threads}
        """ 

rule makeInfoFile:
    input:
        inputFile="src/effectivePop/NEstimatorTemplate.info"
    output:
        outFile=expand("{path}/effectivePopsize/{{pop}}/{{pop}}.info",path=config["inputDir"])
    params:
        outDir=expand("{path}/effectivePopsize/{{pop}}/",path=config["inputDir"])
    shell:
        """sed 's|inputDir|{params.outDir}/|g' {input.inputFile}"""

rule runNEestimator:
    input:
        inFile=expand("{path}/effectivePopsize/{{pop}}/{{pop}}.info",path=config["inputDir"]),
        genepop=expand("{path}/effectivePopsize/{{pop}}/populations.snps.genepop",path=config["inputDir"])
    output:
        EffPop=expand("{path}/effectivePopsize/{{pop}}/{{pop}}_NE.txt",path=config["inputDir"])
    shell:
        """src/effectivePop/NE2L i:{input.inFile}"""
    