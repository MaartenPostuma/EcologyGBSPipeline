configfile: "src/effectivePop/configEP.yaml"
import pandas as pd
import os
import random


min_samples= config["min_samples"]
colnames=["sample","population"]
df = pd.read_csv(os.path.join(config["inputDir"],"stacksFiles/popmapFiltered.tsv"), sep='\t', dtype="object",names=colnames,header=None)
vc = df.population.value_counts()
POPS=vc[vc>=min_samples].index


rule all:
    input:
        EffPop=expand("{path}/effectivePopsize/{population}/{population}_Ne.txt",path=config["inputDir"],population=POPS),
        combined=expand("{path}/effectivePopsize/combinedEffectivePopSize.tsv",path=config["inputDir"])



rule splitPopmap:
    input:
        popmap=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["inputDir"]),
    output:
        popmap=expand("{path}/effectivePopsize/{population}/{population}map.tsv",path=config["inputDir"],population=POPS)
    conda:
        "env/R.yaml"
    params:
        outputDir=expand("{path}/",path=config["inputDir"]),
        min_samples=config["min_samples"]
    resources:
            mem_mb=1000,
            runtime=15,
            cpus_per_task=1        
    shell:
        """
        Rscript src/effectivePop/effectivePopSplit.R {input.popmap} {params.outputDir}/effectivePopsize/ {params.min_samples}
        """

rule perPopStacks:
    input:
        popmap=expand("{path}/effectivePopsize/{{population}}/{{population}}map.tsv",path=config["inputDir"])
    output:
        vcf=expand("{path}/effectivePopsize/{{population}}/populations.snps.genepop",path=config["inputDir"]),
    params:
        inputDir=expand("{path}/",path=config["inputDir"]),
        parDir=expand("{path}/effectivePopsize/{{population}}/",path=config["inputDir"]),
        maf=config["maf"],
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
        populations -M {input.popmap} -P {params.inputDir}/stacks/ -R {params.max_missing} --min-mac 2  --genepop -O {params.parDir} --threads {threads}
        """ 

rule makeInfoFile:
    input:
        inputFile="src/effectivePop/NEstimatorTemplate.info"
    output:
        outFile=expand("{path}/effectivePopsize/{{population}}/{{population}}.info",path=config["inputDir"])
    params:
        outDir=expand("{path}/effectivePopsize/{{population}}/",path=config["inputDir"]),
        population="{population}"
    resources:
        mem_mb=1000,
        runtime=15,
        cpus_per_task=1            
    shell:
        """cat {input.inputFile} | sed 's|inputDir|{params.outDir}/|g'   | sed  's|pop1|{params.population}|g' > {output.outFile}"""

rule runNEestimator:
    input:
        inFile=expand("{path}/effectivePopsize/{{population}}/{{population}}.info",path=config["inputDir"]),
        genepop=expand("{path}/effectivePopsize/{{population}}/populations.snps.genepop",path=config["inputDir"])
    output:
        EffPop=expand("{path}/effectivePopsize/{{population}}/{{population}}_Ne.txt",path=config["inputDir"])
    params:
        NeEstimator=config["NeEstimatorLoc"]
    resources:
        mem_mb=10000,
        runtime=120,
        cpus_per_task=1            
    shell:
        """{params.NeEstimator}/Ne2L i:{input.inFile}"""
    
rule combineResults:
    input:
        EffPop=expand("{path}/effectivePopsize/{population}/{population}_Ne.txt",path=config["inputDir"],population=POPS)
    output:
        combined=expand("{path}/effectivePopsize/combinedEffectivePopSize.tsv",path=config["inputDir"])
    conda:
        "env/R.yaml"
    resources:
            mem_mb=1000,
            runtime=15,
            cpus_per_task=1        
    shell:
        """
        Rscript src/effectivePop/combineEffectivePopSize.R {output.combined} {input.EffPop}
        """