configfile: "config.yaml"
import pandas as pd
import os
import random

df = pd.read_csv(os.path.join(config["inputDir"],config["barcodeFile"]), sep='\t', dtype="object").set_index('sample')
#Read the barcode file and do some management to create all of the info we need (sample file / seprate into runs / make a dictionary of the both of them etc.)
df['run'] = df['rawR1'].str.replace("_R1.fq.gz","",regex=False)
df['sample']=df.index
SAMPLES = df.index
RAWREADSR1 = df.rawR1.str.replace(".fq.gz","",regex=False).unique()
RAWREADSR2 = df.rawR2.str.replace(".fq.gz","",regex=False).unique()
RUN = df.rawR1.str.replace("_R1.fq.gz","",regex=False).unique()
MODE=config["mode"]
paramspace = Paramspace(pd.read_csv("src/filterAndFigures/paramTest.tsv", sep="\t"))

rule all:
    input:
        stacksFiles=expand("{path}/archive/stacksFiles/", path=config["outputDir"]),
        vcf=expand("{path}/archive/populations.snps.vcf",path=config["outputDir"]),
        vcfFilt=expand("{path}/archive/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.instance_patterns),
        log=expand("{path}/archive/logs/{run}/perInd.tsv",path=config["outputDir"],run=RUN),
        report_out=expand("{path}/archive/report{mode}.html",path=config["outputDir"],mode=MODE),
        pcaDataAll=expand("{path}/archive/pcaAll.tsv",path=config["outputDir"]),
        treeLabels=expand("{path}/archive/treeLabelsAll.tsv",path=config["outputDir"]),
        treeSegments=expand("{path}/archive/treeSegmentsAll.tsv",path=config["outputDir"]),
        popStats=expand("{path}/archive/popStatsAll.tsv",path=config["outputDir"]),
        reads1=expand("{input}/raw/{run}_R1.fq.gz",input=config["outputDir"],run=RUN),
        reads2=expand("{input}/raw/{run}_R2.fq.gz",input=config["outputDir"],run=RUN),
        barcodes=expand("{path}/raw/{bar}", path=config["outputDir"], bar=config["barcodeFile"])


rule mvStacksFiles:
    input:
        stacksFiles=expand("{path}/stacksFiles/", path=config["outputDir"]),    
    output:
        stacksFiles=expand("{path}/archive/stacksFiles/", path=config["outputDir"]),
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
    shell:
        """mv {input.stacksFiles} {output.stacksFiles}"""

rule mvVCF:
    input:
        vcf=expand("{path}/stacks/population.snps.vcf", path=config["outputDir"]),    
    output:
        vcf=expand("{path}/archive/populations.snps.vcf",path=config["outputDir"])
 	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
   shell:
        """mv {input.vcf} {output.vcf}"""

rule mvVCFfilt:
    input:
        vcf=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern),
    output:
        vcf=expand("{path}/archive/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern),
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
    shell:
        """mv {input.vcf} {output.vcf}"""

rule mvDemultiLogs:
    input:
        log=expand("{path}/demultiplx/logs/{{run}}/perInd.tsv",path=config["outputDir"]),
    output:
        log=expand("{path}/archive/logs/{{run}}/perInd.tsv",path=config["outputDir"]),
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
    shell:
        """mv {input.log} {output.log}"""

rule mvReport:
    input:
        report_out=expand("{path}/report{mode}.html",path=config["outputDir"],mode=MODE)
    output:
        report_out=expand("{path}/archive/report{mode}.html",path=config["outputDir"],mode=MODE)
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
    shell:
        """mv {input.report_out} {output.report_out}"""

rule mvReportData:
    input:
        pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"]),
        treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
        treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"]),
        popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"])
    output:
        pcaDataAll=expand("{path}/archive/pcaAll.tsv",path=config["outputDir"]),
        treeLabels=expand("{path}/archive/treeLabelsAll.tsv",path=config["outputDir"]),
        treeSegments=expand("{path}/archive/treeSegmentsAll.tsv",path=config["outputDir"]),
        popStats=expand("{path}/archive/popStatsAll.tsv",path=config["outputDir"])
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1
    shell:
        """
        mv {input.pcaDataAll} {output.pcaDataAll}
        mv {input.treeLabels} {output.treeLabels}
        mv {input.treeSegments} {output.treeSegments}
        mv {input.popStats} {output.popStats}
        """

rule mvRawData:
    input:
        reads1=expand("{input}/{{run}}_R1.fq.gz",input=config["inputDir"]),
        reads2=expand("{input}/{{run}}_R2.fq.gz",input=config["inputDir"]),
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"])
    output:
        reads1=expand("{input}/raw/{{run}}_R1.fq.gz",input=config["outputDir"]),
        reads2=expand("{input}/raw/{{run}}_R2.fq.gz",input=config["outputDir"]),
        barcodes=expand("{path}/raw/{bar}", path=config["outputDir"], bar=config["barcodeFile"])
	resources:
		mem_mb= 1000,
		runtime= 1,
		cpus_per_task= 1        
    shell:
        """
        mv {input.pcaDataAll} {output.pcaDataAll}
        mv {input.treeLabels} {output.treeLabels}
        mv {input.treeSegments} {output.treeSegments}
        """
