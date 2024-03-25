configfile: "config.yaml"
import pandas as pd
import os
import random
projectName=random.randint(1,1000000) #To ensure non overlapping tmp directories

#df = pd.read_csv(os.path.join("data/barcodes.txt"), sep='\t', dtype="object").set_index('sample')

df = pd.read_csv(os.path.join(config["inputDir"],config["barcodeFile"]), sep='\t', dtype="object").set_index('sample')
df['run'] = df['rawR1'].str.replace("_R1.fq.gz","",regex=False)
df['sample']=df.index
SAMPLES = df.index
RAWREADSR1 = df.rawR1.str.replace(".fq.gz","",regex=False).unique()
RAWREADSR2 = df.rawR2.str.replace(".fq.gz","",regex=False).unique()
RUN = df.rawR1.str.replace("_R1.fq.gz","",regex=False).unique()
THREADSPERRUN=workflow.cores/RUN.size
MODE=config["mode"]
from snakemake.utils import Paramspace
paramspace = Paramspace(pd.read_csv("src/filterAndFigures/paramTest.tsv", sep="\t"))
grouped = df.groupby("run")["sample"].apply(set)
LANESAMPLE = grouped.to_dict()
DUPES=df['sample'].duplicated().any()

#TODO get this to work, make a seperate unique sample list use that to cat all samples with the same name together?
#For loop? or snakemake weird shit?
#TODO fix one lane!
#
if config["mode"]== "StacksTest":
    rule all:
        input:
            DeDuplR1=expand("{path}/demultiplex/clone_filter/{run}_R1.1.fq.gz",path=config["outputDir"],run=RUN),
            DeDuplR2=expand("{path}/demultiplex/clone_filter/{run}_R2.2.fq.gz",path=config["outputDir"],run=RUN),
            popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"]),
            popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv", path=config["outputDir"]),
            perRUNDemulti=expand("{path}/demultiplex/logs/{run}/process_radtags.clone_filter.log",path=config["outputDir"],run=RUN),
            samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
            samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
            MparameterPNG=expand("{dir}/stacksTest/parameter.png",dir=config["outputDir"])


if config["mode"]== "StacksTest":
    include: "src/demultiplex.smk"
    include: "src/stacksParameterTest.smk"


if config["mode"]== "Denovo":
    include: "src/demultiplex.smk"
    include: "src/stacks.smk"
    include: "src/filterAndFigures.smk"

if config["mode"]== "Denovo":
    rule all:
        input:
            DeDuplR1=expand("{path}/demultiplex/clone_filter/{run}_R1.1.fq.gz",path=config["outputDir"],run=RUN),
            DeDuplR2=expand("{path}/demultiplex/clone_filter/{run}_R2.2.fq.gz",path=config["outputDir"],run=RUN),
            popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"]),
            popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv", path=config["outputDir"]),
            popmapFiltered=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
            perRUNDemulti=expand("{path}/demultiplex/logs/{run}/process_radtags.clone_filter.log",path=config["outputDir"],run=RUN),
            samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
            samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
            vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"]),
            vcfFilt=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.instance_patterns),
            sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
            gds=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.instance_patterns),
            pcaData=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
            pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"]),
            treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
            treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"]),
            popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"]),
            report_out=expand("{path}/report{mode}.html",path=config["outputDir"],mode=MODE)


if config["mode"]== "Reference":
    include: "src/demultiplex.smk"
    include: "src/reference.smk"
    include: "src/filterAndFigures.smk"


if config["mode"]== "Reference":
    rule all:
        input:
            DeDuplR1=expand("{path}/demultiplex/clone_filter/{run}_R1.1.fq.gz",path=config["outputDir"],run=RUN),
            DeDuplR2=expand("{path}/demultiplex/clone_filter/{run}_R2.2.fq.gz",path=config["outputDir"],run=RUN),
            perRUNDemulti=expand("{path}/demultiplex/logs/{run}/process_radtags.clone_filter.log",path=config["outputDir"],run=RUN),
            samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
            samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
            RGBam=expand("{path}/refOut/merged.bam",path=config["outputDir"]),
            vcf=expand("{path}/refOut/populations.vcf.gz",path=config["outputDir"]),
            pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"]),
            treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
            treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"]),
            popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"]),
            report_out=expand("{path}/report{mode}.html",path=config["outputDir"],mode=MODE)


if config["mode"]== "ReportReference":
    rule all:
        input:
            vcf=expand("{path}/refOut/populations.vcf.gz",path=config["outputDir"]),
            pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"]),
            treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
            treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"]),
            popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"]),
            report_out=expand("{path}/report{mode}.html",path=config["outputDir"],mode=MODE)


if config["mode"]== "ReportDenovo":
    rule all:
        input:
            vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"]),
            vcfFilt=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.instance_patterns),
            sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
            gds=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.instance_patterns),
            pcaData=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
            pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"]),
            treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
            treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"]),
            popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"]),
            report_out=expand("{path}/report{mode}.html",path=config["outputDir"],mode=MODE)
