configfile: "config.yaml"
import pandas as pd
import os
import random
projectName=random.randint(1,1000000) #To ensure non overlapping tmp directories


df = pd.read_csv(os.path.join(config["inputDir"],config["barcodeFile"]), sep='\t', dtype="object").set_index('sample')
SAMPLES = df.index
RAWREADSR1 = df.rawR1.str.replace(".fq.gz","",regex=False).unique()
RAWREADSR2 = df.rawR2.str.replace(".fq.gz","",regex=False).unique()
RUN = df.rawR1.str.replace("_R1.fq.gz","",regex=False).unique()

rule all:
    input:
        R1=expand("{path}/demultiplex/clone_filter/{R1}_R1.1.fq.gz",path=config["outputDir"],R1=RUN),
        R2=expand("{path}/demultiplex/clone_filter/{R2}_R2.2.fq.gz",path=config["outputDir"],R2=RUN),
        popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"])

include: "src/demultiplex.smk"