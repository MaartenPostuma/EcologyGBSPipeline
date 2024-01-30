configfile: "config.yaml"
import pandas as pd
import os
import random
df = pd.read_csv(os.path.join("src/demultiplex/barcodeTemplate.csv"), sep=';', dtype="object").set_index('sample')

#df = pd.read_csv(os.path.join(config["input_dir"],config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index
RAWREADSR1 = df.rawR1.str.replace(".fq.gz","").unique()
RAWREADSR2 = df.rawR2.str.replace(".fq.gz","").unique()
RUN = df.rawR1.str.replace("_R1.fq.gz","").unique()

rule all:
    input:
        R1=expand("{path}/demultiplex/{R1}.1_R1.fq.gz",path=config["outputDir"],R1=RUN),
        R2=expand("{path}/demultiplex/{R2}.2_R2.fq.gz",path=config["outputDir"],R2=RUN)