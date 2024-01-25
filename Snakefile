configfile: "configParameter.yaml"
import pandas as pd
import os
import random
df = pd.read_csv(os.path.join("src/demultiplex/barcodeTemplate.csv"), sep=';', dtype="object").set_index('sample')

#df = pd.read_csv(os.path.join(config["input_dir"],config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index
RAWREADSR1 = pd.unique(df.rawR1).astype(str)
RAWREADSR2 = pd.unique(df.rawR2)


RAWREADSR1.replace(".fq.gz","")