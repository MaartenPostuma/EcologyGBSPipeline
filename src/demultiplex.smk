
#TODO MAKE THE CONDA ENVIRONMENTS

param_oligo=config["param_demultiplex"]["nOligo"]

#We always use 3 as the number of oligo for our adapters. (If we would ever change this we can account for this)
def getParam_oligo(param_oligo):
    if param_oligo == "default" or param_oligo == "":
        id = 3
    else:
        id = param_oligo
    return id

#Clone_filter removes PCR duplicates
#It compares reads that are removes reads that are completely identical 
#including the UMI/wobble/oligo's (random bases that are added at the beginning of the reads)
#as if all of these are identical it will most likely be a PCR artifact
rule clone_filter:
    input:
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"]),
        R1=expand("{path}/{{run}}_R1.fq.gz",path=config["inputDir"],R1=RUN),
        R2=expand("{path}/{{run}}_R2.fq.gz",path=config["inputDir"],R2=RUN)
    params:
        outputDir=expand("{path}/{dir}", path=config["outputDir"],dir=projectName),
        outputdir=expand("{path}/demultiplex", path=config["outputDir"]),
        param_oligo=getParam_oligo(param_oligo)
    output:
        R1=expand("{path}/demultiplex/clone_filter/{{run}}_R1.1.fq.gz",path=config["outputDir"]),
        R2=expand("{path}/demultiplex/clone_filter/{{run}}_R2.2.fq.gz",path=config["outputDir"])
    conda:
        "env/stacks.yaml"
    threads: 1
    shell: 
        """
        mkdir -p {params.outputDir}
        clone_filter -1 {input.R1} -2 {input.R2} -o {params.outputdir}/clone_filter/ --oligo_len_1 {params.param_oligo} --oligo_len_2 {params.param_oligo} --inline_inline -i gzfastq
        """

#Stacks and the rest of the pipeline need to have specific files for barcodes, 
#the format is different and we add the control nucleotide and we need to split it per run so we can demultiplex them in parallel
#popmap (to which population do samples belong),
#and an input file for SNPFilter report to add nice colours to the plots based on a priori clustering
rule make_stacks_files:
    input:
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"])
    output:
        popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"], bar=config["barcodeFile"]),
        barcodes=expand("{path}/stacksFiles/barcodeStacks{run}.tsv", path=config["outputDir"], bar=config["barcodeFile"],run=RUN)
    params:
        outputDir=expand("{path}/stacksFiles",path=config["outputDir"])
    conda:
        "env/R.yaml"
    shell:
        "Rscript src/demultiplex/createFilesFromBarcodeFile.R {input.barcodes} {params.outputDir}"

#We then do some snakemake magic to run process radtags for each run.
#This outputs 2 files per sample (the R1 and R2) 
#TODO figure out to have this not be completely dependent on the process_radtags.log file being generated.
#Howver this is complicated due to the samples needing to be split up per file...
#TODO make the output be in the tempdir if it works
rule process_radtags:
    input:
        barcodes=expand("{path}/stacksFiles/barcodeStacks{{run}}.tsv", path=config["outputDir"], bar=config["barcodeFile"]),
        R1=expand("{path}/demultiplex/clone_filter/{{run}}_R1.1.fq.gz",path=config["outputDir"]),
        R2=expand("{path}/demultiplex/clone_filter/{{run}}_R2.2.fq.gz",path=config["outputDir"])
    params:
        outputDir=expand("{path}/demultiplex/logs/{{run}}/",path=config["outputDir"])
    conda:
        "env/stacks.yaml"
    output:
        log=expand("{path}/demultiplex/logs/{{run}}/process_radtags.clone_filter.log",path=config["outputDir"])
    shell:
        "process_radtags -1 {input.R1} -2 {input.R2} -o {params.outputDir} -b {input.barcodes} --renz_1 aseI --renz_2 nsiI -c --inline-inline"

#This moves all samples into the demultiplex file
rule moveDemultiplexFiles:
    input:
        log=expand("{path}/demultiplex/logs/{run}/process_radtags.clone_filter.log",path=config["outputDir"],run=RUN)
    params:
        inputDir=expand("{path}/demultiplex/logs/",path=config["outputDir"]),
        outputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"]),
        log=expand("{path}/logs/",path=config["outputDir"])
    output:
        log=expand("{path}/demultiplex/logs/{run}/process_radtags.log",path=config["outputDir"],run=RUN),
        samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
        samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES)
    shell:
        """
        mv {params.inputDir}/*/*.fq.gz {params.outputDir}/
        """