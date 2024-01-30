
#TODO MAKE THE CONDA ENVIRONMENTS

param_oligo=config["param_demultiplex"]["nOligo"]

def getParam_oligo(param_oligo):
    if param_oligo == "default" or param_oligo == "":
        id = 0.97
    else:
        id = param_oligo
    return id


rule: "clone_filter"
    input:
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"]),
        R1=expand("{path}/{R1}_R1.fq.gz",path=config["inputDir"],R1=RUN),
        R2=expand("{path}/{R2}_R2.fq.gz",path=config["inputDir"],R2=RUN)
    params:
        tmpdir=expand("{path}/{dir}", path=config["tmpDir"],dir=projectName),
        outputdir=expand("{path}/demultiplex", path=config["ouputDir"]),
        param_oligo=getParam_oligo(param_oligo)
    output:
        R1=expand("{path}/demultiplex/clone_filter/{R1}.1_R1.fq.gz",path=config["outputDir"],R1=RUN),
        R2=expand("{path}/demultiplex/clone_filter/{R2}.2_R2.fq.gz",path=config["outputDir"],R2=RUN)
    shell: 
        """
        mkdir -p {params.tmpdir}
        clone_filter -1 {input.R1} -2 {input.R2} -o {params.outputdir}/clone_filter/ --oligo_len_1 {params.param_oligo} --oligo_len_2 {params.param_oligo} --inline_inline -i gzfastq
        """


#rule: "make_stacks_files"
#    input:
#        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodes"])
#    output:
#        popmap=expand("{path}/popmap.tsv", path=config["inputDir"], bar=config["barcodes"]),
#        barcodes=expand("{path}/barcode_stacks.tsv", path=config["inputDir"], bar=config["barcodes"])
#    params:
#        inputDir=config["inputDir"]
#    conda:
#
#    shell:

