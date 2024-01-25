
#TODO MAKE THE CONDA ENVIRONMENTS

param_oligo=config["param_denovo"]["identity"]

def getParam_oligo(param_oligo):
    if param_oligo == "default" or param_oligo == "":
        id = 0.97
    else:
        id = param_oligo
    return id


rule: "clone_filter"
    input:
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodes"]),
        R1=expand("{path}/{R1}",path=config["inputDir"],R1=RAWREADSR1),
        R2=expand("{path}/{R2}",path=config["inputDir"],R2=RAWREADSR2)
    params:
        tmpdir=expand("{path}/{dir}", path=config["tmpdir"],dir=projectName),
        outputdir=expand("{path}/output_demultiplex",  path=config["output_dir"]) 
        param_oligo=getParam_oligo(param_oligo)
    output:
    conda:
    shell: 
        """
        mkdir -p {params.tmpdir}
        clone_filter -1 {input.R1} -2 {input.R2} -o {output.outputdir}/clone_filter/ --oligo_len_1 {params.param_oligo} --oligo_len_2 {params.param_oligo} --inline_inline -i gzfastq
        """


rule: "make_stacks_files"
    input:
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodes"])
    output:
        popmap=expand("{path}/popmap.tsv", path=config["inputDir"], bar=config["barcodes"]),
        barcodes=expand("{path}/barcode_stacks.tsv", path=config["inputDir"], bar=config["barcodes"])
    params:
        inputDir=config["inputDir"]
    conda: 

    shell:

