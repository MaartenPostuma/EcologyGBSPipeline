
rule denovo_map:
    input: 
        samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
        samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
        popmapSub= expand("{dir}/stacksFiles/popmap.tsv",dir=config["outputDir"])
    params:
        M=config["M"],
        outputDir=expand("{path}/stacks/",path=config["outputDir"]),
        inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
    output:
        vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"])
    conda:
        "env/stacks.yaml"
    threads:workflow.cores
    shell:
        "denovo_map.pl --samples {params.inputDir} --popmap {input.popmap} -T  {threads}       -o outputDir -n {params.M} -m {params.M} -X 'populations: --vcf'"