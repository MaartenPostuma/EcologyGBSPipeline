rule map_bwa:
    input:
        samplesR1=expand("{path}/demultiplex/samples/{{samples}}.1.fq.gz",path=config["outputDir"]),
        samplesR2=expand("{path}/demultiplex/samples/{{samples}}.2.fq.gz",path=config["outputDir"]),
        ref=expand("{ref}",ref=config["reference"])
    output:
        bam=expand("{path}/refMapping/{{samples}}.bam",path=config["outputDir"]),
    conda:
        "env/bwa-mem.yaml"
    threads: THREADSPERRUN
    shell:
        """bwa-mem2 mem -t {threads} {input.ref} {input.samplesR1} {input.samplesR2} | samtools view -buS - > {output.bam}"""
