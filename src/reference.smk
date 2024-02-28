rule make_bwa_mem_index:
    input:
        ref=expand("{ref}",ref=config["reference"])
    output:
        refIndex=expand("{ref}.bwt.2bit.64",ref=config["reference"])
    conda:
        "env/bwa-mem.yaml"
    shell:
        "bwa-mem2 index {input.ref}"

rule map_bwa:
    input:
        samplesR1=expand("{path}/demultiplex/samples/{{samples}}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
        samplesR2=expand("{path}/demultiplex/samples/{{samples}}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
        ref=expand("{ref}",ref=config["reference"]),
        refIndex=expand("{ref}.bwt.2bit.64",ref=config["reference"])
    output:
        bam=temp(expand("{path}/refMapping/firstBam/{{samples}}.bam",path=config["tmpDir"])),
    threads:
        workflow.cores/5
    conda:
        "env/bwa-mem.yaml"
    shell:
        "bwa-mem2 mem -t {threads} {input.ref} {input.samplesR1} {input.samplesR2} | samtools view -buS - > {output.bam}"

rule sort_picard:
    input:
        bam=expand("{path}/refMapping/firstBam/{{samples}}.bam",path=config["tmpDir"]),
    output:
        sortedBam=temp(expand("{path}/refMapping/sorted/{{samples}}.bam",path=config["tmpDir"]))
    conda:
        "env/picard.yaml"
    shell:
        "picard SortSam -I {input.bam} -O {output.sortedBam} -SO coordinate"

rule add_RG:
    input:
        sortedBam=expand("{path}/refMapping/sorted/{{samples}}.bam",path=config["tmpDir"])
    output:
        RGBam=expand("{path}/refMapping/RGBams/{{samples}}.bam",path=config["tmpDir"])
    conda:
        "env/picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups \
        I={input.sortedBam}\
        O={output.RGBam} \
        RGLB={wildcards.samples} \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM={wildcards.samples} \
        RGID={wildcards.samples}
        """

rule merge_sort_bam:
    input:
        RGBam=expand("{path}/refMapping/RGBams/{samples}.bam",path=config["tmpDir"],samples=SAMPLES)
    output:
        mergedBam=expand("{path}/refOut/merged.bam",path=config["outputDir"])
    conda:
        "env/samtools.yaml"
    threads: workflow.cores
    shell:
        "samtools merge - {input.RGBam} | samtools sort -@ {threads} > {output.mergedBam}"
  
rule indexRef:
    input:
        ref=expand("{ref}",ref=config["reference"])
    output:
        refIndex=expand("{ref}.fai",ref=config["reference"])
    conda:
        "env/samtools.yaml"
    shell:
        "samtools faidx {input.ref}"

rule indexBam:
    input:
        RGBam=expand("{path}/refOut/merged.bam",path=config["outputDir"]),
    output:
        RGBamIndex=expand("{path}/refOut/merged.bam.bai",path=config["outputDir"]),
    conda:
        "env/samtools.yaml"
    shell:
        "samtools index {input.RGBam}"

rule variantCall:
    input:
        RGBamIndex=expand("{path}/refOut/merged.bam.bai",path=config["outputDir"],samples=SAMPLES),
        refIndex=expand("{ref}.fai",ref=config["reference"]),
        ref=expand("{ref}",ref=config["reference"]),
        bam=expand("{path}/refOut/merged.bam",path=config["outputDir"])
    output:
        vcf=expand("{path}/refVCF/output.vcf.gz",path=config["outputDir"])
    threads: min(workflow.cores,4096/(len(SAMPLES)*2)-4)
    conda:
        "env/freebayes.yaml"
    shell:
        """
        freebayes-parallel <(fasta_generate_regions.py {input.refIndex} 100000) {threads} -f {input.ref} {input.bam} --no-partial-observations --report-genotype-likelihood-max --genotype-qualities --min-coverage 0 --min-base-quality 1 --min-mapping-quality 10 | bgzip -c > {output.vcf}
        """   