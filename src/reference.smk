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
        bam=expand("{path}/refMapping/firstBam/{{samples}}.bam",path=config["outputDir"]),
    threads:
        workflow.cores/5
    conda:
        "env/bwa-mem.yaml"
    shell:
        "bwa-mem2 mem -t {threads} {input.ref} {input.samplesR1} {input.samplesR2} | samtools view -buS - > {output.bam}"

rule sort_picard:
    input:
        bam=expand("{path}/refMapping/firstBam/{{samples}}.bam",path=config["outputDir"]),
    output:
        sortedBam=expand("{path}/refMapping/sorted/{{samples}}.bam",path=config["outputDir"])
    conda:
        "env/picard.yaml"
    shell:
        "picard SortSam -I {input.bam} -O {output.sortedBam} -SO coordinate"

rule add_RG:
    input:
        sortedBam=expand("{path}/refMapping/sorted/{{samples}}.bam",path=config["outputDir"])
    output:
        RGBam=expand("{path}/refMapping/{{samples}}.bam",path=config["outputDir"])
    conda:
        "env/picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups \
        I={input.sortedBam}\
        O={output.RGBam} \
        RGLB={wildcards.sample} \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM={wildcards.sample}
        """

# rule makeBamList:
#     input:
#         RGBam=expand("{path}/refMapping/{samples}.bam",path=config["outputDir"],samples=SAMPLES),
#     output:
#         bamList=expand("{path}/refMapping/bamList.txt",path=config["outputDir"])
#     shell:
#         "ls {input.RGBam} > {output.bamList}"
  
# rule indexRef:
#     input:
#         ref=expand("{ref}",ref=config["reference"])
#     output:
#         refIndex=expand("{ref}.fai",ref=config["reference"])
#     conda:
#         "env/samtools.yaml"
#     shell:
#         "samtools faidx {input.ref}"

# rule indexBam:
#     input:
#         RGBam=expand("{path}/refMapping/{{samples}}.bam",path=config["outputDir"]),
#     output:
#         RGBamIndex=expand("{path}/refMapping/{{samples}}.bam.bai",path=config["outputDir"]),
#     conda:
#         "env/samtools.yaml"
#     shell:
#         "samtools index {input.RGBam}"

# rule variantCall:
#     input:
#         RGBamIndex=expand("{path}/refMapping/{samples}.bam.bai",path=config["outputDir"],samples=SAMPLES),
#         refIndex=expand("{ref}.fai",ref=config["reference"]),
#         ref=expand("{ref}",ref=config["reference"]),
#         bamList=expand("{path}/refMapping/bamList.txt",path=config["outputDir"])
#     output:
#         vcf=expand("{path}/refVCF/output.vcf.gz",path=config["outputDir"])
#     conda:
#         "env/freebayes.html"
#     shell:
#         """
#         freebayes-parallel <(fasta_generate_regions.py {input.refIndex} 100000) 6 -f {input.ref}  \
#             --bam-list {input.bamList} \
#             --no-partial-observations \
#             --report-genotype-likelihood-max \
#             --genotype-qualities \
#             --min-coverage 0 \
#             --min-base-quality 1 \
#             --min-mapping-quality 10 | bgzip -c > {output.vcf}
#         """   