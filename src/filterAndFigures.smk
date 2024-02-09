from snakemake.utils import Paramspace
configfile: "config.yaml"


max_missing_range=config["param_filtering"]["max_missing_range"]
individual_missingness=config["param_filtering"]["individual_missingness"]
maf_range=config["param_filtering"]["maf_range"]
paramspace = Paramspace(pd.read_csv("filterAndFigures/paramTest.tsv", sep="\t"))


rule makeGDS:
	input:
	    vcf="vcf/{filter}/{filter}{Value}.recode.vcf"
	output:
	    gds_out="vcf/{filter}/{filter}{Value}.recode.gds"
	shell:
		"""
		R -e "SNPRelate::snpgdsVCF2GDS("{input.vcf}","{output.gds_out}")"
		"""

rule makeReport:
    input:
        gds_out1=expand("vcf/filter1/filter1_{MAX_MISSING}.recode.gds",MAX_MISSING=max_missing_range),
        gds_out2=expand("vcf/filter2/filter2_{MAF}.recode.gds",MAF=maf_range),
        gds_out5=expand("vcf/filter5/filter5_{MAX_MISSING}.recode.gds",MAX_MISSING=max_missing_range),
        gds_out6=expand("vcf/filter6/filter6_{MAF}.recode.gds",MAF=maf_range)
    output:
        report_out="report.html"
    shell:
        """
        R -e "rmarkdown::render("report.Rmd",output_file="report.html")"
        """

#So this does a lot of things 
rule step_1:
    input:
        vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"]),
    output:
        vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"])
    conda:
        "env/vcftools.yaml"
    params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
        indMissing=individual_missingness,
    shell:
        """vcftools --vcf {input.vcf}  --missing-indv --out {params.outputDir}/missingIndvs
        mawk "$5 > {params.indMissing}" {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/lowDP.step1.indv
        mawk "$5 < {params.indMissing}" {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/highDP.step1.indv
        cat {params.outputDir}/stacksFiles/popmap.tsv | grep -f  {params.outputDir}/filters/highDP.step1.indv > {params.outputDir}/stacksFiles/popmapFiltered.tsv
        mkdir -p {params.outputDir}/finalVCF/
        """

rule filter:
    input:
        vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"])
    output:
        vcf=expand("{path}/filters/{params}/step1.recode.p.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern),
        sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
    params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
        parDir=expand("{path}/filters/{params}",path=config["outputDir"],params=paramspace.wildcard_pattern),
        maf=str(paramspace.wildcard_pattern).split("/")[1].split("~")[1]
    threads:
        THREADSPERRUN//1
    conda:
        "env/stacks.yaml"
    shell:
        "populations -M {params.outputDir}/stacksFiles/popmapFiltered.tsv -P {params.outputDir}/stacks -R {wildcards.max_missing} --min-maf {params.maf} --vcf -O {params.outputDir}/filters/ --threads {threads}"    


rule filter1:
    input:
        vcf=expand("{path}/filters/step1.recode.vcf",path=config["outputDir"])
    output:
        vcf=expand("vcf/filter1/filter1_{{max_missing_range}}.recode.vcf")
    params:
        outputDir=expand("{path}/filter/",path=config["outputDir"]),
        maxMissing="{max_missing_range}"
    shell:
        "vcftools --vcf {input.vcf} --recode --recode-INFO-all --out vcf/filter1/filter1_{params.maxMissing} --max-missing {params.maxMissing}"


rule filter2:
    input:
        vcf=expand("{path}/filters/step1.recode.vcf",path=config["outputDir"])
    output:
        vcf=expand("vcf/filter2/filter2_{{maf_range}}.recode.vcf")
    params:
        outputDir=expand("{path}/filter/",path=config["outputDir"]),
        mafRange="{maf_range}"
    shell:
        "vcftools --vcf {input.vcf} --recode --recode-INFO-all --out vcf/filter2/filter2_{params.mafRange} --maf {params.mafRange}"


rule filter5:
    input:
        vcf=expand("{path}/filters/step1.recode.vcf",path=config["outputDir"])
    output:
        vcf=expand("vcf/filter5/filter5_{{max_missing_range}}.recode.vcf")
    params:
        outputDir=expand("{path}/filter/",path=config["outputDir"]),
        maxMissing="{max_missing_range}",
	    hweVal=config["hwe_val"],
	    mafVal=config["maf_val"],
	    dpVal=config["DP_val"]
    shell:
        "vcftools --vcf {input.vcf} --recode --recode-INFO-all --out vcf/filter5/filter5_{params.maxMissing} --max-missing {params.maxMissing} --hwe {params.hweVal} --maf {params.mafVal} --minDP {params.dpVal} "



rule filter6:
    input:
        vcf=expand("{path}/filters/step1.recode.vcf",path=config["outputDir"])
    output:
        vcf=expand("vcf/filter6/filter6_{{maf_range}}.recode.vcf")
    params:
        outputDir=expand("{path}/filter/",path=config["outputDir"]),
        maxMissing=config["max_missing_val"],
	    hweVal=config["hwe_val"],
	    mafVal="{maf_range}",
	    dpVal=config["DP_val"]
    shell:
        "vcftools --vcf {input.vcf} --recode --recode-INFO-all --out vcf/filter6/filter6_{params.mafVal} --max-missing {params.maxMissing} --hwe {params.hweVal} --maf {params.mafVal} --minDP {params.dpVal} "