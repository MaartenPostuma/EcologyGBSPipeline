configfile: "config.yaml"
max_missing_range=config["max_missing_range"]
maf_range=config["maf_range"]
hwe_range=config["hwe_range"]
DP_range=config["DP_range"]


rule makeGDS:
	input:
	    vcf_out="vcf/{filter}/{filter}{Value}.recode.vcf"
	output:
	    gds_out="vcf/{filter}/{filter}{Value}.recode.gds"
	shell:
		'''
		R -e "SNPRelate::snpgdsVCF2GDS('{input.vcf_out}','{output.gds_out}')"
		'''

rule makeReport:
    input:
        gds_out1=expand("vcf/filter1/filter1_{MAX_MISSING}.recode.gds",MAX_MISSING=config["max_missing_range"]),
        gds_out2=expand("vcf/filter2/filter2_{MAF}.recode.gds",MAF=config["maf_range"]),
	    gds_out3=expand("vcf/filter3/filter3_{HWE}.recode.gds",HWE=config["hwe_range"]),
	    gds_out4=expand("vcf/filter4/filter4_{DP}.recode.gds",DP=config["DP_range"]),
        gds_out5=expand("vcf/filter5/filter5_{MAX_MISSING}.recode.gds",MAX_MISSING=config["max_missing_range"]),
        gds_out6=expand("vcf/filter6/filter6_{MAF}.recode.gds",MAF=config["maf_range"]),
	    gds_out7=expand("vcf/filter7/filter7_{HWE}.recode.gds",HWE=config["hwe_range"]),
	    gds_out8=expand("vcf/filter8/filter8_{DP}.recode.gds",DP=config["DP_range"])
    output:
        report_out="report.html"
    shell:
        """
        R -e "rmarkdown::render('report.Rmd',output_file='report.html')"
        """


rule step_1:
    input:
        vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"]),
    output:
        vcf_out=expand("{path}/filters/step1.recode.vcf",path=config["outputDir"])
    conda:
        "env/vcftools.yaml"
    params:
        indMissing=config["individual_missingness"],
        outputDir=expand("{path}/filters/",path=config["outputDir"])
    shell:
        '''vcftools --vcf {input.vcf_in}  --missing-indv --out {params.outputDir}/missingIndvs
        mawk '$5 < {params.indMissing}' vcf/missingIndvs.imiss | cut -f1 > {params.outputDir}/highDP.step1.indv
        vcftools --vcf {input.vcf_in}  --remove {params.outputDir}/lowDP.step2.indv --recode --out {params.outputDir}/step1 --recode-INFO-all'''

rule filter1:
    input:
        vcf_in='vcf/step1.recode.vcf'
    output:
        vcf_out=expand("vcf/filter1/filter1_{{max_missing_range}}.recode.vcf")
    params:
        maxMissing='{max_missing_range}'
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter1/filter1_{params.maxMissing} --max-missing {params.maxMissing}'


rule filter2:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter2/filter2_{{maf_range}}.recode.vcf")
    params:
        mafRange='{maf_range}'
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter2/filter2_{params.mafRange} --maf {params.mafRange}'



rule filter3:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter3/filter3_{{hwe_range}}.recode.vcf")
    params:
        hweRange='{hwe_range}'
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter3/filter3_{params.hweRange} --hwe {params.hweRange}'

rule filter4:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter4/filter4_{{DP_range}}.recode.vcf")
    params:
        dpRange='{DP_range}'
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter4/filter4_{params.dpRange} --minDP {params.dpRange}'


rule filter5:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter5/filter5_{{max_missing_range}}.recode.vcf")
    params:
        maxMissing='{max_missing_range}',
	    hweVal=config["hwe_val"],
	    mafVal=config["maf_val"],
	    dpVal=config["DP_val"]
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter5/filter5_{params.maxMissing} --max-missing {params.maxMissing} --hwe {params.hweVal} --maf {params.mafVal} --minDP {params.dpVal} '



rule filter6:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter6/filter6_{{maf_range}}.recode.vcf")
    params:
        maxMissing=config["max_missing_val"],
	    hweVal=config["hwe_val"],
	    mafVal='{maf_range}',
	    dpVal=config["DP_val"]
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter6/filter6_{params.mafVal} --max-missing {params.maxMissing} --hwe {params.hweVal} --maf {params.mafVal} --minDP {params.dpVal} '

rule filter7:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter7/filter7_{{hwe_range}}.recode.vcf")
    params:
        maxMissing=config["max_missing_val"],
	    hweVal='{hwe_range}',
	    mafVal=config["maf_val"],
	    dpVal=config["DP_val"]
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter7/filter7_{params.hweVal} --max-missing {params.maxMissing} --hwe {params.hweVal} --maf {params.mafVal} --minDP {params.dpVal} '

	
rule filter8:
    input:
        vcf_in='vcf/step2.recode.vcf'
    output:
        vcf_out=expand("vcf/filter8/filter8_{{DP_range}}.recode.vcf")
    params:
        maxMissing=config["max_missing_val"],
	    hweVal=config["hwe_val"],
	    mafVal=config["maf_val"],
	    dpVal='{DP_range}'
    shell:
        'vcftools --vcf {input.vcf_in} --recode --recode-INFO-all --out vcf/filter8/filter8_{params.dpVal} --max-missing {params.maxMissing} --hwe {params.hweVal} --maf {params.mafVal} --minDP {params.dpVal} '
