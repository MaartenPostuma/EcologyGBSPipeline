configfile: "config.yaml"


max_missing_range=config["param_filtering"]["max_missing_range"]
individual_missingness=config["param_filtering"]["individual_missingness"]
maf_range=config["param_filtering"]["maf_range"]


rule makeGDS:
    input:
        vcf_in=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern)
    output:
        gds_out=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.wildcard_pattern)
    conda:
        "env/R.yaml"
    shell:
        '''
        R -e "SNPRelate::snpgdsVCF2GDS('{input.vcf_in}','{output.gds_out}')"
        '''

rule makeReport:
    input:
        pcaPlot=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.instance_patterns), 
        pcaVar=expand("{path}/filters/{params}/pcaVar.tsv",path=config["outputDir"],params=paramspace.instance_patterns), 
    output:
        report_out=expand("{path}/report.html",path=config["outputDir"])
    params:
        outputDir=expand("{path}/filters/",path=config["outputDir"])
    shell:
        '''
        R -e "rmarkdown::render("report.Rmd",output_file="report.html",params=list(args={params.outputDir}))"
        '''

#So this does a lot of things 
rule step_1:
    input:
        vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"]),
        popmap=expand("{path}/stacksFiles/popmap.tsv", path=config["outputDir"]),
    output:
        vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"])
    conda:
        "env/vcftools.yaml"
    params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
        indMissing=individual_missingness,
    shell:
        """vcftools --vcf {input.vcf}  --missing-indv --out {params.outputDir}/missingIndvs
        mawk '$5 > {params.indMissing}' {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/lowDP.step1.indv
        mawk '$5 < {params.indMissing}' {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/highDP.step1.indv
        cat {input.popmap} | grep -f  {params.outputDir}/highDP.step1.indv > {output.vcf}
        mkdir -p {params.outputDir}/finalVCF/
        """

rule filter:
    input:
        vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"])
    output:
        vcf=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern),
        sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
    params:
        outputDir=expand("{path}",path=config["outputDir"]),
        parDir=expand("{path}/filters/{params}",path=config["outputDir"],params=paramspace.wildcard_pattern),
        maf=str(paramspace.wildcard_pattern).split("/")[1].split("~")[1]
    threads:
        THREADSPERRUN//1
    conda:
        "env/stacks.yaml"
    shell:
        "populations -M {params.outputDir}/stacksFiles/popmapFiltered.tsv -P {params.outputDir}/stacks -R {wildcards.max_missing} --min-maf {params.maf} --vcf -O {params.parDir} --threads {threads}"    

rule makePCAData:
    input:
        gds=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.wildcard_pattern),
        popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv",path=config["outputDir"]),
        popmap=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
    output:
        pcaData=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
    conda:
        "env/R.yaml"
    params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
    shell:
        "Rscript src/filterAndFigures/pca.R {input.popmapSNPFilter} {input.gds} {input.popmap} {output.pcaData}"

rule combinePCAData:
    input:
        pcaData=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.instance_patterns)
    output:
        pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"])
    shell:
        "cat <(cat {input.pcaData} | head -n 1) <(cat {input.pcaData} | grep -v sample.id)  > {output.pcaDataAll}"