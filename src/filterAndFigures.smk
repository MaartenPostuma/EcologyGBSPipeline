configfile: "config.yaml"


individual_missingness=config["param_filtering"]["individual_missingness"]


rule makeGDS:
     input:
        vcf=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern)
     output:
        gds_out=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.wildcard_pattern)
     conda:
        "env/R.yaml"
     resources:
        mem_mb=lambda wc, input: (0.5 * input.vcf.size_mb),
        runtime= 15,
        cpus_per_task= 1
     shell:
        '''
        R -e "SNPRelate::snpgdsVCF2GDS('{input.vcf}','{output.gds_out}')"
        '''

rule makeReport:
     input:
        iMissing=expand('{path}/filters/missingIndvs.imiss',path=config["outputDir"]),
        pcaPlot=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
        treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
        treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"]),
        popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"]),
        fstStats=expand("{path}/filters/fstStatsAll.tsv",path=config["outputDir"]),
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"])
     output:
        report_out=expand("{path}/report{mode}.html",path=config["outputDir"],mode=MODE)
     params:
        outputDir=expand("{path}/",path=config["outputDir"]),
     resources:
        mem_mb=lambda wc, input: (0.5 * input.vcf.size_mb),
        runtime= 15,
        cpus_per_task= 1
     conda:
        "env/R.yaml"
     shell:
        '''
        R -e "rmarkdown::render('src/filterAndFigures/report.Rmd',output_file='report.html',params=list(args=c('{params.outputDir}','{input.barcodes}')))"
        mv src/filterAndFigures/report.html {output.report_out}
        '''

#So this does a lot of things
#It determines the % missingness at 80% max-missing per SNP and uses this to filter!
if config["mode"]== "Denovo": 
     rule step_1:
        input:
             vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"]),
             popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
        output:
             vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
             iMissing=expand('{path}/filters/missingIndvs.imiss',path=config["outputDir"])
        conda:
             "env/vcftools.yaml"
        params:
             outputDir=expand("{path}/filters/",path=config["outputDir"]),
             indMissing=individual_missingness,
        resources:
            mem_mb=lambda wc, input: (0.5 * input.vcf.size_mb),
            runtime= 30,
            cpus_per_task= 1
        shell:
             """vcftools --vcf {input.vcf}  --missing-indv --out {params.outputDir}/missingIndvs --max-missing 0.5
             mawk '$5 > {params.indMissing}' {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/lowDP.step1.indv
             mawk '$5 < {params.indMissing}' {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/highDP.step1.indv
             cat {input.popmap} | grep -f  {params.outputDir}/highDP.step1.indv > {output.vcf}
             """

if config["mode"]== "Denovo":
     rule filter:
        input:
             vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"])
        output:
             vcf=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern),
             sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern),
             FstSumstats=expand("{path}/filters/{params}/populations.fst_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
        params:
             outputDir=expand("{path}",path=config["outputDir"]),
             parDir=expand("{path}/filters/{params}",path=config["outputDir"],params=paramspace.wildcard_pattern),
             maf=str(paramspace.wildcard_pattern).split("/")[1].split("~")[1]
        threads:
             4
        conda:
             "env/stacks.yaml"
        shell:
             "populations -M {params.outputDir}/stacksFiles/popmapFiltered.tsv -P {params.outputDir}/stacks -R {wildcards.max_missing} --min-maf {params.maf} --vcf -O {params.parDir} --fstats --threads {threads}"     

if config["mode"]== "Reference": 
     rule step_1:
        input:
             vcf=expand("{path}/refOut/populations.vcf.gz",path=config["outputDir"]),
             popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
        output:
             vcf=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
             iMissing=expand('{path}/filters/missingIndvs.imiss',path=config["outputDir"])
        conda:
             "env/vcftools.yaml"
        params:
             outputDir=expand("{path}/filters/",path=config["outputDir"]),
             indMissing=individual_missingness,
        resources:
            mem_mb=lambda wc, input: (0.5 * input.vcf.size_mb),
            runtime= 30,
            cpus_per_task= 4
     
        shell:
             """
             vcftools --gzvcf {input.vcf}  --missing-indv --out {params.outputDir}/missingIndvs --max-missing 0.5
             mawk '$5 > {params.indMissing}' {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/lowDP.step1.indv
             mawk '$5 < {params.indMissing}' {params.outputDir}/missingIndvs.imiss | cut -f1 > {params.outputDir}/highDP.step1.indv
             cat {input.popmap} | grep -f  {params.outputDir}/highDP.step1.indv > {output.vcf}
             mkdir -p {params.outputDir}/finalVCF/
             """


if config["mode"]== "Reference":
     rule filter:
        input:
             popmap=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
             vcf=expand("{path}/refOut/populations.vcf.gz",path=config["outputDir"])
        output:
             vcf=expand("{path}/filters/{params}/populations.snps.vcf",path=config["outputDir"],params=paramspace.wildcard_pattern),
             sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern),
             FstSumstats=expand("{path}/filters/{params}/populations.fst_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
        params:
             outputDir=expand("{path}",path=config["outputDir"]),
             parDir=expand("{path}/filters/{params}",path=config["outputDir"],params=paramspace.wildcard_pattern),
             maf=str(paramspace.wildcard_pattern).split("/")[1].split("~")[1]
        threads:
             4
        resources:
            mem_mb=lambda wc, input: (0.5 * input.vcf.size_mb),
            runtime= 30,
            cpus_per_task= 4

        conda:
             "env/stacks.yaml"
        shell:
             """
             populations -M {params.outputDir}/stacksFiles/popmapFiltered.tsv -V {input.vcf} -R {wildcards.max_missing} --min-maf {params.maf} --vcf -O {params.parDir} --threads {threads} --fstats
             mv {params.parDir}/populations.p.snps.vcf {params.parDir}/populations.snps.vcf
             mv {params.parDir}/populations.p.sumstats_summary.tsv {params.parDir}/populations.sumstats_summary.tsv
             """     


rule makePCAData:
     input:
        gds=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.wildcard_pattern),
        popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv",path=config["outputDir"]),
        popmap=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
     output:
        pcaData=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
     conda:
        "env/R.yaml"
     resources:
        mem_mb=lambda wc, input: (0.5 * input.gds.size_mb),
        runtime= 10,
        cpus_per_task= 1     
     params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
     shell:
        "Rscript src/filterAndFigures/pca.R {input.popmapSNPFilter} {input.gds} {input.popmap} {output.pcaData}"

rule combinePCAData:
     input:
        pcaData=expand("{path}/filters/{params}/pcaPlot.tsv",path=config["outputDir"],params=paramspace.instance_patterns)
     output:
        pcaDataAll=expand("{path}/filters/pcaAll.tsv",path=config["outputDir"])
     resources:
        mem_mb= 100,
        runtime= 10,
        cpus_per_task= 1
     
     shell:
        "cat <(cat {input.pcaData} | head -n 1) <(cat {input.pcaData} | grep -v sample.id)  > {output.pcaDataAll}"


rule makeTreeData:
     input:
        gds=expand("{path}/filters/{params}/populations.snps.gds",path=config["outputDir"],params=paramspace.wildcard_pattern),
        popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv",path=config["outputDir"]),
        popmap=expand("{path}/stacksFiles/popmapFiltered.tsv",path=config["outputDir"]),
     output:
        treeLabels=expand("{path}/filters/{params}/treeLabels.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern),
        treeSegments=expand("{path}/filters/{params}/treeSegments.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
     conda:
        "env/R.yaml"
     resources:
        mem_mb=lambda wc, input: (0.5 * input.gds.size_mb),
        runtime= 10,
        cpus_per_task= 1
     params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
     shell:
        "Rscript src/filterAndFigures/cluster.R {input.popmapSNPFilter} {input.gds} {input.popmap} {output.treeLabels} {output.treeSegments}"

rule combineTreeData:
     input:
        treeLabels=expand("{path}/filters/{params}/treeLabels.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
        treeSegments=expand("{path}/filters/{params}/treeSegments.tsv",path=config["outputDir"],params=paramspace.instance_patterns)
     output:
        treeLabels=expand("{path}/filters/treeLabelsAll.tsv",path=config["outputDir"]),
        treeSegments=expand("{path}/filters/treeSegmentsAll.tsv",path=config["outputDir"])
     resources:
        mem_mb= 100,
        runtime= 10,
        cpus_per_task= 1

     shell:
        """
        cat <(cat {input.treeLabels} | head -n 1) <(cat {input.treeLabels} | grep -v max_missing)  > {output.treeLabels}
        cat <(cat {input.treeSegments} | head -n 1) <(cat {input.treeSegments} | grep -v max_missing)  > {output.treeSegments}
        """


rule makePopData:
     input:
        sumstats=expand("{path}/filters/{params}/populations.sumstats_summary.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
     output:
        popStats=expand("{path}/filters/{params}/popStats.tsv",path=config["outputDir"],params=paramspace.wildcard_pattern)
     params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
     resources:
        mem_mb= 100,
        runtime= 10,
        cpus_per_task= 1
     shell:
        """
        paste <(cat {input.sumstats} | grep -v positions | sed '/Poly/q' | grep -v "Poly" | cut -f1,2,9,15,21,24) <(cat {input.sumstats} | sed -n '/Poly/,$p' | cut -f 3,4,5) | sed 's/# Pop ID/pop/' > {output.popStats}
        """

rule combinePopData:
     input:
        popStats=expand("{path}/filters/{params}/popStats.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
        popmapSNPFilter=expand("{path}/stacksFiles/SNPFilterPopMap.tsv",path=config["outputDir"]),
     output:
        popStats=expand("{path}/filters/popStatsAll.tsv",path=config["outputDir"])
     conda:
        "env/R.yaml"
     params:
        outputDir=expand("{path}/filters/",path=config["outputDir"]),
     resources:
        mem_mb= 100,
        runtime= 10,
        cpus_per_task= 1

     shell:
        """
        Rscript src/filterAndFigures/popData.R {output.popStats} {input.popmapSNPFilter} {input.popStats} 
        """

rule combineFSTData:
     input:
        FstSumstats=expand("{path}/filters/{params}/populations.fst_summary.tsv",path=config["outputDir"],params=paramspace.instance_patterns),
        barcodes=expand("{path}/{bar}", path=config["inputDir"], bar=config["barcodeFile"])
     output:
        fstStats=expand("{path}/filters/fstStatsAll.tsv",path=config["outputDir"])
     conda:
        "env/R.yaml"
     resources:
        mem_mb= 100,
        runtime= 10,
        cpus_per_task= 1
     shell:
        """
        Rscript src/filterAndFigures/fst.R {output.fstStats} {input.barcodes} {input.FstSumstats} 
        """

