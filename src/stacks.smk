
if config["mode"]=="Denovo":
    if len(df.index)<101:
        rule denovo_map:
            input: 
                samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
                samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
                popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
            params:
                M=config["M"],
                outputDir=expand("{path}/stacks/",path=config["outputDir"]),
                inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
            output:
                vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"])
            conda:
                "env/stacks.yaml"
            threads:workflow.cores
            resources:
	            mem_mb: 30000,
			    runtime: 24:00:00,
			    cpus_per_task: workflow.cores

            shell:
                "denovo_map.pl --samples {params.inputDir} --popmap {input.popmap} -T  {threads}  --paired     -o {params.outputDir} -n {params.M} -m {params.M} -X 'populations: --vcf'"

    if len(df.index) > 100: 
        rule subset_popmap:
            input:
                popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
            output:
                expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"])
            params:
                90
      		resources:
    	        mem_mb: 100,
	    		runtime: 10:00,
		    	cpus_per_task: 1
            shell:
                "paste <(shuf -n {params} {input.popmap} | cut -f1) <(yes opt | head -n {params}) > {output}"
        
        rule denovo_map:
            input:
                popmap_sub=expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"]),
                samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
                samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
                popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
            params:
                M=config["M"],
                outputDir=expand("{path}/stacks/",path=config["outputDir"]),
                inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
            output:
                vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"])
            conda:
                "env/stacks.yaml"
            threads:workflow.cores
            resources:
	            mem_mb: 30000,
			    runtime: 24:00:00,
			    cpus_per_task: workflow.cores
            shell:
                "denovo_map.pl --samples {params.inputDir} --popmap {input.popmap} -T  {threads} --catalog-popmap {input.popmap_sub} --paired      -o {params.outputDir} -n {params.M} -m {params.M} -X 'populations: --vcf'"