if config["mode"]=="Denovo":
      if len(df.index)<1:
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
               threads:16
               resources:
                  mem_mb= 30000,
                  runtime= 24*60,
                  cpus_per_task= 16

               shell:
                  "denovo_map.pl --samples {params.inputDir} --popmap {input.popmap} -T  {threads}  --paired      -o {params.outputDir} -n {params.M} -m {params.M} -X 'populations: --vcf'"

      if len(df.index) > 10: 
         rule subset_popmap_denovo:
               input:
                  popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
               output:
                  expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"])
               params:
                  90
               resources:
                  mem_mb= 1000,
                  runtime= 5,
                  cpus_per_task= 1
               shell:
                  "paste <(shuf -n {params} {input.popmap} | cut -f1) <(yes opt | head -n {params}) | grep '[0-9]' > {output}"
         
         rule ustacks:
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
                  popmap_sub=expand("{path}/stacks/popmapSub.tsv",path=config["outputDir"])
            resources:
                  mem_mb= 30000,
                  runtime= 6*60,
                  cpus_per_task= 16
            conda:
                  "env/stacks.yaml"
            shell:
               """
               for sample in $(cut -f1 {input.popmap}) 
               do echo 
               ustacks --in-type gzfastq --file {params.inputDir}/$sample.1.fq.gz --out-path {params.outputDir} -m {params.M} -M {params.M} --name $sample --threads 16
               done
               cat {input.popmap_sub} > {output.popmap_sub}
               """

         rule cstacks:
            input:
                  popmap_sub=expand("{path}/stacks/popmapSub.tsv",path=config["outputDir"])
            params:
                  M=config["M"],
                  outputDir=expand("{path}/stacks/",path=config["outputDir"]),
            output:
                  catalog=expand("{path}/stacks/catalog.alleles.tsv.gz",path=config["outputDir"]),
                  cstackslog=expand("{path}/stacks/cstacks.log",path=config["outputDir"])
            resources:
                  mem_mb= 30000,
                  runtime= 6*60,
                  cpus_per_task= 16
            conda:
                  "env/stacks.yaml"
            shell:
               """
                  cstacks --threads 16 -n {params.M} -M {input.popmap_sub} -P {params.outputDir} > {output.cstackslog}
               """

         rule ssstacks:
            input:
                  popmap_sub=expand("{path}/stacks/popmapSub.tsv",path=config["outputDir"]),
                  catalog=expand("{path}/stacks/catalog.fa.gz",path=config["outputDir"]),
                  cstackslog=expand("{path}/stacks/cstacks.log",path=config["outputDir"])
            params:
                  M=config["M"],
                  outputDir=expand("{path}/stacks/",path=config["outputDir"]),
            output:
                  ssstackslog=expand("{path}/stacks/sstacks.log",path=config["outputDir"])
            resources:
                  mem_mb= 30000,
                  runtime= 6*60,
                  cpus_per_task= 16
            conda:
                  "env/stacks.yaml"
            shell:
               """
                  sstacks --threads 16  -M {input.popmap_sub} -P {params.outputDir} > {output.ssstackslog}
               """

         rule tsv2bam:
            input:
                  ssstackslog=expand("{path}/stacks/sstacks.log",path=config["outputDir"]),
                  popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"])
            params:
                  M=config["M"],
                  outputDir=expand("{path}/stacks/",path=config["outputDir"]),
                  inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
            output:
                  tsv2bamlog=expand("{path}/stacks/tsv2bam.log",path=config["outputDir"])
            resources:
                  mem_mb= 30000,
                  runtime= 1*60,
                  cpus_per_task= 16
            conda:
                  "env/stacks.yaml"
            shell:
               """
               tsv2bam -P {params.outputDir}  --threads 16 -M {input.popmap} -R {params.inputDir} > {output.tsv2bamlog}
               """
         rule gstacks:
            input:
                  tsv2bamlog=expand("{path}/stacks/tsv2bam.log",path=config["outputDir"]),
                  popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"])
            params:
                  M=config["M"],
                  outputDir=expand("{path}/stacks/",path=config["outputDir"]),
                  inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
            output:
                  gstackslog=expand("{path}/stacks/gstacks.log",path=config["outputDir"])
            resources:
                  mem_mb= 30000,
                  runtime= 3*60,
                  cpus_per_task= 16
            conda:
                  "env/stacks.yaml"
            shell:
               """
               gstacks -P {params.outputDir} --threads 16 -M {input.popmap} > {output.gstackslog}
               """
         rule population:
            input:
                  gstackslog=expand("{path}/stacks/gstacks.log",path=config["outputDir"]),
                  popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"])
            params:
                  M=config["M"],
                  outputDir=expand("{path}/stacks/",path=config["outputDir"]),
                  inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
            output:
                  vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"])
            resources:
                  mem_mb= 30000,
                  runtime= 1*60,
                  cpus_per_task= 16
            conda:
                  "env/stacks.yaml"
            shell:
               """
               populations -P {params.outputDir} --threads 16  --vcf -M {input.popmap}
               """


   #        rule denovo_map:
   #             input:
   #                popmap_sub=expand("{dir}/stacksTest/popmapSub.tsv",dir=config["outputDir"]),
   #                samplesR1=expand("{path}/demultiplex/samples/{samples}.1.fq.gz",path=config["outputDir"],samples=SAMPLES),
   #                samplesR2=expand("{path}/demultiplex/samples/{samples}.2.fq.gz",path=config["outputDir"],samples=SAMPLES),
   #                popmap=expand("{path}/stacksFiles/popmapFiltDemulti.tsv", path=config["outputDir"]),
   #             params:
   #                M=config["M"],
   #                outputDir=expand("{path}/stacks/",path=config["outputDir"]),
   #                inputDir=expand("{path}/demultiplex/samples/",path=config["outputDir"])
   #             output:
   #                vcf=expand("{path}/stacks/populations.snps.vcf",path=config["outputDir"])
   #             conda:
   #                "env/stacks.yaml"
   #             threads:16
   #             resources:
   #                mem_mb= 30000,
   #                runtime= 24*60,
   #                cpus_per_task= 16
   #             shell:
   #                "denovo_map.pl --samples {params.inputDir} --popmap {input.popmap} -T  {threads} --catalog-popmap {input.popmap_sub} --paired        -o {params.outputDir} -n {params.M} -m {params.M} -X 'populations: --vcf'"