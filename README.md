# EcologyGBSPipeline
Pipeline to analyse GBS data using the inhouse protocol of the radboud Ecology department and the Plant ecology and Nature conservation deparatment at Wageningen University.

This pipeline uses bwa_mem2 and freebayes for species for which a reference genome is available and stacks for species for which no reference exists.
It contains the following steps:

* Trimming (fastp)
* PCR duplicate removal (clone_filter stacks)
* Demultiplexing (process_radtags stacks)
* Reference
    * Mapping (bwa mem)
    * Variant calling (freebayes)
* Denovo
    * denovo_map.pl (stacks)
* Filtering
    * based on max-missing
    * based on minor allele frequency
* Analysis
    * PCA (SNPrelate)
    * Cluster (SNPrelate)
    * Population statistics (stacks)
    * Fst (stacks)
    * Isolation by distance (stacks)

## Prerequisites for running the pipeline

- A basic knowledge of Linux:
	- Knowledge, about how to work with files and directories (cd, ls, nano)
	- Being able to execute commands (git, conda, snakemake)
- Linux server:
	- e.g. Ubuntu 16.04
	- Sudo rights not necessary
	- Miniconda installed or
		- Download with `wget` https://docs.conda.io/en/latest/miniconda.html
		- Installing like described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), can be installed in user's home without sudo
			- `bash Miniconda3-latest-Linux-x86_64.sh`
- Adapters include:
	- Control nucleotide
	- Wobble (running without is also possible)
- Sequencing:
	- paired-end Illumina-sequencing
	- adjust output (number of reads) according to your genome size and expected number of fragments (e.g. based on an *in silico* digest)
- Readfiles:
	- un-demultiplexed but standard Illumina adapter trimmed (usually already done by sequencing agency)
- No reference genome required


# How to run the pipeline

In order to run the pipeline you will first need to download the pipeline using:

* ` git clone https://github.com/MaartenPostuma/EcologyGBSPipeline.git `

Enter the created directory
* `cd EcologyGBSPipeline`

Make a conda environment for snakemake if snakemake is not installed globally on the server. You do not need administrator rights to do this but conda has to be installed (see [Prerequisites for running the pipeline](#prerequisites-for-running-the-pipeline)).
* `conda create -n snakemake`
* `conda install -c bioconda snakemake=8.20.5`

Create the barcode file
* This needs to be done via a very specific format: where each column is tab separated. NOTE: Make sure that the column headers match the ones below other wise the pipeline will not run!
```
sample	pop metaPop	barcode1	barcode2	rawR1	rawR2
POP1_1	POP1	test1	AGGC	AACT	run1_R1.fq.gz	run1_R2.fq.gz
POP1_6	POP1	test1	AGGC	ACTA	run1_R1.fq.gz	run1_R2.fq.gz
POP2_1	POP2	test2	AGGC	AACT	run2_R1.fq.gz	run2_R2.fq.gz
POP3_1	POP3	test2	AGATGC	ATAC	run2_R1.fq.gz	run2_R2.fq.gz
```
* sample is the name of the sample, population is the population. metaPop is a classification of the populations for colouring in plots. barcode1/barcode2 are the barcodes and rawR1/rawR2 are the names of the files in which the raw reads are located.
* This barcode file needs to be saved in the same directory as the raw reads.

Rename the raw reads
* In order to run the pipeline the rawread files need to be renamed as the R1 reads need to end with _R1.fq.gz and the R2 reads need to end with _R2.fq.gz.
Everything befeore _R1.fq.gz will be considered the "sequencing run". So these need to be the same for R1 and R2.

Fill in the config.yaml

* open the config.yaml using `nano config.yaml` and edit so the location of the input/output files match.
* set the pipeline mode to Reference, Denovo or StacksTest. 
* If reference set the location and name of the reference genome.
* If you want the pipeline to run in denovo mode it is recommend to first determine the best parameter using stacksTest.
    * This mode will generate the plots in [paris et al. 2017](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12775) (See also the [stacks manual](https://catchenlab.life.illinois.edu/stacks/param_tut.php))
* Set the parameters for the denovo run.

Check if the pipeline works using `snakemake -n` to do run a dry run
Run the pipeline using `snakemake -j NumberOfJobs --use-conda` to run fully run it

This pipeline also supports slurm via the snakemake slurm functionality.
First you'll need to install the slurm executor plugin in the snakemake environment using: `conda install bioconda::snakemake-executor-plugin-slurm` 
Then you can run `snakemake -j NumberOfJobs --use-conda --executor slurm --latency-wait 60` to run each job as it's own separate slurm command. Resources have been assigned to each job. and should suffice in most cases.
