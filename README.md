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

