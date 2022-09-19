# scWGBS processing with Biscuit and Nextflow

This is a work-in-progress meant to record my current thoughts/process for analyzing (sc)WGBS data with Biscuit and related tools.

# Quick start:
You'll need the following to get started:
* Files:
	* Paired-end fastq files from scWGBS sequencing
	* Biscuit-indexed reference genome
	* Biscuit QC assets for your reference genome
	* A copy of the Biscuit QC script (at least for now...)
	* A 'target' bed file (regions covered in this file will be represented in the output R object)
* Tools:
	* Singularity or similar container engine ([Nextflow Docs](https://www.nextflow.io/docs/latest/container.html))
	* Nextflow

## Reference genome and Biscuit QC config:
Open `conf/nextflow.config` in your favorite text editor and fill in the following values in the `params` section:
* `bqc_path` : path to the Biscuit QC script
* `ref` and `bqc_assets`: if you plan to use the hg38 genome, fill in paths to the reference and Biscuit QC assets for `hg38`. Likewise, if you plan to use to use the mouse genome, fill in the paths for `mm39`. You can also add your own genomes by adding another key, e.g. 
```
	bqc_assets = [
		hg38: '',
		mm39: '',
		yourNewGenome: '/path/to/bqc/assets'
	]
```
## Cluster configuration
Depending on the resources available at your site, you may need to specify a different execution engine, job queue, or other variables. Place these in the `cluster.config` file (see [Nextflow docs](https://www.nextflow.io/docs/latest/config.html#scope-executor) for details).

## Run-specific configuration
Create a copy of the config file `example_profile.config`. Open it in your favorite text editor and update the following fields:
* `in_dir` : path to directory containing input fastq files
* `out_dir` : path where you want your output files (will be created if it doesn't exist)
* `target_bed` : BED file containing regions of the genome you want to analyze (only these regions will be reported in the final R dataframes/files)
* `genome` : Genome to align reads to. Must match one of the keys in the `ref` map mentioned above. E.g. `hg38`, `mm39`, or even `yourNewGenome` if you added a custom genome.
* `bsconv_filter` : reads with more than this fraction of _unconverted_ CpH bases will be discarded. The suggested default (0.1) is probably sensible.

Save your new config file, e.g. as `nathan.config`. 

## Run it!
```
nextflow run analysis.nf -c ./conf/nathan.config
```
Note: on at least some systems, the login or head node limits memory allocation to individual users, preventing the Nextflow JVM from starting. I get around this with the following alias:
```
alias nfr='NXF_OPTS="-Xmx500m" MALLOC_ARENA_MAX=4 nextflow run'
```
