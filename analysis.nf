include { 
	fastqc; 
	multiqc;
	trim;
	align;
	index;
	stats; 
	bsconv;
	index as index_bsconv;
	pileup;
	biscuit_qc;
	bamtobed;
	preseq;
	filter_bam;
	tabix_vcf; 
    vcf2bed; 
    tabix_bed as tabix_b1;
    tabix_bed as tabix_b2;
    mergecg;
    biscuiteer;
} from "$projectDir/procs.nf"

// ------------ Subworkflows ------------- //
workflow initial_qc {
	take: fastqs

	main:
		fastqc(fastqs, "pre-trimming")
		multiqc(fastqc.out.zips.collect(), "pre-trimming")
}

workflow trim_and_qc {
	take: paired_fastqs

	main:
		trim(paired_fastqs)
		fastqc(trim.out.trimmed_fastqs, "trimmed")

		trim.out.reports
			.mix(fastqc.out.zips)
			.collect()
			.set{ multiqc_files }
		
		multiqc(multiqc_files, "trimmed")

	emit: trim.out.trimmed_pairs
}

workflow align_bsconv {
    take: trimmed_pairs

    main:
	reference = Channel.fromPath(params.ref.(params.genome)).first()
    align( trimmed_pairs, reference )
	index(align.out, "raw") 
	bsconv(index.out, reference)
	index_bsconv(bsconv.out, "bsconv-filter")

    emit: index_bsconv.out
}

workflow align_qc {
	take: bams_with_index
	
	main:
	reference = Channel.fromPath(params.ref.(params.genome)).first()
	assets = Channel.fromPath(params.bqc_assets.(params.genome)).first()

	pileup(bams_with_index, reference)
	biscuit_qc(bams_with_index.join(pileup.out, by: 0), reference, assets)
	bamtobed(bams_with_index) | preseq
	stats(bams_with_index, reference)

	biscuit_qc.out
		.mix(preseq.out)
		.mix(stats.out)
		.collect()
		.set{final_multiqc_files}
	
	multiqc(final_multiqc_files, "post-align")
}

workflow bam_to_r {
	take: 
		bams_with_index
		target_bed

	main:
	reference = Channel.fromPath(params.ref.(params.genome)).first()

	filter_bam(bams_with_index, params.target_bed)
	pileup(filter_bam.out, reference) | tabix_vcf | vcf2bed | tabix_b1
	mergecg(tabix_b1.out, reference) | tabix_b2
	biscuiteer(tabix_b2.out.join(tabix_vcf.out))

}
// --------------------------------------- //

fastqs = Channel.fromPath(params.in_dir+"/**.f{astq.gz,q.gz}")
paired_fastqs = Channel.fromFilePairs(params.in_dir+"/**_{1,2}*.f{astq.gz,q.gz}")
	.mix(Channel.fromFilePairs(params.in_dir+"/**_R{1,2}*.f{astq.gz,q.gz}"))

workflow {
	initial_qc(fastqs)
	trim_and_qc(paired_fastqs) | align_bsconv | align_qc 
	bam_to_r(align_bsconv.out, params.target_bed)
}

