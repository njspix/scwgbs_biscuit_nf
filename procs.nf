process fastqc {
    input:
        path(fastq)
        val(tag)

    output:
        path("*.zip"), emit: zips
        path("*.html"), emit: htmls

    script:
    """
        fastqc --nogroup ${fastq}
    """
}

process multiqc {
    input:
        path(files)
        val(tag)

    output:
        path("*.html")

    script:
    """
        multiqc . --title "${tag} Report" --filename "${tag}_multiqc.html"
    """
}

process trim {
    input:
        tuple val(id), path(fastqs)

    output: 
        tuple val(id), path('*_trimmed_R{1,2}.fq.gz'), emit: trimmed_pairs // for subsequent mapping, etc
        path('*.gz'), emit: trimmed_fastqs // for fastqc
        path('*.txt'), emit: reports // for multiqc

    script:
    """
        cutadapt \
            --cores 4 \
            --error-rate 0.1 \
            --overlap 1 \
            --quality-cutoff 20 \
            --trim-n \
            --minimum-length 20 \
            -n 7 \
            -a AGATCGGAAGAGC \
            -A AGATCGGAAGAGC \
            -o ${id}_trimmed_R1.fq.gz -p ${id}_trimmed_R2.fq.gz \
            ${fastqs[0]} ${fastqs[1]} > ${id}_trimming_report.txt
    """
}

process align {
    input:
        tuple val(id), path(fastqs)
        path(ref)

    output:
        tuple val(id), path("*.bam")

    script:
    ref = params.ref.(params.genome)
    """
        biscuit align -@ 14 ${ref} ${fastqs[0]} ${fastqs[1]} | samblaster | samtools sort -O BAM -o ${id}.bam -
    """
}

process index {
    input:
        tuple val(id), path(bam)
        val(tag)

    output:
        tuple val(id), path('*.bam', includeInputs: true), path("*.bai")

    script:
    """
        samtools index -@ 4 ${bam}
    """
}

process bsconv {
    input:
        tuple val(id), path(bam), path(bai)
        path(ref)

    output:
        tuple val(id), path("*_bsconv-cph-filter.bam")

    script:
    ref = params.ref.(params.genome)
    """
        biscuit bsconv -f ${params.bsconv_filter} ${ref} ${bam} ${id}_bsconv-cph-filter.bam
    """
}

process stats {
    input:
        tuple val(id), path(bam), path(index)
        path(ref)

    output:
        path("*.txt")

    script:
    ref = params.ref.(params.genome)
    """
        samtools stats -@ 2 -r ${ref} ${bam} > ${id}_stats.txt
    """
}

process filter_bam {
    input:
        tuple val(id), path(bam), path(bai)
        path(bed_file)

    output:
        tuple val(id), path("*.bam"), path("*.bai")

    script:
    """
        samtools view -hb -L ${bed_file} ${bam} > ${id}_filtered.bam
        samtools index ${id}_filtered.bam
    """
}

process pileup {
	input:
	tuple val(id), path(bam), path(bai)
    path(ref)

	output:
	tuple val(id), path("*_pileup.vcf.gz")

	script:
    ref = params.ref.(params.genome)
	"""
		biscuit pileup -@ 4 ${ref} ${bam} | bgzip -f > ${id}_pileup.vcf.gz
	"""
}

process tabix_vcf {
    input:
        tuple val(id), path(vcf)  

    output:
        tuple val(id), path("*.vcf.gz", includeInputs: true), path("*.tbi")

    script:
    """
        tabix -p vcf ${vcf}
    """
}

process vcf2bed {

    input:
        tuple val(id), path(vcf), path(tbi)

    output:
        tuple val(id), path("*.bed.gz")

    script:
    """
        biscuit vcf2bed -k 1 -t cg -s ALL ${vcf} | bgzip > ${id}_pileup.bed.gz
    """
}

process tabix_bed {
    input:
        tuple val(id), path(bed)

    output:
        tuple val(id), path("*.bed.gz", includeInputs: true), path("*.tbi")

    script:
    """
        tabix -p bed ${bed}
    """
}

process mergecg {
    input:
        tuple val(id), path(bed), path(tbi)
        path(ref)

    output:
        tuple val(id), path("*.bed.gz")

    script:
    ref = params.ref.(params.genome)
    """
        biscuit mergecg ${params.ref} ${bed} | bgzip > ${id}_pileup_mergecg.bed.gz
    """
}

process biscuiteer {
    input:
        tuple val(id), path(bed), path(bed_tbi), path (vcf), path(vcf_tbi)

    output:
        tuple val(id), path("*.rds")

    script:
    """
        #!/usr/bin/env Rscript --vanilla 
        library(biscuiteer)

        bisc <- readBiscuit(BEDfile = "${bed}", VCFfile = "${vcf}", merged = T, genome = "${params.genome}", sparse = T)
        saveRDS(bisc, file = "${id}.rds", compress = F)
    """
}

process biscuit_qc {
	input:
	tuple val(id), path(bam), path(bai), path(vcf)
    path(ref)
    path(assets)

	output:
	path("*.txt")

	script:
    ref = params.ref.(params.genome)
    assets = params.bqc_assets.(params.genome)
	"""
	${params.bqc_path} -o . -v ${vcf} ${assets} ${ref} ${id} ${bam}
	"""
}

process bamtobed {
    input:
        tuple val(id), path(bam), path(bai)

    output:
        tuple val(id), path("*.bed")

    script:
    """
        bedtools bamtobed -i ${bam} | LC_ALL=C sort -k1,1 -k2,2n | cut -f1-3 > ${id}.bed
    """
}

process preseq {
    input:
        tuple val(id), path(bed)

    output:
        path("*.tsv")

    script:
    defects = ["", "-defects", "-defects -Q"][task.attempt - 1]
    // defects = "-defects -Q"
    """
        preseq gc_extrap \
            -o preseq_gc-extrap_${id}.tsv \
            -e 3E9 \
            -s 10E6 \
            -bed \
            ${defects} ${bed}
    """
}