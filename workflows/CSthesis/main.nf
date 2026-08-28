nextflow.enable.dsl = 2

/*
 * SET UP PARAMETERS
 */
params.metadata_url   = "metadata.tsv"
params.classifier_url = "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
params.reads_dir      = null

params.outdir         = "${launchDir}/qiime2_results"
params.barcode_col    = "barcode-sequence"

// Primer sequences (adjust to your actual primers)
params.primer_f       = "YCAGCMGCCGCGGTAA"
params.primer_r       = "TACNVGGGTWTCTAAT"
//mockdata: YCAGCMGCCGCGGTAA, TACNVGGGTWTCTAAT
//colombian: CCTACGGGNGGCWGCAG, GACTACHVGGGTATCTAATCC

// dada2 trimming parameters
params.trim_left      = 0

params.trunc_len_f    = 200
params.trunc_len_r    = 180
// mockdata: 200/180
// colombian: 283/229

// DADA2 quality filtering
params.max_ee              = 2

// Diversity analysis
params.sampling_depth = 5167
// mockdata: 5167
// colombian: 27892


/*
 * RENAME FASTQ FILES TO QIIME2 CASAVA FORMAT
 */
process RENAME_READS {

    publishDir "${params.outdir}/renamed_reads", mode: 'copy'

    input:
    path reads_dir

    output:
    path "*.fastq.gz"

    script:
    """
    mkdir renamed

    for r1 in ${reads_dir}/*_1.fastq*; do

        sample=\$(basename \$r1)
        sample=\${sample%%_1.fastq*}

        # locate R2
        if [[ -f ${reads_dir}/\${sample}_2.fastq.gz ]]; then
            r2=${reads_dir}/\${sample}_2.fastq.gz
        else
            r2=${reads_dir}/\${sample}_2.fastq
        fi

        #
        # R1
        #
        if [[ \$r1 == *.gz ]]; then
            cp \$r1 renamed/\${sample}_S1_L001_R1_001.fastq.gz
        else
            gzip -c \$r1 > renamed/\${sample}_S1_L001_R1_001.fastq.gz
        fi

        #
        # R2
        #
        if [[ \$r2 == *.gz ]]; then
            cp \$r2 renamed/\${sample}_S1_L001_R2_001.fastq.gz
        else
            gzip -c \$r2 > renamed/\${sample}_S1_L001_R2_001.fastq.gz
        fi

    done

    mv renamed/* .
    """
}


/*
 * IMPORT PAIRED-END READS INTO QIIME2
 */
process IMPORT_PAIRED {

    publishDir "${params.outdir}/imported", mode: 'copy'

    input:
    path reads

    output:
    path "paired-end-demux.qza"

    script:
    """
    mkdir input

    cp ${reads} input/

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path input \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt \
        --output-path paired-end-demux.qza
    """
}

//process DEMUX {
//    publishDir "${params.outdir}/demux", mode: 'copy'
//
//    input:
//    path emp_seqs
//    path metadata
//
//    output:
//    path "demux.qza", emit: qza
//    path "demux-details.qza", emit: details
//
//    script:
//    """
//    qiime demux emp-paired \
//        --i-seqs ${emp_seqs} \
//        --m-barcodes-file ${metadata} \
//        --m-barcodes-column ${params.barcode_col} \
//        --o-per-sample-sequences demux.qza \
//        --o-error-correction-details demux-details.qza \
//        --p-chimera-method none
//    """
//}

process VISUALIZE_DEMUX {
    publishDir "${params.outdir}/demux", mode: 'copy'

    input:
    path demux_qza

    output:
    path "demux.qzv"

    script:
    """
    qiime demux summarize \
        --i-data ${demux_qza} \
        --o-visualization demux.qzv
    """
}

process CUTADAPT {
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path demux_qza

    output:
    path "trimmed-seqs.qza", emit: trimmed_seqs

    script:
    """
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${demux_qza} \
        --p-front-f ${params.primer_f} \
        --p-front-r ${params.primer_r} \
        --p-error-rate 0.1 \
        --p-discard-untrimmed \
        --o-trimmed-sequences trimmed-seqs.qza \
        --verbose
    """
}

process DADA2 {
    publishDir "${params.outdir}/dada2", mode: 'copy'

    input:
    path demux_qza

    output:
    path "table.qza",            emit: table
    path "rep-seqs.qza",         emit: rep_seqs
    path "denoising-stats.qza",  emit: stats

    script:
    """
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${demux_qza} \
        --p-trim-left-f ${params.trim_left} \
        --p-trim-left-r ${params.trim_left} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --p-max-ee-f ${params.max_ee} \
        --p-max-ee-r ${params.max_ee} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --p-n-threads ${task.cpus}
    """
}

process CLASSIFY_TAXONOMY {
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path rep_seqs
    path classifier

    output:
    path "taxonomy.qza"

    script:
    """
    qiime feature-classifier classify-sklearn \
        --i-classifier ${classifier} \
        --i-reads ${rep_seqs} \
        --o-classification taxonomy.qza \
        --p-n-jobs ${task.cpus}
    """
}

process BARPLOT {
    publishDir "${params.outdir}/visualization", mode: 'copy'

    input:
    path table
    path taxonomy
    path metadata

    output:
    path "taxa-bar-plots.qzv"

    script:
    """
    qiime taxa barplot \
        --i-table ${table} \
        --i-taxonomy ${taxonomy} \
        --m-metadata-file ${metadata} \
        --o-visualization taxa-bar-plots.qzv
    """
}

process PHYLOGENY {
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path rep_seqs

    output:
    path "rooted-tree.qza", emit: rooted_tree
    path "unrooted-tree.qza", emit: unrooted_tree

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${rep_seqs} \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza \
        --p-n-threads ${task.cpus}
    """
}

process CORE_DIVERSITY {
    publishDir "${params.outdir}/diversity", mode: 'copy'

    input:
    path rooted_tree
    path table
    path metadata

    output:
    path "core-metrics-results", emit: results

    script:
    """
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny ${rooted_tree} \
        --i-table ${table} \
        --p-sampling-depth ${params.sampling_depth} \
        --m-metadata-file ${metadata} \
        --output-dir core-metrics-results
    """
}

process DIVERSITY_VISUALIZATIONS {

    publishDir "${params.outdir}/diversity_visualizations", mode: 'copy'

    input:
    path core_metrics_dir
    path metadata

    output:
    path "*.qzv", optional: true

    script:
    """
    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${core_metrics_dir}/shannon_vector.qza \
        --m-metadata-file ${metadata} \
        --o-visualization shannon-group-significance.qzv || true

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${core_metrics_dir}/observed_features_vector.qza \
        --m-metadata-file ${metadata} \
        --o-visualization observed-features-group-significance.qzv || true

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${core_metrics_dir}/evenness_vector.qza \
        --m-metadata-file ${metadata} \
        --o-visualization evenness-group-significance.qzv || true

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${core_metrics_dir}/faith_pd_vector.qza \
        --m-metadata-file ${metadata} \
        --o-visualization faith-pd-group-significance.qzv || true

    qiime emperor plot \
        --i-pcoa ${core_metrics_dir}/jaccard_pcoa_results.qza \
        --m-metadata-file ${metadata} \
        --o-visualization jaccard-emperor.qzv || true

    qiime emperor plot \
        --i-pcoa ${core_metrics_dir}/bray_curtis_pcoa_results.qza \
        --m-metadata-file ${metadata} \
        --o-visualization bray-curtis-emperor.qzv || true

    qiime emperor plot \
        --i-pcoa ${core_metrics_dir}/unweighted_unifrac_pcoa_results.qza \
        --m-metadata-file ${metadata} \
        --o-visualization unweighted-unifrac-emperor.qzv || true

    qiime emperor plot \
        --i-pcoa ${core_metrics_dir}/weighted_unifrac_pcoa_results.qza \
        --m-metadata-file ${metadata} \
        --o-visualization weighted-unifrac-emperor.qzv || true
    """
}

process CHAO1_ALPHA_DIVERSITY {

    publishDir "${params.outdir}/alpha_diversity", mode: 'copy'

    input:
    path table
    path rooted_tree
    path metadata

    output:
    path "faith-pd-vector.qza", emit: faith_pd
    path "chao1-vector.qza", emit: vector
    path "*.qzv", optional: true, emit: viz

    script:
    """
    qiime diversity alpha-phylogenetic \
        --i-table ${table} \
        --i-phylogeny ${rooted_tree} \
        --p-metric faith_pd \
        --o-alpha-diversity faith-pd-vector.qza

    qiime diversity alpha \
        --i-table ${table} \
        --p-metric chao1 \
        --o-alpha-diversity chao1-vector.qza

    qiime diversity alpha-group-significance \
        --i-alpha-diversity chao1-vector.qza \
        --m-metadata-file ${metadata} \
        --o-visualization chao1-group-significance.qzv || true

    qiime diversity alpha-group-significance \
        --i-alpha-diversity faith-pd-vector.qza \
        --m-metadata-file ${metadata} \
        --o-visualization faith-pd-group-significance.qzv || true
    """
}

process DADA2_STATS_VIZ {
    publishDir "${params.outdir}/dada2", mode: 'copy'

    input:
    path stats

    output:
    path "denoising-stats.qzv"

    script:
    """
    qiime metadata tabulate \
        --m-input-file ${stats} \
        --o-visualization denoising-stats.qzv
    """
}

process NORMALIZE_TABLE {
    publishDir "${params.outdir}/normalized_table", mode: 'copy'

    input:
    path table

    output:
    path "normalized-table.qza", emit: normalized_table

    script:
    """
    qiime feature-table relative-frequency \
        --i-table ${table} \
        --o-relative-frequency-table normalized-table.qza
    """
}

/*
 * WORKFLOW CONTROL
 */
workflow {

    ch_metadata   = file(params.metadata_url)
    ch_classifier = file(params.classifier_url)
    ch_reads_dir  = file(params.reads_dir)

    RENAME_READS(ch_reads_dir)

    IMPORT_PAIRED(RENAME_READS.out)

    VISUALIZE_DEMUX(IMPORT_PAIRED.out)

    CUTADAPT(IMPORT_PAIRED.out)

    DADA2(CUTADAPT.out.trimmed_seqs)

    DADA2_STATS_VIZ(DADA2.out.stats)

    PHYLOGENY(DADA2.out.rep_seqs)

    CORE_DIVERSITY(
        PHYLOGENY.out.rooted_tree,
        DADA2.out.table,
        ch_metadata
    )

    DIVERSITY_VISUALIZATIONS(
        CORE_DIVERSITY.out.results,
        ch_metadata
    )

    CHAO1_ALPHA_DIVERSITY(
        DADA2.out.table,
        PHYLOGENY.out.rooted_tree,
        ch_metadata
    )

    CLASSIFY_TAXONOMY(
        DADA2.out.rep_seqs,
        ch_classifier
    )

    BARPLOT(
        DADA2.out.table,
        CLASSIFY_TAXONOMY.out,
        ch_metadata
    )

    NORMALIZE_TABLE(
        DADA2.out.table
    )
}

