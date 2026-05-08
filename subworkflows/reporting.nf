include {MOSDEPTH_PLOTDIST} from "../modules/local/mosdepth/plotdist/main.nf"
include {REPORT_INDIVIDUAL} from "../modules/local/report/main.nf"
include {REPORT_BOOK} from "../modules/local/report/main.nf"

workflow reporting {
    take:
        mosdepth_report_results   // [meta, global_dist_txt]
        haplotagged_samples       // [meta, bams, bais] — kept for API compatibility, not used here
        whatshap_stats_blocks     // [meta, whatshap_stats_txt, blocks_tsv]
        clustered_reads           // [meta, clustered_tsv, skew_tsv]
        cgi_bed                   // [meta_id, bed_file]

    main:
        // mosdepth coverage plot
        ch_plot_dist_script = Channel.fromPath("https://raw.githubusercontent.com/brentp/mosdepth/v0.3.6/scripts/plot-dist.py")
        (ch_global_dist_bysample, ch_plot_dist_script_rep) = mosdepth_report_results
            .combine(ch_plot_dist_script)
            .multiMap{it ->
                dists: tuple(it[0], it[1])
                scripts: it[2]
            }
        ch_mosdepth_dist_report = MOSDEPTH_PLOTDIST(ch_plot_dist_script_rep, ch_global_dist_bysample)

        // Build clustered-reads channel keyed by (id, sample).
        // clustered_reads has meta.sample as a plain string; recover the grouped sample list
        // from ch_mosdepth_dist_report (which already has meta.sample as a sorted list).
        ch_clustered_keyed = clustered_reads
            .map{meta, clustered, skew -> tuple(meta.id, clustered, skew)}
            .groupTuple()
            .join(ch_mosdepth_dist_report.map{meta, html -> tuple(meta.id, meta.sample)})
            .map{id, clustered_list, skew_list, sample -> tuple(id, sample, clustered_list, skew_list)}

        // Join all report inputs by (id, sample) using by:[0,1] so that the duplicate
        // sample element introduced by plain .join() is avoided.
        // Result tuple after all joins: [id, sample, mosdepth_html, whatshap_txt, whatshap_blocks, clustered_list, skew_list]
        ch_combined_qc_reports = ch_mosdepth_dist_report
            .map{meta, html -> tuple(meta.id, meta.sample, html)}
            .join(
                whatshap_stats_blocks.map{meta, txt, blocks -> tuple(meta.id, meta.sample, txt, blocks)},
                by: [0, 1]
            )
            .join(ch_clustered_keyed, by: [0, 1])
            .map{id, sample, mosdepth_html, whatshap_txt, whatshap_blocks, clustered_list, skew_list ->
                tuple([id: id, sample: sample], [mosdepth_html], whatshap_txt, whatshap_blocks, clustered_list, skew_list)
            }
            .combine(cgi_bed.map{it -> it[1]})
            .combine(channel.fromPath("${projectDir}/assets/report-templates/individual_report.qmd", checkIfExists: true))

        ch_reporting_files = REPORT_INDIVIDUAL(ch_combined_qc_reports)

        ch_book_template_files = channel.fromPath([
            "${projectDir}/assets/report-templates/_quarto_template.yml",
            "${projectDir}/assets/report-templates/index.qmd"
        ], checkIfExists: true).collect()

        book = REPORT_BOOK(
            ch_book_template_files,
            ch_reporting_files.qmds.collect(),
            ch_reporting_files.htmls.collect(),
            ch_reporting_files.whatshap_stats.collect(),
            ch_reporting_files.whatshap_blocks.collect(),
            ch_reporting_files.clustered_reads.collect(),
            ch_reporting_files.skew_tsv.collect(),
            cgi_bed.map{it -> it[1]}
        )

    emit:
        book
}
