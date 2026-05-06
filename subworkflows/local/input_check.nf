//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv //contains the samplesheet_valid.csv
        .splitCsv( header:true, sep:',' )
        .map{row -> tuple(row.individual, row.sample, row.modbam_5mCG)}
        .set{ch_sample}

    emit:
    ch_sample                                     // channel: [ val(individual), val(sample) path(modbam) ]
}

workflow INPUT_CHECK_PHASED {
    take:
    samplesheet // file: /path/to/samplesheet.csv with columns: individual,sample,haplotagged_bam

    main:
    samplesheet
        .splitCsv( header:true, sep:',' )
        .map{row -> tuple([id: row.individual, sample: row.sample], file(row.haplotagged_bam))}
        .set{ch_sample}

    emit:
    ch_sample                                     // channel: [ val(meta), path(haplotagged_bam) ]
}

