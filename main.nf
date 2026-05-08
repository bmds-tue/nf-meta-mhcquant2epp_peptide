params.samplesheet = null
params.mhcquant_outdir = null
params.join_key = "sample"

workflow mhcquant2epp {
    main:
    // --- 1. Parse samplesheet → [join_key_value, meta_map]
    samplesheet_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row[params.join_key], row) }

    // --- 2. Discover TSV files → [sample_id, condition, file]
    files_ch = Channel
        .fromPath("${params.mhcquant_outdir}/*.tsv")
        .map { file ->
            def base      = file.baseName
            def condition = base.substring(base.lastIndexOf('_') + 1)
            def sample_id = base.substring(0, base.lastIndexOf('_'))
            tuple("${sample_id}_${condition}", file)
        }

    // --- 3. Join on the shared key, then merge into one meta map
    csv_ch = samplesheet_ch
        .join(files_ch, remainder: false)
        .map { key, meta, file ->
            meta + [filename: file]
        }

    // --- 4. Collect all rows, sort by original index, write CSV
    csv_ch.collectFile(
        keepHeader: true,
        skip: 1,
        sort: true,
        storeDir: "${params.outdir}")
        { meta ->
            ["samplesheet.csv", "sample,alleles,mhc_calss,filename\n${meta.sample},${meta.alleles},${meta.mhc_class},${meta.filename}\n"]
        }
}

workflow  {
    mhcquant2epp()
}
