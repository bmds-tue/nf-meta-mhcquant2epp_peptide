nextflow.enable.dsl = 2

params.samplesheet     = null   // the SAME file mhcquant's --input uses (or a superset of it)
params.mhcquant_outdir = null
params.outdir          = null
params.sample_col      = "Sample"
params.condition_col   = "Condition"
params.alleles_col     = "alleles"
params.mhc_class_col   = "mhc_class"

workflow mhcquant2epp {
    main:
    // --- 1. Parse samplesheet, deriving the Sample_Condition key ourselves
    //        rather than requiring it as a pre-existing column.
    parsed_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def key = "${row[params.sample_col]}_${row[params.condition_col]}"
            tuple(key, [
                sample:    key,
                alleles:   row[params.alleles_col],
                mhc_class: row[params.mhc_class_col],
            ])
        }

    // --- 2. Report (but don't silently drop) rows that share a derived key
    //        but disagree on alleles/mhc_class -- a real data-entry error.
    conflicts_ch = parsed_ch
        .groupTuple()
        .filter { key, metas -> metas.unique().size() > 1 }
        .map { key, metas ->
            def variants = metas.unique().collect { "${it.alleles}/${it.mhc_class}" }
            "${key}\t${variants.join(' | ')}"
        }

    conflicts_ch
        .toSortedList()
        .map { rows -> rows.isEmpty() ? null : (["sample_condition\tconflicting_alleles_mhc_class"] + rows).join('\n') + '\n' }
        .filter { it != null }
        .collectFile(name: 'conflicts.tsv', storeDir: "${params.outdir}")

    // --- 3. Collapse replicate-level rows down to one row per Sample_Condition.
    samplesheet_ch = parsed_ch.unique { it[0] }

    // --- 4. Discover mhcquant's per-Sample_Condition output TSVs.
    files_ch = Channel
        .fromPath("${params.mhcquant_outdir}/*.tsv")
        .map { file ->
            def base      = file.baseName
            def condition = base.substring(base.lastIndexOf('_') + 1)
            def sample_id = base.substring(0, base.lastIndexOf('_'))
            tuple("${sample_id}_${condition}", file)
        }

    // --- 5. Join on the shared key, then merge into one meta map
    csv_ch = samplesheet_ch
        .join(files_ch, remainder: false)
        .map { key, meta, file ->
            meta + [filename: file]
        }

    // --- 6. Collect all rows, sorted by sample, into one CSV
    csv_ch
        .toSortedList { a, b -> a.sample <=> b.sample }
        .map { rows ->
            def lines = rows.collect { "${it.sample},${it.alleles},${it.mhc_class},${it.filename}" }
            (["sample,alleles,mhc_class,filename"] + lines).join('\n') + '\n'
        }
        .collectFile(name: 'samplesheet.csv', storeDir: "${params.outdir}")

    emit:
    samplesheet = csv_ch
}

workflow {
    mhcquant2epp()
}
