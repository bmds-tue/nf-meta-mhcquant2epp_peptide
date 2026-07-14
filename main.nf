nextflow.enable.dsl = 2

include { validateParameters; paramsHelp } from 'plugin/nf-schema'

params.help            = false
params.samplesheet     = null   // the SAME file mhcquant's --input uses (or a superset of it)
params.mhcquant_outdir = null
params.outdir          = null
params.sample_col      = "Sample"
params.condition_col   = "Condition"
params.alleles_col     = "Alleles"
params.mhc_class_col   = "Mhc_class"

// Delimiter is auto-detected from the samplesheet's extension, since real mhcquant
// samplesheets are commonly tab-separated despite splitCsv defaulting to comma.
def detectSep(path) {
    def name = path.toString().toLowerCase()
    return (name.endsWith('.tsv') || name.endsWith('.tab')) ? '\t' : ','
}

workflow mhcquant2epp {
    main:
    // --- 0. Fail fast if a configured column doesn't actually exist in the
    //        samplesheet header, rather than silently reading null for every row.
    def sheetPath = file(params.samplesheet, checkIfExists: true)
    def sheetSep  = detectSep(sheetPath)
    def sheetCols = sheetPath.readLines()[0].split(sheetSep, -1)*.trim()
    def requiredCols = [params.sample_col, params.condition_col, params.alleles_col, params.mhc_class_col].unique()
    def missingCols = requiredCols.findAll { !(it in sheetCols) }
    if (missingCols) {
        error "samplesheet '${params.samplesheet}' is missing required column(s): ${missingCols.join(', ')} (found: ${sheetCols.join(', ')})"
    }

    // --- 1. Parse samplesheet, deriving the Sample_Condition key ourselves
    //        rather than requiring it as a pre-existing column.
    parsed_ch = Channel
        .fromPath(sheetPath)
        .splitCsv(header: true, sep: sheetSep)
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

    // collectFile must always receive exactly one item -- a channel that emits zero
    // items (e.g. via a .filter{} that drops everything) never signals completion to
    // collectFile's underlying FileCollector, hanging the whole session indefinitely
    // even though every other channel has already finished. toSortedList() always
    // emits exactly one (possibly empty) list, so map it straight to file content --
    // an empty string when there are no conflicts -- instead of filtering to nothing.
    conflicts_ch
        .toSortedList()
        .map { rows -> rows.isEmpty() ? '' : (["sample_condition\tconflicting_alleles_mhc_class"] + rows).join('\n') + '\n' }
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
    if (params.help) {
        log.info paramsHelp(command: "nextflow run main.nf --samplesheet <path> --mhcquant_outdir <dir> --outdir <dir>")
        exit 0
    }
    validateParameters()

    mhcquant2epp()
}
