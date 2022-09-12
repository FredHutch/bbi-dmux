/*
** Notes:
**   o  hashed experiments can use many wells with a single RT
**      barcode so all 'samples' are demultiplexed into a single
**      fastq file, which is likely to become very large. Our
**      solution is to split the demultiplexed fastq file when
**      a single sample occupies more than
**      params.max_wells_per_sample wells. (The split-sample
**      fastq files are re-combined into a single file, which
**      is used to make a single cell_data_set.)
*/
nextflow.enable.dsl=2

DEFAULT = "default"
default_rt2_barcode_file = "$baseDir/bin/barcode_files/rt2.txt"
default_rt3_barcode_file = "$baseDir/bin/barcode_files/rt.txt"
default_p5_barcode_file = "$baseDir/bin/barcode_files/p5.txt"
default_p7_barcode_file = "$baseDir/bin/barcode_files/p7.txt"
default_lig_barcode_file = "$baseDir/bin/barcode_files/ligation.txt"
default_star_file = "$baseDir/bin/star_file.txt"

// Import sub-workflows
include {
    generate_sheets
    check_sample_sheet
    make_sample_sheet
} from './modules/samplesheets'

include {
    bcl2fastq
} from './modules/bcl2fastq'

include {
    seg_sample_fastqs
} from './modules/fastq_splitting'

include {
    demux_dash
} from './modules/dashboard'

include {
    run_recovery
    sum_recovery
} from './modules/recovery'


def helpMessage() {
    log.info"""
        BBI sci-RNA-seq Demultiplexer
        --------------------------------

        For reproducibility, please specify all parameters to a config file
        by specifying -c CONFIG_FILE.config.

        Usage:
            nextflow run bbi-dmux -c CONFIG_FILE

        Help:
            --help                                     Show this message and exit.

        Required parameters (specify in your config file):
            params.run_dir = RUN_DIRECTORY             Path to the sequencer output.
            params.output_dir OUTPUT DIRECTORY         Output directory.
            params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.
            params.level = 3                           Level of run - either 2 or 3.

        Required parameters (one of the pairs below is required - p7_wells and p5_wells or p7_rows and p5_cols or pcr_index_pair_file or mult-exp):
            params.p7_wells = "A1 B1 C1"               Alternative to p7_rows and p5_cols - specify specific PCR wells instead of full rows/columns. Must match order of params.p5_wells.
            params.p5_wells = "A1 A2 A3"               Alternative to p7_rows and p5_cols - specify specific PCR wells instead of full rows/columns. Must match order of params.p7_wells.
            params.p7_rows = "A B C"                   The PCR rows used - must match order of params.p5_cols.
            params.p5_cols = "1 2 3"                   The PCR columns used - must match order of params.p7_rows.
            params.pcr_index_pair_file                 <file_name> The path to a PCR rxn primer pair file.'
            params.multi_exp = "see config"            The PCR columns used for each experiment in map format - see example.config.


        Optional parameters (specify in your config file):
            params.rt_barcode_file = "default"         The path to a custom RT barcode file. If "default", default BBI barcodes will be used.
            params.p7_barcode_file = "default"         The path to a custom p7 barcode file. If "default", default BBI barcodes will be used.
            params.p5_barcode_file = "default"         The path to a custom p5 barcode file. If "default", default BBI barcodes will be used.
            params.lig_barcode_file = "default"        The path to a custom ligation barcode file. If "default", default BBI barcodes will be used.
            params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.
            process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.
            process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted.
            params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.
            params.bcl_max_mem = 40                    The maximum number of GB of RAM to assign for bcl2fastq
            params.max_wells_per_sample = 20           The maximum number of wells per sample - if a sample is in more wells, the fastqs will be split then reassembled.
            --run_recovery true                        Add this to run the recovery script AFTER running the normal pipeline.
            --generate_samplesheets input_csv          Add this to generate the necessary samplesheet from the BBI universal input sheet.

        Leave issue reports at "https://github.com/bbi-lab/bbi-dmux/issues".
    """.stripIndent()
}

workflow {
    // Show help message if the user specifies help flag or if required params are not provided
    if (params.help || !params.run_dir || !params.output_dir || !params.sample_sheet ) {
        helpMessage()
        exit 1
    }

    // Set up input files and directories
    sample_sheet_file = file(params.sample_sheet)
    run_dir = Channel.fromPath(params.run_dir)
    run_parameters_file = file("${params.run_dir}/[Rr]unParameters.xml*")
    if (!run_parameters_file) {
        exit 1, "Run parameters file not found in ${params.run_dir}"
    }

    if (params.rt_barcode_file == DEFAULT) {
        if (params.level == 2) {
            rt_barcode_file = file(default_rt2_barcode_file)
        }
        if (params.level == 3) {
            rt_barcode_file = file(default_rt3_barcode_file)
        }
    } else {
        rt_barcode_file = file(params.rt_barcode_file)
    }

    star_file = params.star_file == DEFAULT ? file(default_star_file) : file(params.star_file)
    lig_barcode_file = params.lig_barcode_file == DEFAULT ? file(default_lig_barcode_file) : file(params.lig_barcode_file)
    p5_barcode_file = params.p5_barcode_file == DEFAULT ? file(default_p5_barcode_file) : file(params.p5_barcode_file)
    p7_barcode_file = params.p7_barcode_file == DEFAULT ? file(default_p7_barcode_file) : file(params.p7_barcode_file)

    if (params.generate_samplesheets) {
        bbi_universal_sheet_file = Channel.fromPath(params.generate_sample_sheets)
        generate_sheets(bbi_universal_sheet_file)

        sample_sheet_file = generate_sheets.out.sample_sheet
    } else {
        check_sample_sheet(sample_sheet_file, star_file, rt_barcode_file, params.level, params.max_wells_per_sample)
        
        sample_sheet_file = check_sample_sheet.out.good_sample_sheet
    }

    if (!params.run_recovery) {
        make_sample_sheet(run_parameters_file, sample_sheet_file)

        bcl_sample_sheet = make_sample_sheet.out.bcl_sample_sheet
    }

    max_cores_bcl = Math.min(16, params.max_cores)

    bcl2fastq(
        run_dir, 
        bcl_sample_sheet, 
        max_cores_bcl, 
        params.bcl_max_mem,
        params.bcl2fastq_barcode_mismatches,
        params.minimum_read_length_after_trim
    )

    // fastqs is a tuple of R1, R2 files separated by each lane
    fastqs = bcl2fastq.out.fastqs.transpose()

    seg_sample_fastqs(
        fastqs,
        run_parameters_file,
        sample_sheet_file,
        rt_barcode_file,
        p5_barcode_file,
        p7_barcode_file,
        lig_barcode_file
    )

    csv_stats = seg_sample_fastqs.out.csv_stats
    json_stats = seg_sample_fastqs.out.json_stats.flatten()

    // Generate dashboard
    skeleton_dash = file("$baseDir/bin/skeleton_dash")
    demux_dash(csv_stats.collect(), json_stats.collect(), sample_sheet_file, skeleton_dash)

    if (params.run_recovery) {
        undetermined_fastqs = Channel.fromPath("${params.demux_out}/Undetermined*")

        run_recovery(
            undetermined_fastqs,
            sample_sheet,
            run_parameters_file,
            rt_barcode_file,
            p5_barcode_file,
            p7_barcode_file,
            lig_barcode_file
        )

        summary = run_recovery.out.summaries.collect()

        sum_recovery(summary)
    }
}

workflow.onComplete {
    println "bbi-dmux completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
