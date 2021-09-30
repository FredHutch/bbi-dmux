
process generate_sheets {
    publishDir path: "${params.output_dir}", pattern: "SampleSheet.csv", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}", pattern: "SampleMap.csv", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}", pattern: "GarnettSheet.csv", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}/sample_id_maps", pattern: "*_SampleIDMap.csv", mode: 'copy', overwrite: true

    input:
        file( sample_sheet )
    output:
        path ("*Sheet.csv"), emit: sample_sheet
        path ("SampleMap.csv") optional true
        path ("*_SampleIDMap.csv") optional true

"""/bin/bash
set -Eeuo pipefail

generate_sample_sheets.py $sample_sheet
"""

}

process check_sample_sheet {
    input:
        file( sample_sheet )
        file( star_file )
        file( rt_barcode_file )
        val( level )
        val( max_wells_per_sample )

    output:
        path ("*.csv" ), emit: good_sample_sheet

"""/bin/bash
set -Eeuo pipefail

echo "test"

check_sample_sheet.py \
    --sample_sheet $sample_sheet \
    --star_file $star_file \
    --level $level \
    --rt_barcode_file $rt_barcode_file \
    --max_wells_per_samp $max_wells_per_sample    
"""

}

/**
 * Generates sample sheet used for bcl2fastq
 */
process make_sample_sheet {
    cache 'lenient'

    input:
        file( run_parameters_file )
        file( good_sample_sheet )

    output:
        path( "SampleSheet.csv" ), emit: bcl_sample_sheet

"""/bin/bash
set -Eeuo pipefail

make_sample_sheet.py --run_directory .

"""    
}
