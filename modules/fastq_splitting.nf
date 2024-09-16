process seg_sample_fastqs {
    container "${params.container__biopython}"

    publishDir path: "${params.output_dir}/", pattern: "demux_out/*fastq.gz", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}/stats", pattern: "demux_out/*.csv", saveAs: { file(it).getName() }, mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}/stats", pattern: "demux_out/*.json", saveAs: { file(it).getName() }, mode: 'copy', overwrite: true

    input:
        tuple file(R1), file(R2)
        file(run_parameters_file)
        file(sample_sheet_file)
        file(rt_barcode_file)
        file(p5_barcode_file)
        file(p7_barcode_file)
        file(lig_barcode_file)

    output:
        path "demux_out/*", emit: seg_output
        path "demux_out/*.fastq.gz", emit: samp_fastqs_check
        path "demux_out/*.stats.json", emit: json_stats
        path "demux_out/*.csv", emit: csv_stats
    script:
"""
set -euo pipefail

mkdir demux_out

make_sample_fastqs.py --run_directory . \
    --read1 <(zcat $R1) --read2 <(zcat $R2) \
    --file_name $R1 --sample_layout $sample_sheet_file \
    --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
    --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
    --pcr_index_pair_file 0 \
    --rt_barcode_file $rt_barcode_file \
    --p5_barcode_file $p5_barcode_file \
    --p7_barcode_file $p7_barcode_file \
    --lig_barcode_file $lig_barcode_file \
    --multi_exp "$params.multi_exp" \
    --output_dir ./demux_out --level $params.level

pigz -p 8 demux_out/*.fastq    
"""
}
