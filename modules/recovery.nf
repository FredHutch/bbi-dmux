save_recovery2 = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz-summary.txt/ + "-recovery_summary.txt"}
save_recovery = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz.txt.gz/ + "-recovery_table.txt.gz"}

process run_recovery {
    container "${params.container__python}"

    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery, pattern: "*.gz.txt.gz", mode: 'link', overwrite: true
    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery2, pattern: "*-summary.txt", mode: 'link', overwrite: true

    input:
        file(input)
        file(sample_sheet)
        file(run_parameters_file)
        file(rt_barcode_file)
        file(p7_barcode_file)
        file(p5_barcode_file)
        file(lig_barcode_file)

    output:
        path("*gz.txt.gz")
        path("*summary.txt"), emit: summaries

    script:
"""
set -euo pipefail

recovery_script.py --input_file <(zcat $input) --output_file ${input}.txt \
    --run_directory . \
    --sample_layout $sample_sheet_file \
    --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
    --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
    --p7_barcode_file $p7_barcode_file \
    --p5_barcode_file $p5_barcode_file \
    --lig_barcode_file $lig_barcode_file \
    --level $params.level \
    --rt_barcodes $rt_barcode_file
    pigz -p 1 *.fastq.gz.txt
"""
}

process sum_recovery {
    container "${params.container__python}"

    publishDir path: "${params.output_dir}/demux_dash/js/", pattern: "recovery_summary.js", mode: 'move', overwrite: true

    input:
        file(summary)

    output: 
        path "*summary.js"

    script:
"""
set -euo pipefail

echo "const log_data = {" > recovery_summary.js
for file in $summary
do
    filename=\$(basename \$file); 
    part=\${filename/Undetermined-L00/};
    lane=\${part/.fastq.gz-summary.txt/};   
    printf "\$lane : \\`" >> recovery_summary.js;
    cat \$file >> recovery_summary.js;
    printf "\\`," >> recovery_summary.js;
done
echo "}" >> recovery_summary.js   
"""
}
