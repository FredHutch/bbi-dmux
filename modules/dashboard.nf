// Creates dashboard of results
process prep {
    container "${params.container__rscript}"

    input:
        file(demux_stats)
        file(jsons)
        file(sample_sheet)
        file(skeleton_dash)
    output:
        path(demux_dash)

    script:
        out_dir_str = params.output_dir.replaceAll("/\\z", "");
        project_name = out_dir_str.substring(out_dir_str.lastIndexOf("/")+1);

"""
set -euo pipefail
mkdir demux_dash
cp -R $skeleton_dash/* demux_dash/
generate_html.R \
    "." --p7_rows "$params.p7_rows" --p5_cols "$params.p5_cols" --p7_wells "$params.p7_wells" \
    --p5_wells "$params.p5_wells" --level "$params.level" --project_name "${project_name}" \
    --sample_sheet "$sample_sheet"
"""
}


process combine {
    publishDir path: "${params.output_dir}/", mode: 'copy', overwrite: true
    container "${params.container__python}"

    input:
        path "*"

    output:
        path "demux_dash.html"

    """#!/bin/bash

set -euo pipefail

# Copy all of the staged files into the working directory
cp -rL demux_dash/* ./

# Copy the log data into the folder with the other JS files
cp log_data.js js/

# Generate a single-page HTML
generate_single_page.py

# Overwrite the template
mv demux_dash_new.html demux_dash.html

"""
}


workflow demux_dash {

    take:
        demux_stats
        jsons
        sample_sheet
        skeleton_dash

    main:

        prep(
            demux_stats,
            jsons,
            sample_sheet,
            skeleton_dash
        )

        combine(
            prep.out
        )

}