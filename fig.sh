#!/usr/bin/env bash

for fig_file in fig_*.py; do
    python3 "$fig_file"
done

if [ -n $( command -v "gs" ) ]; then
    # if ghostscript is installed, optimize PDF figures
    find figs/ -maxdepth 1 -name "*.pdf" \
        -exec gs -q -o {}.opt -dNoOutputFonts -sPAPERSIZE=a4 -sDEVICE=pdfwrite {} \; \
        -exec mv {}.opt {} \;
fi
