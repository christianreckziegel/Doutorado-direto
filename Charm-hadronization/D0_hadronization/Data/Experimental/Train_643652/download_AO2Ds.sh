#!/bin/bash

############################################
# CONFIGURATION
############################################
INPUT_LIST="paths.txt"
OUTPUT_DIR="AO2Ds_non_merged"
LOGFILE="download.log"
N_PARALLEL=1        # number of parallel downloads
N_RETRIES=3         # retries per file

############################################
# PREPARATION
############################################
mkdir -p "$OUTPUT_DIR"
touch "$LOGFILE"

echo "===== AO2D DOWNLOAD START =====" | tee -a "$LOGFILE"
date | tee -a "$LOGFILE"

############################################
# LOAD O2 ENVIRONMENT
############################################
echo "Loading O2Physics environment..."
eval "$(alienv shell-helper)"

alienv enter O2Physics/latest <<EOF

############################################
# FUNCTION TO DOWNLOAD ONE FILE
############################################
download_one() {
    path="\$1"

    # Skip empty lines
    [ -z "\$path" ] && return

    # Ensure /AOD exists
    if [[ "\$path" != */AOD ]]; then
        path="\${path}/AOD"
    fi

    # Extract job ID
    jobid=\$(echo "\$path" | grep -o 'hy_[0-9]*')

    # Fallback if extraction fails
    if [ -z "\$jobid" ]; then
        jobid=\$(echo "\$path" | sed 's#/#_#g')
    fi

    outfile="${OUTPUT_DIR}/AO2D_\${jobid}.root"

    # Skip if already exists
    if [ -f "\$outfile" ]; then
        echo "[SKIP] \$outfile already exists"
        return
    fi

    echo "[INFO] Downloading \$path → \$outfile"

    # Retry loop
    for ((i=1; i<=${N_RETRIES}; i++)); do
        echo "  Attempt \$i..."

        alien_cp alien:\${path}/001/AO2D.root file:\${outfile}

        if [ \$? -eq 0 ]; then
            echo "[OK] \$outfile"
            return
        else
            echo "[WARN] Failed attempt \$i for \$path"
            sleep 2
        fi
    done

    echo "[ERROR] Failed after ${N_RETRIES} attempts: \$path"
}

export -f download_one
export OUTPUT_DIR
export N_RETRIES

############################################
# RUN DOWNLOADS (PARALLEL)
############################################
cat "$INPUT_LIST" | tr ',' '\n' | \
    xargs -I {} -P ${N_PARALLEL} bash -c 'download_one "{}"'

############################################
# FINISH
############################################
echo "===== DOWNLOAD FINISHED ====="
date

EOF