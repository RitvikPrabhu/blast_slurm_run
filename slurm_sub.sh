#!/bin/bash

#set -x

DBFILE=$1
QUERYFILE=$2
NTHREADS=$3
ELAPSE=${4:-30:00}
NCBI_BLAST_PATH=${5:-"/lustre/scratch/rprabhu/ncbi-blast-2.13.0+/bin"}
TIME_LOG_FILE=${6:-"blastn_time.log"}
BLASTP_OUTPUT=${7:-"blastn_output.out"}
DATA_DIR=${8:-"data"}
OUTPUT_FILE=${9:-"slurm-${NTHREADS}_threads"}

USAGE="$0 \${NTHREADS}"
if [ -z ${NTHREADS} ]; then
  echo "NTHREADS not set!"
  echo ${USAGE}
  exit 1
fi
if [ ! -f ${DATA_DIR}/${QUERYFILE} ]; then
    echo "Could not find data/${QUERYFILE}"
    echo ${USAGE}
    exit 1;
fi

SLURM_ARGS=(
 -N 1
 --ntasks=1
 --cpus-per-task=${NTHREADS}
 -p cpu
 --exclusive
 --time ${ELAPSE}
 -o ${OUTPUT_FILE}
 --dependency=afterok:176325
 )


TMPFILE=$(mktemp)
cat > $TMPFILE << EOF
#!/bin/bash
# Record the start time
START_TIME=\$(date +%s)


echo "Running BLASTN..."

/usr/bin/time -v -o "${TIME_LOG_FILE}" ${NCBI_BLAST_PATH}/blastn -query "${DATA_DIR}/${QUERYFILE}" -db "${DATA_DIR}/${DBFILE}" -out "${BLASTP_OUTPUT}" -num_threads ${NTHREADS} -outfmt 6 -evalue 0.0000000001 -max_target_seqs 10 -max_hsps 1 -qcov_hsp_perc 60 -perc_identity 60

# Record the end time
END_TIME=\$(date +%s)

# Calculate the duration
DURATION=\$((END_TIME - START_TIME))

# Output the duration
echo "Job started at: \$(date -d @\$START_TIME)"
echo "Job ended at: \$(date -d @\$END_TIME)"
echo "Job duration: \$DURATION seconds"

echo "BLAST job completed."
EOF

sbatch ${SLURM_ARGS[@]} $TMPFILE 
