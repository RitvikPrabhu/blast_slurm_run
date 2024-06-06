#!/bin/bash

# Arguments for slurm_sub.sh
DBFILE="nr"
QUERYFILE="non-rRNA-reads_sample50.fa"
ELAPSE="12:00:00"
NCBI_BLAST_PATH="/lustre/scratch/rprabhu/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin"
DATA_DIR="/scratch/user/u.rp167879/sparkleblast_data"
THREAD_COUNTS=(1 2 4 8 16 32 48)
NUM_RUNS=10

submit_job() {
    local nthreads=$1
    local run_number=$2
    local time_log_file="blastp_time_${nthreads}_threads_run${run_number}.log"
    local blastp_output="blastp_output_${nthreads}_threads_run${run_number}.out"
    local output_file="slurm_${nthreads}_threads_run${run_number}"

    ./slurm_sub.sh "$DBFILE" "$QUERYFILE" "$nthreads" "$ELAPSE" "$NCBI_BLAST_PATH" "$time_log_file" "$blastp_output" "$DATA_DIR" "$output_file"
}

for nthreads in "${THREAD_COUNTS[@]}"; do
    for run_number in $(seq 1 $NUM_RUNS); do
        submit_job $nthreads $run_number
    done
done

