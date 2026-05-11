#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -pe smp 24
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N KSbuild_allcand_v013

set -euo pipefail

K_VALUES="77 101" \
MAX_PER_SPECIES_PER_K="50000" \

log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO  $*" >&2
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR $*" >&2
}

fail() {
    log_error "$*"
    exit 1
}

run_command() {
    local label="$1"
    shift

    local start_time
    local end_time
    local start_epoch
    local end_epoch
    local runtime_seconds
    local exit_status

    start_time="$(date '+%Y-%m-%d %H:%M:%S')"
    start_epoch="$(date +%s)"

    log_info "Starting command: ${label}"
    log_info "Command: $*"

    set +e
    "$@"
    exit_status="$?"
    set -e

    end_time="$(date '+%Y-%m-%d %H:%M:%S')"
    end_epoch="$(date +%s)"
    runtime_seconds=$((end_epoch - start_epoch))

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${label}" \
        "${start_time}" \
        "${end_time}" \
        "${runtime_seconds}" \
        "${exit_status}" >> "${COMMAND_TIMING_TSV}"

    log_info "Finished command: ${label}; runtime=${runtime_seconds}s; exit_status=${exit_status}"

    if [ "${exit_status}" -ne 0 ]; then
        fail "Command failed: ${label}"
    fi
}

###############################################################################
# User-adjustable settings
###############################################################################

DB_ROOT="${DB_ROOT:-/home/pthorpe001/data/databases/kmersutra_db}"
SOURCE_CONFIG="${SOURCE_CONFIG:-${DB_ROOT}/ncbi_genomes_plasmodium_outgroups_v3/kmersutra_genome_config_targets_supported_roles.tsv}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
BUILD_ROOT="${BUILD_ROOT:-${DB_ROOT}/kmersutra_builds}"
K_VALUES="${K_VALUES:-71}"
SAFE_K_VALUES="$(printf '%s' "${K_VALUES}" | tr ' ' '_')"
OUT_DIR="${OUT_DIR:-${BUILD_ROOT}/kmersutra_plasmodium_outgroups_v3_all_candidate_k${SAFE_K_VALUES}_${RUN_STAMP}}"
TAXONOMY_DIR="${TAXONOMY_DIR:-${DB_ROOT}/ncbi_taxonomy}"
EVIDENCE_RANKS="${EVIDENCE_RANKS:-species genus family order class phylum superkingdom}"
KMERSUTRA_THREADS="${KMERSUTRA_THREADS:-${NSLOTS:-24}}"
SQLITE_BATCH_SIZE="${SQLITE_BATCH_SIZE:-50000}"
MAX_PER_SPECIES_PER_K="${MAX_PER_SPECIES_PER_K:-50000}"
RAM_LOG_INTERVAL_SECONDS="${RAM_LOG_INTERVAL_SECONDS:-30}"
CANDIDATE_ROLES="${CANDIDATE_ROLES:-}"

export MALLOC_ARENA_MAX="${MALLOC_ARENA_MAX:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

###############################################################################
# Output paths
###############################################################################

INPUTS_DIR="${OUT_DIR}/inputs"
METRICS_DIR="${OUT_DIR}/metrics"
TARGET_CONFIG="${INPUTS_DIR}/kmersutra_genome_config_all_candidate.tsv"
RUN_METADATA_TSV="${OUT_DIR}/run_metadata.tsv"
COMMAND_TSV="${OUT_DIR}/kmersutra_build.command.tsv"
COMMAND_TIMING_TSV="${METRICS_DIR}/command_timing.tsv"

mkdir -p "${INPUTS_DIR}" "${METRICS_DIR}"

log_info "Starting KmerSutra v0.13 all-candidate Plasmodium/outgroup build"
log_info "Host: $(hostname)"
log_info "Working directory: $(pwd)"
log_info "DB root: ${DB_ROOT}"
log_info "Source config: ${SOURCE_CONFIG}"
log_info "Output directory: ${OUT_DIR}"
log_info "K values: ${K_VALUES}"
log_info "Threads: ${KMERSUTRA_THREADS}"
log_info "Max per species/evidence/k: ${MAX_PER_SPECIES_PER_K}"
log_info "Candidate roles: ${CANDIDATE_ROLES:-all non-background candidates}"

command -v kmersutra-build-panel >/dev/null 2>&1 || fail "kmersutra-build-panel not found on PATH"
command -v python >/dev/null 2>&1 || fail "python not found on PATH"
[ -s "${SOURCE_CONFIG}" ] || fail "Source config missing or empty: ${SOURCE_CONFIG}"

printf 'label\tstart_time\tend_time\truntime_seconds\texit_status\n' > "${COMMAND_TIMING_TSV}"

###############################################################################
# Copy and normalise config roles conservatively
###############################################################################

awk -F '\t' 'BEGIN {OFS="\t"}
NR == 1 {print; next}
$2 == "Plasmodium falciparum" || $2 == "Plasmodium vivax" || $2 == "Plasmodium knowlesi" {
    $5 = "target_species"
}
$5 == "apicomplexan_outgroup" || $5 == "distant_outgroup" {
    $5 = "outgroup"
}
{print}' "${SOURCE_CONFIG}" > "${TARGET_CONFIG}"

[ -s "${TARGET_CONFIG}" ] || fail "Target config was not created: ${TARGET_CONFIG}"

HEADER_NF="$(awk -F '\t' 'NR == 1 {print NF; exit}' "${TARGET_CONFIG}")"
FIRST_ROW_NF="$(awk -F '\t' 'NR == 2 {print NF; exit}' "${TARGET_CONFIG}")"
if [ "${HEADER_NF}" != "${FIRST_ROW_NF}" ]; then
    fail "Header and first data row have different field counts"
fi

MISSING_FASTA_COUNT="$(
    awk -F '\t' 'NR > 1 {print $1}' "${TARGET_CONFIG}" \
        | while read -r fasta_path; do
            if [ ! -s "${fasta_path}" ]; then
                echo "MISSING_OR_EMPTY\t${fasta_path}"
            fi
        done \
        | tee "${METRICS_DIR}/missing_or_empty_fastas.tsv" \
        | wc -l
)"
if [ "${MISSING_FASTA_COUNT}" -ne 0 ]; then
    fail "Found ${MISSING_FASTA_COUNT} missing or empty FASTA files"
fi

###############################################################################
# Write metadata
###############################################################################

{
    printf 'field\tvalue\n'
    printf 'run_stamp\t%s\n' "${RUN_STAMP}"
    printf 'db_root\t%s\n' "${DB_ROOT}"
    printf 'source_config\t%s\n' "${SOURCE_CONFIG}"
    printf 'target_config\t%s\n' "${TARGET_CONFIG}"
    printf 'out_dir\t%s\n' "${OUT_DIR}"
    printf 'taxonomy_dir\t%s\n' "${TAXONOMY_DIR}"
    printf 'k_values\t%s\n' "${K_VALUES}"
    printf 'kmersutra_threads\t%s\n' "${KMERSUTRA_THREADS}"
    printf 'evidence_ranks\t%s\n' "${EVIDENCE_RANKS}"
    printf 'sqlite_batch_size\t%s\n' "${SQLITE_BATCH_SIZE}"
    printf 'max_per_species_per_k\t%s\n' "${MAX_PER_SPECIES_PER_K}"
    printf 'ram_log_interval_seconds\t%s\n' "${RAM_LOG_INTERVAL_SECONDS}"
    printf 'candidate_roles\t%s\n' "${CANDIDATE_ROLES:-default}"
    printf 'job_id\t%s\n' "${JOB_ID:-NA}"
    printf 'nslots\t%s\n' "${NSLOTS:-NA}"
    printf 'start_time\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
} > "${RUN_METADATA_TSV}"

###############################################################################
# Build command
###############################################################################

read -r -a K_VALUE_ARRAY <<< "${K_VALUES}"
read -r -a EVIDENCE_RANK_ARRAY <<< "${EVIDENCE_RANKS}"

BUILD_COMMAND=(
    kmersutra-build-panel
    --genome_config "${TARGET_CONFIG}"
    --out_dir "${OUT_DIR}"
    --k_values "${K_VALUE_ARRAY[@]}"
    --taxonomy_dir "${TAXONOMY_DIR}"
    --download_taxonomy_if_missing
    --evidence_ranks "${EVIDENCE_RANK_ARRAY[@]}"
    --all_candidate_evidence
    --threads "${KMERSUTRA_THREADS}"
    --sqlite_batch_size "${SQLITE_BATCH_SIZE}"
    --max_per_species_per_k "${MAX_PER_SPECIES_PER_K}"
    --ram_log_interval_seconds "${RAM_LOG_INTERVAL_SECONDS}"
    --profile
    --verbose
)

if [ -n "${CANDIDATE_ROLES}" ]; then
    read -r -a CANDIDATE_ROLE_ARRAY <<< "${CANDIDATE_ROLES}"
    BUILD_COMMAND+=(--candidate_roles "${CANDIDATE_ROLE_ARRAY[@]}")
fi

printf '%q ' "${BUILD_COMMAND[@]}" > "${COMMAND_TSV}"
printf '\n' >> "${COMMAND_TSV}"

run_command "kmersutra_build_panel_all_candidate" "${BUILD_COMMAND[@]}"

###############################################################################
# Basic post-build checks
###############################################################################

PANEL_TSV_GZ="${OUT_DIR}/species_kmer_panel.tsv.gz"
SUMMARY_TSV="${OUT_DIR}/kmer_uniqueness_summary.tsv"
RAM_TSV="${OUT_DIR}/ram_usage.tsv"
ALL_CANDIDATE_SQLITE="${OUT_DIR}/all_candidate_evidence.sqlite"

[ -s "${PANEL_TSV_GZ}" ] || fail "Expected panel file missing or empty: ${PANEL_TSV_GZ}"
[ -s "${SUMMARY_TSV}" ] || fail "Expected summary missing or empty: ${SUMMARY_TSV}"
[ -s "${ALL_CANDIDATE_SQLITE}" ] || fail "Expected all-candidate SQLite missing or empty: ${ALL_CANDIDATE_SQLITE}"

log_info "Panel file size:"
ls -lh "${PANEL_TSV_GZ}" >&2

log_info "Summary head:"
column -t -s $'\t' "${SUMMARY_TSV}" | head -n 40 >&2 || head -n 40 "${SUMMARY_TSV}" >&2

if [ -s "${RAM_TSV}" ]; then
    log_info "RAM usage tail:"
    tail -n 10 "${RAM_TSV}" >&2
fi

{
    printf 'end_time\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'panel_tsv_gz\t%s\n' "${PANEL_TSV_GZ}"
    printf 'summary_tsv\t%s\n' "${SUMMARY_TSV}"
    printf 'ram_tsv\t%s\n' "${RAM_TSV}"
    printf 'all_candidate_sqlite\t%s\n' "${ALL_CANDIDATE_SQLITE}"
    printf 'command_timing_tsv\t%s\n' "${COMMAND_TIMING_TSV}"
} >> "${RUN_METADATA_TSV}"

log_info "Done"
log_info "Output directory: ${OUT_DIR}"

if [ -n "${JOB_ID:-}" ]; then
    log_info "After the job has exited, run:"
    log_info "qacct -j ${JOB_ID} | grep -E 'jobnumber|hostname|exit_status|failed|maxvmem|ru_wallclock|cpu|mem|io'"
fi
