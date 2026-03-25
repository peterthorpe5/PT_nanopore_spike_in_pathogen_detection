#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -jc short
#$ -mods l_hard mfree 16G
#$ -adds l_hard h_vmem 16G
#$ -N make_krona

set -euo pipefail

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

require_exe() {
    command -v "$1" >/dev/null 2>&1 || {
        log_error "Required executable not found on PATH: $1"
        exit 1
    }
}

if [[ -n "${REPO_DIR:-}" ]]; then
    REPO_DIR="${REPO_DIR}"
elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then
    REPO_DIR="${SGE_O_WORKDIR}"
else
    REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"
fi

RUNS_DIR="${RUNS_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs}"
OUT_MANIFEST="${OUT_MANIFEST:-${RUNS_DIR}/krona_manifest.tsv}"
UPDATE_TAXONOMY="${UPDATE_TAXONOMY:-false}"

require_exe ktImportTaxonomy
require_exe cut
require_exe find

if [[ "${UPDATE_TAXONOMY}" == "true" ]]; then
    require_exe ktUpdateTaxonomy.sh
    log_info "Updating Krona taxonomy"
    ktUpdateTaxonomy.sh
fi

mkdir -p "$(dirname "${OUT_MANIFEST}")"
printf 'source_tsv\tkrona_input_tsv\toutput_html\n' > "${OUT_MANIFEST}"

mapfile -t KRKN_FILES < <(
    find "${RUNS_DIR}" -type f \
        \( -name 'kraken.classifications.tsv' -o -name 'assembly.kraken.classifications.tsv' \) \
        | sort
)

if [[ "${#KRKN_FILES[@]}" -eq 0 ]]; then
    log_error "No Kraken classification files found under ${RUNS_DIR}"
    exit 1
fi

log_info "Found ${#KRKN_FILES[@]} Kraken classification files"

for kraken_tsv in "${KRKN_FILES[@]}"; do
    run_dir="$(dirname "${kraken_tsv}")"
    krona_dir="${run_dir}/krona"
    mkdir -p "${krona_dir}"

    base_name="$(basename "${kraken_tsv}" .tsv)"
    krona_input="${krona_dir}/${base_name}.krona_input.tsv"
    out_html="${krona_dir}/${base_name}.krona.html"

    log_info "Preparing Krona input for ${kraken_tsv}"

    cut -f 2,3 "${kraken_tsv}" > "${krona_input}"

    log_info "Creating Krona plot ${out_html}"
    ktImportTaxonomy \
        -o "${out_html}" \
        "${krona_input}"

    printf '%s\t%s\t%s\n' "${kraken_tsv}" "${krona_input}" "${out_html}" >> "${OUT_MANIFEST}"
done

log_info "Done"
log_info "Manifest: ${OUT_MANIFEST}"