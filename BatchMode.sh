#!/usr/bin/env bash

set -euo pipefail

# -------------------------
# Usage function
# -------------------------
usage() {
  echo "Usage: $0 config.txt [--verbose]"
  exit 1
}

# -------------------------
# Argument parsing
# -------------------------
VERBOSE=0
CONFIG_FILE=""

for arg in "$@"; do
  case "$arg" in
  --verbose)
    VERBOSE=1
    ;;
  *)
    if [[ -z "${CONFIG_FILE}" ]]; then
      CONFIG_FILE="$arg"
    else
      usage
    fi
    ;;
  esac
done

[[ -z "${CONFIG_FILE}" ]] && usage
[[ ! -f "${CONFIG_FILE}" ]] && {
  echo "Error: ${CONFIG_FILE} not found"
  exit 1
}

# -------------------------
# Main loop
# -------------------------
# Skip the header (NR > 1)
awk 'NR > 1 { print $1, $2 }' "${CONFIG_FILE}" | while read -r E THETA; do
  echo "Running PBosted with E=${E}, theta=${THETA}"

  CMD=(
    ./Fit
    --input "params.txt"
    --output "output/E=${E}_theta=${THETA}.txt"
    --p "beam_energy=${E}"
    --p "outgoing_lepton_angle=${THETA}"
  )

  if [[ "${VERBOSE}" -eq 1 ]]; then
    "${CMD[@]}"
  else
    "${CMD[@]}" >/dev/null 2>&1
  fi
done
