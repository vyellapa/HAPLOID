# Sourced by the GLIMPSE2 scripts (needs ROOT set): defines glimpse <tool> <args...>, which runs a GLIMPSE2 tool
# through the runner in GLIMPSE_RUNNER:
#   native    : static release binaries in GLIMPSE_BIN (default resources/glimpse_bin, from fetch_glimpse.sh),
#               named GLIMPSE2_<tool>_static or GLIMPSE2_<tool>          (HPC default: no container needed)
#   docker    : the GLIMPSE_IMG image (convenient for local testing)
#   apptainer / singularity : GLIMPSE_SIF image file
# Default: native when the binaries exist, else docker.
GLIMPSE_BIN=${GLIMPSE_BIN:-$ROOT/resources/glimpse_bin}
GLIMPSE_IMG=${GLIMPSE_IMG:-simrub/glimpse:v2.0.0-27-g0919952_20221207}
GLIMPSE_SIF=${GLIMPSE_SIF:-$ROOT/resources/glimpse.sif}
if [ -z "${GLIMPSE_RUNNER:-}" ]; then
  if [ -x "$GLIMPSE_BIN/GLIMPSE2_phase_static" ] || [ -x "$GLIMPSE_BIN/GLIMPSE2_phase" ]; then GLIMPSE_RUNNER=native
  else GLIMPSE_RUNNER=docker; fi
fi
export GLIMPSE_BIN GLIMPSE_IMG GLIMPSE_SIF GLIMPSE_RUNNER

glimpse() {
  local tool=$1; shift
  case $GLIMPSE_RUNNER in
    native)
      if [ -x "$GLIMPSE_BIN/${tool}_static" ]; then "$GLIMPSE_BIN/${tool}_static" "$@"; else "$GLIMPSE_BIN/$tool" "$@"; fi ;;
    docker)
      docker run --rm --platform linux/amd64 -v "$ROOT:$ROOT" -w "$ROOT" "$GLIMPSE_IMG" "$tool" "$@" ;;
    apptainer|singularity)
      "$GLIMPSE_RUNNER" exec -B "$ROOT" --pwd "$ROOT" "$GLIMPSE_SIF" "$tool" "$@" ;;
    *) echo "unknown GLIMPSE_RUNNER=$GLIMPSE_RUNNER" >&2; return 2 ;;
  esac
}
export -f glimpse
