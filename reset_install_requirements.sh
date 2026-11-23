#!/usr/bin/env bash
set -euo pipefail
usage() {
  printf '%s\n' 'Usage: bash reset_install_requirements.sh /absolute/path/to/trusted/basilisk/src' >&2
  exit 2
}
[[ $# -eq 1 && "$1" = /* && -x "$1/qcc" && -f "$1/view.h" ]] || usage
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
src="$(cd -- "$1" && pwd -P)"
if [[ "$src" == "$root/basilisk/src" ]]; then
  printf '%s\n' 'The bundled Basilisk reference tree is not an installation.' >&2
  exit 2
fi
umask 077
config="$(mktemp "$root/.project_config.XXXXXX")"
trap 'rm -f -- "$config"' EXIT
printf 'export BASILISK=%q\nexport PATH="$BASILISK:$PATH"\n' "$src" > "$config"
mv -f -- "$config" "$root/.project_config"
