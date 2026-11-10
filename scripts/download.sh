#!/usr/bin/env bash
# Download every row of metadata/sources.tsv, check its md5, and unpack it.
# Run from the repo root:  pixi run download
set -euo pipefail

SOURCES="metadata/sources.tsv"

tail -n +2 "$SOURCES" | while IFS=$'\t' read -r name url md5 unpack dest; do
  [[ -z "$name" ]] && continue

  if [[ -e "$dest" ]]; then
    echo "SKIP $name: $dest already exists"
    continue
  fi

  echo "GET  $name"
  tmp=$(mktemp)
  curl -fsSL --retry 3 -o "$tmp" "$url"

  got=$(md5sum "$tmp" | cut -d' ' -f1)
  if [[ "$got" != "$md5" ]]; then
    echo "FAIL $name: md5 $got, expected $md5" >&2
    rm -f "$tmp"
    exit 1
  fi
  echo "OK   $name md5"

  case "$unpack" in
    tar)  mkdir -p "$dest"
          tar -xzf "$tmp" -C "$dest" --strip-components=1 ;;  # drop the top-level folder
    zip)  mkdir -p "$dest"
          unzip -q "$tmp" -d "$dest" ;;
    none) mkdir -p "$(dirname "$dest")"
          mv "$tmp" "$dest" ;;
    *)    echo "FAIL $name: unknown unpack '$unpack'" >&2; exit 1 ;;
  esac
  rm -f "$tmp"

  chmod -R a-w "$dest"   # raw data is read-only
done
