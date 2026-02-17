#!/usr/bin/env bash
set -euo pipefail

GRID_WIDTH="${1:-128}"
GRID_HEIGHT="${2:-128}"
EMITTER_COUNT="${3:-8}"

javac src/*.java
java -cp src Main "$GRID_WIDTH" "$GRID_HEIGHT" "$EMITTER_COUNT"
