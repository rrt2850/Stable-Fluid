#!/usr/bin/env bash
set -euo pipefail

GRID_WIDTH="${1:-128}"
GRID_HEIGHT="${2:-128}"
EMITTER_COUNT="${3:-8}"
SIMULATION_STEPS="${4:-180}"

javac src/*.java
java -cp src Main "$GRID_WIDTH" "$GRID_HEIGHT" "$EMITTER_COUNT" "$SIMULATION_STEPS"
