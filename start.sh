#!/usr/bin/env bash
set -euo pipefail

GRID_WIDTH="${1:-128}"
GRID_HEIGHT="${2:-128}"
EMITTER_COUNT="${3:-8}"
SIMULATION_STEPS="${4:-180}"
EXPORT_VIDEO="${5:-false}"
LOG_EVERY_STEP="${6:-false}"

BUILD_DIR="out/classes"
mkdir -p "$BUILD_DIR"

javac -d "$BUILD_DIR" src/*.java
java -cp "$BUILD_DIR" Main "$GRID_WIDTH" "$GRID_HEIGHT" "$EMITTER_COUNT" "$SIMULATION_STEPS" "$EXPORT_VIDEO" "$LOG_EVERY_STEP"
