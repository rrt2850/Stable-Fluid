#!/usr/bin/env bash
set -euo pipefail

EMITTER_COUNT="${1:-8}"

javac src/*.java
java -cp src Main "$EMITTER_COUNT"
