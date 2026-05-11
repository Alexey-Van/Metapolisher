#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

docker pull google/deepvariant:1.9.0
docker pull kishwars/pepper_deepvariant:r0.8

docker build -t metapolisher-align:1.0   "${SCRIPT_DIR}/containers/align"
docker build -t metapolisher-polish:1.0  "${SCRIPT_DIR}/containers/polish"
docker build -t metapolisher-sv:1.0      "${SCRIPT_DIR}/containers/sv"
docker build -t metapolisher-ml:1.0      "${SCRIPT_DIR}/containers/ml"