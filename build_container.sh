#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF_DIR="${SIF_DIR:-${SCRIPT_DIR}/containers/sif}"

docker pull google/deepvariant:1.10.0
docker pull kishwars/pepper_deepvariant:r0.8
docker pull mobinasri/flagger:v1.2.1

docker build -t metapolisher-align:2.0   "${SCRIPT_DIR}/containers/align"
docker build -t metapolisher-polish:2.0  "${SCRIPT_DIR}/containers/polish"
docker build -t metapolisher-sv:2.0      "${SCRIPT_DIR}/containers/sv"
docker build -t metapolisher-ml:2.0      "${SCRIPT_DIR}/containers/ml"
docker build -t metapolisher-truth:1.0   "${SCRIPT_DIR}/containers/truth"

if command -v apptainer >/dev/null 2>&1; then
    SIF_BUILDER=apptainer
elif command -v singularity >/dev/null 2>&1; then
    SIF_BUILDER=singularity
else
    echo "ERROR: install Apptainer or Singularity to create SIF images" >&2
    exit 127
fi

mkdir -p "${SIF_DIR}"

build_sif() {
    local docker_image="$1"
    local sif_name="$2"
    local sif_path="${SIF_DIR}/${sif_name}"

    echo "Building ${sif_path} from ${docker_image}"
    "${SIF_BUILDER}" build --force "${sif_path}" "docker-daemon:${docker_image}"
}

build_sif metapolisher-align:2.0 metapolisher-align-2.0.sif
build_sif metapolisher-polish:2.0 metapolisher-polish-2.0.sif
build_sif metapolisher-sv:2.0 metapolisher-sv-2.0.sif
build_sif metapolisher-ml:2.0 metapolisher-ml-2.0.sif
build_sif metapolisher-truth:1.0 metapolisher-truth-1.0.sif
build_sif google/deepvariant:1.10.0 google-deepvariant-1.10.0.sif
build_sif kishwars/pepper_deepvariant:r0.8 kishwars-pepper-deepvariant-r0.8.sif
build_sif mobinasri/flagger:v1.2.1 mobinasri-flagger-v1.2.1.sif

echo "SIF images created in ${SIF_DIR}"
