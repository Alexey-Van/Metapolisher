process CHECK {

    tag "check_containers"
    executor 'local'

    input:
        val containers

    output:
        val(true)

    script:
    def container_args = containers.collect { image ->
        def escaped_image = image.toString().replace("'", "'\\''")
        "'${escaped_image}'"
    }.join(' ')

    """
    set -euo pipefail

    echo "=== Checking containers ==="

    for image in ${container_args}; do

        if [[ "\$image" == *.sif ]]; then

            # Singularity / Apptainer image
            if [[ -f "\$image" && -r "\$image" ]]; then
                echo "[OK] SIF image exists: \$image"

            else
                echo "[ERROR] Missing or unreadable SIF image: \$image"
                exit 1
            fi

        else

            # Docker image
            if docker image inspect "\$image" >/dev/null 2>&1; then
                echo "[OK] Docker image exists: \$image"

            else
                echo "[ERROR] Missing Docker image: \$image"
                exit 1
            fi

        fi

    done

    echo "All containers available"
    """
}
