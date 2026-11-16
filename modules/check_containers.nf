process CHECK {

    tag "check_containers"
    executor 'local'

    input:
        val containers

    output:
        val(true)

    script:
    def container_args = containers.collect { it.toString() }.join(' ')

    """
    echo "=== Checking containers ==="

    for image in ${container_args}; do

        if docker image inspect \$image >/dev/null 2>&1; then
            echo "[OK] \$image exists"
        else
            echo "[ERROR] Missing image: \$image"
            exit 1
        fi
    done

    echo "All containers available"
    """
}
