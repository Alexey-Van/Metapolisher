process CHECK {

    tag "check_containers"
    executor 'local'

    output:
        val(true)

    script:
    """
    echo "=== Checking containers ==="

    containers=(
        "google/deepvariant:1.9.0"
        "kishwars/pepper_deepvariant:r0.8"
        "metapolisher-align:1.0"
        "metapolisher-polish:1.0"
        "metapolisher-sv:1.0"
        "metapolisher-ml:1.0"
    )

    for image in "\${containers[@]}"; do

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
