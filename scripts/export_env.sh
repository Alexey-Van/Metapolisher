#!/bin/bash

# Путь к окружениям micromamba
ENV_DIR="$HOME/micromamba/envs"

# Папка для экспорта
mkdir -p exported_envs

for env in "$ENV_DIR"/*; do
    if [ -d "$env" ]; then
        env_name=$(basename "$env")
        echo "Экспортирую $env_name..."

        micromamba env export -n "$env_name" > "exported_envs/${env_name}.yml"
    fi
done

echo "Готово!"
