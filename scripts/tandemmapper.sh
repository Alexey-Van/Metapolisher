#!/bin/bash
set -e

### === 1. Создание окружения ===
echo "[1/6] Создаю micromamba-окружение 'tandemtools'..."
micromamba create -y -n tandemtools
micromamba activate tandemtools

### === 2. Скачивание TandemTools ===
echo "[2/6] Клонирую репозиторий TandemTools..."
git clone https://github.com/ablab/TandemTools.git
cd TandemTools

### === 3. Установка зависимостей ===
echo "[3/6] Устанавливаю зависимости..."
micromamba install -y --file requirements.txt
micromamba install -y jellyfish

### === 4. Подготовка входных данных ===
# Укажи свои файлы здесь:
READS="$1"          # ONT или PacBio CLR
ASSEMBLY="$2"    # твоя ETR-сборка
TYPE="$3"
OUTPUT="$4"

# Если сборка одна — используем её дважды (синтаксис требует 2 файла)
ASM1="$ASSEMBLY"
ASM2="$ASSEMBLY"

### === 5. Запуск TandemMapper ===
echo "[4/6] Запускаю TandemMapper..."
python tandemmapper.py \
    --"$TYPE" "$READS" \
    -o "$OUTPUT" \
    "$ASM1" "$ASM2"

### === 6. Запуск TandemQUAST ===
echo "[5/6] Запускаю TandemQUAST..."
python tandemquast.py \
    --"$TYPE" "$READS" \
    -o "$OUTPUT"/quast_out \
    "$ASM1" "$ASM2"

echo "[6/6] Готово!"
echo "Результаты:"
echo "  mapper_out/  — выравнивания (SAM/BED)"
echo "  quast_out/report/ — метрики, k-mer анализ, breakpoint plot"
