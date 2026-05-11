#!/bin/bash

# Проверка аргументов
if [ "$#" -lt 3 ]; then
    echo "Использование: $0 <reference.fasta> <contigs.fasta> <threads>"
    exit 1
fi

REFERENCE=$1
CONTIGS=$2
THREADS=$3

echo "=== RagTag pipeline ==="
echo "Reference: $REFERENCE"
echo "Contigs:   $CONTIGS"
echo "Threads:   $THREADS"
echo

# Создаём директории
mkdir -p ragtag_correct
mkdir -p ragtag_scaffold

echo "=== Шаг 1: Коррекция сборки ==="
ragtag.py correct \
    -t $THREADS \
    $REFERENCE \
    $CONTIGS \
    -o ragtag_correct

CORRECTED="ragtag_correct/ragtag.correct.fasta"

echo
echo "Коррекция завершена. Используем файл:"
echo "$CORRECTED"
echo

echo "=== Шаг 2: Скафолдинг (сшивание контигов) ==="
ragtag.py scaffold \
    -t $THREADS \
    $REFERENCE \
    $CORRECTED \
    -o ragtag_scaffold

echo
echo "=== Готово! ==="
echo "Результаты:"
echo "  ragtag_scaffold/ragtag.scaffold.fasta   — итоговые хромосомы"
echo "  ragtag_scaffold/ragtag.scaffold.agp     — порядок и ориентация"
echo "  ragtag_scaffold/ragtag.scaffold.paf     — выравнивания"
echo
