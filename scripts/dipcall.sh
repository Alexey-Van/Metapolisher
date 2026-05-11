#!/bin/bash
set -e

# Путь к каталогу dipcall.kit
DIPCALL=/mnt/data/vanichkinao/MetaPolisher/dipcall.kit

# Входные файлы
REF="$1"
PAT="$2"
MAT="$3"

# Префикс для выходных файлов
PREFIX=dipcall_output
# Генерация Makefile
$DIPCALL/run-dipcall $PREFIX $REF $PAT $MAT > $PREFIX.mak

# Запуск пайплайна
make -j2 -f $PREFIX.mak
