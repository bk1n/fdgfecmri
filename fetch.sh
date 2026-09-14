#!/usr/bin/env bash
set -euo pipefail

source .env

SRC=$(wslpath "$DATA_PATH")

rsync -avzL "$SRC/" ./data/
