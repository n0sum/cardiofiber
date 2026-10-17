#!/usr/bin/env bash
#
# Ricostruisce data/ e output/ a partire dal codice e dalla
# configurazione versionati.
#
# NOTA: il passo AugmentA non e' automatizzabile (Docker + sessione
# desktop). Lo script si ferma e lo chiede, invece di fingere.
#
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p data/{external,augmenta,reference,clinical} output

# TODO: riempire mano a mano che gli stadi vengono consolidati.
#   1. verifica checksums.sha256 su data/external
#   2. simulator/stage1/*.py
#   3. STOP: AugmentA a mano, vedi README passo 2
#   4. simulator/stage1/opencarp/run_eik.sh
#   5. simulator/convert/*.py
#   6. estimator, su ogni experiments/*.conf
#   7. analysis/analisi.py

echo "non ancora implementato: vedi i passi 1-6 del README"
exit 1
