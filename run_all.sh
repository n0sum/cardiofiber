#!/usr/bin/env bash
#
# Ricostruisce data/reference/ e output/ a partire dall'uscita di AugmentA
# (data/augmenta/). Va lanciato con l'ambiente conda tesi attivo.
#
# I passi prima di AugmentA (make_mean, split_atria, extract_la) e AugmentA
# stesso (Docker + sessione desktop) restano a mano: vedi README passi 1-2.
#
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f data/augmenta/LA_endo_with_fiber.pts ]; then
    echo "manca data/augmenta/LA_endo_with_fiber.*: serve prima AugmentA, a mano (README passo 2)"
    exit 1
fi

python simulator/stage1/opencarp/mm_to_um.py      # -> data/reference/LA_um.*
python simulator/stage1/opencarp/run_sims.py experiments/region_cv/aniso_const.csv   # -> data/reference/sims/aniso_const/run_*/act.igb
python simulator/convert/multiLatToPly.py aniso_const                             # -> data/reference/LA_LAT_aniso_const.ply
python simulator/convert/fibersToPly.py           # -> data/reference/LA_fibers_gt.ply

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build
./estimator/build/Tesi0 experiments/aniso_const.conf
./estimator/build/Tesi0 experiments/aniso_const_nomask.conf

python analysis/analisi.py aniso_const
