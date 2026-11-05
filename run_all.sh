#!/usr/bin/env bash
#
# Ricostruisce data/reference/ e output/ a partire dall'uscita di AugmentA
# (data/augmenta/, nel repository). Va lanciato con l'ambiente conda tesi attivo.
#
# Rigenerare l'uscita di AugmentA dall'SSM e' facoltativo e resta a mano:
# vedi README.
#
set -euo pipefail
cd "$(dirname "$0")"

python simulator/stage1/opencarp/mm_to_um.py      # -> data/reference/LA_um.*

# scenari aniso_zone_N: N zone casuali (seme 0), una delle 25 combinazioni di CV di CV_CSV per zona;
# N = 1: CV uniformi, N = 43270: una zona per faccia
ZONE="1 50 500 5000 25000 43270"
CV_CSV=experiments/region_cv/aniso_random.csv
for n in $ZONE; do
    s=aniso_zone_$n; mesh=data/reference/$s/LA_um
    python simulator/stage1/opencarp/random_tags.py $s $CV_CSV 0 $n   # -> $mesh.*
    python simulator/stage1/opencarp/run_sims.py $s $CV_CSV $mesh     # -> data/reference/sims/$s/run_*/act.igb
    python simulator/convert/toPly.py $s $CV_CSV $mesh                # -> data/reference/LA_{LAT,GT}_$s.ply
done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build

python experiments/zone.py $ZONE                  # -> output/zone/confs/*.conf
./estimator/build/Tesi0 output/zone/confs/*.conf
python analysis/analisi.py zone/aniso_zone_1      # -> output/zone/aniso_zone_1/figure/
python analysis/analisi.py --zone                 # -> output/zone/figure/

python experiments/subsets.py aniso_zone_1 output/subsets   # -> output/subsets/confs/*.conf
./estimator/build/Tesi0 output/subsets/confs/*.conf
python analysis/analisi.py --sottoinsiemi         # -> output/subsets/figure/

# metodo 2 di Roney et al. 2019 contro lo stimatore per faccia, solo terne di mappe
for n in $ZONE; do
    python experiments/subsets.py aniso_zone_$n output/roney/aniso_zone_$n/faccia 3
done
./estimator/build/Tesi0 output/roney/*/faccia/confs/*.conf
python estimator/roney.py experiments/roney.conf $(for n in $ZONE; do echo aniso_zone_$n; done)
python analysis/analisi.py --roney                # -> output/roney/figure/

# mappe degli errori per faccia in PLY, su aniso_zone_1: tutte le mappe e, per braccio, le terne con
# angolo q50 minimo e massimo (da output/roney/figure/terne.csv); tre PLY per caso, uno per errore
python experiments/mappe.py | while read -r csv base; do   # -> output/mappe/confs/*.conf
    ./estimator/build/Tesi0 --ply data/reference/LA_LAT_aniso_zone_1.ply "$csv" "$base"
done
./estimator/build/Tesi0 output/mappe/confs/*.conf  # -> output/mappe/*_{angolo,cv_l,rapporto}.ply
