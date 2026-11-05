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
# d = sqrt(4 A / (pi N)) ~ 123, 44, 15, 5.0 mm, passi uguali in scala log; N = 1: CV uniformi,
# N = 600: zone grandi come il disco di Roney (5 mm), che presuppone CV uniformi al suo interno
ZONE="1 8 70 600"
CV_CSV=experiments/region_cv/aniso_random.csv
for n in $ZONE; do
    s=aniso_zone_$n; mesh=data/reference/$s/LA_um
    python simulator/stage1/opencarp/random_tags.py $s $CV_CSV 0 $n   # -> $mesh.*
    python simulator/stage1/opencarp/run_sims.py $s $CV_CSV $mesh     # -> data/reference/sims/$s/run_*/act.igb
    python simulator/convert/toPly.py $s $CV_CSV $mesh                # -> data/reference/LA_{LAT,GT}_$s.ply
done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build

rm -rf output                                     # niente risultati vecchi nelle figure

# griglia: per ogni scenario tutti i sottoinsiemi di almeno 3 mappe con lo stimatore per faccia,
# e tutte le terne con il metodo 2 di Roney et al. 2019
for n in $ZONE; do
    python experiments/subsets.py aniso_zone_$n   # -> output/aniso_zone_$n/confs/*.conf
done
./estimator/build/Tesi0 output/aniso_zone_*/confs/*.conf   # -> output/SCENARIO/faccia/sub_*/errors.csv
# -> output/SCENARIO/roney/sub_i-j-k/errors.csv
python estimator/roney.py experiments/roney.conf $(for n in $ZONE; do echo aniso_zone_$n; done)
python analysis/analisi.py                        # -> output/figure/

# mappe degli errori in PLY: per faccia con tutte le mappe su aniso_zone_1 e aniso_zone_600, e i due
# metodi sulla terna di aniso_zone_600 con l'area appaiata piu' grande; tre PLY per caso, uno per errore
python experiments/mappe.py | while read -r mesh csv prefisso; do
    ./estimator/build/Tesi0 --ply "$mesh" "$csv" "$prefisso"   # -> output/mesh/*_{angolo,cv_l,rapporto}.ply
done
