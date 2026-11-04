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
python simulator/stage1/opencarp/random_tags.py aniso_random experiments/region_cv/aniso_random.csv 0 43270   # -> data/reference/aniso_random/LA_um.*

# numero di zone degli scenari aniso_zone_N (43270 = una zona per faccia)
ZONE="43270 10000 2500 600 150 40"
for n in $ZONE; do
    python simulator/stage1/opencarp/random_tags.py aniso_zone_$n experiments/region_cv/aniso_random.csv 0 $n   # -> data/reference/aniso_zone_$n/LA_um.*
done

# scenario NOME: velocita' per tag da CSV, simulato sulla mesh MESH_BASE (in um, con i tag dello scenario)
simula() {
    s=$1; csv=$2; mesh=$3
    python simulator/stage1/opencarp/run_sims.py $s $csv $mesh   # -> data/reference/sims/$s/run_*/act.igb
    python simulator/convert/multiLatToPly.py $s                 # -> data/reference/LA_LAT_$s.ply
    python simulator/convert/fibersToPly.py $s $csv $mesh        # -> data/reference/LA_GT_$s.ply
}
simula aniso_const    experiments/region_cv/aniso_const.csv    data/reference/LA_um
simula aniso_regional experiments/region_cv/aniso_regional.csv data/reference/LA_um
simula aniso_random   experiments/region_cv/aniso_random.csv   data/reference/aniso_random/LA_um
for n in $ZONE; do
    simula aniso_zone_$n experiments/region_cv/aniso_random.csv data/reference/aniso_zone_$n/LA_um
done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build
./estimator/build/Tesi0 experiments/aniso_const.conf
./estimator/build/Tesi0 experiments/aniso_const_nomask.conf
./estimator/build/Tesi0 experiments/aniso_regional.conf
./estimator/build/Tesi0 experiments/aniso_random.conf
./estimator/build/Tesi0 experiments/aniso_random_nomask.conf

for e in aniso_const aniso_const_nomask aniso_regional aniso_random aniso_random_nomask; do
    python analysis/analisi.py $e                 # -> output/$e/figure/
done

for e in aniso_random aniso_random_nomask; do
    python analysis/vicini.py $e                  # -> output/$e/figure/vicini.csv
done

python experiments/subsets.py experiments/aniso_const.conf output/subsets   # -> output/subsets/confs/*.conf
./estimator/build/Tesi0 output/subsets/confs/*.conf
python analysis/analisi.py --sottoinsiemi         # -> output/subsets/figure/

python experiments/zone.py $ZONE                  # -> output/zone/confs/*.conf
./estimator/build/Tesi0 output/zone/confs/*.conf
python analysis/analisi.py --zone                 # -> output/zone/figure/

# metodo 2 di Roney et al. 2019 contro lo stimatore per faccia, solo terne di mappe
python experiments/subsets.py experiments/aniso_const.conf output/roney/aniso_const/faccia 3
for n in $ZONE; do
    python experiments/subsets.py output/zone/confs/aniso_zone_$n.conf output/roney/aniso_zone_$n/faccia 3
done
./estimator/build/Tesi0 output/roney/*/faccia/confs/*.conf
python estimator/roney.py experiments/roney.conf aniso_const $(for n in $ZONE; do echo aniso_zone_$n; done)
python analysis/analisi.py --roney                # -> output/roney/figure/
