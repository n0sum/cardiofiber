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

if [ ! -f data/augmenta/LA_endo_with_fiber.pts ]; then
    echo "manca data/augmenta/LA_endo_with_fiber.pts: e' nel repository, il clone e' incompleto"
    exit 1
fi

python simulator/stage1/opencarp/mm_to_um.py      # -> data/reference/LA_um.*

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
simula aniso_const experiments/region_cv/aniso_const.csv data/reference/LA_um
for n in $ZONE; do
    simula aniso_zone_$n experiments/region_cv/aniso_random.csv data/reference/aniso_zone_$n/LA_um
done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build
./estimator/build/Tesi0 experiments/aniso_const.conf
python analysis/analisi.py aniso_const            # -> output/aniso_const/figure/

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

# mappe degli errori per faccia in PLY, su aniso_const: sei mappe e, per braccio, le terne con
# angolo q50 minimo e massimo (da output/roney/figure/terne.csv)
python experiments/mappe.py | while read -r csv ply; do   # -> output/mappe/confs/*.conf
    ./estimator/build/Tesi0 --ply data/reference/LA_LAT_aniso_const.ply "$csv" "$ply"
done
./estimator/build/Tesi0 output/mappe/confs/*.conf  # -> output/mappe/*.ply
