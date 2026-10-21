#!/bin/bash
# uso: ./run_eik.sh MESH_um NOME X_um Y_um Z_um ARGOMENTI_GREGION...   (MESH_um senza estensione)
# Le gregion con le velocita' le costruisce run_sims.py dal CSV dello scenario.
MESH=$1; NAME=$2; X=$3; Y=$4; Z=$5; shift 5
openCARP \
 -meshname $MESH -simID run_$NAME \
 -dream.solve 0 -dream.fim.tol 0.001 -dream.fim.max_iter 10000000 \
 -dream.output.atfile act -tend 500 \
 -num_phys_regions 1 \
 -phys_region[0].name "Eikonal domain" -phys_region[0].ptype 2 -phys_region[0].num_IDs 7 \
 -phys_region[0].ID[0] 1 -phys_region[0].ID[1] 2 -phys_region[0].ID[2] 3 -phys_region[0].ID[3] 4 \
 -phys_region[0].ID[4] 5 -phys_region[0].ID[5] 6 -phys_region[0].ID[6] 7 \
 -num_stim 1 -stim[0].elec.geom_type 1 \
 -stim[0].elec.p0[0] $X -stim[0].elec.p0[1] $Y -stim[0].elec.p0[2] $Z \
 -stim[0].elec.radius 3000 -stim[0].ptcl.start 0.0 \
 "$@"
