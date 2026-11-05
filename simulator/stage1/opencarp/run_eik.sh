#!/bin/bash
# uso: ./run_eik.sh MESH_um NOME X_um Y_um Z_um ALTRI_ARGOMENTI...   (MESH_um senza estensione)
# Dominio eikonale (phys_region) e gregion con le velocita' li passa
# run_sims.py. Il FIM si ferma dopo 1000 iterazioni: le LAT sono ferme da
# circa 120, anche nelle run che non scendono sotto fim.tol.
MESH=$1; NAME=$2; X=$3; Y=$4; Z=$5; shift 5
openCARP \
 -meshname $MESH -simID run_$NAME \
 -dream.solve 0 -dream.fim.tol 0.001 -dream.fim.max_iter 1000 \
 -dream.output.atfile act -tend 500 \
 -num_stim 1 -stim[0].elec.geom_type 1 \
 -stim[0].elec.p0[0] $X -stim[0].elec.p0[1] $Y -stim[0].elec.p0[2] $Z \
 -stim[0].elec.radius 3000 -stim[0].ptcl.start 0.0 \
 "$@"
