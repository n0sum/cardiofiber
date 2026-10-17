#!/bin/bash
# uso: ./run_eik.sh NOME X_um Y_um Z_um
NAME=$1; X=$2; Y=$3; Z=$4
openCARP \
 -meshname LA_um -simID run_$NAME \
 -dream.solve 0 -dream.fim.tol 0.001 -dream.fim.max_iter 10000000 \
 -dream.output.atfile act -tend 500 \
 -num_phys_regions 1 \
 -phys_region[0].name "Eikonal domain" -phys_region[0].ptype 2 -phys_region[0].num_IDs 7 \
 -phys_region[0].ID[0] 1 -phys_region[0].ID[1] 2 -phys_region[0].ID[2] 3 -phys_region[0].ID[3] 4 \
 -phys_region[0].ID[4] 5 -phys_region[0].ID[5] 6 -phys_region[0].ID[6] 7 \
 -num_gregions 1 -gregion[0].num_IDs 7 \
 -gregion[0].ID[0] 1 -gregion[0].ID[1] 2 -gregion[0].ID[2] 3 -gregion[0].ID[3] 4 \
 -gregion[0].ID[4] 5 -gregion[0].ID[5] 6 -gregion[0].ID[6] 7 \
 -gregion[0].dream.vel_l 1050 -gregion[0].dream.vel_t 525 -gregion[0].dream.vel_n 525 \
 -num_stim 1 -stim[0].elec.geom_type 1 \
 -stim[0].elec.p0[0] $X -stim[0].elec.p0[1] $Y -stim[0].elec.p0[2] $Z \
 -stim[0].elec.radius 3000 -stim[0].ptcl.start 0.0
