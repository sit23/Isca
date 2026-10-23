#!/usr/bin/env bash
# Short correctness tests for the Isca-MOM6 coupling, one after another on 8+8 cores
# (the machine has 48 physical cores, shared with the two long runs).
source ~/project_specific_codes/isca_ocean_coupling/env.sh
cd "$(dirname "$0")"
mkdir -p logs
run() { echo "$(date) START $1"; env "${@:2}" python frierson_mom6_aquaplanet.py > logs/$1.log 2>&1; echo "$(date) END $1 exit $?"; }
C="NCORES=8 NCORES_OCEAN=8"
run T2_restart_60d   EXP_NAME=T2_restart_60d   DAYS=60 NMONTHS=1 $C
run T2_restart_2x30d EXP_NAME=T2_restart_2x30d DAYS=30 NMONTHS=2 $C
run T3_layout_8      EXP_NAME=T3_layout_8      DAYS=30 NMONTHS=1 $C
run T3_layout_4      EXP_NAME=T3_layout_4      DAYS=30 NMONTHS=1 NCORES=8 NCORES_OCEAN=4
run T6_regular2deg   EXP_NAME=T6_regular2deg   OCEAN_DLAT=2 DAYS=30 NMONTHS=24 $C
run T5_cpl3h         EXP_NAME=T5_cpl3h  COUPLING_PERIOD=10800 DAYS=30 NMONTHS=24 $C
run T5_cpl6h         EXP_NAME=T5_cpl6h  COUPLING_PERIOD=21600 DAYS=30 NMONTHS=24 $C
run T7_maxlat70      EXP_NAME=T7_maxlat70 MAX_LAT=70 DAYS=30 NMONTHS=24 $C
echo "$(date) START R1_slab20m"
REF_NAME=R1_slab20m_reference DEPTH=20 NMONTHS=120 python frierson_slab_oasis.py reference > logs/R1_slab20m.log 2>&1
echo "$(date) END R1_slab20m exit $?"
