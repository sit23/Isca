#!/usr/bin/env bash
# Short correctness tests for the Isca-MOM6 coupling, run in two lanes (16 cores each).
source ~/project_specific_codes/isca_ocean_coupling/env.sh
cd "$(dirname "$0")"
run() { echo "$(date) START $1"; env "${@:2}" python frierson_mom6_aquaplanet.py > logs/$1.log 2>&1; echo "$(date) END $1 exit $?"; }
mkdir -p logs
lane1() {
  run T2_restart_60d EXP_NAME=T2_restart_60d DAYS=60 NMONTHS=1 NCORES=8 NCORES_OCEAN=8
  run T2_restart_2x30d EXP_NAME=T2_restart_2x30d DAYS=30 NMONTHS=2 NCORES=8 NCORES_OCEAN=8
  run T3_layout_8 EXP_NAME=T3_layout_8 DAYS=30 NMONTHS=1 NCORES=8 NCORES_OCEAN=8
  run T3_layout_16 EXP_NAME=T3_layout_16 DAYS=30 NMONTHS=1 NCORES=8 NCORES_OCEAN=16
  run T5_cpl3h EXP_NAME=T5_cpl3h COUPLING_PERIOD=10800 DAYS=30 NMONTHS=24 NCORES=8 NCORES_OCEAN=8
  run T7_maxlat70 EXP_NAME=T7_maxlat70 MAX_LAT=70 DAYS=30 NMONTHS=24 NCORES=8 NCORES_OCEAN=8
}
lane2() {
  run T6_regular2deg EXP_NAME=T6_regular2deg OCEAN_DLAT=2 DAYS=30 NMONTHS=24 NCORES=8 NCORES_OCEAN=8
  run T5_cpl6h EXP_NAME=T5_cpl6h COUPLING_PERIOD=21600 DAYS=30 NMONTHS=24 NCORES=8 NCORES_OCEAN=8
  echo "$(date) START R1_slab20m"
  REF_NAME=R1_slab20m_reference DEPTH=20 NMONTHS=120 python frierson_slab_oasis.py reference > logs/R1_slab20m.log 2>&1
  echo "$(date) END R1_slab20m exit $?"
}
lane1 > logs/lane1.log 2>&1 &
lane2 > logs/lane2.log 2>&1 &
wait
