#!/usr/bin/env bash
# Start the Earth-like Isca-MOM6 run (E1) once the short-test campaign has finished.
cd "$(dirname "$0")"
until grep -q "END R1_slab20m" logs/campaign.log; do sleep 120; done
source ~/project_specific_codes/isca_ocean_coupling/env.sh
echo "$(date) starting E1"
NCORES=8 NCORES_OCEAN=8 EXP_NAME=E1_realistic_mom6 NMONTHS=240 python realistic_mom6.py > E1_realistic_mom6.log 2>&1
echo "$(date) E1 finished exit $?"
