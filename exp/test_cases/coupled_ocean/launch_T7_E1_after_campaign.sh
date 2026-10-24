#!/usr/bin/env bash
# After the short-test campaign: finish T7 (killed by mistake at month 6), then start E1.
cd "$(dirname "$0")"
until grep -q "END R1_slab20m" logs/campaign.log; do sleep 120; done
source ~/project_specific_codes/isca_ocean_coupling/env.sh
echo "$(date) resuming T7 at month 7"
SKIP_COMPILE=1 FIRST_MONTH=7 EXP_NAME=T7_maxlat70 MAX_LAT=70 DAYS=30 NMONTHS=24 NCORES=8 NCORES_OCEAN=8 \
    python frierson_mom6_aquaplanet.py > logs/T7_maxlat70_resume.log 2>&1
echo "$(date) T7 finished exit $?"
echo "$(date) starting E1"
NCORES=8 NCORES_OCEAN=8 EXP_NAME=E1_realistic_mom6 NMONTHS=240 python realistic_mom6.py > E1_realistic_mom6.log 2>&1
echo "$(date) E1 finished exit $?"
