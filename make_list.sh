#!/bin/bash

#outfile="DsepInputFiles.list"
#outfile="DseAuInputFiles.list"
outfile="DISepInputFiles_q2_1000toINF.list"

get_files () {
  DID=$1
  RSE=$(rucio list-rules "$DID" | grep -E "BNL-XRD|EIC-XRD" | grep -v "JLAB-TAPE-SE" | awk '{print $5}' | head -n1)

  if [ -z "$RSE" ]; then
    echo "No BNL-XRD/EIC-XRD found for $DID" >&2
    return
  fi

  echo "Using RSE: $RSE for $DID" >&2
  rucio replica list file --pfns --rses "$RSE" "$DID"
}

echo "Running Ds filtered sample..."

#get_files "epic:/RECO/26.07.1/epic_craterlake/SIDIS/Ds_ABCONV/HFsim-PYTHIA/pythia8.316-1.1/ep/9x275/q2_1" > $outfile
#get_files "epic:/RECO/26.07.1/epic_craterlake/SIDIS/Ds_ABCONV/HFsim-BeAGLE/BeAGLE1.03.02-1.0/eAu/9x100/q2_1" > $outfile

get_files "epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x275/q2_1000toINF" > $outfile


#epic:/RECO/26.07.1/epic_craterlake/SIDIS/Ds_ABCONV/HFsim-PYTHIA/pythia8.316-1.1/ep/9x275/q2_1
#epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x275/q2_1to10
#epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x275/q2_10to100
#epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x275/q2_100to1000
#epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x275/q2_1000toINF

#epic:/RECO/26.07.1/epic_craterlake/SIDIS/Ds_ABCONV/HFsim-BeAGLE/BeAGLE1.03.02-1.0/eAu/9x100/q2_1
#epic:/RECO/26.07.1/epic_craterlake/DIS/BeAGLE1.03.02-3.1/eAu/9x100/q2_100to10000
#epic:/RECO/26.07.1/epic_craterlake/DIS/BeAGLE1.03.02-3.1/eAu/9x100/q2_10to100
#epic:/RECO/26.07.1/epic_craterlake/DIS/BeAGLE1.03.02-3.1/eAu/9x100/q2_1to10

#epic:/RECO/26.04.1/epic_craterlake/SIDIS/Lc_ABCONV/HFsim-PYTHIA/pythia8.312-2.0/ep/10x250/q2_1
#epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_10to100
