#!/bin/bash


############## COLOR CODE ##############
error(){
    echo -e "${RED}$1${NC}"
}
warn(){
    echo -e "${YELLOW}$1${NC}"
}
info(){
    echo -e "${GREEN}$1${NC}"
}

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color
########################################

if [ $HOSTNAME != 'prime' -a $HOSTNAME != 'tamsa2' ]; then
    error "This script is allowed only in tamsa2 and prime servers!!"
    exit;
fi

##################
# -- Set path -- #
##################
export DM_DATA_PATH=""
export DM_ANALYZER_PATH=""
export DM_DY_TOOLS_PATH=""
export DM_ROOT_PATH=""
export DM_OUTPUT_PATH=""
if [ $HOSTNAME == "prime" ]; then
	DM_DATA_PATH="/scratch/DYntuple"
	DM_ANALYZER_PATH="/home/dmpai/dy_analysis/EventSelection/HEADER"
	DM_DY_TOOLS_PATH="/scratch/DYanalysis"
	DM_OUTPUT_PATH="/home/dmpai/dy_analysis/EventSelection/RESULT"
elif [ $HOSTNAME == "tamsa2.snu.ac.kr" ]; then
	DM_DATA_PATH="/data9/DATA/DYntuple"
	DM_ANALYZER_PATH="/data9/Users/dmpai/DY/Ntuple80X/HEADER"
	DM_DY_TOOLS_PATH="/data5/Users/dmpai"
fi

################################
# -- ROOT setup (in server) -- #
################################
export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

## ROOT thisroot.sh PATH
DM_ROOT_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_20
cd $DM_ROOT_PATH/src

## cmsenv
eval `scramv1 runtime -sh`

cd ~

################
# -- Result -- #
################
echo "================ environment ================"
## DATA
if [ -z $DM_DATA_PATH ]; then
	warn "DATA is N/A in this machine"
else
	info "DM_DATA_PATH: $DM_DATA_PATH"
fi
## ANALYZER
if [ -z $DM_ANALYZER_PATH ]; then
	warn "ANALYZER is N/A in this machine"
else
	info "DM_ANALYZER_PATH: $DM_ANALYZER_PATH"
fi
## TOOLS
if [ -z $DM_DY_TOOLS_PATH ]; then
	warn "TOOLS are N/A in this machine"
else
	info "DM_DY_TOOLS_PATH: $DM_DY_TOOLS_PATH"
fi
## ROOT
if [ -z $DM_ROOT_PATH ]; then
	warn "ROOT is N/A in this machine"
else
	info "DM_ROOT_PATH: $DM_ROOT_PATH"
fi
## OUTPUT
if [ -z $DM_OUTPUT_PATH ]; then
	warn "OUTPUT path is not set in this machine"
else
	info "DM_OUTPUT_PATH: $DM_OUTPUT_PATH"
fi
echo "============================================="
echo "Setup is finished!"
