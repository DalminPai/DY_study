#!/bin/bash


############## COLOR CODE ##############
error(){
		echo -e "${RED}$1${NC}"
}
warn(){
#    echo -e "${YELLOW}$1${NC}"
		echo -e "${MAGENTA}$1${NC}"
}
info(){
		echo -e "${GREEN}$1${NC}"
}

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color
########################################

warn "You should run this script in ROOT available environment!"
warn "Please put your root file into \"setup4fr\", if you didn't."
##########################################################
# -- Usage : ./setup4fr.sh setup4fr/MY_ROOT_FILE.root -- #
##########################################################

rootfile=$(basename $1)
array_est=("denominator_pt" "denominator_pt_barrel" "denominator_pt_endcap" ## vs. pt
           "numerator_pt" "numerator_pt_barrel" "numerator_pt_endcap"       ## vs. pt
           "denominator_eta" "numerator_eta"                                ## vs. eta
           "denominator" "denominator_barrel" "denominator_endcap"          ## vs. iso
           "numerator" "numerator_barrel" "numerator_endcap"                ## vs. iso
           )
array_app=("histDijet1" "histDijet2" "histSameDijet1" "histSameDijet2"
           "fitDijet1" "fitDijet2" "fitSameDijet1" "fitSameDijet2"
           "rapDijet1" "rapDijet2" "rapSameDijet1" "rapSameDijet2"
           "histWJets1" "histWJets2" "histSameWJets1" "histSameWJets2"
           "fitWJets1" "fitWJets2" "fitSameWJets1" "fitSameWJets2"
           "rapWJets1" "rapWJets2" "rapSameWJets1" "rapSameWJets2"
           "histPass" "histFail"
           ## -- before applying FR -- ##
           "histDijet" "histSameDijet" "fitDijet" "fitSameDijet" "rapDijet" "rapSameDijet"
           "histWJets" "histSameWJets" "fitWJets" "fitSameWJets" "rapWJets" "rapSameWJets"
           )
array_var=${array_est[@]}
if [[ "$rootfile" =~ "fake" ]];then
	array_var=${array_app[@]}	
fi
array_name=""

cd setup4fr/

for i_var in ${array_var[@]}; do
	echo
	echo "work for $i_var ..."

	cmd="root -l -b -q 'MakeHist.cc("
	cmd=$cmd"\"$rootfile\", \"$i_var\")"
	cmd=$cmd"'"
#	echo $cmd
	eval $cmd

	array_name=$array_name" "$i_var"__"$rootfile
done

#echo "please wait until rootfiles become complete..."
#sleep 10s

echo
hadd Re_$rootfile $array_name

mv Re_$rootfile ../bin/histograms

rm -f $array_name

cd ../

echo
info "job is finished!"
