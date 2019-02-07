#!/bin/bash



name="ROOTFile_20190117_MuMu_for_SMP_17_010_with_All_SMP_17_010_Corrections_and_Updated_RoccoR_and_DoubleMuon_PD"
rename="ROOTFile_20190117_MuMu_for_SMP_17_010_with_All_SMP_17_010_Corrections_and_Updated_RoccoR_on_DoubleMuon"

types=(
"_1_0"
"_1_1"
"_2_0"
"_3_0"
"_4_0"
"_5_0"
"_6_0"
"_6_1"
"_7_0"
"_7_1"
)


for ii in ${types[@]}; do

	command="mv ${name}${ii}_0.root ${rename}${ii}_0.root"
	echo $command
#	eval $command
done

echo ""
echo "job is completed"


