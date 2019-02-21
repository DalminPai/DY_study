#!/bin/bash

name="ROOTFile_20190208_MuMu_for_validation"

types_data=(
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

types_mc=(
"_10_0"
"_11_0"
"_11_1"
"_12_0"
"_13_0"
"_21_0"
"_22_0"
"_23_0"
"_24_0"
"_25_0"
"_31_0"
"_31_1"
"_32_0"
"_41_0"
)

command=""
chain=""

for ii in ${types_mc[@]}; do
	chain=${chain}"${name}${ii}_0.root "
done

command="hadd ${name}_on_MC.root "${chain}
echo $command
eval $command

echo ""
echo "job is completed"


