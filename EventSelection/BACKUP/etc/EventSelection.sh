#!/bin/bash

array_type=("[1] Run2016B" "[2] Run2016C" "[3] Run2016D" "[4] Run2016E" "[5] Run2016F" "[6] Run2016G" "[7] Run2016H"
			"[11] DYMuMu_M10to50" "[12] DYMuMu_M50toInf"
			"[21] ttbar" "[22] ttbarBackup"
			"[31] DYTauTau_M10to50" "[32] DYTauTau_M50toInf"
			"[41] VVnST"
			"[51] WJetsToLNu")

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

max_rmdr=0
opt=0
iter=0
debug=0
confirm=0
skip_sleep=0

if [ $# -eq 0 ]; then
	error "At least Macro C file name should be provided as the first argument!!!"
	exit 1;
fi

macro=$1
type_all=(1 2 3 4 5 6 7 11 12 21 22 31 32 41 51)
type_data=$(seq 1 7)
type_mc=(11 12 21 22 31 32 41 51)
type_usr=(41)
types=${type_usr[@]}

for arg in $*; do
	case $arg in
		$1)
			echo "Macro path is $macro"
			echo "---------------------------------------------------------------------------"
			;;
		all)
			echo "Running over all the types"
			types=${type_all[@]}
			;;
		mc)
			echo "Running over mc"
			types=${type_mc[@]}
			;;
		data)
			echo "Running over data"
			types=${type_data[@]}
			;;
		usr)
			echo "Running over user defined types"
			types=${type_usr[@]}
			;;
		-h)
			echo "Usage : $0 Macro.C (type) (-d) (-c) (-s)"
			echo ""
			echo "  (type)    it will run the macro for some specific types (ex: mc, data, or some user defined types)"
			echo "  (-d)      it will read small number of entries from ntuple for test and debug"
			echo "  (-c)      it will run jobs"
			echo "  (-s)      it will skip \"sleep\" for compiling time"
			exit 0;
			;;
		-d)
			echo "debug flag is SET and it will read small number of entries from ntuple"
			debug=1
			;;
		-c)
			echo "confirm flag is SET and it will run jobs"
			confirm=1
			;;
		-s)
			echo "skip flag is SET and it will skip \"sleep\" for compiling time"
			skip_sleep=1
			;;
		*)
			error "$arg is not supported option"
			;;
	esac
done

macro_filename=$(basename $macro)
macro_name=$(echo $macro_filename | cut -f1 -d".")
## -- # takes out prefix, % takes out suffix -- ##
pathname=${macro%$macro_filename}
cd $pathname

## -- print option (-c) -- ##
if [ $confirm -eq 0 ]; then
	warn "It is not running. Please add -c to confirm to run!!"
fi

#for i_type in 1 2 3 4 5 6 7 11 12 21 22 31 32 41 51; do
for i_type in ${types[@]}; do
#for i_type in 1 12; do

	## -- print type name -- ##
	for i_array in "${array_type[@]}"; do
		if [[ "$i_array" =~ "[$i_type]" ]]; then
			type_name=${i_array#"[$i_type] "}
			echo ""
			echo For $type_name...
		fi
	done;

	if [ $i_type -eq 1 -o $i_type -eq 6 -o $i_type -eq 7 -o $i_type -eq 11 -o $i_type -eq 31 -o $i_type -eq 51 ]; then
		max_rmdr=1
	fi

	for i_rmdr in $(seq 0 $max_rmdr); do

		## -- submit only if total jobs <= 10 -- ##
		N_jobs=`ps -ef | grep dmpai | grep -v "grep" | grep "root.exe" | wc -l`
		while [ $N_jobs -gt 10 ]; do
			sleep 10m
		done

		result_path=$DM_BASE_PATH/RESULT/$(echo $macro_name | cut -f1 -d"_")
		command="root -l -b -q \"$macro_filename+($debug ,$i_type, $i_rmdr, $opt)\" > $result_path/$macro_name-$i_type-$i_rmdr-$opt.log  2>&1 &"

		## -- option (-c) -- ##
		if [[ $confirm =~ ^0 ]]; then
			echo $command
		else
			if [ ! -d $result_path ]; then
				warn "$result_path directory does not exist!"
				mkdir -p $result_path
				info "$result_path directory has been created!!"
			fi

			eval $command
		fi

		## -- time for compiling -- ##
		if [ $iter -eq 0 -a $confirm -ne 0 -a $skip_sleep -eq 0 ]; then
			warn "This is the first job and let him take a time to compile if needed"
			sleep 30
		fi

		iter=$(echo "$iter + 1" | bc)
		if [ $iter -eq 1 ]; then
			info "$iter job is submitted until now..."
		else
			info "$iter jobs are submitted until now..."
		fi

	done

done

echo ""
echo "job is completed"

