# DY_study
DY analysis with 2016 dataset

Warnning: Currently, the macro doesn't work because there is no ntuple in here. Please set about ntuple.

## exmacro
This is example macros for simple dimuon event selection.

Usage:

	cd DY_study/exmacro
	root -l -b -q 'MuMu.C++((debug), (type))'
	* (debug) : You can choose 0 or 1. If you choose 1, then code will be run using only small number of events.
	* (type)  : You can choose 1 to 7 for Data, and some other numbers for MC. For defult setting, please see L38-54 in 'MuM.C'.
	* ex) root -l -b -q 'MuMu.C++(0, 1)'

## PileUp
This is macros for pile-up re-weighting.

Usage:

	cd DY_study/PileUp
	python mkDataPileupRoot.py MY_JSON_FILE.json
	python mkPUReWeight.py MY_ERA MY_DATE

## EventSelection
This is macros for event selection.

Usage:

	cd DY_study/EventSelection/SCRIPT
	source setup.sh
	./EventSelection.sh MY_MACRO.C -h

## BkgEst
This is macros for background estimation.

Usage:

	cd DY_study/BkgEst/EMuMethod
	root -l -b -q 'emuCheck.cc' # -- Check e-mu distributions -- #
	root -l -b -q 'estimateBkg.cc' # -- Estimate not fake background -- #

