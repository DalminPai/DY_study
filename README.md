# DY_study
DY analysis with 2016 dataset<br>
Warning: Macros have to be used with a proper ntuple. Please check setup of the ntuple first.<br>
Warning: Usage may be outdated. It should be updated later.

	git clone -b v20190207_to_be_tidied https://DalminPai@github.com/DalminPai/DY_study.git

## exmacro
This is example macros for simple dimuon event selection.<br>
Usage:

	cd DY_study/exmacro
	root -l -b -q 'MuMu.C++((debug), (type))'
	* (debug) : You can choose 0 or 1. If you choose 1, then code will be run using only small number of events.
	* (type)  : You can choose 1 to 7 for Data, and some other numbers for MC. For defult setting, please see L38-54 in 'MuM.C'.
	* ex) root -l -b -q 'MuMu.C++(0, 1)'

## PileUp
This is macros for pile-up re-weighting.<br>
Usage:

	cd DY_study/PileUp
	python mkDataPileupRoot.py MY_JSON_FILE.json
	python mkPUReWeight.py MY_ERA MY_DATE

## EventSelection
This is macros for event selection.<br>
Usage:

	cd DY_study/EventSelection/SCRIPT
	source setup.sh
	./EventSelection.sh MY_MACRO.C -h

## BkgEst
This is macros for background estimation.<br>
Usage:

	cd DY_study/BkgEst/EMuMethod
	root -l -b -q 'emuCheck.cc' # -- Check e-mu distributions -- #
	root -l -b -q 'estimateBkg.cc' # -- Estimate not fake background -- #
