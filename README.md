# DY_study

## exmacro

This is example macros for simple dimuon event selection.
Usage:
> git clone https://github.com/DalminPai/DY_study.git
> cd DY_study/exmacro
> root -l -b -q 'MuMu.C++((debug), (type))'
> * (debug) : You can choose 0 or 1. If you choose 1, then code will be run using only small number of events.
> * (type)  : You can choose 1 to 7 for Data, and some other numbers for MC. For defult setting, please see L38-54 in 'MuM.C'.
> ex) root -l -b -q 'MuMu.C++(0, 1)'
