import sys
import os

Era = sys.argv[1] ## for example, B
Date = sys.argv[2] ## for example, 20170817

for i in range(55, 76):
    minBiasXsec = str(i)

    command = "root -l -b -q 'PileUp.C(\""+Era+"\", \""+Date+"_"+minBiasXsec+"mb\")' >&PUReWeight_Run2016"+Era+"_v"+Date+"_"+minBiasXsec+"mb.log&"
    os.system(command)

    print('-' * 50)
    print('PU reweight ROOT file is finished')
    print('\tEra : Run2016'+Era)
    print('\tDate : '+Date)
    print('\tminBiasXsec : '+minBiasXsec)

