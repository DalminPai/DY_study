import sys
import os

MyAnalysisJSON = sys.argv[1]
JSONname = MyAnalysisJSON.split(".")[0]
JSONname = JSONname.split("/")[-1]
#Date = '20180320'
Date = '20180604'

## -- minBiasXsec from 55 to 75 -- ##
#for i in range(55, 76):
#    minBiasXsec = str(i)

#    command = 'pileupCalc.py -i '+MyAnalysisJSON+' --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec '+minBiasXsec+'000 --maxPileupBin 75 --numPileupBins 75  DataPileupHistogram_'+JSONname+'_v'+Date+'_'+minBiasXsec+'mb.root >&'+JSONname+'_v'+Date+'_'+minBiasXsec+'.log&'
#    os.system(command)

## -- minBiasXsec from 64.1 to 64.9 -- ##
for i in range(1, 10):
    x1000 = '64'+str(i)+'00'
    minBiasXsec = '64p'+str(i)

    command = 'pileupCalc.py -i '+MyAnalysisJSON+' --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec '+x1000+' --maxPileupBin 75 --numPileupBins 75  DataPileupHistogram_'+JSONname+'_v'+Date+'_'+minBiasXsec+'mb.root >&'+JSONname+'_v'+Date+'_'+minBiasXsec+'.log&'
#    print(command)
    os.system(command)


    print('-' * 50)
    print('DataPileup ROOT file is finished')
    print('\tMyAnalysisJSON : '+MyAnalysisJSON)
    print('\tminBiasXsec    : '+minBiasXsec)
