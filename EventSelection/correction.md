# Instructions on how to apply corrections (with 2016 dataset)
This is brief instructions on how to apply necessary corrections for 2016 dataset in DY analysis.<br>
The list of corrections is as follows:<br>
* PU reweighting
* Rochester correction
* Energy scale correction
* Efficiency SF
* PVz reweighting
* L1 prefiring reweighting
* Top pt reweighting

## PU reweighting
It should be applied to MC in both channels.<br>
You need a root file which contains weight histogram. (The weight is as a function of number of pileup.)<br>
* How to use in macro

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L116
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L324-L325
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L525

## Rochester correction
It should be applied to both data and MC in moun channel.<br>
You need a proper package which can be downloaded from [1].<br>
* How to use in macro

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L24
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L119-L120
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L467-L484

* Reference

	[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon

## Energy scale correction
It should be applied to both data and MC in electron channel.<br>
It is already applied when ntuple was made.<br>
If you do not want to use this correction, then you can use electron variables which are tagged with "UnCorr".<br>
* Reference

	https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer#!ECAL%20scale%20and%20resolution%20corre

## Efficiency SF
It should be applied to MC in both channels.<br>
You need proper SF root files which contains SF histograms or graphs.<br>
* How to use in macro
  * muon channel

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L123-L129
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L332
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L514-L521
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L525
  * electron channel

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L104
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L323
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L480-L483
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L487

## PVz reweighting
It should be applied to MC in electron channel.<br>
You need a root file which contains weight histogram. (The weight is as a function of PVz.)<br>

	WARN: It is not a standard correction in CMS.
	It was invented by Dalmin and reported to EGM conveners via email, but it was not officially approved yet by EGM POG.
	The PVz was reweighted in MC by multiplying the ratio between data and MC in PVz distribution, which was obtained after event selection.

* How to use in macro

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L101
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L315-L316
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_for_validation.C#L487

* Reference

	https://indico.cern.ch/event/774783/contributions/3219702/attachments/1763588/2862200/DY_working_meeting_20181130.pdf (S21-26)
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/EE/EE_PVz_20181012.C

## L1 prefiring reweighting
It should be applied to MC in both channels.<br>
The weight was calculated and saved in ntuple when it was made.<br>
So, you can easily apply it by calling the weight from ntuple.<br>
* How to use in macro

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L192
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L328-L329
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L525

* Reference

	https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
	http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_418_v10.pdf (section 13)

## Top pt reweighting
It can be applied to ttbar(MC) in both channels.<br>
The weight can be calculated using momentum of top quark, which can be called from ntuple.<br>

	WARN: From early 2018, Dalmin has not used this correction, because an additional study for bias and uncertainty is needed.

* How to use in macro

	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L193-L194
	https://github.com/DalminPai/DY_study/blob/v20190221/EventSelection/MACRO/MuMu/MuMu_for_validation.C#L349-L357

* Reference

	https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
	https://indico.cern.ch/event/680797/contributions/2789439/attachments/1559884/2455119/Preliminary_Result_of_Emu_Method_20171116.pdf (S11-12)


