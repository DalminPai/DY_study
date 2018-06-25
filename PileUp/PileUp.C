void PileUp(TString era, TString Date_minBiasXsec)
{
	Double_t nPileUp[75] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74};

	Double_t puWeight[75] = {1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05};

///////////////////////////////////////////////////////////////////

	TString Baselocation = "./ROOT/Run2016"+era+"/";

	TString JSONname;
	if( era == "B" ) JSONname = "Cert_272007-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_eraB";
	else if( era == "C" ) JSONname = "Cert_275657-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_eraC";
	else if( era == "D" ) JSONname = "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_eraD";
	else if( era == "E" ) JSONname = "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_eraE";
	else if( era == "F" ) JSONname = "Cert_277772-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_eraF";
	else if( era == "G" ) JSONname = "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_eraG";
	else if( era == "H" ) JSONname = "Cert_280919-284044_13TeV_PromptReco_Collisions16_JSON_eraH";
	else if( era == "BtoH" ) JSONname = "284044_13TeV_23Sep2016ReReco_Collisions16_JSON";
	else cout << "You typed wrong JSONname..." << endl;

	TString inputfile = Baselocation+"DataPileupHistogram_"+JSONname+"_v"+Date_minBiasXsec+".root";

	TFile f_input(inputfile, "read");
	TH1D *h_data = (TH1D*)f_input.Get("pileup");
	Double_t norm = 1/h_data->Integral();
	h_data->Scale(norm);

	// https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py
	TH1D *h_mc = new TH1D("h_mc", "", 75, 0, 75);
	h_mc->FillN(75, nPileUp, puWeight, 1);

////////////////////////////////////////////////////////////////////////

	TH1D *h_PUReWeights = new TH1D("h_PUReWeights", "", 75, 0, 75);
	h_PUReWeights->Divide(h_data, h_mc, 1, 1);

///////////////////////////////////////////////////////////////////////////

	TFile f_output("ROOTFile_PUReWeight_80X_Run2016"+era+"_v"+Date_minBiasXsec+".root", "recreate");
	h_data->Write();
	h_mc->Write();
	h_PUReWeights->Write();

	cout << "Job is finished." << endl;
}
