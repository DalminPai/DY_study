#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TMath.h>
#include <THistPainter.h>
#include <TFormula.h>
#include <vector>
// -- Customized Analyzer for Drel-Yan Analysis -- //
//#include "../../HEADER/DYAnalyzer_TightID_PFIso_v20180521.h"
#include "../../HEADER/DYAnalyzer_TightID_PFIso_legacy.h" //legacy

static inline void loadBar(int x, int n, int r, int w);

// -- Electron Channel -- //
// -- Modified from muon channel's macro: 19 Jun. 2018 -- //
// -- Ntuples were re-produced due to energy scale correction problem: 28 Aug. 2018 -- //
// -- Bug was fixed for h_raw histograms: 31 Aug. 2018 -- //
// -- PVz reweighting was introduced: 03 Sep. 2018 -- //
// -- Gen-level distributions were added: 04 Oct. 2018 -- //
// -- Prefiring weight was tested: 05 Dec. 2018 -- //
void EE_legacy_wo_PVzWeight(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12")
{
	gROOT->SetBatch(kTRUE);
	//TString NtupleLocation = gSystem->Getenv("DM_DATA_PATH"); //it doesnt' work in electron channel. why?
	//TString BaseDir = gSystem->Getenv("DM_BASE_PATH"); //it doesn't work in electron channel. why?
	// -- legacy -- //
	//TString NtupleLocation = "/u/user/dmpai/SE_UserHome/_legacy_";
	TString NtupleLocation = "/u/user/dmpai/SE_UserHome/_legacy_with_prefiring_";
	TString BaseDir = "/u/user/dmpai/prime/dy_analysis/EventSelection";

	// -- Run2016 luminosity [/pb] -- //
	Double_t lumi = Lumi; //BtoH

	TString DataLocation, Type;
	// -- Data samples -- //
	/*if( type == 1 ) DataLocation = "DoubleEG_Run2016B";
	else if( type == 2 ) DataLocation = "DoubleEG_Run2016C";
	else if( type == 3 ) DataLocation = "DoubleEG_Run2016D";
	else if( type == 4 ) DataLocation = "DoubleEG_Run2016E";
	else if( type == 5 ) DataLocation = "DoubleEG_Run2016F";
	else if( type == 6 ) DataLocation = "DoubleEG_Run2016G";
	else if( type == 7 ) DataLocation = "DoubleEG_Run2016H";*/
	// -- legacy -- //
	if( type == 1 ) DataLocation = _DoubleEG_Run2016B_;
	else if( type == 2 ) DataLocation = _DoubleEG_Run2016C_;
	else if( type == 3 ) DataLocation = _DoubleEG_Run2016D_;
	else if( type == 4 ) DataLocation = _DoubleEG_Run2016E_;
	else if( type == 5 ) DataLocation = _DoubleEG_Run2016F_;
	else if( type == 6 ) DataLocation = _DoubleEG_Run2016G_;
	else if( type == 7 ) DataLocation = _DoubleEG_Run2016H_;
	// -- Signal MC samples -- //
	else if( type == 11 ) Type = "DYEE_M10to50";
	else if( type == 12 ) Type = "DYEE_M50toInf";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	else if( type == 31 ) Type = "DYTauTau_M10to50";
	else if( type == 32 ) Type = "DYTauTau_M50toInf";
	else if( type == 41 ) Type = "VVnST";
	else if( type == 51 ) Type = "WJetsToLNu";
	// -- Alternative signal MC samples -- //
	else if( type == 10 ) Type = "ZToEE_powheg";

	Bool_t isMC = kTRUE;
	if( type < 10  )
	{
		isMC = kFALSE;
		Type = "Data (" + DataLocation + ")";
		//if( type == 7 ) DataLocation += "ver2"; //for 80X re-miniaod
	}

	Double_t Div = 1;
	if( type == 1 || type == 6 || type == 7 || type == 11 || type == 31 || type == 51 ) Div = 2;

	TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	cout << "Type: " << Type << endl;

	TStopwatch totaltime;
	totaltime.Start();

	DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

	// -- Pile-up setup -- //
	analyzer->SetupPileUpReWeighting_80X( isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

	// -- PVz setup -- //
	analyzer->SetupPVzReWeighting_80X( isMC, "PVz_legacy.root" );

	// -- Efficiency SF setup -- //
	//if( isMC == kTRUE ) analyzer->SetupEfficiencyScaleFactor_electron();
	if( isMC == kTRUE ) analyzer->SetupEfficiencyScaleFactor_electron_legacy(); //legacy

	// -- Output ROOTFile -- //	
	//TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20180911_EE_legacy_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20180911_EE_legacy_PVz_reweighting_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20181004_EE_legacy_PVz_reweighting_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
	TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20181205_EE_legacy_with_PrefiringWeight_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
								+TString::Itoa(isTopPtReweighting,10)+".root";
	if( debug ) Output_ROOTFile = "test.root";
	TFile *f = new TFile(Output_ROOTFile, "recreate");

	// -- Each ntuple directory & corresponding Tags -- //
	vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
	//if( isMC == kTRUE ) analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);
	if( isMC == kTRUE ) analyzer->SetupMCsamples_16Legacy(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents); //legacy
	else
	{
		ntupleDirectory.push_back( "" );
		Tag.push_back( "Data" );
	}

	//Loop for all samples
	const Int_t Ntup = ntupleDirectory.size();
	for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
	{
		TStopwatch looptime;
		looptime.Start();

		cout << "\t<" << Tag[i_tup] << ">" << endl;

		TChain *chain = new TChain("recoTree/DYTree");
		//Set MC chain
		if( isMC == kTRUE )
		{
			//TString version = "v2.1";
			//if( type == 10 ) version = "v2.3";

			//if( remainder == 9999 )
			//	chain->Add(NtupleLocation+"/"+version+"/"+ntupleDirectory[i_tup]+"/*.root");
			//else
			//	for(Double_t ii=1; ii<=1500; ii++)
			//		if(ii - TMath::Floor(ii/Div) * Div == remainder)
			//			chain->Add(NtupleLocation+"/"+version+"/"+ntupleDirectory[i_tup]+"/*_"+TString::Itoa(ii,10)+".root");

			// -- legacy -- //
			if( remainder == 9999 )
				chain->Add(NtupleLocation+"/"+ntupleDirectory[i_tup]+"/*.root");
			else
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/"+ntupleDirectory[i_tup]+"/*_"+TString::Itoa(ii,10)+".root");
		}
		//Set Data chain
		else
		{
			if( remainder == 9999 )
			{
				//chain->Add(NtupleLocation+"/v2.0/"+DataLocation+"/*.root");
				//if(type==7) chain->Add(NtupleLocation+"/v2.0/DoubleEG_Run2016Hver3/*.root");

				// -- legacy -- //
				chain->Add(NtupleLocation+"/"+DataLocation+"/*.root");
				if(type==1) chain->Add(NtupleLocation+"/"+_DoubleEG_Run2016B_2_+"/*.root");
				else if(type==6) chain->Add(NtupleLocation+"/"+_DoubleEG_Run2016G_2_+"/*.root");
				else if(type==7) chain->Add(NtupleLocation+"/"+_DoubleEG_Run2016H_2_+"/*.root");
			}
			else
			{
				//for(Double_t ii=1; ii<=1500; ii++)
				//	if(ii - TMath::Floor(ii/Div) * Div == remainder)
				//		chain->Add(NtupleLocation+"/v2.0/"+DataLocation+"/*_"+TString::Itoa(ii,10)+".root");
				//if(type==7 && remainder==0) chain->Add(NtupleLocation+"/v2.0/DoubleEG_Run2016Hver3/*.root");

				// -- legacy -- //
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
					{
						chain->Add(NtupleLocation+"/"+DataLocation+"/*_"+TString::Itoa(ii,10)+".root");
						if(type==1) chain->Add(NtupleLocation+"/"+_DoubleEG_Run2016B_2_+"/*_"+TString::Itoa(ii,10)+".root");
						else if(type==6) chain->Add(NtupleLocation+"/"+_DoubleEG_Run2016G_2_+"/*_"+TString::Itoa(ii,10)+".root");
						else if(type==7) chain->Add(NtupleLocation+"/"+_DoubleEG_Run2016H_2_+"/*_"+TString::Itoa(ii,10)+".root");
					}
			}
		}

		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Electron();
		if( isMC == kTRUE )
		{
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			//ntuple->TurnOnBranches_GenOthers(); // for quarks
			ntuple->TurnOnBranches_prefiring(); // for prefiring weight
		}

		////////////////////////////
		// -- Making Histogram -- //
		////////////////////////////
		//with full corections
		TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_diPt = new TH1D("h_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_rapi = new TH1D("h_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_pT = new TH1D("h_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_leadPt = (TH1D*)h_pT->Clone("h_leadPt_"+Tag[i_tup]);
		TH1D *h_subPt = (TH1D*)h_pT->Clone("h_subPt_"+Tag[i_tup]);
		TH1D *h_eta = new TH1D("h_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_leadEta = (TH1D*)h_eta->Clone("h_leadEta_"+Tag[i_tup]);
		TH1D *h_subEta = (TH1D*)h_eta->Clone("h_subEta_"+Tag[i_tup]);
		TH1D *h_etaSC = (TH1D*)h_eta->Clone("h_etaSC_"+Tag[i_tup]);
		TH1D *h_leadEtaSC = (TH1D*)h_eta->Clone("h_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_subEtaSC = (TH1D*)h_eta->Clone("h_subEtaSC_"+Tag[i_tup]);
		TH1D *h_phi = new TH1D("h_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_leadPhi = (TH1D*)h_phi->Clone("h_leadPhi_"+Tag[i_tup]);
		TH1D *h_subPhi = (TH1D*)h_phi->Clone("h_subPhi_"+Tag[i_tup]);
		TH1D *h_PVz = new TH1D("h_PVz_"+Tag[i_tup], "", 100, -25, 25);
		TH1D *h_nVTX_noPU = new TH1D("h_nVTX_noPU_"+Tag[i_tup], "", 1000, 0, 1000); //to check pu reweighting
		TH1D *h_nVTX_withPU = (TH1D*)h_nVTX_noPU->Clone("h_nVTX_withPU_"+Tag[i_tup]); //to check pu reweighting

		// Gen-level
		TH1D *h_Gen_mass_fine = (TH1D*)h_mass_fine->Clone("h_Gen_mass_fine_"+Tag[i_tup]);
		TH1D *h_Gen_diPt = (TH1D*)h_diPt->Clone("h_Gen_diPt_"+Tag[i_tup]);
		TH1D *h_Gen_rapi = (TH1D*)h_rapi->Clone("h_Gen_rapi_"+Tag[i_tup]);
		TH1D *h_Gen_pT = (TH1D*)h_pT->Clone("h_Gen_pT_"+Tag[i_tup]);
		TH1D *h_Gen_leadPt = (TH1D*)h_pT->Clone("h_Gen_leadPt_"+Tag[i_tup]);
		TH1D *h_Gen_subPt = (TH1D*)h_pT->Clone("h_Gen_subPt_"+Tag[i_tup]);
		TH1D *h_Gen_eta = (TH1D*)h_eta->Clone("h_Gen_eta_"+Tag[i_tup]);
		TH1D *h_Gen_leadEta = (TH1D*)h_eta->Clone("h_Gen_leadEta_"+Tag[i_tup]);
		TH1D *h_Gen_subEta = (TH1D*)h_eta->Clone("h_Gen_subEta_"+Tag[i_tup]);
		TH1D *h_Gen_phi = (TH1D*)h_phi->Clone("h_Gen_phi_"+Tag[i_tup]);
		TH1D *h_Gen_leadPhi = (TH1D*)h_phi->Clone("h_Gen_leadPhi_"+Tag[i_tup]);
		TH1D *h_Gen_subPhi = (TH1D*)h_phi->Clone("h_Gen_subPhi_"+Tag[i_tup]);

		TH1D *h_restrict_mass_fine = new TH1D("h_restrict_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_restrict_diPt = new TH1D("h_restrict_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_restrict_rapi = new TH1D("h_restrict_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_restrict_pT = new TH1D("h_restrict_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_restrict_leadPt = (TH1D*)h_restrict_pT->Clone("h_restrict_leadPt_"+Tag[i_tup]);
		TH1D *h_restrict_subPt = (TH1D*)h_restrict_pT->Clone("h_restrict_subPt_"+Tag[i_tup]);
		TH1D *h_restrict_eta = new TH1D("h_restrict_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_restrict_leadEta = (TH1D*)h_restrict_eta->Clone("h_restrict_leadEta_"+Tag[i_tup]);
		TH1D *h_restrict_subEta = (TH1D*)h_restrict_eta->Clone("h_restrict_subEta_"+Tag[i_tup]);
		TH1D *h_restrict_etaSC = (TH1D*)h_restrict_eta->Clone("h_restrict_etaSC_"+Tag[i_tup]);
		TH1D *h_restrict_leadEtaSC = (TH1D*)h_restrict_eta->Clone("h_restrict_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_restrict_subEtaSC = (TH1D*)h_restrict_eta->Clone("h_restrict_subEtaSC_"+Tag[i_tup]);
		TH1D *h_restrict_phi = new TH1D("h_restrict_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_restrict_leadPhi = (TH1D*)h_restrict_phi->Clone("h_restrict_leadPhi_"+Tag[i_tup]);
		TH1D *h_restrict_subPhi = (TH1D*)h_restrict_phi->Clone("h_restrict_subPhi_"+Tag[i_tup]);
		TH1D *h_restrict_PVz = new TH1D("h_restrict_PVz_"+Tag[i_tup], "", 100, -25, 25);

		//without any corrections
		TH1D *h_raw_mass_fine = new TH1D("h_raw_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_raw_diPt = new TH1D("h_raw_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_raw_rapi = new TH1D("h_raw_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_raw_pT = new TH1D("h_raw_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_raw_leadPt = (TH1D*)h_raw_pT->Clone("h_raw_leadPt_"+Tag[i_tup]);
		TH1D *h_raw_subPt = (TH1D*)h_raw_pT->Clone("h_raw_subPt_"+Tag[i_tup]);
		TH1D *h_raw_eta = new TH1D("h_raw_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_raw_leadEta = (TH1D*)h_raw_eta->Clone("h_raw_leadEta_"+Tag[i_tup]);
		TH1D *h_raw_subEta = (TH1D*)h_raw_eta->Clone("h_raw_subEta_"+Tag[i_tup]);
		TH1D *h_raw_etaSC = (TH1D*)h_raw_eta->Clone("h_raw_etaSC_"+Tag[i_tup]);
		TH1D *h_raw_leadEtaSC = (TH1D*)h_raw_eta->Clone("h_raw_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_raw_subEtaSC = (TH1D*)h_raw_eta->Clone("h_raw_subEtaSC_"+Tag[i_tup]);
		TH1D *h_raw_phi = new TH1D("h_raw_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_raw_leadPhi = (TH1D*)h_raw_phi->Clone("h_raw_leadPhi_"+Tag[i_tup]);
		TH1D *h_raw_subPhi = (TH1D*)h_raw_phi->Clone("h_raw_subPhi_"+Tag[i_tup]);
		TH1D *h_raw_PVz = new TH1D("h_raw_PVz_"+Tag[i_tup], "", 100, -25, 25);
		TH1D *h_raw_nVTX_noPU = (TH1D*)h_nVTX_noPU->Clone("h_raw_nVTX_noPU_"+Tag[i_tup]); //to check pu reweighting
		TH1D *h_raw_nVTX_withPU = (TH1D*)h_nVTX_noPU->Clone("h_raw_nVTX_withPU_"+Tag[i_tup]); //to check pu reweighting

		TH1D *h_raw_restrict_mass_fine = new TH1D("h_raw_restrict_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_raw_restrict_diPt = new TH1D("h_raw_restrict_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_raw_restrict_rapi = new TH1D("h_raw_restrict_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_raw_restrict_pT = new TH1D("h_raw_restrict_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_raw_restrict_leadPt = (TH1D*)h_raw_restrict_pT->Clone("h_raw_restrict_leadPt_"+Tag[i_tup]);
		TH1D *h_raw_restrict_subPt = (TH1D*)h_raw_restrict_pT->Clone("h_raw_restrict_subPt_"+Tag[i_tup]);
		TH1D *h_raw_restrict_eta = new TH1D("h_raw_restrict_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_raw_restrict_leadEta = (TH1D*)h_raw_restrict_eta->Clone("h_raw_restrict_leadEta_"+Tag[i_tup]);
		TH1D *h_raw_restrict_subEta = (TH1D*)h_raw_restrict_eta->Clone("h_raw_restrict_subEta_"+Tag[i_tup]);
		TH1D *h_raw_restrict_etaSC = (TH1D*)h_raw_restrict_eta->Clone("h_raw_restrict_etaSC_"+Tag[i_tup]);
		TH1D *h_raw_restrict_leadEtaSC = (TH1D*)h_raw_restrict_eta->Clone("h_raw_restrict_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_raw_restrict_subEtaSC = (TH1D*)h_raw_restrict_eta->Clone("h_raw_restrict_subEtaSC_"+Tag[i_tup]);
		TH1D *h_raw_restrict_phi = new TH1D("h_raw_restrict_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_raw_restrict_leadPhi = (TH1D*)h_raw_restrict_phi->Clone("h_raw_restrict_leadPhi_"+Tag[i_tup]);
		TH1D *h_raw_restrict_subPhi = (TH1D*)h_raw_restrict_phi->Clone("h_raw_restrict_subPhi_"+Tag[i_tup]);
		TH1D *h_raw_restrict_PVz = new TH1D("h_raw_restrict_PVz_"+Tag[i_tup], "", 100, -25, 25);

		//add pileup reweighting
		TH1D *h_pu_mass_fine = new TH1D("h_pu_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_pu_diPt = new TH1D("h_pu_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_pu_rapi = new TH1D("h_pu_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_pu_pT = new TH1D("h_pu_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_pu_leadPt = (TH1D*)h_pu_pT->Clone("h_pu_leadPt_"+Tag[i_tup]);
		TH1D *h_pu_subPt = (TH1D*)h_pu_pT->Clone("h_pu_subPt_"+Tag[i_tup]);
		TH1D *h_pu_eta = new TH1D("h_pu_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_pu_leadEta = (TH1D*)h_pu_eta->Clone("h_pu_leadEta_"+Tag[i_tup]);
		TH1D *h_pu_subEta = (TH1D*)h_pu_eta->Clone("h_pu_subEta_"+Tag[i_tup]);
		TH1D *h_pu_etaSC = (TH1D*)h_pu_eta->Clone("h_pu_etaSC_"+Tag[i_tup]);
		TH1D *h_pu_leadEtaSC = (TH1D*)h_pu_eta->Clone("h_pu_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_pu_subEtaSC = (TH1D*)h_pu_eta->Clone("h_pu_subEtaSC_"+Tag[i_tup]);
		TH1D *h_pu_phi = new TH1D("h_pu_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_pu_leadPhi = (TH1D*)h_pu_phi->Clone("h_pu_leadPhi_"+Tag[i_tup]);
		TH1D *h_pu_subPhi = (TH1D*)h_pu_phi->Clone("h_pu_subPhi_"+Tag[i_tup]);
		TH1D *h_pu_PVz = new TH1D("h_pu_PVz_"+Tag[i_tup], "", 100, -25, 25);

		TH1D *h_pu_restrict_mass_fine = new TH1D("h_pu_restrict_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_pu_restrict_diPt = new TH1D("h_pu_restrict_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_pu_restrict_rapi = new TH1D("h_pu_restrict_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_pu_restrict_pT = new TH1D("h_pu_restrict_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_pu_restrict_leadPt = (TH1D*)h_pu_restrict_pT->Clone("h_pu_restrict_leadPt_"+Tag[i_tup]);
		TH1D *h_pu_restrict_subPt = (TH1D*)h_pu_restrict_pT->Clone("h_pu_restrict_subPt_"+Tag[i_tup]);
		TH1D *h_pu_restrict_eta = new TH1D("h_pu_restrict_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_pu_restrict_leadEta = (TH1D*)h_pu_restrict_eta->Clone("h_pu_restrict_leadEta_"+Tag[i_tup]);
		TH1D *h_pu_restrict_subEta = (TH1D*)h_pu_restrict_eta->Clone("h_pu_restrict_subEta_"+Tag[i_tup]);
		TH1D *h_pu_restrict_etaSC = (TH1D*)h_pu_restrict_eta->Clone("h_pu_restrict_etaSC_"+Tag[i_tup]);
		TH1D *h_pu_restrict_leadEtaSC = (TH1D*)h_pu_restrict_eta->Clone("h_pu_restrict_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_pu_restrict_subEtaSC = (TH1D*)h_pu_restrict_eta->Clone("h_pu_restrict_subEtaSC_"+Tag[i_tup]);
		TH1D *h_pu_restrict_phi = new TH1D("h_pu_restrict_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_pu_restrict_leadPhi = (TH1D*)h_pu_restrict_phi->Clone("h_pu_restrict_leadPhi_"+Tag[i_tup]);
		TH1D *h_pu_restrict_subPhi = (TH1D*)h_pu_restrict_phi->Clone("h_pu_restrict_subPhi_"+Tag[i_tup]);
		TH1D *h_pu_restrict_PVz = new TH1D("h_pu_restrict_PVz_"+Tag[i_tup], "", 100, -25, 25);

		//add energy scale correction
		TH1D *h_es_mass_fine = new TH1D("h_es_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_es_diPt = new TH1D("h_es_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_es_rapi = new TH1D("h_es_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_es_pT = new TH1D("h_es_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_es_leadPt = (TH1D*)h_es_pT->Clone("h_es_leadPt_"+Tag[i_tup]);
		TH1D *h_es_subPt = (TH1D*)h_es_pT->Clone("h_es_subPt_"+Tag[i_tup]);
		TH1D *h_es_eta = new TH1D("h_es_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_es_leadEta = (TH1D*)h_es_eta->Clone("h_es_leadEta_"+Tag[i_tup]);
		TH1D *h_es_subEta = (TH1D*)h_es_eta->Clone("h_es_subEta_"+Tag[i_tup]);
		TH1D *h_es_etaSC = (TH1D*)h_es_eta->Clone("h_es_etaSC_"+Tag[i_tup]);
		TH1D *h_es_leadEtaSC = (TH1D*)h_es_eta->Clone("h_es_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_es_subEtaSC = (TH1D*)h_es_eta->Clone("h_es_subEtaSC_"+Tag[i_tup]);
		TH1D *h_es_phi = new TH1D("h_es_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_es_leadPhi = (TH1D*)h_es_phi->Clone("h_es_leadPhi_"+Tag[i_tup]);
		TH1D *h_es_subPhi = (TH1D*)h_es_phi->Clone("h_es_subPhi_"+Tag[i_tup]);
		TH1D *h_es_PVz = new TH1D("h_es_PVz_"+Tag[i_tup], "", 100, -25, 25);

		TH1D *h_es_restrict_mass_fine = new TH1D("h_es_restrict_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_es_restrict_diPt = new TH1D("h_es_restrict_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_es_restrict_rapi = new TH1D("h_es_restrict_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_es_restrict_pT = new TH1D("h_es_restrict_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_es_restrict_leadPt = (TH1D*)h_es_restrict_pT->Clone("h_es_restrict_leadPt_"+Tag[i_tup]);
		TH1D *h_es_restrict_subPt = (TH1D*)h_es_restrict_pT->Clone("h_es_restrict_subPt_"+Tag[i_tup]);
		TH1D *h_es_restrict_eta = new TH1D("h_es_restrict_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_es_restrict_leadEta = (TH1D*)h_es_restrict_eta->Clone("h_es_restrict_leadEta_"+Tag[i_tup]);
		TH1D *h_es_restrict_subEta = (TH1D*)h_es_restrict_eta->Clone("h_es_restrict_subEta_"+Tag[i_tup]);
		TH1D *h_es_restrict_etaSC = (TH1D*)h_es_restrict_eta->Clone("h_es_restrict_etaSC_"+Tag[i_tup]);
		TH1D *h_es_restrict_leadEtaSC = (TH1D*)h_es_restrict_eta->Clone("h_es_restrict_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_es_restrict_subEtaSC = (TH1D*)h_es_restrict_eta->Clone("h_es_restrict_subEtaSC_"+Tag[i_tup]);
		TH1D *h_es_restrict_phi = new TH1D("h_es_restrict_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_es_restrict_leadPhi = (TH1D*)h_es_restrict_phi->Clone("h_es_restrict_leadPhi_"+Tag[i_tup]);
		TH1D *h_es_restrict_subPhi = (TH1D*)h_es_restrict_phi->Clone("h_es_restrict_subPhi_"+Tag[i_tup]);
		TH1D *h_es_restrict_PVz = new TH1D("h_es_restrict_PVz_"+Tag[i_tup], "", 100, -25, 25);

		Double_t SumWeight = 0, SumWeight_Separated = 0;

		Int_t NEvents = 10000; // test using small events
		if( !debug )  NEvents = chain->GetEntries();
		cout << "\t[Total Events: " << NEvents << "]" << endl;
		for(Int_t i=0; i<NEvents; i++)		
		{	
			loadBar(i+1, NEvents, 100, 100);
		
			ntuple->GetEvent(i);

			/////////////////////////////
			// -- Bring the weights -- //
			/////////////////////////////
			// -- Positive/Negative Gen-weights -- //
			Double_t GenWeight;
			ntuple->GENEvt_weight < 0 ? GenWeight = -1 : GenWeight = 1;
			SumWeight += GenWeight;

			// -- Pileup-Reweighting -- //
			Double_t PUWeight = 1;
			if( isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( ntuple->nPileUp );

			// -- PVz-Reweighting -- //
			Double_t PVzWeight = 1;
			//if( isMC == kTRUE ) PVzWeight = analyzer->PVzWeightValue_80X( ntuple->PVz );

			// -- Prefiring weight -- //
			Double_t PrefiringWeight = 1;
			if( isMC == kTRUE ) PrefiringWeight = ntuple->_prefiringweight;

			// -- efficiency weights -- //
			Double_t effweight = 1;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			// -- Separate ttbar samples -- //
			Bool_t GenFlag_top = kTRUE;
			//Bool_t GenFlag_top = kFALSE;
			//vector<GenOthers> GenTopCollection;
			//GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

			if( GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				SumWeight_Separated += GenWeight;

				// -- Top Pt Reweighting -- //
				/*if( isTopPtReweighting == 1 && Tag[i_tup].Contains("ttbar") )
				{
					GenOthers t1 = GenTopCollection[0];
					GenOthers t2 = GenTopCollection[1];

					Double_t SF1 = exp(0.0615 - 0.0005*(t1.Pt));
					Double_t SF2 = exp(0.0615 - 0.0005*(t2.Pt));
					GenWeight = GenWeight*sqrt(SF1*SF2);
				}*/
			}

			// -- Normalization -- //
			Double_t TotWeight = GenWeight;
			if( isMC == kTRUE ) TotWeight = (lumi*Xsec[i_tup]/nEvents[i_tup])*GenWeight;

			/////////////////////////////////////
			// -- Generator level selection -- //
			/////////////////////////////////////
			if( isMC == kTRUE && type < 20 && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				vector< GenLepton > GenLeptonCollection;
				Int_t NGenLeptons = ntuple->gnpair;
				for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
				{
					GenLepton genlep;
					genlep.FillFromNtuple(ntuple, i_gen);
					if( genlep.isElectron() && genlep.fromHardProcessFinalState )
						GenLeptonCollection.push_back( genlep );
				}

				if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 electrons from hard-process -- //
				{
					GenLepton genlep1, genlep2; // lead: 1, sub: 2
					if( GenLeptonCollection[0].Pt > GenLeptonCollection[1].Pt )
					{
						genlep1 = GenLeptonCollection[0];
						genlep2 = GenLeptonCollection[1];
					}
					else
					{
						genlep1 = GenLeptonCollection[1];
						genlep2 = GenLeptonCollection[0];
					}

					Bool_t isPassAcc_GenLepton = kFALSE;
					//isPassAcc_GenLepton = analyzer->isPassAccCondition_GenLepton(genlep1, genlep2);
					isPassAcc_GenLepton = analyzer->isPassAccCondition_GenLepton_ECALGAP(genlep1, genlep2);
					if( isPassAcc_GenLepton == kTRUE )
					{
						Double_t gen_M = (genlep1.Momentum + genlep2.Momentum).M();
						Double_t gen_Pt = (genlep1.Momentum + genlep2.Momentum).Pt();
						Double_t gen_rapi = (genlep1.Momentum + genlep2.Momentum).Rapidity();

						if( 60 < gen_M && gen_M < 120 ) // Z-peak
						{
							h_Gen_mass_fine->Fill( gen_M, TotWeight );
							h_Gen_diPt->Fill( gen_Pt, TotWeight );
							h_Gen_rapi->Fill( gen_rapi, TotWeight );
							h_Gen_pT->Fill( genlep1.Pt, TotWeight );
							h_Gen_pT->Fill( genlep2.Pt, TotWeight );
							h_Gen_eta->Fill( genlep1.eta, TotWeight );
							h_Gen_eta->Fill( genlep2.eta, TotWeight );
							h_Gen_phi->Fill( genlep1.phi, TotWeight );
							h_Gen_phi->Fill( genlep2.phi, TotWeight );

							h_Gen_leadPt->Fill( genlep1.Pt, TotWeight );
							h_Gen_leadEta->Fill( genlep1.eta, TotWeight );
							h_Gen_leadPhi->Fill( genlep1.phi, TotWeight );
							h_Gen_subPt->Fill( genlep2.Pt, TotWeight );
							h_Gen_subEta->Fill( genlep2.eta, TotWeight );
							h_Gen_subPhi->Fill( genlep2.phi, TotWeight );

						} //Z-peak mass cut

					} //isPassAcc_GenLepton

				} //the events containing 2 electrons from hard-process

			} //End of Gen-level selection


			TotWeight *= PrefiringWeight; // multiply prefiring weight


			////////////////////////////////
			// -- Reco level selection -- //
			////////////////////////////////
			Bool_t TriggerFlag = kFALSE;
			TriggerFlag = ntuple->isTriggered( analyzer->HLT );

			if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				vector< Electron > ElectronCollection;
				Int_t NLeptons = ntuple->Nelectrons;
				for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
				{
					Electron ele;
					ele.FillFromNtuple(ntuple, i_reco);

					ElectronCollection.push_back( ele );
				}	

				// -- Event Selection -- //
				vector< Electron > SelectedElectronCollection;
				Bool_t isPassEventSelection = kFALSE;
				isPassEventSelection = analyzer->EventSelection_ElectronChannel_Zpeak(ElectronCollection, ntuple, &SelectedElectronCollection);

				if( isPassEventSelection == kTRUE )
				{
					Electron ele1, ele2; // lead: 1, sub: 2
					if( SelectedElectronCollection[0].Pt > SelectedElectronCollection[1].Pt )
					{
						ele1 = SelectedElectronCollection[0];
						ele2 = SelectedElectronCollection[1];
					}
					else
					{
						ele1 = SelectedElectronCollection[1];
						ele2 = SelectedElectronCollection[0];
					}

					Double_t reco_M = (ele1.Momentum + ele2.Momentum).M();
					Double_t reco_Pt = (ele1.Momentum + ele2.Momentum).Pt();
					Double_t reco_rapi = (ele1.Momentum + ele2.Momentum).Rapidity();

					// -- Efficiency scale factor -- //
					if( isMC == kTRUE )
					{
						//effweight = analyzer->EfficiencySF_EventWeight_electron( ele1, ele2 );
						effweight = analyzer->EfficiencySF_EventWeight_electron_legacy( ele1, ele2 ); //legacy
					}

					h_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight * effweight );
					h_nVTX_noPU->Fill( ntuple->nVertices, TotWeight * PVzWeight *  effweight ); //to check pu reweighting
					h_nVTX_withPU->Fill( ntuple->nVertices, TotWeight * PVzWeight * PUWeight * effweight ); //to check pu reweighting
					h_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
					h_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight * effweight );
					h_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
					h_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
					h_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
					h_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
					h_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
					h_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
					h_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
					h_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
					h_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );

					h_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
					h_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
					h_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
					h_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
					h_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
					h_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
					h_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
					h_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );

					h_es_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight );
					h_es_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight );
					h_es_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight );
					h_es_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight );
					h_es_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight );
					h_es_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight );
					h_es_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight );
					h_es_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight );
					h_es_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight );
					h_es_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight );
					h_es_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight );
					h_es_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight );

					h_es_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight );
					h_es_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight );
					h_es_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight );
					h_es_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight );
					h_es_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight );
					h_es_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight );
					h_es_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight );
					h_es_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight );

					if( fabs(ntuple->PVz) < 2 )
					{
						h_restrict_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );

						h_restrict_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_restrict_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );

						h_es_restrict_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight );

						h_es_restrict_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight );
						h_es_restrict_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight );
					}
				} //End of event selection

				// -- Event Selection without electron energy scale corrections -- //
				vector< Electron > SelectedElectronCollection0;
				Bool_t isPassEventSelection0 = kFALSE;
				isPassEventSelection0 = analyzer->EventSelection_ElectronChannel0(ElectronCollection, ntuple, &SelectedElectronCollection0);

				if( isPassEventSelection0 == kTRUE )
				{
					Electron ele1, ele2; // lead: 1, sub: 2
					if( SelectedElectronCollection0[0].pTUnCorr > SelectedElectronCollection0[1].pTUnCorr )
					{
						ele1 = SelectedElectronCollection0[0];
						ele2 = SelectedElectronCollection0[1];
					}
					else
					{
						ele1 = SelectedElectronCollection0[1];
						ele2 = SelectedElectronCollection0[0];
					}

					Double_t reco_M = (ele1.Momentum_UnCorr + ele2.Momentum_UnCorr).M();
					Double_t reco_Pt = (ele1.Momentum_UnCorr + ele2.Momentum_UnCorr).Pt();
					Double_t reco_rapi = (ele1.Momentum_UnCorr + ele2.Momentum_UnCorr).Rapidity();

					h_raw_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight );
					h_raw_nVTX_noPU->Fill( ntuple->nVertices, TotWeight * PVzWeight ); //to check pu reweighting
					h_raw_nVTX_withPU->Fill( ntuple->nVertices, TotWeight * PVzWeight * PUWeight ); //to check pu reweighting
					h_raw_mass_fine->Fill( reco_M, TotWeight * PVzWeight );
					h_raw_diPt->Fill( reco_Pt, TotWeight * PVzWeight );
					h_raw_rapi->Fill( reco_rapi, TotWeight * PVzWeight );
					h_raw_pT->Fill( ele1.pTUnCorr, TotWeight * PVzWeight ); //
					h_raw_pT->Fill( ele2.pTUnCorr, TotWeight * PVzWeight );
					h_raw_eta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight );
					h_raw_eta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight );
					h_raw_etaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight );
					h_raw_etaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight );
					h_raw_phi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight );
					h_raw_phi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight );

					h_raw_leadPt->Fill( ele1.pTUnCorr, TotWeight * PVzWeight );
					h_raw_leadEta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight );
					h_raw_leadEtaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight );
					h_raw_leadPhi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight );
					h_raw_subPt->Fill( ele2.pTUnCorr, TotWeight * PVzWeight );
					h_raw_subEta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight );
					h_raw_subEtaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight );
					h_raw_subPhi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight );

					h_pu_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight );
					h_pu_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight );
					h_pu_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight );
					h_pu_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight );
					h_pu_pT->Fill( ele1.pTUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_pT->Fill( ele2.pTUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_eta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_eta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_etaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_etaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_phi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_phi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight * PUWeight );

					h_pu_leadPt->Fill( ele1.pTUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_leadEta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_leadEtaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_leadPhi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_subPt->Fill( ele2.pTUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_subEta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_subEtaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
					h_pu_subPhi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight * PUWeight );

					if( fabs(ntuple->PVz) < 2 )
					{
						h_raw_restrict_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight );
						h_raw_restrict_mass_fine->Fill( reco_M, TotWeight * PVzWeight );
						h_raw_restrict_diPt->Fill( reco_Pt, TotWeight * PVzWeight );
						h_raw_restrict_rapi->Fill( reco_rapi, TotWeight * PVzWeight );
						h_raw_restrict_pT->Fill( ele1.pTUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_pT->Fill( ele2.pTUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_eta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_eta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_etaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_etaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_phi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_phi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight );

						h_raw_restrict_leadPt->Fill( ele1.pTUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_leadEta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_leadEtaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_leadPhi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_subPt->Fill( ele2.pTUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_subEta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_subEtaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight );
						h_raw_restrict_subPhi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight );

						h_pu_restrict_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_pT->Fill( ele1.pTUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_pT->Fill( ele2.pTUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_eta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_eta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_etaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_etaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_phi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_phi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight * PUWeight );

						h_pu_restrict_leadPt->Fill( ele1.pTUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_leadEta->Fill( ele1.etaUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_leadEtaSC->Fill( ele1.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_leadPhi->Fill( ele1.phiUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_subPt->Fill( ele2.pTUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_subEta->Fill( ele2.etaUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_subEtaSC->Fill( ele2.etaSCUnCorr, TotWeight * PVzWeight * PUWeight );
						h_pu_restrict_subPhi->Fill( ele2.phiUnCorr, TotWeight * PVzWeight * PUWeight );
					}
				} //End of event selection without electron energy scale correction

			} //End of if( isTriggered )

		} //End of event iteration

		//with full corrections
		h_PVz->Write();
		h_nVTX_noPU->Write(); //to check pu reweighting
		h_nVTX_withPU->Write(); //to check pu reweighting
		h_mass_fine->Write();
		h_diPt->Write();
		h_rapi->Write();
		h_pT->Write();
		h_leadPt->Write();
		h_subPt->Write();
		h_eta->Write();
		h_leadEta->Write();
		h_subEta->Write();
		h_etaSC->Write();
		h_leadEtaSC->Write();
		h_subEtaSC->Write();
		h_phi->Write();
		h_leadPhi->Write();
		h_subPhi->Write();

		h_restrict_PVz->Write();
		h_restrict_mass_fine->Write();
		h_restrict_diPt->Write();
		h_restrict_rapi->Write();
		h_restrict_pT->Write();
		h_restrict_leadPt->Write();
		h_restrict_subPt->Write();
		h_restrict_eta->Write();
		h_restrict_leadEta->Write();
		h_restrict_subEta->Write();
		h_restrict_etaSC->Write();
		h_restrict_leadEtaSC->Write();
		h_restrict_subEtaSC->Write();
		h_restrict_phi->Write();
		h_restrict_leadPhi->Write();
		h_restrict_subPhi->Write();

		//without any corrections
		h_raw_PVz->Write();
		h_raw_nVTX_noPU->Write(); //to check pu reweighting
		h_raw_nVTX_withPU->Write(); //to check pu reweighting
		h_raw_mass_fine->Write();
		h_raw_diPt->Write();
		h_raw_rapi->Write();
		h_raw_pT->Write();
		h_raw_leadPt->Write();
		h_raw_subPt->Write();
		h_raw_eta->Write();
		h_raw_leadEta->Write();
		h_raw_subEta->Write();
		h_raw_etaSC->Write();
		h_raw_leadEtaSC->Write();
		h_raw_subEtaSC->Write();
		h_raw_phi->Write();
		h_raw_leadPhi->Write();
		h_raw_subPhi->Write();

		h_raw_restrict_PVz->Write();
		h_raw_restrict_mass_fine->Write();
		h_raw_restrict_diPt->Write();
		h_raw_restrict_rapi->Write();
		h_raw_restrict_pT->Write();
		h_raw_restrict_leadPt->Write();
		h_raw_restrict_subPt->Write();
		h_raw_restrict_eta->Write();
		h_raw_restrict_leadEta->Write();
		h_raw_restrict_subEta->Write();
		h_raw_restrict_etaSC->Write();
		h_raw_restrict_leadEtaSC->Write();
		h_raw_restrict_subEtaSC->Write();
		h_raw_restrict_phi->Write();
		h_raw_restrict_leadPhi->Write();
		h_raw_restrict_subPhi->Write();

		//add pile-up re-weighting
		h_pu_PVz->Write();
		h_pu_mass_fine->Write();
		h_pu_diPt->Write();
		h_pu_rapi->Write();
		h_pu_pT->Write();
		h_pu_leadPt->Write();
		h_pu_subPt->Write();
		h_pu_eta->Write();
		h_pu_leadEta->Write();
		h_pu_subEta->Write();
		h_pu_etaSC->Write();
		h_pu_leadEtaSC->Write();
		h_pu_subEtaSC->Write();
		h_pu_phi->Write();
		h_pu_leadPhi->Write();
		h_pu_subPhi->Write();

		h_pu_restrict_PVz->Write();
		h_pu_restrict_mass_fine->Write();
		h_pu_restrict_diPt->Write();
		h_pu_restrict_rapi->Write();
		h_pu_restrict_pT->Write();
		h_pu_restrict_leadPt->Write();
		h_pu_restrict_subPt->Write();
		h_pu_restrict_eta->Write();
		h_pu_restrict_leadEta->Write();
		h_pu_restrict_subEta->Write();
		h_pu_restrict_etaSC->Write();
		h_pu_restrict_leadEtaSC->Write();
		h_pu_restrict_subEtaSC->Write();
		h_pu_restrict_phi->Write();
		h_pu_restrict_leadPhi->Write();
		h_pu_restrict_subPhi->Write();

		//add energy scale correction
		h_es_PVz->Write();
		h_es_mass_fine->Write();
		h_es_diPt->Write();
		h_es_rapi->Write();
		h_es_pT->Write();
		h_es_leadPt->Write();
		h_es_subPt->Write();
		h_es_eta->Write();
		h_es_leadEta->Write();
		h_es_subEta->Write();
		h_es_etaSC->Write();
		h_es_leadEtaSC->Write();
		h_es_subEtaSC->Write();
		h_es_phi->Write();
		h_es_leadPhi->Write();
		h_es_subPhi->Write();

		h_es_restrict_PVz->Write();
		h_es_restrict_mass_fine->Write();
		h_es_restrict_diPt->Write();
		h_es_restrict_rapi->Write();
		h_es_restrict_pT->Write();
		h_es_restrict_leadPt->Write();
		h_es_restrict_subPt->Write();
		h_es_restrict_eta->Write();
		h_es_restrict_leadEta->Write();
		h_es_restrict_subEta->Write();
		h_es_restrict_etaSC->Write();
		h_es_restrict_leadEtaSC->Write();
		h_es_restrict_subEtaSC->Write();
		h_es_restrict_phi->Write();
		h_es_restrict_leadPhi->Write();
		h_es_restrict_subPhi->Write();

		if( isMC == kTRUE && type < 20 )
		{
			h_Gen_mass_fine->Write();
			h_Gen_diPt->Write();
			h_Gen_rapi->Write();
			h_Gen_pT->Write();
			h_Gen_eta->Write();
			h_Gen_phi->Write();
			h_Gen_leadPt->Write();
			h_Gen_subPt->Write();
			h_Gen_leadEta->Write();
			h_Gen_subEta->Write();
			h_Gen_leadPhi->Write();
			h_Gen_subPhi->Write();
		}

		printf("\tTotal sum of weights: %.1lf\n", SumWeight);
		printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
		if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", lumi*Xsec[i_tup]/nEvents[i_tup]);

		Double_t LoopRunTime = looptime.CpuTime();
		cout << "\tLoop RunTime(" << Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

	} //end of i_tup iteration

	Double_t TotalRunTime = totaltime.CpuTime();
	cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

	TTimeStamp ts_end;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
    	cout << endl;

    if ( x % (n/r +1) != 0 ) return;

 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";
 
    for (int x=c; x<w; x++) cout << " ";
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
	cout << "]\r" << flush;
}

