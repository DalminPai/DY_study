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
#include "../../HEADER/DYAnalyzer_TightID_PFIso_v20180521.h"

static inline void loadBar(int x, int n, int r, int w);

// -- Electron Channel -- //
// -- Modified from muon channel's macro : 19 Jun. 2018 -- //
// -- Prefiring weight was tested: 05 Dec. 2018 -- //
// -- Test _prefiringweightup(down): 18 Dec. 2018 -- //
void EE_PVz_20181012(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "Ele23Ele12")
{
	gROOT->SetBatch(kTRUE);
	//TString NtupleLocation = gSystem->Getenv("DM_DATA_PATH");
	//TString BaseDir = gSystem->Getenv("DM_BASE_PATH");
	//TString NtupleLocation = "/scratch/DYntuple";
	//TString BaseDir = "/home/dmpai/dy_analysis/EventSelection";
	TString NtupleLocation = "/u/user/dmpai/SE_UserHome/_prime_/DYntuple";
	TString BaseDir = "/u/user/dmpai/prime/dy_analysis/EventSelection";
	
	// -- Choose prefiringweight -- //
	TString type_prefiringweight = "";
	type_prefiringweight = "up";
	//type_prefiringweight = "down";


	// -- Run2016 luminosity [/pb] -- //
	Double_t lumi = Lumi; //BtoH

	TString DataLocation, Type;
	// -- Data samples -- //
	if( type == 1 ) DataLocation = "DoubleEG_Run2016B";
	else if( type == 2 ) DataLocation = "DoubleEG_Run2016C";
	else if( type == 3 ) DataLocation = "DoubleEG_Run2016D";
	else if( type == 4 ) DataLocation = "DoubleEG_Run2016E";
	else if( type == 5 ) DataLocation = "DoubleEG_Run2016F";
	else if( type == 6 ) DataLocation = "DoubleEG_Run2016G";
	else if( type == 7 ) DataLocation = "DoubleEG_Run2016H";
	// -- Signal MC samples -- //
	else if( type == 11 ) Type = "DYEE_M10to50";
	else if( type == 12 ) Type = "DYEE_M50toInf";
	//else if( type == 12 ) Type = "DYEE_M50to200";
	//else if( type == 13 ) Type = "DYEE_M200toInf";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	//else if( type == 21 ) Type = "ttbar_Mto700";
	//else if( type == 22 ) Type = "ttbarBackup_Mto700";
	//else if( type == 23 ) Type = "ttbar_M700toInf";
	else if( type == 31 ) Type = "DYTauTau_M10to50";
	else if( type == 32 ) Type = "DYTauTau_M50toInf";
	else if( type == 41 ) Type = "VVnST";
	else if( type == 51 ) Type = "WJetsToLNu";

	Bool_t isMC = kTRUE;
	if( type < 10  )
	{
		isMC = kFALSE;
		Type = "Data (" + DataLocation + ")";
		if( type == 7 ) DataLocation += "ver2";
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
	analyzer->SetupPVzReWeighting_80X( isMC, "PVz.root" );

	// -- Efficiency SF setup -- //
	if( isMC == kTRUE ) analyzer->SetupEfficiencyScaleFactor_electron();

	// -- Output ROOTFile -- //	
	//TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20181012_EE_PVz_reweighting_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20181205_EE_PVz_reweighting_with_PrefiringWeight_"
	TString Output_ROOTFile = BaseDir+"/RESULT/EE/ROOTFile_20181218_EE_PVz_reweighting_with_PrefiringWeight"+type_prefiringweight+"_"
					+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"+TString::Itoa(isTopPtReweighting,10)+".root";
	if( debug ) Output_ROOTFile = "test.root";
	TFile *f = new TFile(Output_ROOTFile, "recreate");

	// -- Each ntuple directory & corresponding Tags -- //
	vector<TString> ntupleDirectory; vector<TString> Tag; vector<Double_t> Xsec; vector<Double_t> nEvents;
	if( isMC == kTRUE ) analyzer->SetupMCsamples_Moriond17(Type, &ntupleDirectory, &Tag, &Xsec, &nEvents);
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
			TString version = "v2.5"; // for prefiring weight

			if( remainder == 9999 )
				chain->Add(NtupleLocation+"/"+version+"/"+ntupleDirectory[i_tup]+"/*.root");
			else
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/"+version+"/"+ntupleDirectory[i_tup]+"/*_"+TString::Itoa(ii,10)+".root");
		}
		//Set Data chain
		else
		{
			if( remainder == 9999 )
			{
				chain->Add(NtupleLocation+"/v2.0/"+DataLocation+"/*.root");
				if(type==7) chain->Add(NtupleLocation+"/v2.0/DoubleEG_Run2016Hver3/*.root");
			}
			else
			{
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/v2.0/"+DataLocation+"/*_"+TString::Itoa(ii,10)+".root");
				if(type==7 && remainder==0) chain->Add(NtupleLocation+"/v2.0/DoubleEG_Run2016Hver3/*.root");
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
		// -- Reco-level -- //
		TH1D *h_mass = new TH1D("h_mass_"+Tag[i_tup], "", 43, massbins);
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

		// Dilepton pT cut: 30GeV
		TH1D *h_PtCut_mass = (TH1D*)h_mass->Clone("h_PtCut_mass_"+Tag[i_tup]);
		TH1D *h_PtCut_mass_fine = (TH1D*)h_mass_fine->Clone("h_PtCut_mass_fine_"+Tag[i_tup]);
		TH1D *h_PtCut_diPt = (TH1D*)h_diPt->Clone("h_PtCut_diPt_"+Tag[i_tup]);
		TH1D *h_PtCut_rapi = (TH1D*)h_rapi->Clone("h_PtCut_rapi_"+Tag[i_tup]);
		TH1D *h_PtCut_pT = (TH1D*)h_pT->Clone("h_PtCut_pT_"+Tag[i_tup]);
		TH1D *h_PtCut_leadPt = (TH1D*)h_pT->Clone("h_PtCut_leadPt_"+Tag[i_tup]);
		TH1D *h_PtCut_subPt = (TH1D*)h_pT->Clone("h_PtCut_subPt_"+Tag[i_tup]);
		TH1D *h_PtCut_eta = (TH1D*)h_eta->Clone("h_PtCut_eta_"+Tag[i_tup]);
		TH1D *h_PtCut_leadEta = (TH1D*)h_eta->Clone("h_PtCut_leadEta_"+Tag[i_tup]);
		TH1D *h_PtCut_subEta = (TH1D*)h_eta->Clone("h_PtCut_subEta_"+Tag[i_tup]);
		TH1D *h_PtCut_etaSC = (TH1D*)h_eta->Clone("h_PtCut_etaSC_"+Tag[i_tup]);
		TH1D *h_PtCut_leadEtaSC = (TH1D*)h_eta->Clone("h_PtCut_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_PtCut_subEtaSC = (TH1D*)h_eta->Clone("h_PtCut_subEtaSC_"+Tag[i_tup]);
		TH1D *h_PtCut_phi = (TH1D*)h_phi->Clone("h_PtCut_phi_"+Tag[i_tup]);
		TH1D *h_PtCut_leadPhi = (TH1D*)h_phi->Clone("h_PtCut_leadPhi_"+Tag[i_tup]);
		TH1D *h_PtCut_subPhi = (TH1D*)h_phi->Clone("h_PtCut_subPhi_"+Tag[i_tup]);
		TH1D *h_PtCut_PVz = (TH1D*)h_PVz->Clone("h_PtCut_PVz_"+Tag[i_tup]);

		// 2D measurement
		TH1D *h_M20to30_rapi = (TH1D*)h_rapi->Clone("h_M20to30_rapi_"+Tag[i_tup]);
		TH1D *h_M30to45_rapi = (TH1D*)h_rapi->Clone("h_M30to45_rapi_"+Tag[i_tup]);
		TH1D *h_M45to60_rapi = (TH1D*)h_rapi->Clone("h_M45to60_rapi_"+Tag[i_tup]);
		TH1D *h_M60to120_rapi = (TH1D*)h_rapi->Clone("h_M60to120_rapi_"+Tag[i_tup]);
		TH1D *h_M120to200_rapi = (TH1D*)h_rapi->Clone("h_M120to200_rapi_"+Tag[i_tup]);
		TH1D *h_M200to1500_rapi = (TH1D*)h_rapi->Clone("h_M200to1500_rapi_"+Tag[i_tup]);

		// -- Reco-level in Z-peak -- //
		TH1D *h_Zpeak_mass = (TH1D*)h_mass->Clone("h_Zpeak_mass_"+Tag[i_tup]);
		TH1D *h_Zpeak_mass_fine = (TH1D*)h_mass_fine->Clone("h_Zpeak_mass_fine_"+Tag[i_tup]);
		TH1D *h_Zpeak_diPt = (TH1D*)h_diPt->Clone("h_Zpeak_diPt_"+Tag[i_tup]);
		TH1D *h_Zpeak_rapi = (TH1D*)h_rapi->Clone("h_Zpeak_rapi_"+Tag[i_tup]);
		TH1D *h_Zpeak_pT = (TH1D*)h_pT->Clone("h_Zpeak_pT_"+Tag[i_tup]);
		TH1D *h_Zpeak_leadPt = (TH1D*)h_pT->Clone("h_Zpeak_leadPt_"+Tag[i_tup]);
		TH1D *h_Zpeak_subPt = (TH1D*)h_pT->Clone("h_Zpeak_subPt_"+Tag[i_tup]);
		TH1D *h_Zpeak_eta = (TH1D*)h_eta->Clone("h_Zpeak_eta_"+Tag[i_tup]);
		TH1D *h_Zpeak_leadEta = (TH1D*)h_eta->Clone("h_Zpeak_leadEta_"+Tag[i_tup]);
		TH1D *h_Zpeak_subEta = (TH1D*)h_eta->Clone("h_Zpeak_subEta_"+Tag[i_tup]);
		TH1D *h_Zpeak_etaSC = (TH1D*)h_eta->Clone("h_Zpeak_etaSC_"+Tag[i_tup]);
		TH1D *h_Zpeak_leadEtaSC = (TH1D*)h_eta->Clone("h_Zpeak_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_Zpeak_subEtaSC = (TH1D*)h_eta->Clone("h_Zpeak_subEtaSC_"+Tag[i_tup]);
		TH1D *h_Zpeak_phi = (TH1D*)h_phi->Clone("h_Zpeak_phi_"+Tag[i_tup]);
		TH1D *h_Zpeak_leadPhi = (TH1D*)h_phi->Clone("h_Zpeak_leadPhi_"+Tag[i_tup]);
		TH1D *h_Zpeak_subPhi = (TH1D*)h_phi->Clone("h_Zpeak_subPhi_"+Tag[i_tup]);
		TH1D *h_Zpeak_PVz = (TH1D*)h_PVz->Clone("h_Zpeak_PVz_"+Tag[i_tup]);

		// Dilepton pT cut: 30GeV
		TH1D *h_Zpeak_PtCut_mass = (TH1D*)h_mass->Clone("h_Zpeak_PtCut_mass_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_mass_fine = (TH1D*)h_mass_fine->Clone("h_Zpeak_PtCut_mass_fine_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_diPt = (TH1D*)h_diPt->Clone("h_Zpeak_PtCut_diPt_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_rapi = (TH1D*)h_rapi->Clone("h_Zpeak_PtCut_rapi_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_pT = (TH1D*)h_pT->Clone("h_Zpeak_PtCut_pT_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_leadPt = (TH1D*)h_pT->Clone("h_Zpeak_PtCut_leadPt_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_subPt = (TH1D*)h_pT->Clone("h_Zpeak_PtCut_subPt_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_eta = (TH1D*)h_eta->Clone("h_Zpeak_PtCut_eta_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_leadEta = (TH1D*)h_eta->Clone("h_Zpeak_PtCut_leadEta_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_subEta = (TH1D*)h_eta->Clone("h_Zpeak_PtCut_subEta_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_etaSC = (TH1D*)h_eta->Clone("h_Zpeak_PtCut_etaSC_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_leadEtaSC = (TH1D*)h_eta->Clone("h_Zpeak_PtCut_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_subEtaSC = (TH1D*)h_eta->Clone("h_Zpeak_PtCut_subEtaSC_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_phi = (TH1D*)h_phi->Clone("h_Zpeak_PtCut_phi_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_leadPhi = (TH1D*)h_phi->Clone("h_Zpeak_PtCut_leadPhi_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_subPhi = (TH1D*)h_phi->Clone("h_Zpeak_PtCut_subPhi_"+Tag[i_tup]);
		TH1D *h_Zpeak_PtCut_PVz = (TH1D*)h_PVz->Clone("h_Zpeak_PtCut_PVz_"+Tag[i_tup]);

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

		//add efficiency scale factor
		TH1D *h_eff_mass_fine = new TH1D("h_eff_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_eff_diPt = new TH1D("h_eff_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_eff_rapi = new TH1D("h_eff_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_eff_pT = new TH1D("h_eff_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_eff_leadPt = (TH1D*)h_eff_pT->Clone("h_eff_leadPt_"+Tag[i_tup]);
		TH1D *h_eff_subPt = (TH1D*)h_eff_pT->Clone("h_eff_subPt_"+Tag[i_tup]);
		TH1D *h_eff_eta = new TH1D("h_eff_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_eff_leadEta = (TH1D*)h_eff_eta->Clone("h_eff_leadEta_"+Tag[i_tup]);
		TH1D *h_eff_subEta = (TH1D*)h_eff_eta->Clone("h_eff_subEta_"+Tag[i_tup]);
		TH1D *h_eff_etaSC = (TH1D*)h_eff_eta->Clone("h_eff_etaSC_"+Tag[i_tup]);
		TH1D *h_eff_leadEtaSC = (TH1D*)h_eff_eta->Clone("h_eff_leadEtaSC_"+Tag[i_tup]);
		TH1D *h_eff_subEtaSC = (TH1D*)h_eff_eta->Clone("h_eff_subEtaSC_"+Tag[i_tup]);
		TH1D *h_eff_phi = new TH1D("h_eff_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_eff_leadPhi = (TH1D*)h_eff_phi->Clone("h_eff_leadPhi_"+Tag[i_tup]);
		TH1D *h_eff_subPhi = (TH1D*)h_eff_phi->Clone("h_eff_subPhi_"+Tag[i_tup]);
		TH1D *h_eff_PVz = (TH1D*)h_PVz->Clone("h_eff_PVz_"+Tag[i_tup]);


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
			if( isMC == kTRUE ) PVzWeight = analyzer->PVzWeightValue_80X( ntuple->PVz );

			// -- Prefiring weight -- //
			Double_t PrefiringWeight = 1;
			if( isMC == kTRUE )
			{
				if( type_prefiringweight == "" ) PrefiringWeight = ntuple->_prefiringweight;
				else if( type_prefiringweight == "up" ) PrefiringWeight = ntuple->_prefiringweightup;
				else if( type_prefiringweight == "down" ) PrefiringWeight = ntuple->_prefiringweightdown;
			}

			// -- efficiency weights -- //
			Double_t effweight = 1;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			// -- Separate ttbar samples -- //
			//Bool_t GenFlag_top = kTRUE;
			Bool_t GenFlag_top = kFALSE;
			vector<GenOthers> GenTopCollection;
			GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

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
				isPassEventSelection = analyzer->EventSelection_ElectronChannel(ElectronCollection, ntuple, &SelectedElectronCollection);

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
						effweight = analyzer->EfficiencySF_EventWeight_electron( ele1, ele2 );
					}

					h_mass->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
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
					h_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight * effweight );

					// Dilepton pT cut : 30GeV
					if( 30 < reco_Pt )
					{
						h_PtCut_mass->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_PtCut_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight * effweight );
					}

					// 2D measurement
					if( 20 < reco_M && reco_M < 30 ) h_M20to30_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
					else if( 30 < reco_M && reco_M < 45 ) h_M30to45_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
					else if( 45 < reco_M && reco_M < 60 ) h_M45to60_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
					else if( 60 < reco_M && reco_M < 120 ) h_M60to120_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
					else if( 120 < reco_M && reco_M < 200 ) h_M120to200_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
					else if( 200 < reco_M && reco_M < 1500 ) h_M200to1500_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );

					if( 60 < reco_M && reco_M < 120 )
					{
						h_Zpeak_mass->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );
						h_Zpeak_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight * effweight );

						// Dilepton pT cut : 30GeV
						if( 30 < reco_Pt )
						{
							h_Zpeak_PtCut_mass->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_mass_fine->Fill( reco_M, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_diPt->Fill( reco_Pt, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_rapi->Fill( reco_rapi, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_pT->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_pT->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_eta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_eta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_etaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_etaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_phi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_phi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_leadPt->Fill( ele1.Pt, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_leadEta->Fill( ele1.eta, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_leadEtaSC->Fill( ele1.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_leadPhi->Fill( ele1.phi, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_subPt->Fill( ele2.Pt, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_subEta->Fill( ele2.eta, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_subEtaSC->Fill( ele2.etaSC, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_subPhi->Fill( ele2.phi, TotWeight * PVzWeight * PUWeight * effweight );
							h_Zpeak_PtCut_PVz->Fill( ntuple->PVz, TotWeight * PVzWeight * PUWeight * effweight );
						}

						h_eff_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_eff_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
						h_eff_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
						h_eff_pT->Fill( ele1.Pt, TotWeight * PUWeight * effweight );
						h_eff_pT->Fill( ele2.Pt, TotWeight * PUWeight * effweight );
						h_eff_eta->Fill( ele1.eta, TotWeight * PUWeight * effweight );
						h_eff_eta->Fill( ele2.eta, TotWeight * PUWeight * effweight );
						h_eff_etaSC->Fill( ele1.etaSC, TotWeight * PUWeight * effweight );
						h_eff_etaSC->Fill( ele2.etaSC, TotWeight * PUWeight * effweight );
						h_eff_phi->Fill( ele1.phi, TotWeight * PUWeight * effweight );
						h_eff_phi->Fill( ele2.phi, TotWeight * PUWeight * effweight );
						h_eff_leadPt->Fill( ele1.Pt, TotWeight * PUWeight * effweight );
						h_eff_leadEta->Fill( ele1.eta, TotWeight * PUWeight * effweight );
						h_eff_leadEtaSC->Fill( ele1.etaSC, TotWeight * PUWeight * effweight );
						h_eff_leadPhi->Fill( ele1.phi, TotWeight * PUWeight * effweight );
						h_eff_subPt->Fill( ele2.Pt, TotWeight * PUWeight * effweight );
						h_eff_subEta->Fill( ele2.eta, TotWeight * PUWeight * effweight );
						h_eff_subEtaSC->Fill( ele2.etaSC, TotWeight * PUWeight * effweight );
						h_eff_subPhi->Fill( ele2.phi, TotWeight * PUWeight * effweight );
						h_eff_PVz->Fill( ntuple->PVz, TotWeight * PUWeight * effweight );

						h_es_mass_fine->Fill( reco_M, TotWeight * PUWeight );
						h_es_diPt->Fill( reco_Pt, TotWeight * PUWeight );
						h_es_rapi->Fill( reco_rapi, TotWeight * PUWeight );
						h_es_pT->Fill( ele1.Pt, TotWeight * PUWeight );
						h_es_pT->Fill( ele2.Pt, TotWeight * PUWeight );
						h_es_eta->Fill( ele1.eta, TotWeight * PUWeight );
						h_es_eta->Fill( ele2.eta, TotWeight * PUWeight );
						h_es_etaSC->Fill( ele1.etaSC, TotWeight * PUWeight );
						h_es_etaSC->Fill( ele2.etaSC, TotWeight * PUWeight );
						h_es_phi->Fill( ele1.phi, TotWeight * PUWeight );
						h_es_phi->Fill( ele2.phi, TotWeight * PUWeight );
						h_es_leadPt->Fill( ele1.Pt, TotWeight * PUWeight );
						h_es_leadEta->Fill( ele1.eta, TotWeight * PUWeight );
						h_es_leadEtaSC->Fill( ele1.etaSC, TotWeight * PUWeight );
						h_es_leadPhi->Fill( ele1.phi, TotWeight * PUWeight );
						h_es_subPt->Fill( ele2.Pt, TotWeight * PUWeight );
						h_es_subEta->Fill( ele2.eta, TotWeight * PUWeight );
						h_es_subEtaSC->Fill( ele2.etaSC, TotWeight * PUWeight );
						h_es_subPhi->Fill( ele2.phi, TotWeight * PUWeight );
					} // End of Z-peak

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

					if( 60 < reco_M && reco_M < 120 )
					{
						h_raw_mass_fine->Fill( reco_M, TotWeight );
						h_raw_diPt->Fill( reco_Pt, TotWeight );
						h_raw_rapi->Fill( reco_rapi, TotWeight );
						h_raw_pT->Fill( ele1.Pt, TotWeight );
						h_raw_pT->Fill( ele2.Pt, TotWeight );
						h_raw_eta->Fill( ele1.eta, TotWeight );
						h_raw_eta->Fill( ele2.eta, TotWeight );
						h_raw_etaSC->Fill( ele1.etaSC, TotWeight );
						h_raw_etaSC->Fill( ele2.etaSC, TotWeight );
						h_raw_phi->Fill( ele1.phi, TotWeight );
						h_raw_phi->Fill( ele2.phi, TotWeight );
						h_raw_leadPt->Fill( ele1.Pt, TotWeight );
						h_raw_leadEta->Fill( ele1.eta, TotWeight );
						h_raw_leadEtaSC->Fill( ele1.etaSC, TotWeight );
						h_raw_leadPhi->Fill( ele1.phi, TotWeight );
						h_raw_subPt->Fill( ele2.Pt, TotWeight );
						h_raw_subEta->Fill( ele2.eta, TotWeight );
						h_raw_subEtaSC->Fill( ele2.etaSC, TotWeight );
						h_raw_subPhi->Fill( ele2.phi, TotWeight );

						h_pu_mass_fine->Fill( reco_M, TotWeight * PUWeight );
						h_pu_diPt->Fill( reco_Pt, TotWeight * PUWeight );
						h_pu_rapi->Fill( reco_rapi, TotWeight * PUWeight );
						h_pu_pT->Fill( ele1.Pt, TotWeight * PUWeight );
						h_pu_pT->Fill( ele2.Pt, TotWeight * PUWeight );
						h_pu_eta->Fill( ele1.eta, TotWeight * PUWeight );
						h_pu_eta->Fill( ele2.eta, TotWeight * PUWeight );
						h_pu_etaSC->Fill( ele1.etaSC, TotWeight * PUWeight );
						h_pu_etaSC->Fill( ele2.etaSC, TotWeight * PUWeight );
						h_pu_phi->Fill( ele1.phi, TotWeight * PUWeight );
						h_pu_phi->Fill( ele2.phi, TotWeight * PUWeight );
						h_pu_leadPt->Fill( ele1.Pt, TotWeight * PUWeight );
						h_pu_leadEta->Fill( ele1.eta, TotWeight * PUWeight );
						h_pu_leadEtaSC->Fill( ele1.etaSC, TotWeight * PUWeight );
						h_pu_leadPhi->Fill( ele1.phi, TotWeight * PUWeight );
						h_pu_subPt->Fill( ele2.Pt, TotWeight * PUWeight );
						h_pu_subEta->Fill( ele2.eta, TotWeight * PUWeight );
						h_pu_subEtaSC->Fill( ele2.etaSC, TotWeight * PUWeight );
						h_pu_subPhi->Fill( ele2.phi, TotWeight * PUWeight );
					} //End of Z-peak

				} //End of event selection without electron energy scale correction

			} //End of if( isTriggered )

		} //End of event iteration

		h_mass->Write();
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
		h_PVz->Write();

		h_PtCut_mass->Write();
		h_PtCut_mass_fine->Write();
		h_PtCut_diPt->Write();
		h_PtCut_rapi->Write();
		h_PtCut_pT->Write();
		h_PtCut_leadPt->Write();
		h_PtCut_subPt->Write();
		h_PtCut_eta->Write();
		h_PtCut_leadEta->Write();
		h_PtCut_subEta->Write();
		h_PtCut_etaSC->Write();
		h_PtCut_leadEtaSC->Write();
		h_PtCut_subEtaSC->Write();
		h_PtCut_phi->Write();
		h_PtCut_leadPhi->Write();
		h_PtCut_subPhi->Write();
		h_PtCut_PVz->Write();

		h_Zpeak_mass->Write();
		h_Zpeak_mass_fine->Write();
		h_Zpeak_diPt->Write();
		h_Zpeak_rapi->Write();
		h_Zpeak_pT->Write();
		h_Zpeak_leadPt->Write();
		h_Zpeak_subPt->Write();
		h_Zpeak_eta->Write();
		h_Zpeak_leadEta->Write();
		h_Zpeak_subEta->Write();
		h_Zpeak_etaSC->Write();
		h_Zpeak_leadEtaSC->Write();
		h_Zpeak_subEtaSC->Write();
		h_Zpeak_phi->Write();
		h_Zpeak_leadPhi->Write();
		h_Zpeak_subPhi->Write();
		h_Zpeak_PVz->Write();

		h_Zpeak_PtCut_mass->Write();
		h_Zpeak_PtCut_mass_fine->Write();
		h_Zpeak_PtCut_diPt->Write();
		h_Zpeak_PtCut_rapi->Write();
		h_Zpeak_PtCut_pT->Write();
		h_Zpeak_PtCut_leadPt->Write();
		h_Zpeak_PtCut_subPt->Write();
		h_Zpeak_PtCut_eta->Write();
		h_Zpeak_PtCut_leadEta->Write();
		h_Zpeak_PtCut_subEta->Write();
		h_Zpeak_PtCut_etaSC->Write();
		h_Zpeak_PtCut_leadEtaSC->Write();
		h_Zpeak_PtCut_subEtaSC->Write();
		h_Zpeak_PtCut_phi->Write();
		h_Zpeak_PtCut_leadPhi->Write();
		h_Zpeak_PtCut_subPhi->Write();
		h_Zpeak_PtCut_PVz->Write();

		h_M20to30_rapi->Write();
		h_M30to45_rapi->Write();
		h_M45to60_rapi->Write();
		h_M60to120_rapi->Write();
		h_M120to200_rapi->Write();
		h_M200to1500_rapi->Write();

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

		h_eff_mass_fine->Write();
		h_eff_diPt->Write();
		h_eff_rapi->Write();
		h_eff_pT->Write();
		h_eff_leadPt->Write();
		h_eff_subPt->Write();
		h_eff_eta->Write();
		h_eff_leadEta->Write();
		h_eff_subEta->Write();
		h_eff_etaSC->Write();
		h_eff_leadEtaSC->Write();
		h_eff_subEtaSC->Write();
		h_eff_phi->Write();
		h_eff_leadPhi->Write();
		h_eff_subPhi->Write();
		h_eff_PVz->Write();

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

