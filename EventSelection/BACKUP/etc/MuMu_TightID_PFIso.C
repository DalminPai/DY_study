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
// -- for Rochester Muon momentum correction -- //
#include "../../TOOL/RoccoR/RoccoR.cc"
// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "../../HEADER/DYAnalyzer_TightID_PFIso_v20180521.h"

static inline void loadBar(int x, int n, int r, int w);

// -- Muon Channel -- //
// -- Off high mass ttbar samples, Top Pt Reweighting, and TuneP variables : 12? Mar 2018 -- //
// -- "TRandom3()" -> "TRandom3(0)" : 13 Mar 2018 -- //
// -- Use "PileUpWeightValue_80X" only for MC : 13 Mar 2018 -- //
// -- Add "h_eta_before_PUCorr" and "h_eta_before_RoccoR" : 13 Mar 2018 -- // 
// -- Move "massbins" into DYAnalyzer : 20 Mar 2018 -- //
// -- Correct Rochester correction about updating momentum : 20 Mar 2018 -- //
// -- Add "LeadingMuonEtaSF" : 17 Apr. 2018 -- //
// -- Modifying "LeadingMuonEtaSF" using TH2D (not finished) : 21 Apr. 2018 -- //
// -- Modify the code to fit in with prime server : 17 May. 2018 -- //
// -- Update "LeadingMuonEtaSF" as "LeadEtaCorr" : 21 May. 2018 -- //
void MuMu_TightID_PFIso(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24")
{
	gROOT->SetBatch(kTRUE);
	TString NtupleLocation = gSystem->Getenv("DM_DATA_PATH");
	TString BaseDir= gSystem->Getenv("DM_BASE_PATH");

	// -- Run2016 luminosity [/pb] -- //
	Double_t lumi = Lumi; //BtoH

	TString DataLocation, Type;
	// -- Data samples -- //
	if( type == 1 ) DataLocation = "SingleMuon_Run2016B";
	else if( type == 2 ) DataLocation = "SingleMuon_Run2016C";
	else if( type == 3 ) DataLocation = "SingleMuon_Run2016D";
	else if( type == 4 ) DataLocation = "SingleMuon_Run2016E";
	else if( type == 5 ) DataLocation = "SingleMuon_Run2016F";
	else if( type == 6 ) DataLocation = "SingleMuon_Run2016G";
	else if( type == 7 ) DataLocation = "SingleMuon_Run2016H";
	// -- Signal MC samples -- //
	else if( type == 11 ) Type = "DYMuMu_M10to50";
	else if( type == 12 ) Type = "DYMuMu_M50toInf";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
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

	// -- Rochester setup -- //
	TRandom3 *r1 = new TRandom3(0);
	RoccoR rc("../../TOOL/RoccoR/rcdata.2016.v3");

	// -- Efficiency SF setup -- //
	if( isMC == kTRUE )
	{
		analyzer->SetupEfficiencyScaleFactor_BtoF();
		analyzer->SetupEfficiencyScaleFactor_GtoH();
	}

	// -- Output ROOTFile -- //	
	TFile *f = new TFile(BaseDir+"/RESULT/MuMu/ROOTFile_MuMu_TightID_PFIso_v20180521_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
						+TString::Itoa(isTopPtReweighting,10)+".root", "RECREATE");
	//TFile *f = new TFile("test.root", "recreate");

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
			if( remainder == 9999 )
				chain->Add(NtupleLocation+"/v2.1/"+ntupleDirectory[i_tup]+"/*.root");
			else
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/v2.1/"+ntupleDirectory[i_tup]+"/*_"+TString::Itoa(ii,10)+".root");
		}
		//Set Data chain
		else
		{
			if( remainder == 9999 )
			{
				chain->Add(NtupleLocation+"/v2.0/"+DataLocation+"/*.root");
				if(type==7) chain->Add(NtupleLocation+"/v2.0/SingleMuon_Run2016Hver3/*.root");
			}
			else
			{
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/v2.0/"+DataLocation+"/*_"+TString::Itoa(ii,10)+".root");
				if(type==7 && remainder==0) chain->Add(NtupleLocation+"/v2.0/SingleMuon_Run2016Hver3/*.root");
			}
		}

		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		if( isMC == kTRUE )
		{
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			//ntuple->TurnOnBranches_GenOthers(); // for quarks
		}

		////////////////////////////
		// -- Making Histogram -- //
		////////////////////////////
		TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_diPt = new TH1D("h_diPt_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_rapi = new TH1D("h_rapi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_pT = new TH1D("h_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_leadPt = (TH1D*)h_pT->Clone("h_leadPt_"+Tag[i_tup]);
		TH1D *h_subPt = (TH1D*)h_pT->Clone("h_subPt_"+Tag[i_tup]);
		TH1D *h_eta = new TH1D("h_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_leadEta = (TH1D*)h_eta->Clone("h_leadEta_"+Tag[i_tup]);
		TH1D *h_subEta = (TH1D*)h_eta->Clone("h_subEta_"+Tag[i_tup]);
		TH1D *h_phi = new TH1D("h_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_leadPhi = (TH1D*)h_phi->Clone("h_leadPhi_"+Tag[i_tup]);
		TH1D *h_subPhi = (TH1D*)h_phi->Clone("h_subPhi_"+Tag[i_tup]);
		TH2D *h_leadEta_subEta = new TH2D("h_leadEta_subEta_"+Tag[i_tup], "", 1000, -5, 5, 1000, -5, 5);

		TH1D *h_BtoF_mass_fine = (TH1D*)h_mass_fine->Clone("h_BtoF_mass_fine_"+Tag[i_tup]);
		TH1D *h_BtoF_diPt = (TH1D*)h_diPt->Clone("h_BtoF_diPt_"+Tag[i_tup]);
		TH1D *h_BtoF_rapi = (TH1D*)h_rapi->Clone("h_BtoF_rapi_"+Tag[i_tup]);
		TH1D *h_BtoF_pT = (TH1D*)h_pT->Clone("h_BtoF_pT_"+Tag[i_tup]);
		TH1D *h_BtoF_leadPt = (TH1D*)h_pT->Clone("h_BtoF_leadPt_"+Tag[i_tup]);
		TH1D *h_BtoF_subPt = (TH1D*)h_pT->Clone("h_BtoF_subPt_"+Tag[i_tup]);
		TH1D *h_BtoF_eta = (TH1D*)h_eta->Clone("h_BtoF_eta_"+Tag[i_tup]);
		TH1D *h_BtoF_leadEta = (TH1D*)h_eta->Clone("h_BtoF_leadEta_"+Tag[i_tup]);
		TH1D *h_BtoF_subEta = (TH1D*)h_eta->Clone("h_BtoF_subEta_"+Tag[i_tup]);
		TH1D *h_BtoF_phi = (TH1D*)h_phi->Clone("h_BtoF_phi_"+Tag[i_tup]);
		TH1D *h_BtoF_leadPhi = (TH1D*)h_phi->Clone("h_BtoF_leadPhi_"+Tag[i_tup]);
		TH1D *h_BtoF_subPhi = (TH1D*)h_phi->Clone("h_BtoF_subPhi_"+Tag[i_tup]);

		TH1D *h_GtoH_mass_fine = (TH1D*)h_mass_fine->Clone("h_GtoH_mass_fine_"+Tag[i_tup]);
		TH1D *h_GtoH_diPt = (TH1D*)h_diPt->Clone("h_GtoH_diPt_"+Tag[i_tup]);
		TH1D *h_GtoH_rapi = (TH1D*)h_rapi->Clone("h_GtoH_rapi_"+Tag[i_tup]);
		TH1D *h_GtoH_pT = (TH1D*)h_pT->Clone("h_GtoH_pT_"+Tag[i_tup]);
		TH1D *h_GtoH_leadPt = (TH1D*)h_pT->Clone("h_GtoH_leadPt_"+Tag[i_tup]);
		TH1D *h_GtoH_subPt = (TH1D*)h_pT->Clone("h_GtoH_subPt_"+Tag[i_tup]);
		TH1D *h_GtoH_eta = (TH1D*)h_eta->Clone("h_GtoH_eta_"+Tag[i_tup]);
		TH1D *h_GtoH_leadEta = (TH1D*)h_eta->Clone("h_GtoH_leadEta_"+Tag[i_tup]);
		TH1D *h_GtoH_subEta = (TH1D*)h_eta->Clone("h_GtoH_subEta_"+Tag[i_tup]);
		TH1D *h_GtoH_phi = (TH1D*)h_phi->Clone("h_GtoH_phi_"+Tag[i_tup]);
		TH1D *h_GtoH_leadPhi = (TH1D*)h_phi->Clone("h_GtoH_leadPhi_"+Tag[i_tup]);
		TH1D *h_GtoH_subPhi = (TH1D*)h_phi->Clone("h_GtoH_subPhi_"+Tag[i_tup]);

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

		TH1D *h_Gen_M50to60_mass_fine = (TH1D*)h_mass_fine->Clone("h_Gen_M50to60_mass_fine_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_diPt = (TH1D*)h_diPt->Clone("h_Gen_M50to60_diPt_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_rapi = (TH1D*)h_rapi->Clone("h_Gen_M50to60_rapi_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_pT = (TH1D*)h_pT->Clone("h_Gen_M50to60_pT_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_leadPt = (TH1D*)h_pT->Clone("h_Gen_M50to60_leadPt_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_subPt = (TH1D*)h_pT->Clone("h_Gen_M50to60_subPt_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_eta = (TH1D*)h_eta->Clone("h_Gen_M50to60_eta_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_leadEta = (TH1D*)h_eta->Clone("h_Gen_M50to60_leadEta_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_subEta = (TH1D*)h_eta->Clone("h_Gen_M50to60_subEta_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_phi = (TH1D*)h_phi->Clone("h_Gen_M50to60_phi_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_leadPhi = (TH1D*)h_phi->Clone("h_Gen_M50to60_leadPhi_"+Tag[i_tup]);
		TH1D *h_Gen_M50to60_subPhi = (TH1D*)h_phi->Clone("h_Gen_M50to60_subPhi_"+Tag[i_tup]);

		TH1D *h_Gen_M60to120_mass_fine = (TH1D*)h_mass_fine->Clone("h_Gen_M60to120_mass_fine_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_diPt = (TH1D*)h_diPt->Clone("h_Gen_M60to120_diPt_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_rapi = (TH1D*)h_rapi->Clone("h_Gen_M60to120_rapi_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_pT = (TH1D*)h_pT->Clone("h_Gen_M60to120_pT_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_leadPt = (TH1D*)h_pT->Clone("h_Gen_M60to120_leadPt_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_subPt = (TH1D*)h_pT->Clone("h_Gen_M60to120_subPt_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_eta = (TH1D*)h_eta->Clone("h_Gen_M60to120_eta_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_leadEta = (TH1D*)h_eta->Clone("h_Gen_M60to120_leadEta_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_subEta = (TH1D*)h_eta->Clone("h_Gen_M60to120_subEta_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_phi = (TH1D*)h_phi->Clone("h_Gen_M60to120_phi_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_leadPhi = (TH1D*)h_phi->Clone("h_Gen_M60to120_leadPhi_"+Tag[i_tup]);
		TH1D *h_Gen_M60to120_subPhi = (TH1D*)h_phi->Clone("h_Gen_M60to120_subPhi_"+Tag[i_tup]);
		TH2D *h_Gen_M60to120_leadEta_subEta = new TH2D("h_Gen_M60to120_leadEta_subEta_"+Tag[i_tup], "", 1000, -5, 5, 1000, -5, 5);

		TH1D *h_Gen_M120to200_mass_fine = (TH1D*)h_mass_fine->Clone("h_Gen_M120to200_mass_fine_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_diPt = (TH1D*)h_diPt->Clone("h_Gen_M120to200_diPt_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_rapi = (TH1D*)h_rapi->Clone("h_Gen_M120to200_rapi_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_pT = (TH1D*)h_pT->Clone("h_Gen_M120to200_pT_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_leadPt = (TH1D*)h_pT->Clone("h_Gen_M120to200_leadPt_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_subPt = (TH1D*)h_pT->Clone("h_Gen_M120to200_subPt_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_eta = (TH1D*)h_eta->Clone("h_Gen_M120to200_eta_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_leadEta = (TH1D*)h_eta->Clone("h_Gen_M120to200_leadEta_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_subEta = (TH1D*)h_eta->Clone("h_Gen_M120to200_subEta_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_phi = (TH1D*)h_phi->Clone("h_Gen_M120to200_phi_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_leadPhi = (TH1D*)h_phi->Clone("h_Gen_M120to200_leadPhi_"+Tag[i_tup]);
		TH1D *h_Gen_M120to200_subPhi = (TH1D*)h_phi->Clone("h_Gen_M120to200_subPhi_"+Tag[i_tup]);

		TH1D *h_Gen_M200to400_mass_fine = (TH1D*)h_mass_fine->Clone("h_Gen_M200to400_mass_fine_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_diPt = (TH1D*)h_diPt->Clone("h_Gen_M200to400_diPt_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_rapi = (TH1D*)h_rapi->Clone("h_Gen_M200to400_rapi_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_pT = (TH1D*)h_pT->Clone("h_Gen_M200to400_pT_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_leadPt = (TH1D*)h_pT->Clone("h_Gen_M200to400_leadPt_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_subPt = (TH1D*)h_pT->Clone("h_Gen_M200to400_subPt_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_eta = (TH1D*)h_eta->Clone("h_Gen_M200to400_eta_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_leadEta = (TH1D*)h_eta->Clone("h_Gen_M200to400_leadEta_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_subEta = (TH1D*)h_eta->Clone("h_Gen_M200to400_subEta_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_phi = (TH1D*)h_phi->Clone("h_Gen_M200to400_phi_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_leadPhi = (TH1D*)h_phi->Clone("h_Gen_M200to400_leadPhi_"+Tag[i_tup]);
		TH1D *h_Gen_M200to400_subPhi = (TH1D*)h_phi->Clone("h_Gen_M200to400_subPhi_"+Tag[i_tup]);

		TH1D *h_Gen_M400to500_mass_fine = (TH1D*)h_mass_fine->Clone("h_Gen_M400to500_mass_fine_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_diPt = (TH1D*)h_diPt->Clone("h_Gen_M400to500_diPt_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_rapi = (TH1D*)h_rapi->Clone("h_Gen_M400to500_rapi_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_pT = (TH1D*)h_pT->Clone("h_Gen_M400to500_pT_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_leadPt = (TH1D*)h_pT->Clone("h_Gen_M400to500_leadPt_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_subPt = (TH1D*)h_pT->Clone("h_Gen_M400to500_subPt_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_eta = (TH1D*)h_eta->Clone("h_Gen_M400to500_eta_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_leadEta = (TH1D*)h_eta->Clone("h_Gen_M400to500_leadEta_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_subEta = (TH1D*)h_eta->Clone("h_Gen_M400to500_subEta_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_phi = (TH1D*)h_phi->Clone("h_Gen_M400to500_phi_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_leadPhi = (TH1D*)h_phi->Clone("h_Gen_M400to500_leadPhi_"+Tag[i_tup]);
		TH1D *h_Gen_M400to500_subPhi = (TH1D*)h_phi->Clone("h_Gen_M400to500_subPhi_"+Tag[i_tup]);

		TH1D *h_BEC_mass_fine = (TH1D*)h_mass_fine->Clone("h_BEC_mass_fine_"+Tag[i_tup]);
		TH1D *h_BEC_diPt = (TH1D*)h_diPt->Clone("h_BEC_diPt_"+Tag[i_tup]);
		TH1D *h_BEC_rapi = (TH1D*)h_rapi->Clone("h_BEC_rapi_"+Tag[i_tup]);
		TH1D *h_BEC_pT = (TH1D*)h_pT->Clone("h_BEC_pT_"+Tag[i_tup]);
		TH1D *h_BEC_leadPt = (TH1D*)h_pT->Clone("h_BEC_leadPt_"+Tag[i_tup]);
		TH1D *h_BEC_subPt = (TH1D*)h_pT->Clone("h_BEC_subPt_"+Tag[i_tup]);
		TH1D *h_BEC_eta = (TH1D*)h_eta->Clone("h_BEC_eta_"+Tag[i_tup]);
		TH1D *h_BEC_leadEta = (TH1D*)h_eta->Clone("h_BEC_leadEta_"+Tag[i_tup]);
		TH1D *h_BEC_subEta = (TH1D*)h_eta->Clone("h_BEC_subEta_"+Tag[i_tup]);
		TH1D *h_BEC_phi = (TH1D*)h_phi->Clone("h_BEC_phi_"+Tag[i_tup]);
		TH1D *h_BEC_leadPhi = (TH1D*)h_phi->Clone("h_BEC_leadPhi_"+Tag[i_tup]);
		TH1D *h_BEC_subPhi = (TH1D*)h_phi->Clone("h_BEC_subPhi_"+Tag[i_tup]);

		TH1D *h_Reco_mass_fine = (TH1D*)h_mass_fine->Clone("h_Reco_mass_fine_"+Tag[i_tup]);
		TH1D *h_Reco_diPt = (TH1D*)h_diPt->Clone("h_Reco_diPt_"+Tag[i_tup]);
		TH1D *h_Reco_rapi = (TH1D*)h_rapi->Clone("h_Reco_rapi_"+Tag[i_tup]);
		TH1D *h_Reco_pT = (TH1D*)h_pT->Clone("h_Reco_pT_"+Tag[i_tup]);
		TH1D *h_Reco_leadPt = (TH1D*)h_pT->Clone("h_Reco_leadPt_"+Tag[i_tup]);
		TH1D *h_Reco_subPt = (TH1D*)h_pT->Clone("h_Reco_subPt_"+Tag[i_tup]);
		TH1D *h_Reco_eta = (TH1D*)h_eta->Clone("h_Reco_eta_"+Tag[i_tup]);
		TH1D *h_Reco_leadEta = (TH1D*)h_eta->Clone("h_Reco_leadEta_"+Tag[i_tup]);
		TH1D *h_Reco_subEta = (TH1D*)h_eta->Clone("h_Reco_subEta_"+Tag[i_tup]);
		TH1D *h_Reco_phi = (TH1D*)h_phi->Clone("h_Reco_phi_"+Tag[i_tup]);
		TH1D *h_Reco_leadPhi = (TH1D*)h_phi->Clone("h_Reco_leadPhi_"+Tag[i_tup]);
		TH1D *h_Reco_subPhi = (TH1D*)h_phi->Clone("h_Reco_subPhi_"+Tag[i_tup]);

		TH1D *h_BtoF_Reco_mass_fine = (TH1D*)h_mass_fine->Clone("h_BtoF_Reco_mass_fine_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_diPt = (TH1D*)h_diPt->Clone("h_BtoF_Reco_diPt_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_rapi = (TH1D*)h_rapi->Clone("h_BtoF_Reco_rapi_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_pT = (TH1D*)h_pT->Clone("h_BtoF_Reco_pT_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_leadPt = (TH1D*)h_pT->Clone("h_BtoF_Reco_leadPt_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_subPt = (TH1D*)h_pT->Clone("h_BtoF_Reco_subPt_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_eta = (TH1D*)h_eta->Clone("h_BtoF_Reco_eta_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_leadEta = (TH1D*)h_eta->Clone("h_BtoF_Reco_leadEta_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_subEta = (TH1D*)h_eta->Clone("h_BtoF_Reco_subEta_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_phi = (TH1D*)h_phi->Clone("h_BtoF_Reco_phi_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_leadPhi = (TH1D*)h_phi->Clone("h_BtoF_Reco_leadPhi_"+Tag[i_tup]);
		TH1D *h_BtoF_Reco_subPhi = (TH1D*)h_phi->Clone("h_BtoF_Reco_subPhi_"+Tag[i_tup]);

		TH1D *h_GtoH_Reco_mass_fine = (TH1D*)h_mass_fine->Clone("h_GtoH_Reco_mass_fine_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_diPt = (TH1D*)h_diPt->Clone("h_GtoH_Reco_diPt_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_rapi = (TH1D*)h_rapi->Clone("h_GtoH_Reco_rapi_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_pT = (TH1D*)h_pT->Clone("h_GtoH_Reco_pT_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_leadPt = (TH1D*)h_pT->Clone("h_GtoH_Reco_leadPt_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_subPt = (TH1D*)h_pT->Clone("h_GtoH_Reco_subPt_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_eta = (TH1D*)h_eta->Clone("h_GtoH_Reco_eta_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_leadEta = (TH1D*)h_eta->Clone("h_GtoH_Reco_leadEta_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_subEta = (TH1D*)h_eta->Clone("h_GtoH_Reco_subEta_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_phi = (TH1D*)h_phi->Clone("h_GtoH_Reco_phi_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_leadPhi = (TH1D*)h_phi->Clone("h_GtoH_Reco_leadPhi_"+Tag[i_tup]);
		TH1D *h_GtoH_Reco_subPhi = (TH1D*)h_phi->Clone("h_GtoH_Reco_subPhi_"+Tag[i_tup]);

		TH1D *h_RecoID_mass_fine = (TH1D*)h_mass_fine->Clone("h_RecoID_mass_fine_"+Tag[i_tup]);
		TH1D *h_RecoID_diPt = (TH1D*)h_diPt->Clone("h_RecoID_diPt_"+Tag[i_tup]);
		TH1D *h_RecoID_rapi = (TH1D*)h_rapi->Clone("h_RecoID_rapi_"+Tag[i_tup]);
		TH1D *h_RecoID_pT = (TH1D*)h_pT->Clone("h_RecoID_pT_"+Tag[i_tup]);
		TH1D *h_RecoID_leadPt = (TH1D*)h_pT->Clone("h_RecoID_leadPt_"+Tag[i_tup]);
		TH1D *h_RecoID_subPt = (TH1D*)h_pT->Clone("h_RecoID_subPt_"+Tag[i_tup]);
		TH1D *h_RecoID_eta = (TH1D*)h_eta->Clone("h_RecoID_eta_"+Tag[i_tup]);
		TH1D *h_RecoID_leadEta = (TH1D*)h_eta->Clone("h_RecoID_leadEta_"+Tag[i_tup]);
		TH1D *h_RecoID_subEta = (TH1D*)h_eta->Clone("h_RecoID_subEta_"+Tag[i_tup]);
		TH1D *h_RecoID_phi = (TH1D*)h_phi->Clone("h_RecoID_phi_"+Tag[i_tup]);
		TH1D *h_RecoID_leadPhi = (TH1D*)h_phi->Clone("h_RecoID_leadPhi_"+Tag[i_tup]);
		TH1D *h_RecoID_subPhi = (TH1D*)h_phi->Clone("h_RecoID_subPhi_"+Tag[i_tup]);

		TH1D *h_BtoF_RecoID_mass_fine = (TH1D*)h_mass_fine->Clone("h_BtoF_RecoID_mass_fine_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_diPt = (TH1D*)h_diPt->Clone("h_BtoF_RecoID_diPt_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_rapi = (TH1D*)h_rapi->Clone("h_BtoF_RecoID_rapi_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_pT = (TH1D*)h_pT->Clone("h_BtoF_RecoID_pT_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_leadPt = (TH1D*)h_pT->Clone("h_BtoF_RecoID_leadPt_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_subPt = (TH1D*)h_pT->Clone("h_BtoF_RecoID_subPt_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_eta = (TH1D*)h_eta->Clone("h_BtoF_RecoID_eta_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_leadEta = (TH1D*)h_eta->Clone("h_BtoF_RecoID_leadEta_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_subEta = (TH1D*)h_eta->Clone("h_BtoF_RecoID_subEta_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_phi = (TH1D*)h_phi->Clone("h_BtoF_RecoID_phi_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_leadPhi = (TH1D*)h_phi->Clone("h_BtoF_RecoID_leadPhi_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoID_subPhi = (TH1D*)h_phi->Clone("h_BtoF_RecoID_subPhi_"+Tag[i_tup]);

		TH1D *h_GtoH_RecoID_mass_fine = (TH1D*)h_mass_fine->Clone("h_GtoH_RecoID_mass_fine_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_diPt = (TH1D*)h_diPt->Clone("h_GtoH_RecoID_diPt_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_rapi = (TH1D*)h_rapi->Clone("h_GtoH_RecoID_rapi_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_pT = (TH1D*)h_pT->Clone("h_GtoH_RecoID_pT_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_leadPt = (TH1D*)h_pT->Clone("h_GtoH_RecoID_leadPt_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_subPt = (TH1D*)h_pT->Clone("h_GtoH_RecoID_subPt_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_eta = (TH1D*)h_eta->Clone("h_GtoH_RecoID_eta_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_leadEta = (TH1D*)h_eta->Clone("h_GtoH_RecoID_leadEta_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_subEta = (TH1D*)h_eta->Clone("h_GtoH_RecoID_subEta_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_phi = (TH1D*)h_phi->Clone("h_GtoH_RecoID_phi_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_leadPhi = (TH1D*)h_phi->Clone("h_GtoH_RecoID_leadPhi_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoID_subPhi = (TH1D*)h_phi->Clone("h_GtoH_RecoID_subPhi_"+Tag[i_tup]);

		TH1D *h_RecoIDIso_mass_fine = (TH1D*)h_mass_fine->Clone("h_RecoIDIso_mass_fine_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_diPt = (TH1D*)h_diPt->Clone("h_RecoIDIso_diPt_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_rapi = (TH1D*)h_rapi->Clone("h_RecoIDIso_rapi_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_pT = (TH1D*)h_pT->Clone("h_RecoIDIso_pT_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_leadPt = (TH1D*)h_pT->Clone("h_RecoIDIso_leadPt_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_subPt = (TH1D*)h_pT->Clone("h_RecoIDIso_subPt_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_eta = (TH1D*)h_eta->Clone("h_RecoIDIso_eta_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_leadEta = (TH1D*)h_eta->Clone("h_RecoIDIso_leadEta_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_subEta = (TH1D*)h_eta->Clone("h_RecoIDIso_subEta_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_phi = (TH1D*)h_phi->Clone("h_RecoIDIso_phi_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_leadPhi = (TH1D*)h_phi->Clone("h_RecoIDIso_leadPhi_"+Tag[i_tup]);
		TH1D *h_RecoIDIso_subPhi = (TH1D*)h_phi->Clone("h_RecoIDIso_subPhi_"+Tag[i_tup]);

		TH1D *h_BtoF_RecoIDIso_mass_fine = (TH1D*)h_mass_fine->Clone("h_BtoF_RecoIDIso_mass_fine_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_diPt = (TH1D*)h_diPt->Clone("h_BtoF_RecoIDIso_diPt_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_rapi = (TH1D*)h_rapi->Clone("h_BtoF_RecoIDIso_rapi_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_pT = (TH1D*)h_pT->Clone("h_BtoF_RecoIDIso_pT_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_leadPt = (TH1D*)h_pT->Clone("h_BtoF_RecoIDIso_leadPt_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_subPt = (TH1D*)h_pT->Clone("h_BtoF_RecoIDIso_subPt_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_eta = (TH1D*)h_eta->Clone("h_BtoF_RecoIDIso_eta_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_leadEta = (TH1D*)h_eta->Clone("h_BtoF_RecoIDIso_leadEta_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_subEta = (TH1D*)h_eta->Clone("h_BtoF_RecoIDIso_subEta_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_phi = (TH1D*)h_phi->Clone("h_BtoF_RecoIDIso_phi_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_leadPhi = (TH1D*)h_phi->Clone("h_BtoF_RecoIDIso_leadPhi_"+Tag[i_tup]);
		TH1D *h_BtoF_RecoIDIso_subPhi = (TH1D*)h_phi->Clone("h_BtoF_RecoIDIso_subPhi_"+Tag[i_tup]);

		TH1D *h_GtoH_RecoIDIso_mass_fine = (TH1D*)h_mass_fine->Clone("h_GtoH_RecoIDIso_mass_fine_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_diPt = (TH1D*)h_diPt->Clone("h_GtoH_RecoIDIso_diPt_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_rapi = (TH1D*)h_rapi->Clone("h_GtoH_RecoIDIso_rapi_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_pT = (TH1D*)h_pT->Clone("h_GtoH_RecoIDIso_pT_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_leadPt = (TH1D*)h_pT->Clone("h_GtoH_RecoIDIso_leadPt_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_subPt = (TH1D*)h_pT->Clone("h_GtoH_RecoIDIso_subPt_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_eta = (TH1D*)h_eta->Clone("h_GtoH_RecoIDIso_eta_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_leadEta = (TH1D*)h_eta->Clone("h_GtoH_RecoIDIso_leadEta_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_subEta = (TH1D*)h_eta->Clone("h_GtoH_RecoIDIso_subEta_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_phi = (TH1D*)h_phi->Clone("h_GtoH_RecoIDIso_phi_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_leadPhi = (TH1D*)h_phi->Clone("h_GtoH_RecoIDIso_leadPhi_"+Tag[i_tup]);
		TH1D *h_GtoH_RecoIDIso_subPhi = (TH1D*)h_phi->Clone("h_GtoH_RecoIDIso_subPhi_"+Tag[i_tup]);

		TH1D *h_before_PUCorr_mass_fine = (TH1D*)h_mass_fine->Clone("h_before_PUCorr_mass_fine_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_diPt = (TH1D*)h_diPt->Clone("h_before_PUCorr_diPt_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_rapi = (TH1D*)h_rapi->Clone("h_before_PUCorr_rapi_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_pT = (TH1D*)h_pT->Clone("h_before_PUCorr_pT_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_leadPt = (TH1D*)h_pT->Clone("h_before_PUCorr_leadPt_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_subPt = (TH1D*)h_pT->Clone("h_before_PUCorr_subPt_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_eta = (TH1D*)h_eta->Clone("h_before_PUCorr_eta_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_leadEta = (TH1D*)h_eta->Clone("h_before_PUCorr_leadEta_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_subEta = (TH1D*)h_eta->Clone("h_before_PUCorr_subEta_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_phi = (TH1D*)h_phi->Clone("h_before_PUCorr_phi_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_leadPhi = (TH1D*)h_phi->Clone("h_before_PUCorr_leadPhi_"+Tag[i_tup]);
		TH1D *h_before_PUCorr_subPhi = (TH1D*)h_phi->Clone("h_before_PUCorr_subPhi_"+Tag[i_tup]);

		TH1D *h_before_RoccoR_mass_fine = (TH1D*)h_mass_fine->Clone("h_before_RoccoR_mass_fine_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_diPt = (TH1D*)h_diPt->Clone("h_before_RoccoR_diPt_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_rapi = (TH1D*)h_rapi->Clone("h_before_RoccoR_rapi_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_pT = (TH1D*)h_pT->Clone("h_before_RoccoR_pT_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_leadPt = (TH1D*)h_pT->Clone("h_before_RoccoR_leadPt_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_subPt = (TH1D*)h_pT->Clone("h_before_RoccoR_subPt_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_eta = (TH1D*)h_eta->Clone("h_before_RoccoR_eta_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_leadEta = (TH1D*)h_eta->Clone("h_before_RoccoR_leadEta_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_subEta = (TH1D*)h_eta->Clone("h_before_RoccoR_subEta_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_phi = (TH1D*)h_phi->Clone("h_before_RoccoR_phi_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_leadPhi = (TH1D*)h_phi->Clone("h_before_RoccoR_leadPhi_"+Tag[i_tup]);
		TH1D *h_before_RoccoR_subPhi = (TH1D*)h_phi->Clone("h_before_RoccoR_subPhi_"+Tag[i_tup]);

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

			// -- efficiency weights -- //
			Double_t effweight = 1, effweight_Reco = 1, effweight_RecoID = 1, effweight_RecoIDIso = 1,
					 effweight_BtoF = 1, effweight_BtoF_Reco = 1, effweight_BtoF_RecoID = 1, effweight_BtoF_RecoIDIso = 1,
					 effweight_GtoH = 1, effweight_GtoH_Reco = 1, effweight_GtoH_RecoID = 1, effweight_GtoH_RecoIDIso = 1;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			// -- SF of LeadEtaCorr -- //
			Double_t SF_LeadEtaCorr = 1, SF_BtoF_LeadEtaCorr = 1, SF_GtoH_LeadEtaCorr = 1;

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
			if( isMC == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				vector< GenLepton > GenLeptonCollection;
				Int_t NGenLeptons = ntuple->gnpair;
				for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
				{
					GenLepton genlep;
					genlep.FillFromNtuple(ntuple, i_gen);
					if( genlep.isMuon() && genlep.fromHardProcessFinalState )
						GenLeptonCollection.push_back( genlep );
				}

				if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
				{
					GenLepton genlep1 = GenLeptonCollection[0];
					GenLepton genlep2 = GenLeptonCollection[1];

					Bool_t isPassAcc_GenLepton = kFALSE;
					isPassAcc_GenLepton = analyzer->isPassAccCondition_GenLepton(genlep1, genlep2);
					if( isPassAcc_GenLepton == kTRUE )
					{
						Double_t gen_M = (genlep1.Momentum + genlep2.Momentum).M();
						Double_t gen_Pt = (genlep1.Momentum + genlep2.Momentum).Pt();
						Double_t gen_rapi = (genlep1.Momentum + genlep2.Momentum).Rapidity();

						h_Gen_mass_fine->Fill( gen_M, TotWeight );
						h_Gen_diPt->Fill( gen_Pt, TotWeight );
						h_Gen_rapi->Fill( gen_rapi, TotWeight );
						h_Gen_pT->Fill( genlep1.Pt, TotWeight );
						h_Gen_pT->Fill( genlep2.Pt, TotWeight );
						h_Gen_eta->Fill( genlep1.eta, TotWeight );
						h_Gen_eta->Fill( genlep2.eta, TotWeight );
						h_Gen_phi->Fill( genlep1.phi, TotWeight );
						h_Gen_phi->Fill( genlep2.phi, TotWeight );

						if( genlep1.Pt > genlep2.Pt )
						{
							h_Gen_leadPt->Fill( genlep1.Pt, TotWeight );
							h_Gen_leadEta->Fill( genlep1.eta, TotWeight );
							h_Gen_leadPhi->Fill( genlep1.phi, TotWeight );
							h_Gen_subPt->Fill( genlep2.Pt, TotWeight );
							h_Gen_subEta->Fill( genlep2.eta, TotWeight );
							h_Gen_subPhi->Fill( genlep2.phi, TotWeight );
						}
						else
						{
							h_Gen_leadPt->Fill( genlep2.Pt, TotWeight );
							h_Gen_leadEta->Fill( genlep2.eta, TotWeight );
							h_Gen_leadPhi->Fill( genlep2.phi, TotWeight );
							h_Gen_subPt->Fill( genlep1.Pt, TotWeight );
							h_Gen_subEta->Fill( genlep1.eta, TotWeight );
							h_Gen_subPhi->Fill( genlep1.phi, TotWeight );
						}

						if( type == 12 && 50 < gen_M && gen_M < 60 )
						{
                            
                           //for (int i=M50to60; i<=M200to500; ++i){
                           //if ( min_mom[i] < gen_M && gen_M < max_mom[i] ) {
                             //  h_Gen_mass_fine_Mom[i]->Fill( gen_M, TotWeight );
                           //}	                                
                           //}
							h_Gen_M50to60_mass_fine->Fill( gen_M, TotWeight );
							h_Gen_M50to60_diPt->Fill( gen_Pt, TotWeight );
							h_Gen_M50to60_rapi->Fill( gen_rapi, TotWeight );
							h_Gen_M50to60_pT->Fill( genlep1.Pt, TotWeight );
							h_Gen_M50to60_pT->Fill( genlep2.Pt, TotWeight );
							h_Gen_M50to60_eta->Fill( genlep1.eta, TotWeight );
							h_Gen_M50to60_eta->Fill( genlep2.eta, TotWeight );
							h_Gen_M50to60_phi->Fill( genlep1.phi, TotWeight );
							h_Gen_M50to60_phi->Fill( genlep2.phi, TotWeight );

							if( genlep1.Pt > genlep2.Pt )
							{
								h_Gen_M50to60_leadPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M50to60_leadEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M50to60_leadPhi->Fill( genlep1.phi, TotWeight );
								h_Gen_M50to60_subPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M50to60_subEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M50to60_subPhi->Fill( genlep2.phi, TotWeight );
							}
							else
							{
								h_Gen_M50to60_leadPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M50to60_leadEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M50to60_leadPhi->Fill( genlep2.phi, TotWeight );
								h_Gen_M50to60_subPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M50to60_subEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M50to60_subPhi->Fill( genlep1.phi, TotWeight );
							}
						}
						else if( type == 12 && 60 < gen_M && gen_M < 120 )
						{
							h_Gen_M60to120_mass_fine->Fill( gen_M, TotWeight );
							h_Gen_M60to120_diPt->Fill( gen_Pt, TotWeight );
							h_Gen_M60to120_rapi->Fill( gen_rapi, TotWeight );
							h_Gen_M60to120_pT->Fill( genlep1.Pt, TotWeight );
							h_Gen_M60to120_pT->Fill( genlep2.Pt, TotWeight );
							h_Gen_M60to120_eta->Fill( genlep1.eta, TotWeight );
							h_Gen_M60to120_eta->Fill( genlep2.eta, TotWeight );
							h_Gen_M60to120_phi->Fill( genlep1.phi, TotWeight );
							h_Gen_M60to120_phi->Fill( genlep2.phi, TotWeight );

							if( genlep1.Pt > genlep2.Pt )
							{
								h_Gen_M60to120_leadPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M60to120_leadEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M60to120_leadPhi->Fill( genlep1.phi, TotWeight );
								h_Gen_M60to120_subPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M60to120_subEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M60to120_subPhi->Fill( genlep2.phi, TotWeight );
								h_Gen_M60to120_leadEta_subEta->Fill( genlep1.eta, genlep2.eta, TotWeight );
							}
							else
							{
								h_Gen_M60to120_leadPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M60to120_leadEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M60to120_leadPhi->Fill( genlep2.phi, TotWeight );
								h_Gen_M60to120_subPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M60to120_subEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M60to120_subPhi->Fill( genlep1.phi, TotWeight );
								h_Gen_M60to120_leadEta_subEta->Fill( genlep2.eta, genlep1.eta, TotWeight );
							}
						}
						else if( type == 12 && 120 < gen_M && gen_M < 200 )
						{
							h_Gen_M120to200_mass_fine->Fill( gen_M, TotWeight );
							h_Gen_M120to200_diPt->Fill( gen_Pt, TotWeight );
							h_Gen_M120to200_rapi->Fill( gen_rapi, TotWeight );
							h_Gen_M120to200_pT->Fill( genlep1.Pt, TotWeight );
							h_Gen_M120to200_pT->Fill( genlep2.Pt, TotWeight );
							h_Gen_M120to200_eta->Fill( genlep1.eta, TotWeight );
							h_Gen_M120to200_eta->Fill( genlep2.eta, TotWeight );
							h_Gen_M120to200_phi->Fill( genlep1.phi, TotWeight );
							h_Gen_M120to200_phi->Fill( genlep2.phi, TotWeight );

							if( genlep1.Pt > genlep2.Pt )
							{
								h_Gen_M120to200_leadPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M120to200_leadEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M120to200_leadPhi->Fill( genlep1.phi, TotWeight );
								h_Gen_M120to200_subPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M120to200_subEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M120to200_subPhi->Fill( genlep2.phi, TotWeight );
							}
							else
							{
								h_Gen_M120to200_leadPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M120to200_leadEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M120to200_leadPhi->Fill( genlep2.phi, TotWeight );
								h_Gen_M120to200_subPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M120to200_subEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M120to200_subPhi->Fill( genlep1.phi, TotWeight );
							}
						}
						else if( type == 12 && 200 < gen_M && gen_M < 400 )
						{
							h_Gen_M200to400_mass_fine->Fill( gen_M, TotWeight );
							h_Gen_M200to400_diPt->Fill( gen_Pt, TotWeight );
							h_Gen_M200to400_rapi->Fill( gen_rapi, TotWeight );
							h_Gen_M200to400_pT->Fill( genlep1.Pt, TotWeight );
							h_Gen_M200to400_pT->Fill( genlep2.Pt, TotWeight );
							h_Gen_M200to400_eta->Fill( genlep1.eta, TotWeight );
							h_Gen_M200to400_eta->Fill( genlep2.eta, TotWeight );
							h_Gen_M200to400_phi->Fill( genlep1.phi, TotWeight );
							h_Gen_M200to400_phi->Fill( genlep2.phi, TotWeight );

							if( genlep1.Pt > genlep2.Pt )
							{
								h_Gen_M200to400_leadPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M200to400_leadEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M200to400_leadPhi->Fill( genlep1.phi, TotWeight );
								h_Gen_M200to400_subPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M200to400_subEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M200to400_subPhi->Fill( genlep2.phi, TotWeight );
							}
							else
							{
								h_Gen_M200to400_leadPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M200to400_leadEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M200to400_leadPhi->Fill( genlep2.phi, TotWeight );
								h_Gen_M200to400_subPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M200to400_subEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M200to400_subPhi->Fill( genlep1.phi, TotWeight );
							}
						}
						else if( type == 12 && 400 < gen_M && gen_M < 500 )
						{
							h_Gen_M400to500_mass_fine->Fill( gen_M, TotWeight );
							h_Gen_M400to500_diPt->Fill( gen_Pt, TotWeight );
							h_Gen_M400to500_rapi->Fill( gen_rapi, TotWeight );
							h_Gen_M400to500_pT->Fill( genlep1.Pt, TotWeight );
							h_Gen_M400to500_pT->Fill( genlep2.Pt, TotWeight );
							h_Gen_M400to500_eta->Fill( genlep1.eta, TotWeight );
							h_Gen_M400to500_eta->Fill( genlep2.eta, TotWeight );
							h_Gen_M400to500_phi->Fill( genlep1.phi, TotWeight );
							h_Gen_M400to500_phi->Fill( genlep2.phi, TotWeight );

							if( genlep1.Pt > genlep2.Pt )
							{
								h_Gen_M400to500_leadPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M400to500_leadEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M400to500_leadPhi->Fill( genlep1.phi, TotWeight );
								h_Gen_M400to500_subPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M400to500_subEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M400to500_subPhi->Fill( genlep2.phi, TotWeight );
							}
							else
							{
								h_Gen_M400to500_leadPt->Fill( genlep2.Pt, TotWeight );
								h_Gen_M400to500_leadEta->Fill( genlep2.eta, TotWeight );
								h_Gen_M400to500_leadPhi->Fill( genlep2.phi, TotWeight );
								h_Gen_M400to500_subPt->Fill( genlep1.Pt, TotWeight );
								h_Gen_M400to500_subEta->Fill( genlep1.eta, TotWeight );
								h_Gen_M400to500_subPhi->Fill( genlep1.phi, TotWeight );
							}
						}
					} //isPassAcc_GenLepton
				} //the events containing 2 muons from hard-process
			} //End of Gen-level selection

			Bool_t TriggerFlag = kFALSE;
			TriggerFlag = ntuple->isTriggered( analyzer->HLT );

			////////////////////////////////
			// -- Reco level selection -- //
			////////////////////////////////
			if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				vector< Muon > MuonCollection; vector< Muon > MuonCollection_noRoccoR;
				Int_t NLeptons = ntuple->nMuon;
				for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
				{
					Muon mu;
					mu.FillFromNtuple(ntuple, i_reco);

					// -- Convert to TuneP variables -- //
					//analyzer->ConvertToTunePInfo( mu );

					MuonCollection_noRoccoR.push_back( mu );

					////////////////////////////////
					// -- Rochester correction -- //
					////////////////////////////////
					Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
					Int_t s, m;
						
					if( isMC == kFALSE )
						//SF = rc.kScaleDT(mu.charge, mu.TuneP_pT, mu.TuneP_eta, mu.TuneP_phi, s=0, m=0);
						SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
					else
						//SF = rc.kScaleAndSmearMC(mu.charge, mu.TuneP_pT, mu.TuneP_eta, mu.TuneP_phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
						SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);

					mu.Pt = SF*mu.Pt;
					mu.Momentum.SetPtEtaPhiM( mu.Pt, mu.eta, mu.phi, M_Mu );
				
					// -- Convert to TuneP variables -- //
					//mu.TuneP_pT = SF*mu.TuneP_pT;
					//analyzer->ConvertToTunePInfo( mu );

					MuonCollection.push_back( mu );
				}	

				// -- Event Selection -- //
				vector< Muon > SelectedMuonCollection;
				Bool_t isPassEventSelection = kFALSE;
				isPassEventSelection = analyzer->EventSelection(MuonCollection, ntuple, &SelectedMuonCollection);

				if( isPassEventSelection == kTRUE )
				{
					Muon mu1 = SelectedMuonCollection[0];
					Muon mu2 = SelectedMuonCollection[1];

					Double_t reco_M = (mu1.Momentum + mu2.Momentum).M();
					Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
					Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

					//////////////////////////////
					// -- Apply scale factor -- //
					//////////////////////////////
					if( isMC == kTRUE )
					{
						effweight_BtoF = analyzer->EfficiencySF_EventWeight_HLT_BtoF( mu1, mu2 );
						effweight_BtoF_Reco = analyzer->EfficiencySF_EventWeight_HLT_BtoF_Reco( mu1, mu2 );
						effweight_BtoF_RecoID = analyzer->EfficiencySF_EventWeight_HLT_BtoF_RecoID( mu1, mu2 );
						effweight_BtoF_RecoIDIso = analyzer->EfficiencySF_EventWeight_HLT_BtoF_RecoIDIso( mu1, mu2 );

						effweight_GtoH = analyzer->EfficiencySF_EventWeight_HLT_GtoH( mu1, mu2 );
						effweight_GtoH_Reco = analyzer->EfficiencySF_EventWeight_HLT_GtoH_Reco( mu1, mu2 );
						effweight_GtoH_RecoID = analyzer->EfficiencySF_EventWeight_HLT_GtoH_RecoID( mu1, mu2 );
						effweight_GtoH_RecoIDIso = analyzer->EfficiencySF_EventWeight_HLT_GtoH_RecoIDIso( mu1, mu2 );

						effweight = (Lumi_BtoF * effweight_BtoF + Lumi_GtoH * effweight_GtoH) / Lumi;
						effweight_Reco = (Lumi_BtoF * effweight_BtoF_Reco + Lumi_GtoH * effweight_GtoH_Reco) / Lumi;
						effweight_RecoID = (Lumi_BtoF * effweight_BtoF_RecoID + Lumi_GtoH * effweight_GtoH_RecoID) / Lumi;
						effweight_RecoIDIso = (Lumi_BtoF * effweight_BtoF_RecoIDIso + Lumi_GtoH * effweight_GtoH_RecoIDIso) / Lumi;
							
						if( 60 < reco_M && reco_M < 120 )
						{
							SF_BtoF_LeadEtaCorr = analyzer->LeadEtaCorr( type, "BtoF", mu1, mu2 );
							SF_GtoH_LeadEtaCorr = analyzer->LeadEtaCorr( type, "GtoH", mu1, mu2 );
							SF_LeadEtaCorr = analyzer->LeadEtaCorr( type, "BtoH", mu1, mu2 );
						}
					}

					h_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
					h_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );

					if( mu1.Pt > mu2.Pt ) //leading : mu1
					{
						h_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_leadEta_subEta->Fill( mu1.eta, mu2.eta, TotWeight * PUWeight * effweight );
					}
					else
					{
						h_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight * SF_LeadEtaCorr );
						h_leadEta_subEta->Fill( mu2.eta, mu1.eta, TotWeight * PUWeight * effweight );
					}

					if( isMC == kTRUE )
					{
						h_BtoF_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
						h_BtoF_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );

						h_GtoH_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
						h_GtoH_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );

						h_Reco_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_Reco );
						h_Reco_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_Reco );
						h_Reco_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_Reco );
						h_Reco_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_Reco );
						h_Reco_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_Reco );
						h_Reco_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_Reco );
						h_Reco_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_Reco );
						h_Reco_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_Reco );
						h_Reco_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_Reco );

						h_BtoF_Reco_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_Reco );
						h_BtoF_Reco_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_Reco );

						h_GtoH_Reco_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_Reco );
						h_GtoH_Reco_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_Reco );

						h_RecoID_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_RecoID );
						h_RecoID_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_RecoID );

						h_BtoF_RecoID_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_RecoID );
						h_BtoF_RecoID_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_RecoID );

						h_GtoH_RecoID_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_RecoID );
						h_GtoH_RecoID_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_RecoID );

						h_RecoIDIso_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_RecoIDIso );
						h_RecoIDIso_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_RecoIDIso );

						h_BtoF_RecoIDIso_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						h_BtoF_RecoIDIso_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );

						h_GtoH_RecoIDIso_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						h_GtoH_RecoIDIso_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );

						if( mu1.Pt > mu2.Pt ) //leading : mu1
						{
							h_BtoF_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );

							h_GtoH_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );

							h_Reco_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_Reco );
							h_Reco_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_Reco );
							h_Reco_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_Reco );
							h_Reco_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_Reco );
							h_Reco_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_Reco );
							h_Reco_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_Reco );

							h_BtoF_Reco_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_Reco );

							h_GtoH_Reco_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_Reco );

							h_RecoID_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_RecoID );

							h_BtoF_RecoID_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_RecoID );

							h_GtoH_RecoID_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_RecoID );

							h_RecoIDIso_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_RecoIDIso );

							h_BtoF_RecoIDIso_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						
							h_GtoH_RecoIDIso_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						}
						else //leading : mu2
						{
							h_BtoF_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );
							h_BtoF_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF * SF_BtoF_LeadEtaCorr );

							h_GtoH_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );
							h_GtoH_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH * SF_GtoH_LeadEtaCorr );

							h_Reco_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_Reco );
							h_Reco_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_Reco );
							h_Reco_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_Reco );
							h_Reco_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_Reco );
							h_Reco_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_Reco );
							h_Reco_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_Reco );

							h_BtoF_Reco_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_Reco );
							h_BtoF_Reco_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_Reco );

							h_GtoH_Reco_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_Reco );
							h_GtoH_Reco_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_Reco );

							h_RecoID_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_RecoID );
							h_RecoID_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_RecoID );

							h_BtoF_RecoID_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_RecoID );
							h_BtoF_RecoID_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_RecoID );

							h_GtoH_RecoID_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_RecoID );
							h_GtoH_RecoID_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_RecoID );

							h_RecoIDIso_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_RecoIDIso );
							h_RecoIDIso_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_RecoIDIso );

							h_BtoF_RecoIDIso_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
							h_BtoF_RecoIDIso_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_BtoF_RecoIDIso );
						
							h_GtoH_RecoIDIso_leadPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_leadEta->Fill( mu2.eta, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_leadPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_subPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_subEta->Fill( mu1.eta, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
							h_GtoH_RecoIDIso_subPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight_GtoH_RecoIDIso );
						}

						/////////////////////////////////////////
						// -- Before applying efficiency SF -- //
						/////////////////////////////////////////
						h_BEC_mass_fine->Fill( reco_M, TotWeight * PUWeight );
						h_BEC_diPt->Fill( reco_Pt, TotWeight * PUWeight );
						h_BEC_rapi->Fill( reco_rapi, TotWeight * PUWeight );
						h_BEC_pT->Fill( mu1.Pt, TotWeight * PUWeight );
						h_BEC_pT->Fill( mu2.Pt, TotWeight * PUWeight );
						h_BEC_eta->Fill( mu1.eta, TotWeight * PUWeight );
						h_BEC_eta->Fill( mu2.eta, TotWeight * PUWeight );
						h_BEC_phi->Fill( mu1.phi, TotWeight * PUWeight );
						h_BEC_phi->Fill( mu2.phi, TotWeight * PUWeight );

						if( mu1.Pt > mu2.Pt )
						{
							h_BEC_leadPt->Fill( mu1.Pt, TotWeight * PUWeight );
							h_BEC_leadEta->Fill( mu1.eta, TotWeight * PUWeight );
							h_BEC_leadPhi->Fill( mu1.phi, TotWeight * PUWeight );
							h_BEC_subPt->Fill( mu2.Pt, TotWeight * PUWeight );
							h_BEC_subEta->Fill( mu2.eta, TotWeight * PUWeight );
							h_BEC_subPhi->Fill( mu2.phi, TotWeight * PUWeight );
						}
						else
						{
							h_BEC_leadPt->Fill( mu2.Pt, TotWeight * PUWeight );
							h_BEC_leadEta->Fill( mu2.eta, TotWeight * PUWeight );
							h_BEC_leadPhi->Fill( mu2.phi, TotWeight * PUWeight );
							h_BEC_subPt->Fill( mu1.Pt, TotWeight * PUWeight );
							h_BEC_subEta->Fill( mu1.eta, TotWeight * PUWeight );
							h_BEC_subPhi->Fill( mu1.phi, TotWeight * PUWeight );
						}
					}
				} // End of event selection

				// -- Event Selection without Rochester correction -- //
				vector< Muon > SelectedMuonCollection0;
				Bool_t isPassEventSelection0 = kFALSE;
				isPassEventSelection0 = analyzer->EventSelection(MuonCollection_noRoccoR, ntuple, &SelectedMuonCollection0);

				if( isMC == kTRUE && isPassEventSelection0 == kTRUE )
				{
					Muon mu1 = SelectedMuonCollection0[0];
					Muon mu2 = SelectedMuonCollection0[1];

					Double_t reco_M = (mu1.Momentum + mu2.Momentum).M();
					Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
					Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

					h_before_PUCorr_mass_fine->Fill( reco_M, TotWeight );
					h_before_PUCorr_diPt->Fill( reco_Pt, TotWeight );
					h_before_PUCorr_rapi->Fill( reco_rapi, TotWeight );
					h_before_PUCorr_pT->Fill( mu1.Pt, TotWeight );
					h_before_PUCorr_pT->Fill( mu2.Pt, TotWeight );
					h_before_PUCorr_eta->Fill( mu1.eta, TotWeight );
					h_before_PUCorr_eta->Fill( mu2.eta, TotWeight );
					h_before_PUCorr_phi->Fill( mu1.phi, TotWeight );
					h_before_PUCorr_phi->Fill( mu2.phi, TotWeight );

					h_before_RoccoR_mass_fine->Fill( reco_M, TotWeight * PUWeight );
					h_before_RoccoR_diPt->Fill( reco_Pt, TotWeight * PUWeight );
					h_before_RoccoR_rapi->Fill( reco_rapi, TotWeight * PUWeight );
					h_before_RoccoR_pT->Fill( mu1.Pt, TotWeight * PUWeight );
					h_before_RoccoR_pT->Fill( mu2.Pt, TotWeight * PUWeight );
					h_before_RoccoR_eta->Fill( mu1.eta, TotWeight * PUWeight );
					h_before_RoccoR_eta->Fill( mu2.eta, TotWeight * PUWeight );
					h_before_RoccoR_phi->Fill( mu1.phi, TotWeight * PUWeight );
					h_before_RoccoR_phi->Fill( mu2.phi, TotWeight * PUWeight );

					if( mu1.Pt > mu2.Pt ) //leading : mu1
					{
						h_before_PUCorr_leadPt->Fill( mu1.Pt, TotWeight );
						h_before_PUCorr_leadEta->Fill( mu1.eta, TotWeight );
						h_before_PUCorr_leadPhi->Fill( mu1.phi, TotWeight );
						h_before_PUCorr_subPt->Fill( mu2.Pt, TotWeight );
						h_before_PUCorr_subEta->Fill( mu2.eta, TotWeight );
						h_before_PUCorr_subPhi->Fill( mu2.phi, TotWeight );

						h_before_RoccoR_leadPt->Fill( mu1.Pt, TotWeight * PUWeight );
						h_before_RoccoR_leadEta->Fill( mu1.eta, TotWeight * PUWeight );
						h_before_RoccoR_leadPhi->Fill( mu1.phi, TotWeight * PUWeight );
						h_before_RoccoR_subPt->Fill( mu2.Pt, TotWeight * PUWeight );
						h_before_RoccoR_subEta->Fill( mu2.eta, TotWeight * PUWeight );
						h_before_RoccoR_subPhi->Fill( mu2.phi, TotWeight * PUWeight );
					}
					else //leading : mu2
					{
						h_before_PUCorr_leadPt->Fill( mu2.Pt, TotWeight );
						h_before_PUCorr_leadEta->Fill( mu2.eta, TotWeight );
						h_before_PUCorr_leadPhi->Fill( mu2.phi, TotWeight );
						h_before_PUCorr_subPt->Fill( mu1.Pt, TotWeight );
						h_before_PUCorr_subEta->Fill( mu1.eta, TotWeight );
						h_before_PUCorr_subPhi->Fill( mu1.phi, TotWeight );

						h_before_RoccoR_leadPt->Fill( mu2.Pt, TotWeight * PUWeight );
						h_before_RoccoR_leadEta->Fill( mu2.eta, TotWeight * PUWeight );
						h_before_RoccoR_leadPhi->Fill( mu2.phi, TotWeight * PUWeight );
						h_before_RoccoR_subPt->Fill( mu1.Pt, TotWeight * PUWeight );
						h_before_RoccoR_subEta->Fill( mu1.eta, TotWeight * PUWeight );
						h_before_RoccoR_subPhi->Fill( mu1.phi, TotWeight * PUWeight );
					}
				}

			} //End of if( isTriggered )

		} //End of event iteration

		h_mass_fine->Write();
		h_diPt->Write();
		h_rapi->Write();
		h_pT->Write();
		h_leadPt->Write();
		h_subPt->Write();
		h_eta->Write();
		h_leadEta->Write();
		h_subEta->Write();
		h_phi->Write();
		h_leadPhi->Write();
		h_subPhi->Write();
		h_leadEta_subEta->Write();

		if( isMC == kTRUE )
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

			h_Gen_M50to60_mass_fine->Write();
			h_Gen_M50to60_diPt->Write();
			h_Gen_M50to60_rapi->Write();
			h_Gen_M50to60_pT->Write();
			h_Gen_M50to60_eta->Write();
			h_Gen_M50to60_phi->Write();
			h_Gen_M50to60_leadPt->Write();
			h_Gen_M50to60_subPt->Write();
			h_Gen_M50to60_leadEta->Write();
			h_Gen_M50to60_subEta->Write();
			h_Gen_M50to60_leadPhi->Write();
			h_Gen_M50to60_subPhi->Write();

			h_Gen_M60to120_mass_fine->Write();
			h_Gen_M60to120_diPt->Write();
			h_Gen_M60to120_rapi->Write();
			h_Gen_M60to120_pT->Write();
			h_Gen_M60to120_eta->Write();
			h_Gen_M60to120_phi->Write();
			h_Gen_M60to120_leadPt->Write();
			h_Gen_M60to120_subPt->Write();
			h_Gen_M60to120_leadEta->Write();
			h_Gen_M60to120_subEta->Write();
			h_Gen_M60to120_leadPhi->Write();
			h_Gen_M60to120_subPhi->Write();
			h_Gen_M60to120_leadEta_subEta->Write();

			h_Gen_M120to200_mass_fine->Write();
			h_Gen_M120to200_diPt->Write();
			h_Gen_M120to200_rapi->Write();
			h_Gen_M120to200_pT->Write();
			h_Gen_M120to200_eta->Write();
			h_Gen_M120to200_phi->Write();
			h_Gen_M120to200_leadPt->Write();
			h_Gen_M120to200_subPt->Write();
			h_Gen_M120to200_leadEta->Write();
			h_Gen_M120to200_subEta->Write();
			h_Gen_M120to200_leadPhi->Write();
			h_Gen_M120to200_subPhi->Write();

			h_Gen_M200to400_mass_fine->Write();
			h_Gen_M200to400_diPt->Write();
			h_Gen_M200to400_rapi->Write();
			h_Gen_M200to400_pT->Write();
			h_Gen_M200to400_eta->Write();
			h_Gen_M200to400_phi->Write();
			h_Gen_M200to400_leadPt->Write();
			h_Gen_M200to400_subPt->Write();
			h_Gen_M200to400_leadEta->Write();
			h_Gen_M200to400_subEta->Write();
			h_Gen_M200to400_leadPhi->Write();
			h_Gen_M200to400_subPhi->Write();

			h_Gen_M400to500_mass_fine->Write();
			h_Gen_M400to500_diPt->Write();
			h_Gen_M400to500_rapi->Write();
			h_Gen_M400to500_pT->Write();
			h_Gen_M400to500_eta->Write();
			h_Gen_M400to500_phi->Write();
			h_Gen_M400to500_leadPt->Write();
			h_Gen_M400to500_subPt->Write();
			h_Gen_M400to500_leadEta->Write();
			h_Gen_M400to500_subEta->Write();
			h_Gen_M400to500_leadPhi->Write();
			h_Gen_M400to500_subPhi->Write();

			h_Reco_mass_fine->Write();
			h_Reco_diPt->Write();
			h_Reco_rapi->Write();
			h_Reco_pT->Write();
			h_Reco_eta->Write();
			h_Reco_phi->Write();
			h_Reco_leadPt->Write();
			h_Reco_subPt->Write();
			h_Reco_leadEta->Write();
			h_Reco_subEta->Write();
			h_Reco_leadPhi->Write();
			h_Reco_subPhi->Write();

			h_RecoID_mass_fine->Write();
			h_RecoID_diPt->Write();
			h_RecoID_rapi->Write();
			h_RecoID_pT->Write();
			h_RecoID_eta->Write();
			h_RecoID_phi->Write();
			h_RecoID_leadPt->Write();
			h_RecoID_subPt->Write();
			h_RecoID_leadEta->Write();
			h_RecoID_subEta->Write();
			h_RecoID_leadPhi->Write();
			h_RecoID_subPhi->Write();

			h_RecoIDIso_mass_fine->Write();
			h_RecoIDIso_diPt->Write();
			h_RecoIDIso_rapi->Write();
			h_RecoIDIso_pT->Write();
			h_RecoIDIso_eta->Write();
			h_RecoIDIso_phi->Write();
			h_RecoIDIso_leadPt->Write();
			h_RecoIDIso_subPt->Write();
			h_RecoIDIso_leadEta->Write();
			h_RecoIDIso_subEta->Write();
			h_RecoIDIso_leadPhi->Write();
			h_RecoIDIso_subPhi->Write();

			h_before_PUCorr_mass_fine->Write();
			h_before_PUCorr_diPt->Write();
			h_before_PUCorr_rapi->Write();
			h_before_PUCorr_pT->Write();
			h_before_PUCorr_leadPt->Write();
			h_before_PUCorr_subPt->Write();
			h_before_PUCorr_eta->Write();
			h_before_PUCorr_leadEta->Write();
			h_before_PUCorr_subEta->Write();
			h_before_PUCorr_phi->Write();
			h_before_PUCorr_leadPhi->Write();
			h_before_PUCorr_subPhi->Write();

			h_before_RoccoR_mass_fine->Write();
			h_before_RoccoR_diPt->Write();
			h_before_RoccoR_rapi->Write();
			h_before_RoccoR_pT->Write();
			h_before_RoccoR_leadPt->Write();
			h_before_RoccoR_subPt->Write();
			h_before_RoccoR_eta->Write();
			h_before_RoccoR_leadEta->Write();
			h_before_RoccoR_subEta->Write();
			h_before_RoccoR_phi->Write();
			h_before_RoccoR_leadPhi->Write();
			h_before_RoccoR_subPhi->Write();

			h_BEC_mass_fine->Write();
			h_BEC_diPt->Write();
			h_BEC_rapi->Write();
			h_BEC_pT->Write();
			h_BEC_eta->Write();
			h_BEC_phi->Write();
			h_BEC_leadPt->Write();
			h_BEC_subPt->Write();
			h_BEC_leadEta->Write();
			h_BEC_subEta->Write();
			h_BEC_leadPhi->Write();
			h_BEC_subPhi->Write();

			h_BtoF_mass_fine->Write();
			h_BtoF_diPt->Write();
			h_BtoF_rapi->Write();
			h_BtoF_pT->Write();
			h_BtoF_leadPt->Write();
			h_BtoF_subPt->Write();
			h_BtoF_eta->Write();
			h_BtoF_leadEta->Write();
			h_BtoF_subEta->Write();
			h_BtoF_phi->Write();
			h_BtoF_leadPhi->Write();
			h_BtoF_subPhi->Write();

			h_BtoF_Reco_mass_fine->Write();
			h_BtoF_Reco_diPt->Write();
			h_BtoF_Reco_rapi->Write();
			h_BtoF_Reco_pT->Write();
			h_BtoF_Reco_eta->Write();
			h_BtoF_Reco_phi->Write();
			h_BtoF_Reco_leadPt->Write();
			h_BtoF_Reco_subPt->Write();
			h_BtoF_Reco_leadEta->Write();
			h_BtoF_Reco_subEta->Write();
			h_BtoF_Reco_leadPhi->Write();
			h_BtoF_Reco_subPhi->Write();

			h_BtoF_RecoID_mass_fine->Write();
			h_BtoF_RecoID_diPt->Write();
			h_BtoF_RecoID_rapi->Write();
			h_BtoF_RecoID_pT->Write();
			h_BtoF_RecoID_eta->Write();
			h_BtoF_RecoID_phi->Write();
			h_BtoF_RecoID_leadPt->Write();
			h_BtoF_RecoID_subPt->Write();
			h_BtoF_RecoID_leadEta->Write();
			h_BtoF_RecoID_subEta->Write();
			h_BtoF_RecoID_leadPhi->Write();
			h_BtoF_RecoID_subPhi->Write();

			h_BtoF_RecoIDIso_mass_fine->Write();
			h_BtoF_RecoIDIso_diPt->Write();
			h_BtoF_RecoIDIso_rapi->Write();
			h_BtoF_RecoIDIso_pT->Write();
			h_BtoF_RecoIDIso_eta->Write();
			h_BtoF_RecoIDIso_phi->Write();
			h_BtoF_RecoIDIso_leadPt->Write();
			h_BtoF_RecoIDIso_subPt->Write();
			h_BtoF_RecoIDIso_leadEta->Write();
			h_BtoF_RecoIDIso_subEta->Write();
			h_BtoF_RecoIDIso_leadPhi->Write();
			h_BtoF_RecoIDIso_subPhi->Write();

			h_GtoH_mass_fine->Write();
			h_GtoH_diPt->Write();
			h_GtoH_rapi->Write();
			h_GtoH_pT->Write();
			h_GtoH_leadPt->Write();
			h_GtoH_subPt->Write();
			h_GtoH_eta->Write();
			h_GtoH_leadEta->Write();
			h_GtoH_subEta->Write();
			h_GtoH_phi->Write();
			h_GtoH_leadPhi->Write();
			h_GtoH_subPhi->Write();

			h_GtoH_Reco_mass_fine->Write();
			h_GtoH_Reco_diPt->Write();
			h_GtoH_Reco_rapi->Write();
			h_GtoH_Reco_pT->Write();
			h_GtoH_Reco_eta->Write();
			h_GtoH_Reco_phi->Write();
			h_GtoH_Reco_leadPt->Write();
			h_GtoH_Reco_subPt->Write();
			h_GtoH_Reco_leadEta->Write();
			h_GtoH_Reco_subEta->Write();
			h_GtoH_Reco_leadPhi->Write();
			h_GtoH_Reco_subPhi->Write();

			h_GtoH_RecoID_mass_fine->Write();
			h_GtoH_RecoID_diPt->Write();
			h_GtoH_RecoID_rapi->Write();
			h_GtoH_RecoID_pT->Write();
			h_GtoH_RecoID_eta->Write();
			h_GtoH_RecoID_phi->Write();
			h_GtoH_RecoID_leadPt->Write();
			h_GtoH_RecoID_subPt->Write();
			h_GtoH_RecoID_leadEta->Write();
			h_GtoH_RecoID_subEta->Write();
			h_GtoH_RecoID_leadPhi->Write();
			h_GtoH_RecoID_subPhi->Write();

			h_GtoH_RecoIDIso_mass_fine->Write();
			h_GtoH_RecoIDIso_diPt->Write();
			h_GtoH_RecoIDIso_rapi->Write();
			h_GtoH_RecoIDIso_pT->Write();
			h_GtoH_RecoIDIso_eta->Write();
			h_GtoH_RecoIDIso_phi->Write();
			h_GtoH_RecoIDIso_leadPt->Write();
			h_GtoH_RecoIDIso_subPt->Write();
			h_GtoH_RecoIDIso_leadEta->Write();
			h_GtoH_RecoIDIso_subEta->Write();
			h_GtoH_RecoIDIso_leadPhi->Write();
			h_GtoH_RecoIDIso_subPhi->Write();
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

