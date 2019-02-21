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
#include "../../HEADER/DYAnalyzer_for_SMP_17_010.h"

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
// -- Update "LeadingMuonEtaSF" to "LeadEtaCorr" : 21 May. 2018 -- //
// -- Remove "LeadEtaCorr" : 18 Jun. 2018 -- //
// -- Leading and sub-leading are modified : 18 Jun. 2018 -- //
// -- Update code with "PtCut" and "2D distribution of mass-rapdity" : 19 Jun. 2018 -- //
// -- Add "2D distribution of pT-eta" : 19 Jun. 2018 -- //
// -- Introduce Z-peak mass cut : 19 Jun. 2018 -- //
// -- Add "ZToMuMu_powheg" : 19 Jun. 2018 -- //
// -- Use WJetsToLNu_amcatnlo samples: 14 Nov. 2018 -- //
// -- Add type 52: 15 Nov. 2018 -- //
//--------------------------------------------------------------------------
//    For cross check with SMP-17-010
//--------------------------------------------------------------------------
// -- 2018.11.20: Copied from "MuMu_TightID_PFIso_NewMuonSF_20181012.C"
// -- 2018.11.21: Use Muon POG's official Medium ID SF
// -- 2018.11.28: Introduce SMP-17-010's efficiency SF (except trigger efficiency)
// -- 2018.11.30: Change the trigger path (Add Mu50!)
//                Trigger efficiency is updated using SMP-17-010's result
// -- 2018.12.10: Add "IsoMu22_OR_IsoTkMu22_OR_IsoMu24_OR_IsoTkMu24_OR_Mu50_for_SMP_17_010"
// -- 2018.12.12: Update PU reweighting by using SMP-17-010's values
// -- 2018.12.13: Prefiring weight was tested
// -- 2018.12.19: Add prefiring unc.
//                Remove gen-level and unused histograms
// -- 2018.12.20: PVz reweighting is tested
// -- 2019.01.07: PVz reweighting is removed
//                Rochester correction is updated by using gen-level muon pt
// -- 2019.01.14: DoubleMuon PD & trigger are tested
// -- 2019.01.17: Convention for naming output root file is slightly changed
// -- 2019.01.23: Vetoing using loose muons was tested.(=> No effect, so keep applying)
// -- 2019.01.24: Range of Z peak was modified using Z mass instead of 90 GeV
// -- 2019.01.31: Change muon IDs using ID flags
//
//
//void MuMu_for_SMP_17_010(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24_for_SMP_17_010")
//void MuMu_for_SMP_17_010(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24_OR_Mu50_for_SMP_17_010")
void MuMu_for_SMP_17_010_step3(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu22_OR_IsoTkMu22_OR_IsoMu24_OR_IsoTkMu24_OR_Mu50_for_SMP_17_010")
{
	gROOT->SetBatch(kTRUE);
	TString NtupleLocation = gSystem->Getenv("DM_DATA_PATH");
	TString BaseDir= gSystem->Getenv("DM_BASE_PATH");

	// -- Run2016 luminosity [/pb] -- //
	Double_t lumi = Lumi; //BtoH

	TString DataLocation, Type, DataType;
	//DataType = "DoubleMuon_";
	DataType = "SingleMuon_";

	// -- Data samples -- //
	if( type == 1 ) DataLocation = DataType+"Run2016B";
	else if( type == 2 ) DataLocation = DataType+"Run2016C";
	else if( type == 3 ) DataLocation = DataType+"Run2016D";
	else if( type == 4 ) DataLocation = DataType+"Run2016E";
	else if( type == 5 ) DataLocation = DataType+"Run2016F";
	else if( type == 6 ) DataLocation = DataType+"Run2016G";
	else if( type == 7 ) DataLocation = DataType+"Run2016H";
	// -- Signal MC samples -- //
	else if( type == 11 ) Type = "DYMuMu_M10to50";
	else if( type == 12 ) Type = "DYMuMu_M50toInf";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	else if( type == 31 ) Type = "DYTauTau_M10to50";
	else if( type == 32 ) Type = "DYTauTau_M50toInf";
	else if( type == 41 ) Type = "VVnST";
	else if( type == 51 ) Type = "WJetsToLNu_amcatnlo";
	else if( type == 52 ) Type = "WJetsToLNu_amcatnlo_ext2v5";

	Bool_t isMC;
	if( type < 10  )
	{
		isMC = kFALSE;
		Type = "Data (" + DataLocation + ")";
		if( type == 7 ) DataLocation += "ver2";
	}
	else
	{
		isMC = kTRUE;
		DataType = "MC_";
	}

	Double_t Div = 1;
	if( type == 1 || type == 6 || type == 7 || type == 11 || type == 31 || type == 51 || type == 52 ) Div = 2;

	TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	cout << "Type: " << Type << endl;

	TStopwatch totaltime;
	totaltime.Start();

	DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

	// -- Pile-up setup -- //
	//analyzer->SetupPileUpReWeighting_80X( isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );
	analyzer->SetupPileUpReWeighting_for_SMP_17_010( isMC );

	// -- Rochester setup -- //
	TRandom3 *r1 = new TRandom3(0);
	RoccoR rc("../../TOOL/RoccoR/rcdata.2016.v3");

	// -- Efficiency SF setup -- //
	if( isMC == kTRUE )
	{
	//	analyzer->SetupEfficiencyScaleFactor_BtoF();
	//	analyzer->SetupEfficiencyScaleFactor_GtoH();
		analyzer->SetupEfficiencyScaleFactor_BtoF_SMP_17_010();
		analyzer->SetupEfficiencyScaleFactor_GtoH_SMP_17_010();
	}

	// -- PVz setup -- //
	//analyzer->SetupPVzReWeighting_80X( isMC, "PVz.root" );
	//analyzer->SetupPVzReWeighting_80X( isMC, "PVz_for_muon.root" );



	// -- Output ROOTFile -- //	
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181120_MuMu_for_SMP_17_010_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181120_MuMu_for_SMP_17_010_with_Official_SF_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181128_MuMu_for_SMP_17_010_with_SMP_17_010_SF_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181130_MuMu_for_SMP_17_010_with_All_SMP_17_010_Eff_SF_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181210_MuMu_for_SMP_17_010_with_All_SMP_17_010_Eff_SF_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181212_MuMu_for_SMP_17_010_with_All_SMP_17_010_Eff_SF_and_PU_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181213_MuMu_for_SMP_17_010_with_All_SMP_17_010_Eff_SF_and_PU_and_PrefiringWeight_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181219_MuMu_for_SMP_17_010_with_All_SMP_17_010_Eff_SF_and_PU_and_PrefiringWeight_with_unc_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20181220_MuMu_for_SMP_17_010_with_All_SMP_17_010_Eff_SF_and_PU_and_PrefiringWeight_with_unc_and_PVz_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190107_MuMu_for_SMP_17_010_with_All_SMP_17_010_Corrections_and_Updated_RoccoR_"
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190117_MuMu_for_SMP_17_010_with_All_SMP_17_010_Corrections_and_Updated_RoccoR_on_"+DataType
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190123_MuMu_for_SMP_17_010_without_veto_on_"+DataType
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190124_MuMu_for_SMP_17_010_with_Zmass_on_"+DataType
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190131_MuMu_IDtest_for_SMP_17_010_on_"+DataType
	TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190213_MuMu_for_SMP_17_010_step3_on_"+DataType
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
			//TString version = "v2.4";
			//TString version = "v2.5"; // for prefiring weight
			TString version = "v2.6"; // for Muon ID flags

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
			//TString version = "v2.4";
			TString version = "v2.6"; // for Muon ID flags

			if( remainder == 9999 )
			{
				chain->Add(NtupleLocation+"/"+version+"/"+DataLocation+"/*.root");
				if(type==7) chain->Add(NtupleLocation+"/"+version+"/"+DataType+"Run2016Hver3/*.root");
			}
			else
			{
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/"+version+"/"+DataLocation+"/*_"+TString::Itoa(ii,10)+".root");
				if(type==7 && remainder==0) chain->Add(NtupleLocation+"/"+version+"/"+DataType+"Run2016Hver3/*.root");
			}
		}

		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		ntuple->TurnOnBranches_MuonIDFlags(); // for Muon ID flags
		if( isMC == kTRUE )
		{
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			//ntuple->TurnOnBranches_GenOthers(); // for quarks
			ntuple->TurnOnBranches_prefiring(); // for prefiring weight
		}

		////////////////////////////
		// -- Making Histogram -- //
		////////////////////////////
		vector<TString> wtype; wtype.push_back(""); wtype.push_back("up_"); wtype.push_back("down_");

		// -- Reco-level -- //
		TH1D *h_mass[3];
		TH1D *h_mass_fine[3];
		TH1D *h_diPt[3];
		TH1D *h_rapi[3];
		TH1D *h_pT[3];
		TH1D *h_leadPt[3];
		TH1D *h_subPt[3];
		TH1D *h_eta[3];
		TH1D *h_leadEta[3];
		TH1D *h_subEta[3];
		TH1D *h_phi[3];
		TH1D *h_leadPhi[3];
		TH1D *h_subPhi[3];
		//TH1D *h_PVz[3];

		// -- Reco-level in Z-peak -- //
		TH1D *h_Zpeak_mass[3];
		TH1D *h_Zpeak_mass_fine[3];
		TH1D *h_Zpeak_diPt[3];
		TH1D *h_Zpeak_rapi[3];
		TH1D *h_Zpeak_pT[3];
		TH1D *h_Zpeak_leadPt[3];
		TH1D *h_Zpeak_subPt[3];
		TH1D *h_Zpeak_eta[3];
		TH1D *h_Zpeak_leadEta[3];
		TH1D *h_Zpeak_subEta[3];
		TH1D *h_Zpeak_phi[3];
		TH1D *h_Zpeak_leadPhi[3];
		TH1D *h_Zpeak_subPhi[3];
		//TH1D *h_Zpeak_PVz[3];

		for(int i=0;i<3;i++)
		{	
			h_mass[i] = new TH1D("h_mass_"+wtype[i]+Tag[i_tup], "", 43, massbins);
			h_mass_fine[i] = new TH1D("h_mass_fine_"+wtype[i]+Tag[i_tup], "", 10000, 0, 10000);
			h_diPt[i] = new TH1D("h_diPt_"+wtype[i]+Tag[i_tup], "", 10000, 0, 10000);
			h_rapi[i] = new TH1D("h_rapi_"+wtype[i]+Tag[i_tup], "", 1000, -5, 5);
			h_pT[i] = new TH1D("h_pT_"+wtype[i]+Tag[i_tup], "", 10000, 0, 10000);
			h_leadPt[i] = (TH1D*)h_pT[i]->Clone("h_leadPt_"+wtype[i]+Tag[i_tup]);
			h_subPt[i] = (TH1D*)h_pT[i]->Clone("h_subPt_"+wtype[i]+Tag[i_tup]);
			h_eta[i] = new TH1D("h_eta_"+wtype[i]+Tag[i_tup], "", 1000, -5, 5);
			h_leadEta[i] = (TH1D*)h_eta[i]->Clone("h_leadEta_"+wtype[i]+Tag[i_tup]);
			h_subEta[i] = (TH1D*)h_eta[i]->Clone("h_subEta_"+wtype[i]+Tag[i_tup]);
			h_phi[i] = new TH1D("h_phi_"+wtype[i]+Tag[i_tup], "", 100, -5, 5);
			h_leadPhi[i] = (TH1D*)h_phi[i]->Clone("h_leadPhi_"+wtype[i]+Tag[i_tup]);
			h_subPhi[i] = (TH1D*)h_phi[i]->Clone("h_subPhi_"+wtype[i]+Tag[i_tup]);
			//h_PVz[i] = new TH1D("h_PVz_"+wtype[i]+Tag[i_tup], "", 100, -25, 25);

			h_Zpeak_mass[i] = (TH1D*)h_mass[i]->Clone("h_Zpeak_mass_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_mass_fine[i] = (TH1D*)h_mass_fine[i]->Clone("h_Zpeak_mass_fine_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_diPt[i] = (TH1D*)h_diPt[i]->Clone("h_Zpeak_diPt_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_rapi[i] = (TH1D*)h_rapi[i]->Clone("h_Zpeak_rapi_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_pT[i] = (TH1D*)h_pT[i]->Clone("h_Zpeak_pT_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_leadPt[i] = (TH1D*)h_pT[i]->Clone("h_Zpeak_leadPt_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_subPt[i] = (TH1D*)h_pT[i]->Clone("h_Zpeak_subPt_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_eta[i] = (TH1D*)h_eta[i]->Clone("h_Zpeak_eta_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_leadEta[i] = (TH1D*)h_eta[i]->Clone("h_Zpeak_leadEta_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_subEta[i] = (TH1D*)h_eta[i]->Clone("h_Zpeak_subEta_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_phi[i] = (TH1D*)h_phi[i]->Clone("h_Zpeak_phi_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_leadPhi[i] = (TH1D*)h_phi[i]->Clone("h_Zpeak_leadPhi_"+wtype[i]+Tag[i_tup]);
			h_Zpeak_subPhi[i] = (TH1D*)h_phi[i]->Clone("h_Zpeak_subPhi_"+wtype[i]+Tag[i_tup]);
			//h_Zpeak_PVz[i] = (TH1D*)h_PVz[i]->Clone("h_Zpeak_PVz_"+wtype[i]+Tag[i_tup]);
		}

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
			//if( isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( ntuple->nPileUp );
			if( isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_for_SMP_17_010( ntuple->nPileUp );

			// -- Prefiring weight -- //
			Double_t PrefiringWeight[3];

			for(int i=0;i<3;i++)
				PrefiringWeight[i] = 1;

/*			if( isMC == kTRUE )
			{
				PrefiringWeight[0] = ntuple->_prefiringweight;
				PrefiringWeight[1] = ntuple->_prefiringweightup;
				PrefiringWeight[2] = ntuple->_prefiringweightdown;
			}*/

			// -- efficiency weights -- //
			Double_t effweight = 1, effweight_BtoF = 1, effweight_GtoH = 1;

			// -- PVz-Reweighting -- //
			Double_t PVzWeight = 1;
			//if( isMC == kTRUE ) PVzWeight = analyzer->PVzWeightValue_80X( ntuple->PVz );


			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			Bool_t GenFlag_top = kTRUE;
			// -- Separate ttbar samples -- //
			//GenFlag_top = kFALSE;
			//vector<GenOthers> GenTopCollection;
			//GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

			if( GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				SumWeight_Separated += GenWeight; // <- This part must exist!! 

				// -- Top Pt Reweighting -- //
				//if( isTopPtReweighting == 1 && Tag[i_tup].Contains("ttbar") )
				//{
				//	GenOthers t1 = GenTopCollection[0];
				//	GenOthers t2 = GenTopCollection[1];
				//
				//	Double_t SF1 = exp(0.0615 - 0.0005*(t1.Pt));
				//	Double_t SF2 = exp(0.0615 - 0.0005*(t2.Pt));
				//	GenWeight = GenWeight*sqrt(SF1*SF2);
				//}
			}

			// -- Normalization -- //
			Double_t TotWeight = GenWeight;
			if( isMC == kTRUE ) TotWeight = (lumi*Xsec[i_tup]/nEvents[i_tup])*GenWeight;


			////////////////////////////////
			// -- Reco level selection -- //
			////////////////////////////////
			Bool_t TriggerFlag = kFALSE;
			TriggerFlag = ntuple->isTriggered( analyzer->HLT );

			if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				vector< Muon > MuonCollection;
				Int_t NLeptons = ntuple->nMuon;
				for(Int_t i_reco=0; i_reco<NLeptons; i_reco++)
				{
					Muon mu;
					mu.FillFromNtuple(ntuple, i_reco);

					// -- Rochester correction ----------------------------------------------------------------------
					Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
					Int_t s, m;
						
					if( isMC == kFALSE ) // Data
						SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
					else // MC
					{
						Double_t genPt = analyzer->GenMuonPt("fromHardProcessFinalState", ntuple, mu);

						if( genPt > 0 )
							SF = rc.kScaleFromGenMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, genPt, rndm[0], s=0, m=0);
						else
							SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
					}

					mu.Pt = SF*mu.Pt;

					mu.Momentum.SetPtEtaPhiM( mu.Pt, mu.eta, mu.phi, M_Mu );
					// ----------------------------------------------------------------------------------------------

					MuonCollection.push_back( mu );
				}	

				// -- Event Selection -- //
				vector< Muon > SelectedMuonCollection;
				Bool_t isPassEventSelection = kFALSE;
				//isPassEventSelection = analyzer->EventSelection(MuonCollection, ntuple, &SelectedMuonCollection);
				isPassEventSelection = analyzer->EventSelection_for_SMP_17_010(MuonCollection, ntuple, &SelectedMuonCollection);

				if( isPassEventSelection == kTRUE )
				{
					Muon mu1, mu2; // lead: 1, sub: 2
					if( SelectedMuonCollection[0].Pt > SelectedMuonCollection[1].Pt )
					{
						mu1 = SelectedMuonCollection[0];
						mu2 = SelectedMuonCollection[1];
					}
					else
					{
						mu1 = SelectedMuonCollection[1];
						mu2 = SelectedMuonCollection[0];
					}

					Double_t reco_M = (mu1.Momentum + mu2.Momentum).M();
					Double_t reco_Pt = (mu1.Momentum + mu2.Momentum).Pt();
					Double_t reco_rapi = (mu1.Momentum + mu2.Momentum).Rapidity();

					// -- Efficiency scale factor -- //
					if( isMC == kTRUE )
					{
					//	effweight_BtoF = analyzer->EfficiencySF_EventWeight_HLT_BtoF( mu1, mu2 );
					//	effweight_GtoH = analyzer->EfficiencySF_EventWeight_HLT_GtoH( mu1, mu2 );
						effweight_BtoF = analyzer->EfficiencySF_BtoF_SMP_17_010( mu1, mu2 );
						effweight_GtoH = analyzer->EfficiencySF_GtoH_SMP_17_010( mu1, mu2 );
						effweight = (Lumi_BtoF * effweight_BtoF + Lumi_GtoH * effweight_GtoH) / Lumi;
					}

					for(int i=0;i<3;i++)
					{
						h_mass[i]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_mass_fine[i]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_diPt[i]->Fill( reco_Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_rapi[i]->Fill( reco_rapi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_pT[i]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_pT[i]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_eta[i]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_eta[i]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_phi[i]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_phi[i]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );

						h_leadPt[i]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_leadEta[i]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_leadPhi[i]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_subPt[i]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_subEta[i]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						h_subPhi[i]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						//h_PVz[i]->Fill( ntuple->PVz, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );

						//if( 60 < reco_M && reco_M < 120 )
						//if( 76 < reco_M && reco_M < 106 )
						if( fabs( reco_M - 91.1876 ) < 15 )
						{
							h_Zpeak_mass[i]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_mass_fine[i]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_diPt[i]->Fill( reco_Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_rapi[i]->Fill( reco_rapi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_pT[i]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_pT[i]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_eta[i]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_eta[i]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_phi[i]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_phi[i]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );

							h_Zpeak_leadPt[i]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_leadEta[i]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_leadPhi[i]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_subPt[i]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_subEta[i]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							h_Zpeak_subPhi[i]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
							//h_Zpeak_PVz[i]->Fill( ntuple->PVz, TotWeight * PUWeight * effweight * PrefiringWeight[i] * PVzWeight );
						} //End of Z-peak
					}

				} //End of event selection

			} //End of if( isTriggered )

		} //End of event iteration

		for(int i=0;i<3;i++)
		{
			h_mass[i]->Write();
			h_mass_fine[i]->Write();
			h_diPt[i]->Write();
			h_rapi[i]->Write();
			h_pT[i]->Write();
			h_leadPt[i]->Write();
			h_subPt[i]->Write();
			h_eta[i]->Write();
			h_leadEta[i]->Write();
			h_subEta[i]->Write();
			h_phi[i]->Write();
			h_leadPhi[i]->Write();
			h_subPhi[i]->Write();
			//h_PVz[i]->Write();

			h_Zpeak_mass[i]->Write();
			h_Zpeak_mass_fine[i]->Write();
			h_Zpeak_diPt[i]->Write();
			h_Zpeak_rapi[i]->Write();
			h_Zpeak_pT[i]->Write();
			h_Zpeak_leadPt[i]->Write();
			h_Zpeak_subPt[i]->Write();
			h_Zpeak_eta[i]->Write();
			h_Zpeak_leadEta[i]->Write();
			h_Zpeak_subEta[i]->Write();
			h_Zpeak_phi[i]->Write();
			h_Zpeak_leadPhi[i]->Write();
			h_Zpeak_subPhi[i]->Write();
			//h_Zpeak_PVz[i]->Write();
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

