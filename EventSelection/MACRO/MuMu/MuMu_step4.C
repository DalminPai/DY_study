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
//#include "../../HEADER/DYAnalyzer_TightID_PFIso_v20180521.h"
//#include "../../HEADER/DYAnalyzer_TightID_PFIso_NewMuonSF.h"
#include "../../HEADER/DYAnalyzer_TightID_PFIso_v20190212.h"

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
// -- Prefiring weight was tested: 13 Dec. 2018 -- //
//--------------------------------------------------------------------------
//    For KP's validation task
//--------------------------------------------------------------------------
// -- 2019.02.08: Copied from "MuMu_TightID_PFIso_NewMuonSF_20181012.C"
// -- 2019.02.13: Copied from "MuMu_for_validation.C"
// -- 2019.02.14: Use PU reweighting of SMP-17-010
//                Copied from "MuMu_step3.C"
//
//
void MuMu_step4(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24")
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
	else if( type == 10 ) Type = "DYMuMu_M50toInf";
	else if( type == 11 ) Type = "DYMuMu_M10to50";
	else if( type == 12 ) Type = "DYMuMu_M50to100";
	else if( type == 13 ) Type = "DYMuMu_M100toInf";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	else if( type == 23 ) Type = "ttbar_Mto700";
	else if( type == 24 ) Type = "ttbarBackup_Mto700";
	else if( type == 25 ) Type = "ttbar_M700toInf";
	else if( type == 31 ) Type = "DYTauTau_M10to50";
	else if( type == 32 ) Type = "DYTauTau_M50toInf";
	else if( type == 41 ) Type = "VVnST";
	//else if( type == 51 ) Type = "WJetsToLNu";
	else if( type == 51 ) Type = "WJetsToLNu_amcatnlo";
	else if( type == 52 ) Type = "WJetsToLNu_amcatnlo_ext2v5";

	Bool_t isMC = kTRUE;
	if( type < 10  )
	{
		isMC = kFALSE;
		Type = "Data (" + DataLocation + ")";
		if( type == 7 ) DataLocation += "ver2";
	}

	Double_t Div = 1;
	//if( type == 1 || type == 6 || type == 7 || type == 11 || type == 31 || type == 51 ) Div = 2;
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
		//analyzer->SetupEfficiencyScaleFactor_BtoF();
		//analyzer->SetupEfficiencyScaleFactor_GtoH();
		analyzer->SetupEfficiencyScaleFactor_BtoF_new();
		analyzer->SetupEfficiencyScaleFactor_GtoH_new();
	}

	// -- Output ROOTFile -- //	
	TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20190213_MuMu_step4_MITpu_"
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
			TString version = "v2.6";

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
			TString version = "v2.6";

			if( remainder == 9999 )
			{
				chain->Add(NtupleLocation+"/"+version+"/"+DataLocation+"/*.root");
				if(type==7) chain->Add(NtupleLocation+"/"+version+"/SingleMuon_Run2016Hver3/*.root");
			}
			else
			{
				for(Double_t ii=1; ii<=1500; ii++)
					if(ii - TMath::Floor(ii/Div) * Div == remainder)
						chain->Add(NtupleLocation+"/"+version+"/"+DataLocation+"/*_"+TString::Itoa(ii,10)+".root");
				if(type==7 && remainder==0) chain->Add(NtupleLocation+"/"+version+"/SingleMuon_Run2016Hver3/*.root");
			}
		}

		NtupleHandle *ntuple = new NtupleHandle( chain );
		ntuple->TurnOnBranches_Muon();
		if( isMC == kTRUE )
		{
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			ntuple->TurnOnBranches_prefiring(); // for prefiring weight
			if( Tag[i_tup].Contains("ttbar") && Tag[i_tup].Contains("M") )
				ntuple->TurnOnBranches_GenOthers(); // for quarks
		}

		////////////////////////////
		// -- Making Histogram -- //
		////////////////////////////
		vector<TString> corrtype;
		corrtype.push_back("");
		corrtype.push_back("raw_");
		const Int_t ncorr = (Int_t)corrtype.size();

		// -- Reco-level -- //
		TH1D *h_mass[ncorr];
		TH1D *h_mass_fine[ncorr];
		TH1D *h_diPt[ncorr];
		TH1D *h_rapi[ncorr];
		TH1D *h_pT[ncorr];
		TH1D *h_leadPt[ncorr];
		TH1D *h_subPt[ncorr];
		TH1D *h_eta[ncorr];
		TH1D *h_leadEta[ncorr];
		TH1D *h_subEta[ncorr];
		TH1D *h_phi[ncorr];
		TH1D *h_leadPhi[ncorr];
		TH1D *h_subPhi[ncorr];
		TH1D *h_nVTX_woPU[ncorr];
		TH1D *h_nVTX_withPU[ncorr];

		// -- Reco-level in Z-peak -- //
		TH1D *h_Zpeak_mass[ncorr];
		TH1D *h_Zpeak_mass_fine[ncorr];
		TH1D *h_Zpeak_diPt[ncorr];
		TH1D *h_Zpeak_rapi[ncorr];
		TH1D *h_Zpeak_pT[ncorr];
		TH1D *h_Zpeak_leadPt[ncorr];
		TH1D *h_Zpeak_subPt[ncorr];
		TH1D *h_Zpeak_eta[ncorr];
		TH1D *h_Zpeak_leadEta[ncorr];
		TH1D *h_Zpeak_subEta[ncorr];
		TH1D *h_Zpeak_phi[ncorr];
		TH1D *h_Zpeak_leadPhi[ncorr];
		TH1D *h_Zpeak_subPhi[ncorr];
		TH1D *h_Zpeak_nVTX_woPU[ncorr];
		TH1D *h_Zpeak_nVTX_withPU[ncorr];

		for(int i=0;i<ncorr;i++)
		{	
			h_mass[i] = new TH1D("h_"+corrtype[i]+"mass_"+Tag[i_tup], "", 43, massbins);
			h_mass_fine[i] = new TH1D("h_"+corrtype[i]+"mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
			h_diPt[i] = new TH1D("h_"+corrtype[i]+"diPt_"+Tag[i_tup], "", 10000, 0, 10000);
			h_rapi[i] = new TH1D("h_"+corrtype[i]+"rapi_"+Tag[i_tup], "", 1000, -5, 5);
			h_pT[i] = new TH1D("h_"+corrtype[i]+"pT_"+Tag[i_tup], "", 10000, 0, 10000);
			h_leadPt[i] = (TH1D*)h_pT[i]->Clone("h_"+corrtype[i]+"leadPt_"+Tag[i_tup]);
			h_subPt[i] = (TH1D*)h_pT[i]->Clone("h_"+corrtype[i]+"subPt_"+Tag[i_tup]);
			h_eta[i] = new TH1D("h_"+corrtype[i]+"eta_"+Tag[i_tup], "", 1000, -5, 5);
			h_leadEta[i] = (TH1D*)h_eta[i]->Clone("h_"+corrtype[i]+"leadEta_"+Tag[i_tup]);
			h_subEta[i] = (TH1D*)h_eta[i]->Clone("h_"+corrtype[i]+"subEta_"+Tag[i_tup]);
			h_phi[i] = new TH1D("h_"+corrtype[i]+"phi_"+Tag[i_tup], "", 100, -5, 5);
			h_leadPhi[i] = (TH1D*)h_phi[i]->Clone("h_"+corrtype[i]+"leadPhi_"+Tag[i_tup]);
			h_subPhi[i] = (TH1D*)h_phi[i]->Clone("h_"+corrtype[i]+"subPhi_"+Tag[i_tup]);
			h_nVTX_woPU[i] = new TH1D("h_"+corrtype[i]+"nVTX_woPU_"+Tag[i_tup], "", 1000, 0, 1000);
			h_nVTX_withPU[i] = (TH1D*)h_nVTX_woPU[i]->Clone("h_"+corrtype[i]+"nVTX_withPU_"+Tag[i_tup]);

			h_Zpeak_mass[i] = (TH1D*)h_mass[i]->Clone("h_"+corrtype[i]+"Zpeak_mass_"+Tag[i_tup]);
			h_Zpeak_mass_fine[i] = (TH1D*)h_mass_fine[i]->Clone("h_"+corrtype[i]+"Zpeak_mass_fine_"+Tag[i_tup]);
			h_Zpeak_diPt[i] = (TH1D*)h_diPt[i]->Clone("h_"+corrtype[i]+"Zpeak_diPt_"+Tag[i_tup]);
			h_Zpeak_rapi[i] = (TH1D*)h_rapi[i]->Clone("h_"+corrtype[i]+"Zpeak_rapi_"+Tag[i_tup]);
			h_Zpeak_pT[i] = (TH1D*)h_pT[i]->Clone("h_"+corrtype[i]+"Zpeak_pT_"+Tag[i_tup]);
			h_Zpeak_leadPt[i] = (TH1D*)h_pT[i]->Clone("h_"+corrtype[i]+"Zpeak_leadPt_"+Tag[i_tup]);
			h_Zpeak_subPt[i] = (TH1D*)h_pT[i]->Clone("h_"+corrtype[i]+"Zpeak_subPt_"+Tag[i_tup]);
			h_Zpeak_eta[i] = (TH1D*)h_eta[i]->Clone("h_"+corrtype[i]+"Zpeak_eta_"+Tag[i_tup]);
			h_Zpeak_leadEta[i] = (TH1D*)h_eta[i]->Clone("h_"+corrtype[i]+"Zpeak_leadEta_"+Tag[i_tup]);
			h_Zpeak_subEta[i] = (TH1D*)h_eta[i]->Clone("h_"+corrtype[i]+"Zpeak_subEta_"+Tag[i_tup]);
			h_Zpeak_phi[i] = (TH1D*)h_phi[i]->Clone("h_"+corrtype[i]+"Zpeak_phi_"+Tag[i_tup]);
			h_Zpeak_leadPhi[i] = (TH1D*)h_phi[i]->Clone("h_"+corrtype[i]+"Zpeak_leadPhi_"+Tag[i_tup]);
			h_Zpeak_subPhi[i] = (TH1D*)h_phi[i]->Clone("h_"+corrtype[i]+"Zpeak_subPhi_"+Tag[i_tup]);
			h_Zpeak_nVTX_woPU[i] = (TH1D*)h_nVTX_woPU[i]->Clone("h_"+corrtype[i]+"Zpeak_nVTX_woPU_"+Tag[i_tup]);
			h_Zpeak_nVTX_withPU[i] = (TH1D*)h_nVTX_woPU[i]->Clone("h_"+corrtype[i]+"Zpeak_nVTX_withPU_"+Tag[i_tup]);
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
			Double_t PrefiringWeight = 1;
			if( isMC == kTRUE ) PrefiringWeight = ntuple->_prefiringweight;

			// -- efficiency weights -- //
			Double_t effweight = 1, effweight_BtoF = 1, effweight_GtoH = 1;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);

			//Bool_t GenFlag_top = kTRUE;
			// -- Separate ttbar samples -- //
			Bool_t GenFlag_top = kFALSE;
			vector<GenOthers> GenTopCollection;
			GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);

			if( GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				SumWeight_Separated += GenWeight; // <- This part must exist!! 

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
				isPassEventSelection = analyzer->EventSelection(MuonCollection, ntuple, &SelectedMuonCollection);

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
						//effweight_BtoF = analyzer->EfficiencySF_EventWeight_HLT_BtoF( mu1, mu2 );
						//effweight_GtoH = analyzer->EfficiencySF_EventWeight_HLT_GtoH( mu1, mu2 );
						effweight_BtoF = analyzer->EfficiencySF_EventWeight_HLT_BtoF_new( mu1, mu2 );
						effweight_GtoH = analyzer->EfficiencySF_EventWeight_HLT_GtoH_new( mu1, mu2 );
						effweight = (Lumi_BtoF * effweight_BtoF + Lumi_GtoH * effweight_GtoH) / Lumi;
					}

					Int_t icorr = 0; // full corrected

					h_mass[icorr]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_mass_fine[icorr]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_diPt[icorr]->Fill( reco_Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_rapi[icorr]->Fill( reco_rapi, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_pT[icorr]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_pT[icorr]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_eta[icorr]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_eta[icorr]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_phi[icorr]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_phi[icorr]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_nVTX_woPU[icorr]->Fill( ntuple->nVertices, TotWeight * effweight * PrefiringWeight );
					h_nVTX_withPU[icorr]->Fill( ntuple->nVertices, TotWeight * PUWeight * effweight * PrefiringWeight );

					h_leadPt[icorr]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_leadEta[icorr]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_leadPhi[icorr]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_subPt[icorr]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_subEta[icorr]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
					h_subPhi[icorr]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight );

					if( 60 < reco_M && reco_M < 120 )
					{
						h_Zpeak_mass[icorr]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_mass_fine[icorr]->Fill( reco_M, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_diPt[icorr]->Fill( reco_Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_rapi[icorr]->Fill( reco_rapi, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_pT[icorr]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_pT[icorr]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_eta[icorr]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_eta[icorr]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_phi[icorr]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_phi[icorr]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_nVTX_woPU[icorr]->Fill( ntuple->nVertices, TotWeight * effweight * PrefiringWeight );
						h_Zpeak_nVTX_withPU[icorr]->Fill( ntuple->nVertices, TotWeight * PUWeight * effweight * PrefiringWeight );

						h_Zpeak_leadPt[icorr]->Fill( mu1.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_leadEta[icorr]->Fill( mu1.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_leadPhi[icorr]->Fill( mu1.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_subPt[icorr]->Fill( mu2.Pt, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_subEta[icorr]->Fill( mu2.eta, TotWeight * PUWeight * effweight * PrefiringWeight );
						h_Zpeak_subPhi[icorr]->Fill( mu2.phi, TotWeight * PUWeight * effweight * PrefiringWeight );
					} //End of Z-peak

				} //End of event selection


			} //End of if( isTriggered )

		} //End of event iteration

		for(Int_t i=0;i<ncorr;i++)
		{
			if( i != 0 ) continue;

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
			h_nVTX_woPU[i]->Write();
			h_nVTX_withPU[i]->Write();

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
			h_Zpeak_nVTX_woPU[i]->Write();
			h_Zpeak_nVTX_withPU[i]->Write();
		}


		printf("\tTotal sum of weights: %.1lf\n", SumWeight);
		printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
		if( isMC == kTRUE )
		{
			if( nEvents[i_tup] == SumWeight_Separated ) cout << "\t\t>> nEvents is correct!" << endl;
			else cout << "\t\t>> WARNING: nEvents is different" << endl;
			printf("\tNormalization factor: %.8f\n", lumi*Xsec[i_tup]/nEvents[i_tup]);
		}

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

