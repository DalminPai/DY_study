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
#include "../../HEADER/DYAnalyzer_TightID_PFIso_NewMuonSF.h"

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
void MuMu_TightID_PFIso_NewMuonSF(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24")
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
	// -- Alternative signal MC samples -- //
	else if( type == 19 ) Type = "ZToMuMu_powheg";
	else if( type == 18 ) Type = "DYMuMu_M50toInf_madgraph";

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
		//analyzer->SetupEfficiencyScaleFactor_BtoF();
		//analyzer->SetupEfficiencyScaleFactor_GtoH();
		analyzer->SetupEfficiencyScaleFactor_BtoF_new();
		analyzer->SetupEfficiencyScaleFactor_GtoH_new();
	}

	// -- Output ROOTFile -- //	
	//TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20180802_MuMu_TightID_PFIso_Zpeak_wo_trackingSF_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
	TString Output_ROOTFile = BaseDir+"/RESULT/MuMu/ROOTFile_20180807_MuMu_TightID_PFIso_Zpeak_NewMuonSF_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
								+TString::Itoa(isTopPtReweighting,10)+".root";
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
			TString version = "v2.1";
			//if( type == 10 ) version = "v2.3";
			if( type == 19 || type == 18 ) version = "v2.3";

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
		// Reco-level
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
		TH1D *h_phi = new TH1D("h_phi_"+Tag[i_tup], "", 100, -5, 5);
		TH1D *h_leadPhi = (TH1D*)h_phi->Clone("h_leadPhi_"+Tag[i_tup]);
		TH1D *h_subPhi = (TH1D*)h_phi->Clone("h_subPhi_"+Tag[i_tup]);

		TH1D *h_plus_pT = (TH1D*)h_pT->Clone("h_plus_pT_"+Tag[i_tup]);
		TH1D *h_minus_pT = (TH1D*)h_pT->Clone("h_minus_pT_"+Tag[i_tup]);
		TH1D *h_plus_eta = (TH1D*)h_eta->Clone("h_plus_eta_"+Tag[i_tup]);
		TH1D *h_minus_eta = (TH1D*)h_eta->Clone("h_minus_eta_"+Tag[i_tup]);

		// Dilepton pT cut : 30GeV
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
		TH1D *h_PtCut_phi = (TH1D*)h_phi->Clone("h_PtCut_phi_"+Tag[i_tup]);
		TH1D *h_PtCut_leadPhi = (TH1D*)h_phi->Clone("h_PtCut_leadPhi_"+Tag[i_tup]);
		TH1D *h_PtCut_subPhi = (TH1D*)h_phi->Clone("h_PtCut_subPhi_"+Tag[i_tup]);

		// 2D measurement
		TH1D *h_M20to30_rapi = (TH1D*)h_rapi->Clone("h_M20to30_rapi_"+Tag[i_tup]);
		TH1D *h_M30to45_rapi = (TH1D*)h_rapi->Clone("h_M30to45_rapi_"+Tag[i_tup]);
		TH1D *h_M45to60_rapi = (TH1D*)h_rapi->Clone("h_M45to60_rapi_"+Tag[i_tup]);
		TH1D *h_M60to120_rapi = (TH1D*)h_rapi->Clone("h_M60to120_rapi_"+Tag[i_tup]);
		TH1D *h_M120to200_rapi = (TH1D*)h_rapi->Clone("h_M120to200_rapi_"+Tag[i_tup]);
		TH1D *h_M200to1500_rapi = (TH1D*)h_rapi->Clone("h_M200to1500_rapi_"+Tag[i_tup]);

		// Gen-level
		TH1D *h_Gen_mass = (TH1D*)h_mass->Clone("h_Gen_mass_"+Tag[i_tup]);
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

		TH1D *h_Gen_plus_pT = (TH1D*)h_pT->Clone("h_Gen_plus_pT_"+Tag[i_tup]);
		TH1D *h_Gen_minus_pT = (TH1D*)h_pT->Clone("h_Gen_minus_pT_"+Tag[i_tup]);
		TH1D *h_Gen_plus_eta = (TH1D*)h_eta->Clone("h_Gen_plus_eta_"+Tag[i_tup]);
		TH1D *h_Gen_minus_eta = (TH1D*)h_eta->Clone("h_Gen_minus_eta_"+Tag[i_tup]);

		// LeadEtaStudy
		TH1D *h_pt_leadEta[6];
		TH1D *h_pt_subEta[6];
		TH1D *h_Gen_pt_leadEta[6];
		TH1D *h_Gen_pt_subEta[6];
		for(Int_t i=0; i<6; i++)
		{
			h_pt_leadEta[i] = (TH1D*)h_eta->Clone("h_pt_"+TString::Itoa(i+1,10)+"_leadEta_"+Tag[i_tup]);
			h_pt_subEta[i] = (TH1D*)h_eta->Clone("h_pt_"+TString::Itoa(i+1,10)+"_subEta_"+Tag[i_tup]);
			h_Gen_pt_leadEta[i] = (TH1D*)h_eta->Clone("h_Gen_pt_"+TString::Itoa(i+1,10)+"_leadEta_"+Tag[i_tup]);
			h_Gen_pt_subEta[i] = (TH1D*)h_eta->Clone("h_Gen_pt_"+TString::Itoa(i+1,10)+"_subEta_"+Tag[i_tup]);
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
			if( isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( ntuple->nPileUp );

			// -- efficiency weights -- //
			Double_t effweight = 1, effweight_BtoF = 1, effweight_GtoH = 1;

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
					if( genlep.isMuon() && genlep.fromHardProcessFinalState )
						GenLeptonCollection.push_back( genlep );
				}

				if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
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
					isPassAcc_GenLepton = analyzer->isPassAccCondition_GenLepton(genlep1, genlep2);
					if( isPassAcc_GenLepton == kTRUE )
					{
						Double_t gen_M = (genlep1.Momentum + genlep2.Momentum).M();
						Double_t gen_Pt = (genlep1.Momentum + genlep2.Momentum).Pt();
						Double_t gen_rapi = (genlep1.Momentum + genlep2.Momentum).Rapidity();

						if( 60 < gen_M && gen_M < 120 ) // Z-peak
						{
							h_Gen_mass->Fill( gen_M, TotWeight );
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

							if( genlep1.charge == 1 )
							{
								h_Gen_plus_pT->Fill( genlep1.Pt, TotWeight );
								h_Gen_plus_eta->Fill( genlep1.eta, TotWeight );
							}
							else if( genlep1.charge == -1 )
							{
								h_Gen_minus_pT->Fill( genlep1.Pt, TotWeight );
								h_Gen_minus_eta->Fill( genlep1.eta, TotWeight );
							}

							if( genlep2.charge == 1 )
							{
								h_Gen_plus_pT->Fill( genlep2.Pt, TotWeight );
								h_Gen_plus_eta->Fill( genlep2.eta, TotWeight );
							}
							else if( genlep2.charge == -1 )
							{
								h_Gen_minus_pT->Fill( genlep2.Pt, TotWeight );
								h_Gen_minus_eta->Fill( genlep2.eta, TotWeight );
							}

							// LeadEtaStudy
							for(Int_t i=0; i<6; i++)
							{
								Int_t pt_min = 30 + i*5;
								Int_t pt_max = 35 + i*5;
								if( pt_min < genlep1.Pt && genlep1.Pt < pt_max && pt_min < genlep2.Pt && genlep2.Pt < pt_max )
								{
										h_Gen_pt_subEta[i]->Fill( genlep2.eta, TotWeight );
										h_Gen_pt_leadEta[i]->Fill( genlep1.eta, TotWeight );
								}
							}

						} //Z-peak mass cut

					} //isPassAcc_GenLepton

				} //the events containing 2 muons from hard-process

			} //End of Gen-level selection

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

					// -- Rochester correction -- //
					Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
					Int_t s, m;
						
					if( isMC == kFALSE )
						SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
						//SF = rc.kScaleDT(mu.charge, mu.TuneP_pT, mu.TuneP_eta, mu.TuneP_phi, s=0, m=0);
					else
						SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);
						//SF = rc.kScaleAndSmearMC(mu.charge, mu.TuneP_pT, mu.TuneP_eta, mu.TuneP_phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);

					mu.Pt = SF*mu.Pt;
					//mu.TuneP_pT = SF*mu.TuneP_pT;

					mu.Momentum.SetPtEtaPhiM( mu.Pt, mu.eta, mu.phi, M_Mu );
					//analyzer->ConvertToTunePInfo( mu );

					MuonCollection.push_back( mu );
				}	

				// -- Event Selection -- //
				vector< Muon > SelectedMuonCollection;
				Bool_t isPassEventSelection = kFALSE;
				isPassEventSelection = analyzer->EventSelection_Zpeak(MuonCollection, ntuple, &SelectedMuonCollection);

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

					h_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
					h_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );
					h_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
					h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
					h_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight );
					h_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight );
					h_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight );
					h_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight );
					h_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight );
					h_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight );

					h_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight );
					h_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight );
					h_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight );
					h_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight );
					h_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight );
					h_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight );

					if( mu1.charge == 1 )
					{
						h_plus_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight );
						h_plus_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight );
					}
					else if( mu1.charge == -1 )
					{
						h_minus_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight );
						h_minus_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight );
					}

					if( mu2.charge == 1 )
					{
						h_plus_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight );
						h_plus_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight );
					}
					else if( mu2.charge == -1 )
					{
						h_minus_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight );
						h_minus_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight );
					}

					// Dilepton pT cut : 30GeV
					if( 30 < reco_Pt )
					{
						h_PtCut_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_PtCut_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_PtCut_diPt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
						h_PtCut_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
						h_PtCut_pT->Fill( mu1.Pt, TotWeight * PUWeight * effweight );
						h_PtCut_pT->Fill( mu2.Pt, TotWeight * PUWeight * effweight );
						h_PtCut_eta->Fill( mu1.eta, TotWeight * PUWeight * effweight );
						h_PtCut_eta->Fill( mu2.eta, TotWeight * PUWeight * effweight );
						h_PtCut_phi->Fill( mu1.phi, TotWeight * PUWeight * effweight );
						h_PtCut_phi->Fill( mu2.phi, TotWeight * PUWeight * effweight );

						h_PtCut_leadPt->Fill( mu1.Pt, TotWeight * PUWeight * effweight );
						h_PtCut_leadEta->Fill( mu1.eta, TotWeight * PUWeight * effweight );
						h_PtCut_leadPhi->Fill( mu1.phi, TotWeight * PUWeight * effweight );
						h_PtCut_subPt->Fill( mu2.Pt, TotWeight * PUWeight * effweight );
						h_PtCut_subEta->Fill( mu2.eta, TotWeight * PUWeight * effweight );
						h_PtCut_subPhi->Fill( mu2.phi, TotWeight * PUWeight * effweight );
					}

					// 2D measurement
					if( 20 < reco_M && reco_M < 30 ) h_M20to30_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
					else if( 30 < reco_M && reco_M < 45 ) h_M30to45_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
					else if( 45 < reco_M && reco_M < 60 ) h_M45to60_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
					else if( 60 < reco_M && reco_M < 120 ) h_M60to120_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
					else if( 120 < reco_M && reco_M < 200 ) h_M120to200_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );
					else if( 200 < reco_M && reco_M < 1500 ) h_M200to1500_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );

					// LeadEtaStudy
					for(Int_t i=0; i<6; i++)
					{
						Int_t pt_min = 30 + i*5;
						Int_t pt_max = 35 + i*5;
						if( pt_min < mu1.Pt && mu1.Pt < pt_max && pt_min < mu2.Pt && mu2.Pt < pt_max )
						{
							h_pt_leadEta[i]->Fill( mu1.eta, TotWeight * PUWeight * effweight );
							h_pt_subEta[i]->Fill( mu2.eta, TotWeight* PUWeight * effweight );
						}
					}
				} //End of event selection

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
		h_phi->Write();
		h_leadPhi->Write();
		h_subPhi->Write();

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
		h_PtCut_phi->Write();
		h_PtCut_leadPhi->Write();
		h_PtCut_subPhi->Write();

		h_M20to30_rapi->Write();
		h_M30to45_rapi->Write();
		h_M45to60_rapi->Write();
		h_M60to120_rapi->Write();
		h_M120to200_rapi->Write();
		h_M200to1500_rapi->Write();

		h_plus_pT->Write();
		h_minus_pT->Write();
		h_plus_eta->Write();
		h_minus_eta->Write();

		for(Int_t i=0; i<6; i++)
		{
			h_pt_leadEta[i]->Write();
			h_pt_subEta[i]->Write();
		}

		if( isMC == kTRUE && type < 20 )
		{
			h_Gen_mass->Write();
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

			h_Gen_plus_pT->Write();
			h_Gen_minus_pT->Write();
			h_Gen_plus_eta->Write();
			h_Gen_minus_eta->Write();

			for(Int_t i=0; i<6; i++)
			{
				h_Gen_pt_leadEta[i]->Write();
				h_Gen_pt_subEta[i]->Write();
			}
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

