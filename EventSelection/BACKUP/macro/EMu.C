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

// -- Electron-Muon Event Selection -- //
void EMu(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24")
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
	//else if( type == 12 ) Type = "DYMuMu_M50toInf";
	else if( type == 12 ) Type = "DYMuMu_M50to200";
	else if( type == 13 ) Type = "DYMuMu_M200toInf";
	// -- Background MC samples -- //
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	else if( type == 23 ) Type = "ttbar_M700toInf";
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
		analyzer->SetupEfficiencyScaleFactor_electron();
	}

	// -- Output ROOTFile -- //	
	TString Output_ROOTFile = BaseDir+"/RESULT/EMu/ROOTFile_20180625_EMu_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
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
		ntuple->TurnOnBranches_Electron();
		if( isMC == kTRUE )
		{
			ntuple->TurnOnBranches_GenLepton(); // for all leptons
			ntuple->TurnOnBranches_GenOthers(); // for quarks
		}

		////////////////////////////
		// -- Making Histogram -- //
		////////////////////////////
		TH1D *h_emu_mass = new TH1D("h_emu_mass_"+Tag[i_tup], "", 43, massbins);
		TH1D *h_emu_mass_fine = new TH1D("h_emu_mass_fine_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_mu_pT = new TH1D("h_mu_pT_"+Tag[i_tup], "", 10000, 0, 10000);
		TH1D *h_mu_eta = new TH1D("h_mu_eta_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_mu_phi = new TH1D("h_mu_phi_"+Tag[i_tup], "", 1000, -5, 5);
		TH1D *h_ele_pT = new (TH1D*)h_mu_pT->Clone("h_ele_pT_"+Tag[i_tup]);
		TH1D *h_ele_eta = new (TH1D*)h_mu_eta->Clone("h_ele_eta_"+Tag[i_tup]);
		TH1D *h_ele_etaSC = new (TH1D*)h_mu_eta->Clone("h_ele_etaSC_"+Tag[i_tup]);
		TH1D *h_ele_phi = new (TH1D*)h_mu_phi->Clone("h_ele_phi_"+Tag[i_tup]);

		TH1D *h_emuSS_mass = new (TH1D*)h_emu_mass->Clone("h_emuSS_mass_"+Tag[i_tup]);
		TH1D *h_emuSS_mass_fine = new (TH1D*)h_emu_mass_fine->Clone("h_emuSS_mass_fine_"+Tag[i_tup]);
		TH1D *h_muSS_pT = new (TH1D*)h_mu_pT->Clone("h_muSS_pT_"+Tag[i_tup]);
		TH1D *h_muSS_eta = new (TH1D*)h_mu_eta->Clone("h_muSS_eta_"+Tag[i_tup]);
		TH1D *h_muSS_phi = new (TH1D*)h_mu_phi->Clone("h_muSS_phi_"+Tag[i_tup]);
		TH1D *h_eleSS_pT = new (TH1D*)h_mu_pT->Clone("h_eleSS_pT_"+Tag[i_tup]);
		TH1D *h_eleSS_eta = new (TH1D*)h_mu_eta->Clone("h_eleSS_eta_"+Tag[i_tup]);
		TH1D *h_eleSS_etaSC = new (TH1D*)h_mu_eta->Clone("h_eleSS_etaSC_"+Tag[i_tup]);
		TH1D *h_eleSS_phi = new (TH1D*)h_mu_phi->Clone("h_eleSS_phi_"+Tag[i_tup]);

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

			////////////////////////////////
			// -- Reco level selection -- //
			////////////////////////////////
			Bool_t TriggerFlag = kFALSE;
			TriggerFlag = ntuple->isTriggered( analyzer->HLT );

			if( TriggerFlag == kTRUE && GenFlag == kTRUE && GenFlag_top == kTRUE )
			{
				// -- Muon -- //
				vector< Muon > MuonCollection;
				Int_t N_Muons = ntuple->nMuon;
				for(Int_t i_reco=0; i_reco<N_Muons; i_reco++)
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

				// -- Electron -- //
				vector< Electron > ElectronCollection;
				Int_t N_Electrons = ntuple->Nelectrons;
				for(Int_t i_reco=0; i_reco<N_Electrons; i_reco++)
				{
					Electron ele;
					ele.FillFromNtuple(ntuple, i_reco);

					ElectronCollection.push_back( ele );
				}	

				// -- Event Selection -- //
				vector< Muon > SelectedMuonCollection;
				vector< Electron > SelectedElectronCollection;
				Bool_t isPassEventSelection = kFALSE;
				isPassEventSelection = analyzer->EventSelection_emu_method(MuonCollection, ElectronCollection, ntuple,
																			&SelectedMuonCollection, &SelectedElectronCollection);

				if( isPassEventSelection == kTRUE )
				{
					Muon mu = SelectedMuonCollection[0];
					Electron ele = SelectedElectronCollection[0];

					Double_t reco_M = (mu.Momentum + ele.Momentum).M();

					// -- Apply efficiency correcion -- //
					if( isMC == kTRUE )
					{
						effweight_BtoF = analyzer->EfficiencySF_EventWeight_emu_BtoF( mu, ele );
						effweight_GtoH = analyzer->EfficiencySF_EventWeight_emu_GtoH( mu, ele );
						effweight = (Lumi_BtoF * effweight_BtoF + Lumi_GtoH * effweight_GtoH) / Lumi;
					}

					if( mu.charge != ele.charge ) // opposite charge
					{
						h_emu_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_emu_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_mu_pT->Fill( mu.Pt, TotWeight * PUWeight * effweight );
						h_mu_eta->Fill( mu.eta, TotWeight * PUWeight * effweight );
						h_mu_phi->Fill( mu.phi, TotWeight * PUWeight * effweight );
						h_ele_pT->Fill( ele.Pt, TotWeight * PUWeight * effweight );
						h_ele_eta->Fill( ele.eta, TotWeight * PUWeight * effweight );
						h_ele_etaSC->Fill( ele.etaSC, TotWeight * PUWeight * effweight );
						h_ele_phi->Fill( ele.phi, TotWeight * PUWeight * effweight );
					}
					else // same charge
					{
						h_emuSS_mass->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_emuSS_mass_fine->Fill( reco_M, TotWeight * PUWeight * effweight );
						h_muSS_pT->Fill( mu.Pt, TotWeight * PUWeight * effweight );
						h_muSS_eta->Fill( mu.eta, TotWeight * PUWeight * effweight );
						h_muSS_phi->Fill( mu.phi, TotWeight * PUWeight * effweight );
						h_eleSS_pT->Fill( ele.Pt, TotWeight * PUWeight * effweight );
						h_eleSS_eta->Fill( ele.eta, TotWeight * PUWeight * effweight );
						h_eleSS_etaSC->Fill( ele.etaSC, TotWeight * PUWeight * effweight );
						h_eleSS_phi->Fill( ele.phi, TotWeight * PUWeight * effweight );
					}
				} // End of event selection

			} //End of if( isTriggered )

		} //End of event iteration

		h_emu_mass->Write();
		h_emu_mass_fine->Write();
		h_mu_pT->Write();
		h_mu_eta->Write();
		h_mu_phi->Write();
		h_ele_pT->Write();
		h_ele_eta->Write();
		h_ele_etaSC->Write();
		h_ele_phi->Write();

		h_emuSS_mass->Write();
		h_emuSS_mass_fine->Write();
		h_muSS_pT->Write();
		h_muSS_eta->Write();
		h_muSS_phi->Write();
		h_eleSS_pT->Write();
		h_eleSS_eta->Write();
		h_eleSS_etaSC->Write();
		h_eleSS_phi->Write();

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

