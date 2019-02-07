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
// -- Started at 27 Jul 2018 -- //
void FR_applyFR(Int_t debug, Int_t type, Int_t remainder = 9999, Int_t isTopPtReweighting = 0, TString HLTname = "IsoMu24_OR_IsoTkMu24")
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
	// -- MC samples -- //
	else if( type == 11 ) Type = "DY_M10to50";
	else if( type == 12 ) Type = "DY_M50toInf";
	else if( type == 21 ) Type = "ttbar";
	else if( type == 22 ) Type = "ttbarBackup";
	else if( type == 41 ) Type = "VVnST";
	else if( type == 51 ) Type = "WJetsToLNu";
	//else if( type == 61 ) Type = "QCDMuEnriched_Pt15to170";
	//else if( type == 62 ) Type = "QCDMuEnriched_Pt170to600";
	//else if( type == 63 ) Type = "QCDMuEnriched_Pt600toInf";

	Bool_t isMC = kTRUE;
	if( type < 10  )
	{
		isMC = kFALSE;
		Type = "Data (" + DataLocation + ")";
		if( type == 7 ) DataLocation += "ver2";
	}

	Double_t Div = 1;
	//if( type == 1 || type == 6 || type == 7 || type == 11 || type == 31 || type == 51 ) Div = 2;
	if( type == 1 || type == 6 || type == 7 || type == 11 || type == 31 || type == 51 || type > 60 ) Div = 2;

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
	/*if( isMC == kTRUE )
	{
		analyzer->SetupEfficiencyScaleFactor_BtoF();
		analyzer->SetupEfficiencyScaleFactor_GtoH();
	}*/

	// -- Output ROOTFile -- //	
	TString Output_ROOTFile = BaseDir+"/RESULT/FR/fake_test_"+TString::Itoa(type,10)+"_"+TString::Itoa(remainder,10)+"_"
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
			if( Type.Contains("QCD") ) version = "v2.3";

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
		TH1D* histDijet1 = new TH1D("histDijet1_"+Tag[i_tup],"",43,massbins);
		TH1D* histDijet2 = new TH1D("histDijet2_"+Tag[i_tup],"",43,massbins);
		TH1D* histSameDijet1 = new TH1D("histSameDijet1_"+Tag[i_tup],"",43,massbins);
		TH1D* histSameDijet2 = new TH1D("histSameDijet2_"+Tag[i_tup],"",43,massbins);

		TH1D* fitDijet1 = new TH1D("fitDijet1_"+Tag[i_tup],"",37,15,200);
		TH1D* fitDijet2 = new TH1D("fitDijet2_"+Tag[i_tup],"",37,15,200);
		TH1D* fitSameDijet1 = new TH1D("fitSameDijet1_"+Tag[i_tup],"",37,15,200);
		TH1D* fitSameDijet2 = new TH1D("fitSameDijet2_"+Tag[i_tup],"",37,15,200);

		TH1D* rapDijet1 = new TH1D("rapDijet1_"+Tag[i_tup],"",48,-2.4,2.4);
		TH1D* rapDijet2 = new TH1D("rapDijet2_"+Tag[i_tup],"",48,-2.4,2.4);
		TH1D* rapSameDijet1 = new TH1D("rapSameDijet1_"+Tag[i_tup],"",48,-2.4,2.4);
		TH1D* rapSameDijet2 = new TH1D("rapSameDijet2_"+Tag[i_tup],"",48,-2.4,2.4);

		histDijet1->Sumw2();
		histDijet2->Sumw2();
		histSameDijet1->Sumw2();
		histSameDijet2->Sumw2();

		fitDijet1->Sumw2();
		fitDijet2->Sumw2();
		fitSameDijet1->Sumw2();
		fitSameDijet2->Sumw2();

		rapDijet1->Sumw2();
		rapDijet2->Sumw2();
		rapSameDijet1->Sumw2();
		rapSameDijet2->Sumw2();

		TH1D* histWJets1 = new TH1D("histWJets1_"+Tag[i_tup],"",43,massbins);
		TH1D* histWJets2 = new TH1D("histWJets2_"+Tag[i_tup],"",43,massbins);
		TH1D* histSameWJets1 = new TH1D("histSameWJets1_"+Tag[i_tup],"",43,massbins);
		TH1D* histSameWJets2 = new TH1D("histSameWJets2_"+Tag[i_tup],"",43,massbins);

		TH1D* fitWJets1 = new TH1D("fitWJets1_"+Tag[i_tup],"",37,15,200);
		TH1D* fitWJets2 = new TH1D("fitWJets2_"+Tag[i_tup],"",37,15,200);
		TH1D* fitSameWJets1 = new TH1D("fitSameWJets1_"+Tag[i_tup],"",37,15,200);
		TH1D* fitSameWJets2 = new TH1D("fitSameWJets2_"+Tag[i_tup],"",37,15,200);

		TH1D* rapWJets1 = new TH1D("rapWJets1_"+Tag[i_tup],"",48,-2.4,2.4);
		TH1D* rapWJets2 = new TH1D("rapWJets2_"+Tag[i_tup],"",48,-2.4,2.4);
		TH1D* rapSameWJets1 = new TH1D("rapSameWJets1_"+Tag[i_tup],"",48,-2.4,2.4);
		TH1D* rapSameWJets2 = new TH1D("rapSameWJets2_"+Tag[i_tup],"",48,-2.4,2.4);

		TH1D* histPass = new TH1D("histPass_"+Tag[i_tup],"",10,0,10);
		TH1D* histFail = new TH1D("histFail_"+Tag[i_tup],"",10,0,10);

		histWJets1->Sumw2();
		histWJets2->Sumw2();
		histSameWJets1->Sumw2();
		histSameWJets2->Sumw2();

		fitWJets1->Sumw2();
		fitWJets2->Sumw2();
		fitSameWJets1->Sumw2();
		fitSameWJets2->Sumw2();

		rapWJets1->Sumw2();
		rapWJets2->Sumw2();
		rapSameWJets1->Sumw2();
		rapSameWJets2->Sumw2();

		Double_t SumWeight = 0, SumWeight_Separated = 0;
		Double_t nPass = 0, nFail = 0; // for fake-rate

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
			//Double_t effweight = 1, effweight_BtoF = 1, effweight_GtoH = 1;

			// -- Separate DYLL samples -- //
			Bool_t GenFlag = kTRUE;
			/*Bool_t GenFlag = kFALSE;
			GenFlag = analyzer->SeparateDYLLSample_isHardProcess(Tag[i_tup], ntuple);*/

			// -- Separate ttbar samples -- //
			Bool_t GenFlag_top = kTRUE;
			/*Bool_t GenFlag_top = kFALSE;
			vector<GenOthers> GenTopCollection;
			GenFlag_top = analyzer->Separate_ttbarSample(Tag[i_tup], ntuple, &GenTopCollection);*/

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
			//if( isMC == kTRUE ) TotWeight = (lumi*Xsec[i_tup]/nEvents[i_tup])*GenWeight;

			////////////////////////////////
			// -- Reco level selection -- //
			////////////////////////////////
			Bool_t TriggerFlag = kTRUE;
			//Bool_t TriggerFlag = kFALSE;
			//TriggerFlag = ntuple->isTriggered( analyzer->HLT );

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

				// -- Check # of passing and failing muons -- //
				Bool_t flag_lead = kFALSE;
				vector< Muon > passMuons; vector< Muon > failMuons;
				for(Int_t j=0; j<(int)MuonCollection.size(); j++)
				{
					Muon mu_ = MuonCollection[j];

					if( mu_.isHighPtMuon() && mu_.Pt > 17 && mu_.eta < 2.4 )
					{
						if( mu_.Pt > 28 ) flag_lead = kTRUE;

						if( mu_.trkiso < 0.1 ) passMuons.push_back(MuonCollection[j]);
						else failMuons.push_back(MuonCollection[j]);
					}
				}

				if( flag_lead == kFALSE ) continue;

				histPass->Fill(passMuons.size(), TotWeight * PUWeight);
				histFail->Fill(failMuons.size(), TotWeight * PUWeight);

				nPass += TotWeight * PUWeight * passMuons.size();
				nFail += TotWeight * PUWeight * failMuons.size();

				// -- Event Selection -- //
				//Dijet selection
				vector< Muon > SelectedMuonCollection_Dijet;
				Bool_t isPassEventSelection_Dijet = kFALSE;
				isPassEventSelection_Dijet = analyzer->EventSelection_Dijet(MuonCollection, ntuple, &SelectedMuonCollection_Dijet);
				
				if( isPassEventSelection_Dijet == kTRUE )
				{
					Muon mu1, mu2;
					mu1 = SelectedMuonCollection_Dijet[0];
					mu2 = SelectedMuonCollection_Dijet[1];
					
					Double_t mass = (mu1.Momentum + mu2.Momentum).M();
					Double_t rap = (mu1.Momentum + mu2.Momentum).Rapidity();

					Double_t FR1_template = analyzer->FR_template(mu1);
					Double_t FR2_template = analyzer->FR_template(mu2);
					Double_t FR1_ratio = analyzer->FR_ratio(mu1);
					Double_t FR2_ratio = analyzer->FR_ratio(mu2);

					Double_t weight_template = TotWeight * PUWeight * FR1_template * FR2_template / ((1-FR1_template) * (1-FR2_template));
					Double_t weight_ratio = TotWeight * PUWeight * FR1_ratio * FR2_ratio / ((1-FR1_ratio) * (1-FR2_ratio));

					if( mu1.charge != mu2.charge )
					{
						if( mass > 15 && mass < 3000)
						{
							histDijet1->Fill(mass, weight_template);
							histDijet2->Fill(mass, weight_ratio);
							fitDijet1->Fill(mass, weight_template);
							fitDijet2->Fill(mass, weight_ratio);
							rapDijet1->Fill(rap, weight_template);
							rapDijet2->Fill(rap, weight_ratio);
						}
					}
					else
					{
						if( mass > 15 && mass < 3000)
						{
							histSameDijet1->Fill(mass, weight_template);
							histSameDijet2->Fill(mass, weight_ratio);
							fitSameDijet1->Fill(mass, weight_template);
							fitSameDijet2->Fill(mass, weight_ratio);
							rapSameDijet1->Fill(rap, weight_template);
							rapSameDijet2->Fill(rap, weight_ratio);
						}
					}
				} //End of Dijet selection

				//Wjet selection
				vector< Muon > SelectedMuonCollection_Wjet;
				Bool_t isPassEventSelection_Wjet = kFALSE;
				isPassEventSelection_Wjet = analyzer->EventSelection_Wjet(MuonCollection, ntuple, &SelectedMuonCollection_Wjet);

				if( isPassEventSelection_Wjet == kFALSE )
					cout << "fail Wjet" << endl;

				if( isPassEventSelection_Wjet == kTRUE )
				{
					cout << "pass Wjet" << endl;

					Muon mu1, mu2;
					mu1 = SelectedMuonCollection_Wjet[0]; //passing muon
					mu2 = SelectedMuonCollection_Wjet[1]; //failing muon

					Double_t mass = (mu1.Momentum + mu2.Momentum).M();
					Double_t rap = (mu1.Momentum + mu2.Momentum).Rapidity();
				
					Double_t FR2_template = analyzer->FR_template(mu2);
					Double_t FR2_ratio = analyzer->FR_ratio(mu2);

					Double_t weight_template = TotWeight * PUWeight * FR2_template / (1-FR2_template);
					Double_t weight_ratio = TotWeight * PUWeight * FR2_ratio / (1-FR2_ratio);

					if( mu1.charge != mu2.charge )
					{
							if( mass > 15 && mass < 3000)
							{
								histWJets1->Fill(mass, weight_template);
								histWJets2->Fill(mass, weight_ratio);
								fitWJets1->Fill(mass, weight_template);
								fitWJets2->Fill(mass, weight_ratio);
								rapWJets1->Fill(rap, weight_template);
								rapWJets2->Fill(rap, weight_ratio);
							}
					}
					else
					{
							if( mass > 15 && mass < 3000)
							{
								histSameWJets1->Fill(mass, weight_template);
								histSameWJets2->Fill(mass, weight_ratio);
								fitSameWJets1->Fill(mass, weight_template);
								fitSameWJets2->Fill(mass, weight_ratio);
								rapSameWJets1->Fill(rap, weight_template);
								rapSameWJets2->Fill(rap, weight_ratio);
							}
					}
				} //End of Wjet selection

			} //End of if( isTriggered )

		} //End of event iteration

		histDijet1->Write();
		histDijet2->Write();
		histWJets1->Write();
		histWJets2->Write();
		fitDijet1->Write();
		fitDijet2->Write();
		fitWJets1->Write();
		fitWJets2->Write();
		rapDijet1->Write();
		rapDijet2->Write();
		rapWJets1->Write();
		rapWJets2->Write();

		histSameDijet1->Write();
		histSameDijet2->Write();
		histSameWJets1->Write();
		histSameWJets2->Write();
		fitSameDijet1->Write();
		fitSameDijet2->Write();
		fitSameWJets1->Write();
		fitSameWJets2->Write();
		rapSameDijet1->Write();
		rapSameDijet2->Write();
		rapSameWJets1->Write();
		rapSameWJets2->Write();

		histPass->Write();
		histFail->Write();

		cout<<"# of passing muons = "<<nPass<<endl;
		cout<<"# of failing muons = "<<nFail<<endl;
		cout<<endl;
		cout<<"# of passing muons per event = "<<nPass/SumWeight_Separated<<endl;
		cout<<"# of failing muons per event = "<<nFail/SumWeight_Separated<<endl;
		cout<<endl;

		printf("\tTotal sum of weights: %.1lf\n", SumWeight);
		printf("\tSum of weights of Seperated events: %.1lf\n", SumWeight_Separated);
		//if( isMC == kTRUE ) printf("\tNormalization factor: %.8f\n", lumi*Xsec[i_tup]/nEvents[i_tup]);

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
	int   c	 = ratio * w;
 
	// Show the percentage complete.
	printf("%3d%% [", (int)(ratio*100) );
 
	// Show the load bar.
	for (int x=0; x<c; x++) cout << "=";
 
	for (int x=c; x<w; x++) cout << " ";
 
	// ANSI Control codes to go back to the
	// previous line and clear it.
	cout << "]\r" << flush;

}

