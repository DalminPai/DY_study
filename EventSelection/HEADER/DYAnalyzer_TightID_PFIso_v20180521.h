// -- Class for common functions used in DY differential cross section measurement analysis @ 13 TeV -- //
// -- Author: KyoengPil Lee, 05 Dec. 2015 -- //
// -- Author: Dalmin Pai, 11 Sep. 2017 -- //
// -- Synchronize acceptance cuts in dilepton channels : 19 Jan. 2018 -- //
// -- Add N-1 selection of trkiso : 26 Jan. 2018 -- //
// -- Change isolation (trkiso->RelPFIso_dBeta) in muon channel : 30 Jan. 2018 -- //
// -- Change muon ID (HighPt->Tight) : 07 Feb. 2018 -- //
// -- Change Xsec of DYMuMu_M50toInf from NLO to NNLO : 20 Feb. 2018 -- //
// -- Add "Mu50_OR_TkMu50" trigger : 20 Feb. 2018 -- //
// -- Add trigger matching in "EventSelection_Zdiff_13TeV" : 22 Feb. 2018 -- //
// -- Add test functions to check efficiency SFs separately : 23 Feb. 2018 -- //
// -- Add Tracking SF : 27 Feb. 2018 -- //
// -- Modify "EventSelection" : 02 Mar. 2018 -- //
// -- Add "EventSelection_ElectronChannel0" for without energy scale correction : 08 Mar. 2018 -- //
// -- Add "EfficiencySF_EventWeight_electron_Reco" : 08 Mar. 2018 -- //
// -- Modify wrong NNLO Xsec value of DYLL_M50toInf ( 6025.2/3 -> 1921.8 ) : 09 Mar. 2018 -- //
// -- Add trigger part into "EfficiencySF_EventWeight_electron" : 14 Mar. 2018 -- //
// -- Add "EfficiencySF_EventWeight_electron_RecoID" : 14 Mar. 2018 -- //
// -- Add "massbins" which is used in DY mass distribution : 20 Mar. 2018 -- //
// -- Add "DYMuMu_M200toInf" into SetupMCsamples_Moriond17 : 11 Apr. 2018 -- //
// -- Modify EventSelections of dimuon, dielectron, and emu channels to use full mass range, correct Xsec values of DY samples with KP's values, and correct # of events of WJets : 13 Apr. 2018 -- //
// -- Update "EventSelection_emu_method" : 13 Apr. 2018 -- //
// -- Introduce "Leading muon eta SF" : 17 Apr. 2018 -- //
// -- Modify "EventSelection" to use only Z peak : 17 Apr. 2018 -- //
// -- Update "Leading muon eta SF" as "LeadEtaCorr" : 21 May. 2018 -- //
// -- Insert nEvents of DYEE_M50to200 (27206082.0) : 19 Jun. 2018 -- //
// -- Introduce "EventSelection_Zpeak" and modify "EventSelection" to use full mass : 19 Jun. 2018 -- //
// -- Add "ZToMuMu(EE)_powheg" : 19 Jun. 2018 -- //
// -- Introduce "EventSelection_ElectronChannel_Zpeak" : 19 Jun. 2018 -- //
#pragma once

#include "Object_v01Dec17.h"
#include "NtupleHandle_v31Oct17.h"
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#define Lumi 35867 // -- from Run2016B to Run2016H, JSON. unit: /pb, Updated at 2017.07.30 -- //
#define Lumi_BtoF 19721 // -- from Run2016B to Run2016F, JSON. unit: /pb, Updated at 2018.05.17 -- //
#define Lumi_GtoH 16146 // -- from Run2016G to Run2016H, JSON. unit: /pb, Updated at 2018.05.17 -- //
#define nMassBin 43

const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

class DYAnalyzer
{
public:

	TString HLT;
	Double_t LeadPtCut;
	Double_t SubPtCut;
	Double_t LeadEtaCut;
	Double_t SubEtaCut;

//	Double_t PileUpWeight[52];
	Double_t PileUpWeight[75];

	/////////////////////////
	// -- Efficiency SF -- //
	/////////////////////////
	// -- For efficiency SF of BtoF -- //
	Double_t Eff_Reco_data_BtoF[15][1]; //Tracking SF
	Double_t Eff_Reco_MC_BtoF[15][1];

//	Double_t Eff_ID_data_BtoF[4][7]; //HighPt
//	Double_t Eff_ID_MC_BtoF[4][7];
	Double_t Eff_ID_data_BtoF[4][6]; //Tight
	Double_t Eff_ID_MC_BtoF[4][6];

//	Double_t Eff_Iso_data_BtoF[4][7]; //trkiso
//	Double_t Eff_Iso_MC_BtoF[4][7];
	Double_t Eff_Iso_data_BtoF[4][6]; //pfiso
	Double_t Eff_Iso_MC_BtoF[4][6];

	Double_t Eff_HLT_data_BtoF[4][7];
	Double_t Eff_HLT_MC_BtoF[4][7];

	// -- For efficiency SF of GtoH -- //
	Double_t Eff_Reco_data_GtoH[15][1]; //Tracking SF
	Double_t Eff_Reco_MC_GtoH[15][1];

//	Double_t Eff_ID_data_GtoH[4][7]; //HighPt
//	Double_t Eff_ID_MC_GtoH[4][7];
	Double_t Eff_ID_data_GtoH[4][6]; //Tight
	Double_t Eff_ID_MC_GtoH[4][6];

//	Double_t Eff_Iso_data_GtoH[4][7]; //trkiso
//	Double_t Eff_Iso_MC_GtoH[4][7];
	Double_t Eff_Iso_data_GtoH[4][6]; //pfiso
	Double_t Eff_Iso_MC_GtoH[4][6];

	Double_t Eff_HLT_data_GtoH[4][7];
	Double_t Eff_HLT_MC_GtoH[4][7];

	// -- For efficiency SF of electron -- //
	Double_t Eff_Reco_data[30][1];
	Double_t Eff_Reco_MC[30][1];

	Double_t Eff_ID_data[10][5];
	Double_t Eff_ID_MC[10][5];

	Double_t Eff_HLT_Leg1_data[10][8];
	Double_t Eff_HLT_Leg1_MC[10][8];
	Double_t Eff_HLT_Leg2_data[10][8];
	Double_t Eff_HLT_Leg2_MC[10][8];

	// -- outdated -- //
	Double_t Eff_RecoID_data[5][4];
	Double_t Eff_RecoID_MC[5][4];

	Double_t Eff_Iso_data[5][4];
	Double_t Eff_Iso_MC[5][4];

	Double_t Eff_HLTv4p2_data[5][4];
	Double_t Eff_HLTv4p2_MC[5][4];

	Double_t Eff_HLTv4p3_data[5][4];
	Double_t Eff_HLTv4p3_MC[5][4];

	// -- Constructor -- //
	DYAnalyzer(TString HLTname);

	// -- Setup accetpance cuts -- //
	void AssignAccThreshold(TString HLTname, TString *HLT, Double_t *LeadPtCut, Double_t *SubPtCut, Double_t *LeadEtaCut, Double_t *SubEtaCut);

	////////////////////////////
	// -- Setup MC samples -- //
	////////////////////////////
	void SetupMCsamples_v20160412_76X_MINIAODv2_CheckPremix( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_Moriond17( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_v20160309_76X_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_v20160131_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	void SetupMCsamples_v20160117_MiniAOD_JetMET( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents );
	Bool_t SeparateDYLLSample_isHardProcess(TString Tag, NtupleHandle *ntuple);
	Bool_t Separate_ttbarSample(TString Tag, NtupleHandle *ntuple, vector<GenOthers> *GenTopCollection);

	// -- outdated -- //
	Bool_t SeparateDYLLSample(TString Tag, NtupleHandle *ntuple);

	//////////////////////////////////
	// -- Setup pileup weighting -- //
	//////////////////////////////////
	void SetupPileUpReWeighting( Bool_t isMC );
	Double_t PileUpWeightValue(Int_t PileUp_MC);

	// -- for 76X -- //
	void SetupPileUpReWeighting_76X( Bool_t isMC );
	Double_t PileUpWeightValue_76X(Int_t PileUp_MC);

	// -- for 80X -- //
	void SetupPileUpReWeighting_80X( Bool_t isMC, TString ROOTFileName );
	Double_t PileUpWeightValue_80X(Int_t PileUp_MC);

	/////////////////////////////////////////
	// -- Setup Efficiency scale factor -- //
	/////////////////////////////////////////
	// -- for muon -- //
	void SetupEfficiencyScaleFactor_BtoF();
	void SetupEfficiencyScaleFactor_GtoH();
	Double_t EfficiencySF_EventWeight_HLT_BtoF(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_GtoH(Muon mu1, Muon mu2);
	Int_t Find_muon_PtBin_Reco(Double_t Pt);
	Int_t Find_muon_PtBin_ID(Double_t Pt);
	Int_t Find_muon_PtBin_Iso(Double_t Pt);
	Int_t Find_muon_PtBin_Trig(Double_t Pt);
	Int_t Find_muon_EtaBin_Reco(Double_t eta);
	Int_t Find_muon_EtaBin_ID(Double_t eta);
	Int_t Find_muon_EtaBin_Iso(Double_t eta);
	Int_t Find_muon_EtaBin_Trig(Double_t eta);
	// Test SFs separately
	Double_t EfficiencySF_EventWeight_HLT_BtoF_Reco(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_BtoF_RecoID(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_BtoF_RecoIDIso(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_GtoH_Reco(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_GtoH_RecoID(Muon mu1, Muon mu2);
	Double_t EfficiencySF_EventWeight_HLT_GtoH_RecoIDIso(Muon mu1, Muon mu2);

	// -- for electron -- //
	void SetupEfficiencyScaleFactor_electron();
	Double_t EfficiencySF_EventWeight_electron(Electron ele1, Electron ele2);
	Int_t Find_electron_PtBin_Reco(Double_t Pt);
	Int_t Find_electron_PtBin_ID(Double_t Pt);
	Int_t Find_electron_PtBin_Trig(Double_t Pt);
	Int_t Find_electron_EtaBin_Reco(Double_t eta);
	Int_t Find_electron_EtaBin_ID(Double_t eta);
	Int_t Find_electron_EtaBin_Trig(Double_t eta);
	// Test SFs separately
	Double_t EfficiencySF_EventWeight_electron_Reco(Electron ele1, Electron ele2);
	Double_t EfficiencySF_EventWeight_electron_RecoID(Electron ele1, Electron ele2);

	// -- outdated -- //
	void SetupEfficiencyScaleFactor();
	void SetupEfficiencyScaleFactor(TString ROOTFileName);
	//Double_t EfficiencySF_EventWeight(Muon mu1, Muon mu2, NtupleHandle *ntuple);
	//Double_t EfficiencySF_EventWeight_RecoIdIso(Muon mu1, Muon mu2, NtupleHandle *ntuple);
	//Double_t EfficiencySF_EventWeight_HLTv4p2(Muon mu1, Muon mu2);
	//Double_t EfficiencySF_EventWeight_HLTv4p3(Muon mu1, Muon mu2);
	
	////////////////////////////
	// -- Event Selections -- //
	////////////////////////////
	Bool_t EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Zpeak(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Mu50(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_minusDimuonVtxCut(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Zdiff_13TeV(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Dijet(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_Wjet(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection); // -- output: 2 muons passing event selection conditions -- //
	Bool_t EventSelection_CheckMoreThanOneDimuonCand(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, Bool_t& isMoreThanOneCand); // -- output: 2 muons passing event selection conditions -- //

	// -- for N-1 cuts of muon channel -- //
	Bool_t EventSelection_Zdiff_13TeV_HighPt1(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt2(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt3(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt4(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt5(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt6(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt7(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt8(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt9(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt10(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt11(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);
	Bool_t EventSelection_Zdiff_13TeV_HighPt12(vector< Muon > MuonCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection);

	Bool_t isPassAccCondition_Muon(Muon Mu1, Muon Mu2);
	Bool_t isPassAccCondition_GenLepton(GenLepton genlep1, GenLepton genlep2);
	void CompareMuon(Muon *Mu1, Muon *Mu2, Muon *leadMu, Muon *subMu);
	void CompareGenLepton(GenLepton *genlep1, GenLepton *genlep2, GenLepton *leadgenlep, GenLepton *subgenlep);
	void DimuonVertexProbNormChi2(NtupleHandle *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2);

	// -- for electron channel -- //
	Bool_t EventSelection_Electron(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel_NminusPFIso(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel_Zpeak(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel0(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);

	// -- for N-1 cuts of electron channel -- //
	Bool_t EventSelection_ElectronChannel1(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel2(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel3(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel4(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel5(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel6(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel7(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t EventSelection_ElectronChannel8(vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Electron >* SelectedElectronCollection);
	Bool_t isPassAccCondition_Electron(Electron Elec1, Electron Elec2);
	Bool_t isPassAccCondition_GenLepton_ECALGAP(GenLepton genlep1, GenLepton genlep2);
	void CompareElectron(Electron *Elec1, Electron *Elec2, Electron *leadElec, Electron *subElec);

	// -- pre-FSR functions -- //
	void PostToPreFSR_byDressedLepton(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection);
	void PostToPreFSR_byDressedLepton_AllPhotons(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection);
	TString DecideFSRType(GenLepton preFSR1, GenLepton preFSR2, GenLepton postFSR1, GenLepton postFSR2);
	Double_t Calc_dR_GenLeptons( GenLepton genlep1, GenLepton genlep2 );
	Double_t Calc_dR_GenLepton_GenOthers( GenLepton genlep1, GenOthers genlep2 );

	// -- miscellaneous -- //
	void GenMatching(TString MuonType, NtupleHandle* ntuple, vector<Muon>* MuonCollection);
	void ConvertToTunePInfo( Muon &mu );
	void PrintOutDoubleMuInfo( Muon mu1, Muon mu2 );
	Int_t Count_QuarkNum(TString Tag, NtupleHandle *ntuple);
	Double_t LeadEtaCorr( Int_t type, TString era, Muon mu1, Muon mu2 );

	// -- emu method -- //
	Bool_t EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection); // -- output: 1 muon and 1 electron passing event selection conditions -- //
	Bool_t EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple, vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection); // -- output: 1 muon and 1 electron passing event selection conditions -- //
	void emuVertexProbNormChi2(NtupleHandle *ntuple, Double_t ele_Pt, Double_t mu_Pt, Double_t *VtxProb, Double_t *VtxNormChi2);
	Double_t EfficiencySF_EventWeight_emu_BtoF(Muon mu1, Electron ele2);
	Double_t EfficiencySF_EventWeight_emu_GtoH(Muon mu1, Electron ele2);
};

DYAnalyzer::DYAnalyzer(TString HLTname)
{
	if( HLTname == "None" )
	{
		cout << "===================================================" << endl;
		cout << "[No specific trigger setting ... basic constructor]" << endl;
		cout << "===================================================" << endl;
		
		HLT = "None";
		LeadPtCut = 9999;
		SubPtCut = 9999;
		LeadEtaCut = 9999;
		SubEtaCut = 9999;
	}
	else
	{
		this->AssignAccThreshold(HLTname, &HLT, &LeadPtCut, &SubPtCut, &LeadEtaCut, &SubEtaCut);
		cout << "===========================================================" << endl;
		cout << "Trigger: " << HLT << endl;
		cout << "leading lepton pT Cut: " << LeadPtCut << endl;
		cout << "Sub-leading lepton pT Cut: " << SubPtCut << endl;
		cout << "leading lepton Eta Cut: " << LeadEtaCut << endl;
		cout << "sub-leading lepton Eta Cut: " << SubEtaCut << endl;
		cout << "===========================================================" << endl;
	}

}

void DYAnalyzer::AssignAccThreshold(TString HLTname, TString *HLT, Double_t *LeadPtCut, Double_t *SubPtCut, Double_t *LeadEtaCut, Double_t *SubEtaCut)
{
	if( HLTname == "IsoMu20" )
	{
		*HLT = "HLT_IsoMu20_v*";
		*LeadPtCut = 22;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu20_OR_IsoTkMu20" )
	{
		*HLT = "HLT_IsoMu20_v* || HLT_IsoTkMu20_v*";
		*LeadPtCut = 22;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu20_OR_IsoTkMu20_hw" ) // modified from "hw" at 2017.02.21
	{
		*HLT = "HLT_IsoMu20_v* || HLT_IsoTkMu20_v*";
		*LeadPtCut = 25;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu24" ) // modified from "hw2" at 2017.02.21
	{
		*HLT = "HLT_IsoMu24_v*";
		*LeadPtCut = 26;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu24_OR_IsoTkMu24" ) // added at 2017.08.01
	{
		*HLT = "HLT_IsoMu24_v* || HLT_IsoTkMu24_v*";
		//*LeadPtCut = 26;
		//*SubPtCut = 10;
		*LeadPtCut = 28;
		*SubPtCut = 17;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Mu45_eta2p1" )
	{
		*HLT = "HLT_Mu45_eta2p1_v*";
		*LeadPtCut = 46;
		*SubPtCut = 10;
		*LeadEtaCut = 2.1;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Mu50" )
	{
		*HLT = "HLT_Mu50_v*";
		*LeadPtCut = 53;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Mu50_OR_TkMu50" ) // added at 2018.02.20
	{
		*HLT = "HLT_Mu50_v* || HLT_TkMu50_v*";
		*LeadPtCut = 53;
		*SubPtCut = 10;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "IsoMu20_SymmetricPt25" )
	{
		*HLT = "HLT_IsoMu20_v*";
		*LeadPtCut = 25;
		*SubPtCut = 25;
		*LeadEtaCut = 2.4;
		*SubEtaCut = 2.4;
	}
	else if( HLTname == "Ele17Ele12" )
	{
		*HLT = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Ele22_eta2p1" )
	{
		*HLT = "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.1;
		*SubEtaCut = 2.1;
	}
	else if( HLTname == "Ele22_eta2p1_NoEtaCut" )
	{
		*HLT = "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 25;
		*SubPtCut = 15;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Pt_30_10_eta_2p5" )
	{
		*HLT = "None"; // -- just for acceptance test -- //
		*LeadPtCut = 30;
		*SubPtCut = 10;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Ele23_WPLoose" )
	{
		*HLT = "HLT_Ele23_WPLoose_Gsf_v*"; // -- Exist only for the data; "HLT_Ele22_eta2p1_WP75_Gsf_v*" should be used for MC
		*LeadPtCut = 30;
		*SubPtCut = 10;
		*LeadEtaCut = 2.5;
		*SubEtaCut = 2.5;
	}
	else if( HLTname == "Ele23Ele12" ) // updated at 2017.02.21 by Dalmin Pai
	{
		*HLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
		*LeadPtCut = 28;
		*SubPtCut = 17;
		//*LeadEtaCut = 2.5; // -- Later it should exclude ECAL gap
		//*SubEtaCut = 2.5; // -- Later it should exclude ECAL gap
		*LeadEtaCut = 2.4; // -- Later it should exclude ECAL gap
		*SubEtaCut = 2.4; // -- Later it should exclude ECAL gap
	}
	else
	{ 
		cout << "Wrong HLT name!: " << HLTname << endl;
		return; 
	}

}

void DYAnalyzer::SetupMCsamples_v20160412_76X_MINIAODv2_CheckPremix( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "DYMuMu_PU25" )
	{
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_Classic_PU25" ); Tag->push_back( "DYMuMu_M50_PU25_Classic" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8426438.0 );
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_NonDeterministic_PU25" ); Tag->push_back( "DYMuMu_M50_PU25_NonDet" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8286714.0 );
	}
	else if( Type == "DYEE_PU25" )
	{
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_Classic_PU25" ); Tag->push_back( "DYEE_M50_PU25_Classic" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8439605.0 );
		ntupleDirectory->push_back( "Premix/v20160412_76X_MINIAODv2_CheckPremix_CorrectDataSetName_DYJets_NonDeterministic_PU25" ); Tag->push_back( "DYEE_M50_PU25_NonDet" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 8300442.0 );
	}
}

void DYAnalyzer::SetupMCsamples_Moriond17( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "DYMuMu_M10to50" ) // In case of full production : nEvents = 33344559.0
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		//ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYMuMu_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		//ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYMuMu_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		//ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYMuMu_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYMuMu_M10to50_v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYMuMu_M10to50_v2" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYMuMu_M10to50_ext1v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYMuMu_M50toInf" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		//ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 1921.8 ); nEvents->push_back( 27257791.0 ); //nEvents: sum of DYMuMu weights, NNLO Xsec
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 1952.68432327 ); nEvents->push_back( 27257791.0 ); //nEvents: sum of DYMuMu weights, NNLO Xsec
	}
	else if( Type == "DYMuMu_M50to200" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		//ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 27218076.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 1873.52 + 76.2401 ); nEvents->push_back( 27218076.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYMuMu_M200toInf" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		/*ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56340.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48188.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44984.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40110.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 33360.0 ); //nEvents: sum of DYMuMu weights*/
		ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.67606 ); nEvents->push_back( 56340.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.139728 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.0792496 ); nEvents->push_back( 48188.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.0123176 ); nEvents->push_back( 44984.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.01042 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.00552772 ); nEvents->push_back( 40110.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.000741613 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.000178737 ); nEvents->push_back( 33360.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYMuMu_aMCNLO" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		/*ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYMuMu_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYMuMu_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYMuMu_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 26175605.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M100to200" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3433295.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M100to200_ext" ); Tag->push_back( "DYMuMu_M100to200_ext" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3433295.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56340.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48188.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44984.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40110.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 33360.0 ); //nEvents: sum of DYMuMu weights*/
		ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYMuMu_M10to50_v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYMuMu_M10to50_v2" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYMuMu_M10to50_ext1v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33278866.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 1873.52 + 76.2401 ); nEvents->push_back( 27218076.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.67606 ); nEvents->push_back( 56340.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.139728 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.0792496 ); nEvents->push_back( 48188.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.0123176 ); nEvents->push_back( 44984.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.01042 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.00552772 ); nEvents->push_back( 40110.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.000741613 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.000178737 ); nEvents->push_back( 33360.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "DYEE_M10to50" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYEE_M10to50_v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYEE_M10to50_v2" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYEE_M10to50_ext1v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "DYEE_M50toInf" )
	{
		 cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYEE_M50toInf" ); Xsec->push_back( 1952.68432327 ); nEvents->push_back( 27245327.0 ); //nEvents: sum of DYEE weights, NNLO Xsec
	}
	else if( Type == "DYEE_M50to200" )
	{
		cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYEE_M50to200" ); Xsec->push_back( 1873.52 + 76.2401 ); nEvents->push_back( 27206082.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "DYEE_M200toInf" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 2.67606 ); nEvents->push_back( 56144.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.139728 ); nEvents->push_back( 50420.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.0792496 ); nEvents->push_back( 48039.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.0123176 ); nEvents->push_back( 46114.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.01042 ); nEvents->push_back( 44256.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.00552772 ); nEvents->push_back( 39712.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.000741613 ); nEvents->push_back( 37287.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.000178737 ); nEvents->push_back( 34031.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "DYEE_aMCNLO" )
	{
		cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		/*ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYEE_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYEE_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYEE_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYEE_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 26166194.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M100to200" ); Tag->push_back( "DYEE_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3437885.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M100to200_ext" ); Tag->push_back( "DYEE_M100to200_ext" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 3437885.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 56144.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50420.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 48039.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 46114.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 44256.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 39712.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37287.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 34031.0 ); //nEvents: sum of DYEE weights*/
		ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYEE_M10to50_v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYEE_M10to50_v2" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYEE_M10to50_ext1v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33275218.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYEE_M50to200" ); Xsec->push_back( 1873.52 + 76.2401 ); nEvents->push_back( 27206082.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M200to400" ); Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 2.67606 ); nEvents->push_back( 56144.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M400to500" ); Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.139728 ); nEvents->push_back( 50420.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M500to700" ); Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.0792496 ); nEvents->push_back( 48039.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M700to800" ); Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.0123176 ); nEvents->push_back( 46114.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M800to1000" ); Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.01042 ); nEvents->push_back( 44256.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M1000to1500" ); Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.00552772 ); nEvents->push_back( 39712.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M1500to2000" ); Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.000741613 ); nEvents->push_back( 37287.0 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "DYLL_M2000to3000" ); Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.000178737 ); nEvents->push_back( 34031.0 ); //nEvents: sum of DYEE weights
	}
	else if( Type == "ttbar" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		//ntupleDirectory->push_back( "ttbar" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 154948878.0 ); //ttbar + ttbarBackup
		ntupleDirectory->push_back( "ttbar" ); Tag->push_back( "ttbar" ); Xsec->push_back( 734.577 ); nEvents->push_back( 135949780.0 ); //M(ttbar) < 700GeV, ttbar + ttbarBackup
	}
	else if( Type == "ttbarBackup" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		//ntupleDirectory->push_back( "ttbarBackup" ); Tag->push_back( "ttbarBackup" ); Xsec->push_back( 831.76 ); nEvents->push_back( 154948878.0 ); //ttbar + ttbarBackup
		ntupleDirectory->push_back( "ttbarBackup" ); Tag->push_back( "ttbarBackup" ); Xsec->push_back( 734.577 ); nEvents->push_back( 135949780.0 ); //M(ttbar) < 700GeV, ttbar + ttbarBackup
	}
	else if( Type == "ttbar_M700toInf" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "ttbar_M700to1000" ); Tag->push_back( "ttbar_M700to1000" ); Xsec->push_back( 76.605 ); nEvents->push_back( 38422582.0 ); //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
		ntupleDirectory->push_back( "ttbar_M1000toInf" ); Tag->push_back( "ttbar_M1000toInf" ); Xsec->push_back( 20.578 ); nEvents->push_back( 24561630.0 ); //It is not sure. (https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMonteCarlo)
	}
	else if( Type == "DYTauTau_M10to50" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		/*ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYTauTau_M10to50_v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYTauTau_M10to50_v2" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYTauTau_M10to50_ext1v1" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights*/
		ntupleDirectory->push_back( "DYLL_M10to50_v1" ); Tag->push_back( "DYTauTau_M10to50_v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "DYLL_M10to50_v2" ); Tag->push_back( "DYTauTau_M10to50_v2" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "DYLL_M10to50_ext1v1" ); Tag->push_back( "DYTauTau_M10to50_ext1v1" ); Xsec->push_back( 6016.88 ); nEvents->push_back( 33080379.0 ); //nEvents: sum of DYTauTau weights
	}
	else if( Type == "DYTauTau_M50toInf" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Signal binned samples -- //
		//ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 1921.8 ); nEvents->push_back( 27277866.0 ); //nEvents: sum of DYTauTau weights, NNLO Xsec
		ntupleDirectory->push_back( "DYLL_M50toInf" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 1952.68432327 ); nEvents->push_back( 27277866.0 ); //nEvents: sum of DYTauTau weights, NNLO Xsec
	}
	else if( Type == "VVnST" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "ST_tW" ); Tag->push_back( "tW" ); Xsec->push_back( 35.85 ); nEvents->push_back( 6952830.0 );
		ntupleDirectory->push_back( "ST_tbarW" ); Tag->push_back( "tbarW" ); Xsec->push_back( 35.85 ); nEvents->push_back( 6933093.0 );
		ntupleDirectory->push_back( "ZZ" ); Tag->push_back( "ZZ" ); Xsec->push_back( 16.523 ); nEvents->push_back( 998034.0 );
		ntupleDirectory->push_back( "WZ" ); Tag->push_back( "WZ" ); Xsec->push_back( 47.13 ); nEvents->push_back( 2995828.0 );
		ntupleDirectory->push_back( "WW" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 6987123.0 );
	}
	else if( Type == "WJetsToLNu" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "WJetsToLNu" ); Tag->push_back( "WJetsToLNu" ); Xsec->push_back( 61526.7 ); nEvents->push_back( 86731698.0 );
		ntupleDirectory->push_back( "WJetsToLNu_ext" ); Tag->push_back( "WJetsToLNu_ext" ); Xsec->push_back( 61526.7 ); nEvents->push_back( 86731698.0 );
	}
	else if( Type == "ZToMuMu_powheg" )
	{
		ntupleDirectory->push_back( "ZToMuMu_M50to120" ); Tag->push_back( "ZToMuMu_M50to120" ); Xsec->push_back(1.0); nEvents->push_back(1.0);
		//ntupleDirectory->push_back( "ZToMuMu_M120to200" ); Tag->push_back( "ZToMuMu_M120to200" ); Xsec->push_back(1.0); nEvents->push_back(1.0);
		//ntupleDirectory->push_back( "ZToMuMu_M200to400" ); Tag->push_back( "ZToMuMu_M200to400" ); Xsec->push_back(1.0); nEvents->push_back(1.0);
		//ntupleDirectory->push_back( "ZToMuMu_M400to800" ); Tag->push_back( "ZToMuMu_M400to800" ); Xsec->push_back(1.0); nEvents->push_back(1.0);
	}
	else if( Type == "ZToEE_powheg" )
	{
		ntupleDirectory->push_back( "ZToEE_M50to120" ); Tag->push_back( "ZToEE_M50to120" ); Xsec->push_back(1.0); nEvents->push_back(1.0);
	}
	else if( Type == "QCDMuEnriched" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "QCDMuEnriched_Pt15to20" ); Tag->push_back( "QCDMuEnriched_Pt15to20" ); Xsec->push_back( 720648000*0.00042 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt20to30" ); Tag->push_back( "QCDMuEnriched_Pt20to30" ); Xsec->push_back( 1273190000*0.003 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt30to50" ); Tag->push_back( "QCDMuEnriched_Pt30to50" ); Xsec->push_back( 139803000*0.01182 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt50to80" ); Tag->push_back( "QCDMuEnriched_Pt50to80" ); Xsec->push_back( 19222500*0.02276 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt80to120" ); Tag->push_back( "QCDMuEnriched_Pt80to120" ); Xsec->push_back( 2758420*0.03844 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt80to120_ext1" ); Tag->push_back( "QCDMuEnriched_Pt80to120_ext1" ); Xsec->push_back( 2758420*0.03844  ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt120to170" ); Tag->push_back( "QCDMuEnriched_Pt120to170" ); Xsec->push_back( 469797*0.05362 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt120to170_backup" ); Tag->push_back( "QCDMuEnriched_Pt120to170_backup" ); Xsec->push_back( 469797*0.05362  ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt170to300" ); Tag->push_back( "QCDMuEnriched_Pt170to300" ); Xsec->push_back( 117989*0.07335 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt170to300_ext1" ); Tag->push_back( "QCDMuEnriched_Pt170to300_ext1" ); Xsec->push_back( 117989*0.07335 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt170to300_backup" ); Tag->push_back( "QCDMuEnriched_Pt170to300_backup" ); Xsec->push_back( 117989*0.07335 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt300to470" ); Tag->push_back( "QCDMuEnriched_Pt300to470" ); Xsec->push_back( 7820.25*0.10196 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt300to470_ext1" ); Tag->push_back( "QCDMuEnriched_Pt300to470_ext1" ); Xsec->push_back( 7820.25*0.10196 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt300to470_ext2" ); Tag->push_back( "QCDMuEnriched_Pt300to470_ext2" ); Xsec->push_back( 7820.25*0.10196 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt470to600" ); Tag->push_back( "QCDMuEnriched_Pt470to600" ); Xsec->push_back( 645.528*0.12242 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt470to600_ext1" ); Tag->push_back( "QCDMuEnriched_Pt470to600_ext1" ); Xsec->push_back( 645.528*0.12242 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt470to600_ext2" ); Tag->push_back( "QCDMuEnriched_Pt470to600_ext2" ); Xsec->push_back( 645.528*0.12242 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt600to800" ); Tag->push_back( "QCDMuEnriched_Pt600to800" ); Xsec->push_back( 187.109*0.13412 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt600to800_ext1" ); Tag->push_back( "QCDMuEnriched_Pt600to800_ext1" ); Xsec->push_back( 187.109*0.13412 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt600to800_backup" ); Tag->push_back( "QCDMuEnriched_Pt600to800_backup" ); Xsec->push_back( 187.109*0.13412 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt800to1000" ); Tag->push_back( "QCDMuEnriched_Pt800to1000" ); Xsec->push_back( 32.3486*0.14552 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Tag->push_back( "QCDMuEnriched_Pt800to1000_ext1" ); Xsec->push_back( 32.3486*0.14552 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Tag->push_back( "QCDMuEnriched_Pt800to1000_ext2" ); Xsec->push_back( 32.3486*0.14552 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt1000toInf" ); Tag->push_back( "QCDMuEnriched_Pt1000toInf" ); Xsec->push_back( 10.4305*0.15544 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Tag->push_back( "QCDMuEnriched_Pt1000toInf_ext1" ); Xsec->push_back( 10.4305*0.15544 ); nEvents->push_back( 1.0 );
	}
	else if( Type == "QCDEMEnriched" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "QCDEMEnriched_Pt20to30" ); Tag->push_back( "QCDEMEnriched_Pt20to30" ); Xsec->push_back( 557600000*0.0096 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt30to50" ); Tag->push_back( "QCDEMEnriched_Pt30to50" ); Xsec->push_back( 136000000*0.073 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt30to50_ext1" ); Tag->push_back( "QCDEMEnriched_Pt30to50_ext1" ); Xsec->push_back( 136000000*0.073 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt50to80" ); Tag->push_back( "QCDEMEnriched_Pt50to80" ); Xsec->push_back( 19800000*0.146 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt50to80_ext1" ); Tag->push_back( "QCDEMEnriched_Pt50to80_ext1" ); Xsec->push_back( 19800000*0.146 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt80to120" ); Tag->push_back( "QCDEMEnriched_Pt80to120" ); Xsec->push_back( 2800000*0.125 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt80to120_ext1" ); Tag->push_back( "QCDEMEnriched_Pt80to120_ext1" ); Xsec->push_back( 2800000*0.125 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt120to170" ); Tag->push_back( "QCDEMEnriched_Pt120to170" ); Xsec->push_back( 477000*0.132 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt120to170_ext1" ); Tag->push_back( "QCDEMEnriched_Pt120to170_ext1" ); Xsec->push_back( 477000*0.132 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt170to300" ); Tag->push_back( "QCDEMEnriched_Pt170to300" ); Xsec->push_back( 114000*0.165 ); nEvents->push_back( 1.0 );
		ntupleDirectory->push_back( "QCDEMEnriched_Pt300toInf" ); Tag->push_back( "QCDEMEnriched_Pt300toInf" ); Xsec->push_back( 9000*0.15 ); nEvents->push_back( 1.0 );
	}
	else
		cout << "Wrong Type!" << endl;
}

void DYAnalyzer::SetupMCsamples_v20160309_76X_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "Full" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6302525.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6302525.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "M100to200" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 227522.0 ); //nEvents: sum of weights within 10<M<50
	}
	else if( Type == "M50to200" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6302525.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "M50toInf" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6311695.0 );
	}
	else if( Type == "M10to50_M50toInf" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6311695.0 );
	}
	else if( Type == "aMCNLO_M120Cut" )
	{
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to120" ); Xsec->push_back( 1975 ); nEvents->push_back( 6243307.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M120to200" ); Xsec->push_back( 19.32 ); nEvents->push_back( 55554.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.731 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights
	}

	else if( Type == "Full_Include_M100to200" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 6061181.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 227522.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Full_NoHighMass" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6311695.0 ); // -- sum of weight should be updated! -- //
	}
	else if( Type == "Full_Powheg" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160519_76X_MINIAODv2_Resubmit4_AdjustRunTime_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2971982.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(99600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(97600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99200.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(99200.0);
	}
	else if( Type == "Full_M120Cut" )
	{
		// cout << "# events should be adjusted later" << endl;
		// -- Background Samples -- //
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 985598 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 999996 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 988416 );
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16520811.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7467514.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6309713.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "76X/v20160303_76X_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 97994304 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to120" ); Xsec->push_back( 1975 ); nEvents->push_back( 6243307.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M120to200" ); Xsec->push_back( 19.32 ); nEvents->push_back( 55554.0 );
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.731 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_Include_M100to200")
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7506956.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to100" ); Xsec->push_back( 5869.58346/3.0 ); nEvents->push_back( 6061181.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M100to200_25ns" ); Tag->push_back( "DYMuMu_M100to200" ); Xsec->push_back( 226/3.0 ); nEvents->push_back( 227522.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 170955.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 50136.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 47833.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 44740.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 43496.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 40783.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 37176.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "76X/v20160304_76X_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 23078.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Madgraph" )
	{
		ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 631905.0 );
		ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 6014/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M5to50" ); Xsec->push_back( 7160/3.0 ); nEvents->push_back( 2782834.0 );
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 4895/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50to150" ); Xsec->push_back( (4895 - 6.58)/3.0 ); nEvents->push_back( 3003455.0 ); 
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M150toInf_25ns" ); Tag->push_back( "Madgraph_M150toInf" ); Xsec->push_back( 6.58/3.0 ); nEvents->push_back( 1.0 );

	}
	else if( Type == "MadgraphPowheg" ) // -- for estimation of syst. from unfolding -- //
	{
		ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 631905.0 );
		ntupleDirectory->push_back( "76X/v20160519_76X_MINIAODv2_Resubmit4_AdjustRunTime_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2971982.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(99600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(97600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99200.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(99200.0);


		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M5to50" ); Xsec->push_back( 7160/3.0 ); nEvents->push_back( 2782834.0 );
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 4895/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50to150" ); Xsec->push_back( (4895 - 6.58)/3.0 ); nEvents->push_back( 3003455.0 ); 
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M150toInf_25ns" ); Tag->push_back( "Madgraph_M150toInf" ); Xsec->push_back( 6.58/3.0 ); nEvents->push_back( 1.0 );

	}
	else if( Type == "Powheg" ) // -- for estimation of syst. from unfolding -- //
	{
		ntupleDirectory->push_back( "76X/v20160519_76X_MINIAODv2_Resubmit4_AdjustRunTime_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2971982.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(99999.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(99600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(97600.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99200.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160525_76X_MINIAODv2_Resubmit_HighMass_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000.0);
		ntupleDirectory->push_back( "76X/v20160404_76X_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(99200.0);


		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M5to50_25ns" ); Tag->push_back( "Madgraph_M5to50" ); Xsec->push_back( 7160/3.0 ); nEvents->push_back( 2782834.0 );
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50toInf" ); Xsec->push_back( 4895/3.0 ); nEvents->push_back( 3003455.0 );

		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M50toInf_25ns" ); Tag->push_back( "Madgraph_M50to150" ); Xsec->push_back( (4895 - 6.58)/3.0 ); nEvents->push_back( 3003455.0 ); 
		// ntupleDirectory->push_back( "76X/v20160520_76X_MINIAODv2_Madgraph_LO_M150toInf_25ns" ); Tag->push_back( "Madgraph_M150toInf" ); Xsec->push_back( 6.58/3.0 ); nEvents->push_back( 1.0 );

	}
	else
		cout << "Wrong Type!" << endl;
}

void DYAnalyzer::SetupMCsamples_v20160131_MiniAODv2( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "Full" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Full_NoHighMass" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6422093.0 ); //nEvents: sum of DYMuMu weights
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "M50_M200to400" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6422093.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "Powheg" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2836871);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99600);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(100000);
	}
	else if( Type == "Full_withoutM200to400" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to400" ); Xsec->push_back( 6103.25346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_withoutM200to400" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to400" ); Xsec->push_back( 6103.25346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		// ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_M50toInf" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50toInf" ); Xsec->push_back( 2008.4 ); nEvents->push_back( 6422093.0 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "aMCNLO_M200to400" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO_DYEE" )
	{
		// cout << "Warning: # events should be adjusted using Sum weights of DYEE events (current one: DYMuMu SumWeights)" << endl;
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYEE_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7.29361e+06 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYEE_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6.40938e+06 ); //nEvents: sum of DYEE weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYEE_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18348 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYEE_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17410 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYEE_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17245 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYEE_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 16120 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYEE_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14397 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYEE_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 13857 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYEE_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13495 ); //nEvents: sum of DYEE weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYEE_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12859 ); //nEvents: sum of DYEE weights 
	}
	else if( Type == "aMCNLO_FEWZxSec" )
	{
		// xSec of M10-50 and M50 sample: aMC@NLO -- //
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 2.59583 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.136235 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.0775862 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.0121251 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.010281 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.00546713 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.000735022 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.000176089 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
	}
	else
		cout << "Wrong Type!" << endl;

	return;
}

void DYAnalyzer::SetupMCsamples_v20160117_MiniAOD_JetMET( TString Type, vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	if( Type == "Full" )
	{
		// -- Background Samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996168 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 991232 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 994416 );
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16518173 ); //nEvents: sum of weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7418362 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6430407 ); //nEvents: sum of DYTauTau weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19899492 );
		
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7418362 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6430407 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18339 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 6951 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376 ); //nEvents: sum of DYMuMu weights 
	}
	else if( Type == "aMCNLO" )
	{
		// -- Signal binned samples -- //
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7418362 ); //nEvents: sum of weights within 10<M<50
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6430407 ); //nEvents: sum of DYMuMu weights
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18339 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 6951 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1000to1500_25ns" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855 ); //nEvents: sum of DYMuMu weights 
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376 ); //nEvents: sum of DYMuMu weights
	}
	else if( Type == "Powheg" )
	{
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M50to120_25ns" ); Tag->push_back( "ZMuMu_M50to120" );  Xsec->push_back(1975);  nEvents->push_back(2848071);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M120to200_25ns" ); Tag->push_back( "ZMuMu_M120to200" );  Xsec->push_back(19.32);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M200to400_25ns" ); Tag->push_back( "ZMuMu_M200to400" );  Xsec->push_back(2.731);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M400to800_25ns" ); Tag->push_back( "ZMuMu_M400to800" );  Xsec->push_back(0.241);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M800to1400_25ns" );  Tag->push_back( "ZMuMu_M800to1400" );  Xsec->push_back(0.01678);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M1400to2300_25ns" );  Tag->push_back( "ZMuMu_M1400to2300" );  Xsec->push_back(0.00139);  nEvents->push_back(99600);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M2300to3500_25ns" );  Tag->push_back( "ZMuMu_M2300to3500" );  Xsec->push_back(0.00008948);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M3500to4500_25ns" );  Tag->push_back( "ZMuMu_M3500to4500" );  Xsec->push_back(0.000004135);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M4500to6000_25ns" );  Tag->push_back( "ZMuMu_M4500to6000" );  Xsec->push_back(4.56E-07);  nEvents->push_back(100000);
		ntupleDirectory->push_back( "Spring15DR/25ns/v20160102_MINIAOD_AddJetMET_ZMuMuPowheg_M6000toInf_25ns" );  Tag->push_back( "ZMuMu_M6000toInf" );  Xsec->push_back(2.066E-08);  nEvents->push_back(100000);
	}
	else
		cout << "Wrong Type!" << endl;

	return;
}

Int_t DYAnalyzer::Count_QuarkNum(TString Tag, NtupleHandle *ntuple)
{
	Int_t nQuarks = -9999;

	vector<GenOthers> GenOthersCollection;
	Int_t NGenOthers = ntuple->nGenOthers;
	for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
	{
		GenOthers genothers;
		genothers.FillFromNtuple(ntuple, i_gen);
		if( 0 < abs(genothers.ID) && abs(genothers.ID) < 7 && genothers.isHardProcess )
			GenOthersCollection.push_back( genothers );
	}

	nQuarks = GenOthersCollection.size();
	if( nQuarks > 5 || nQuarks == -9999 ) printf("Number of Quarks = %d\n", nQuarks);

	return nQuarks;
}

Bool_t DYAnalyzer::Separate_ttbarSample(TString Tag, NtupleHandle *ntuple, vector<GenOthers> *GenTopCollection)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate ttbar events -- //
	if( Tag.Contains("ttbar") )
	{
		vector<GenOthers> GenOthersCollection;
		Int_t NGenOthers = ntuple->nGenOthers;
		for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
		{
			GenOthers genothers;
			genothers.FillFromNtuple(ntuple, i_gen);
			if( abs(genothers.ID) == 6 && genothers.isHardProcess )
				GenOthersCollection.push_back( genothers );
		}

		if( GenOthersCollection.size() == 2 ) // -- Select the ttbar events from hard-process -- //
		{
			// -- Check top & anti-top pair -- //
			if( GenOthersCollection[0].ID == GenOthersCollection[1].ID )
				printf("%d %d\n", GenOthersCollection[0].ID, GenOthersCollection[1].ID);

			//if( Tag == "ttbar" ) // -- Select only evetns withtin M < 700 -- //
			if( Tag == "ttbar" || Tag == "ttbarBackup" ) // -- Select only evetns withtin M < 700 -- //
			{
				TLorentzVector v1 = GenOthersCollection[0].Momentum;
				TLorentzVector v2 = GenOthersCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 700 )
				//if( reco_M > -999 )
				{
					GenFlag = kTRUE;
					GenTopCollection->push_back( GenOthersCollection[0] );
					GenTopCollection->push_back( GenOthersCollection[1] );
				}
			}
			else // ex: ttbar_M700to1000, ttbar_M1000toInf
			{
				GenFlag = kTRUE;
				GenTopCollection->push_back( GenOthersCollection[0] );
				GenTopCollection->push_back( GenOthersCollection[1] );
			}
		}
		else
		{
			printf("Wrong? : more than two!!\n");
			printf("%d %d %d\n", GenOthersCollection[0].ID, GenOthersCollection[1].ID, GenOthersCollection[2].ID); //Check upto 3rd one
		}
	}
	// -- other cases(e.g. DY, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

Bool_t DYAnalyzer::SeparateDYLLSample_isHardProcess(TString Tag, NtupleHandle *ntuple)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate DYMuMu events from DYTauTau  -- //
	if( Tag.Contains("DYMuMu") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
		{
			if( Tag == "DYMuMu_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to400" ) // -- Select only evetns withtin 50 < M < 400 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 400 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to100" || Tag == "DYMuMu_Photos_M50to100" ) // -- Select only evetns withtin 50 < M < 100 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 100 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to120" )
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 120 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M120to200" )
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M > 120 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	else if( Tag.Contains("DYEE") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isElectron() && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 electrons from hard-process -- //
		{
			if( Tag == "DYEE_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYEE_M50to100" ) // -- Select only evetns withtin 50 < M < 100 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 100 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- Separate DYTauTau events from MuMu events -- //
	else if( Tag.Contains("DYTauTau") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( abs(genlep.ID) == 15 && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 taus from hard-process -- //
		{
			GenFlag = kTRUE;
		}
	}
	// -- Madgraph sample -- //
	else if( Tag.Contains("Madgraph") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isHardProcess )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 muons from hard-process -- //
		{
			if( Tag == "Madgraph_M50to150" ) // -- Select only evetns withtin 50 < M < 150 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 150 )
					GenFlag = kTRUE;
			}
			else if( Tag == "Madgraph_M10to50" )
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M > 10 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- other cases(e.g. ttbar, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

Bool_t DYAnalyzer::SeparateDYLLSample(TString Tag, NtupleHandle *ntuple)
{
	Bool_t GenFlag = kFALSE;

	// -- Seperate DYMuMu events from DYTauTau  -- //
	if( Tag.Contains("DYMuMu") )
	{
		vector<GenLepton> GenLeptonCollection;
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
			if( Tag == "DYMuMu_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else if( Tag == "DYMuMu_M50to400" ) // -- Select only evetns withtin 50 < M < 400 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 400 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	else if( Tag.Contains("DYEE") )
	{
		vector<GenLepton> GenLeptonCollection;
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
			if( Tag == "DYEE_M50to200" ) // -- Select only evetns withtin 50 < M < 200 -- //
			{
				TLorentzVector v1 = GenLeptonCollection[0].Momentum;
				TLorentzVector v2 = GenLeptonCollection[1].Momentum;
				Double_t reco_M = (v1 + v2).M();
				if( reco_M < 200 )
					GenFlag = kTRUE;
			}
			else
				GenFlag = kTRUE;
		}
	}
	// -- Separate DYTauTau events from MuMu events -- //
	else if( Tag.Contains("DYTauTau") )
	{
		vector<GenLepton> GenLeptonCollection;
		Int_t NGenLeptons = ntuple->gnpair;
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( abs(genlep.ID) == 15 && genlep.fromHardProcessDecayed )
				GenLeptonCollection.push_back( genlep );
		}

		if( GenLeptonCollection.size() == 2 ) // -- Select the events containing 2 taus from hard-process -- //
		{
			GenFlag = kTRUE;
		}
	}
	// -- other cases(e.g. ttbar, WJets, Diboson...): pass
	else
		GenFlag = kTRUE; 

	return GenFlag;
}

void DYAnalyzer::SetupPileUpReWeighting( Bool_t isMC )
{
	if( isMC == kFALSE ) // -- for data -- //
	{
		for(Int_t i=0; i<52; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TFile *f = new TFile("/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_PUReWeight_v20160208_2nd_71mb.root");
	f->cd();
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
	if( h_weight == NULL )
	{
		cout << "ERROR! ... No Weight histogram!"<< endl;
		return;
	}

	for(Int_t i=0; i<52; i++)
	{
		Int_t i_bin = i+1;
		PileUpWeight[i] = h_weight->GetBinContent(i_bin);
	}
}

Double_t DYAnalyzer::PileUpWeightValue(Int_t PileUp_MC)
{
	if( PileUp_MC < 0 || PileUp_MC > 51 )
	{
		cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

void DYAnalyzer::SetupPileUpReWeighting_76X( Bool_t isMC )
{
	if( isMC == kFALSE ) // -- for data -- //
	{
		for(Int_t i=0; i<50; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TString FileLocation = "/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_PUReWeight_76X_v20160404_71mb.root";
	TFile *f = new TFile(FileLocation);
	f->cd();
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
	if( h_weight == NULL )
	{
		cout << "ERROR! ... No Weight histogram!"<< endl;
		return;
	}

	for(Int_t i=0; i<50; i++)
	{
		Int_t i_bin = i+1;
		PileUpWeight[i] = h_weight->GetBinContent(i_bin);
	}
}

Double_t DYAnalyzer::PileUpWeightValue_76X(Int_t PileUp_MC)
{
	if( PileUp_MC < 0 || PileUp_MC > 49 )
	{
		cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

void DYAnalyzer::SetupPileUpReWeighting_80X( Bool_t isMC, TString ROOTFileName )
{
	if( isMC == kFALSE ) // -- for data -- //
	{
		for(Int_t i=0; i<75; i++)
			PileUpWeight[i] = 1;

		return;
	}
	
	// -- Only for the MC -- //
	TString FileLocation = "../../TOOL/PileUp/80X/"+ROOTFileName;
	TFile *f = new TFile(FileLocation);
	TH1D *h_weight = (TH1D*)f->Get("h_PUReWeights");
	if( h_weight == NULL )
	{
		cout << "ERROR! ... No Weight histogram!"<< endl;
		return;
	}

	for(Int_t i=0; i<75; i++)
	{
		Int_t i_bin = i+1;
		PileUpWeight[i] = h_weight->GetBinContent(i_bin);
	}
}

Double_t DYAnalyzer::PileUpWeightValue_80X(Int_t PileUp_MC)
{
	if( PileUp_MC < 0 || PileUp_MC > 74 )
	{
		cout << "[PileUp_MC = " << PileUp_MC << "]: NO CORRESPONDING PU Weight! ... it returns 0" << endl;
		return 0;
	}
	return PileUpWeight[PileUp_MC];
}

Double_t DYAnalyzer::LeadEtaCorr( Int_t type, TString era, Muon mu1, Muon mu2 )
{
	Double_t Event_SF = 1;

	if( type == 12 )
	{
		Event_SF = 0;

		Muon leadMu, subMu;
		CompareMuon(&mu1, &mu2, &leadMu, &subMu);

		// Parameters of 8th order polynomial, starting from lower order
		Double_t BtoF[9] = {1.020775e+00, 1.141423e-02, 2.538435e-02, -1.809580e-02, -3.598146e-02, 7.254615e-03, 9.862204e-03, -8.151668e-04, -8.768882e-04};
		Double_t GtoH[9] = {1.013401e+00, 2.817957e-03, 9.194193e-03, -1.000880e-02, -7.600872e-03, 4.595650e-03, -2.718031e-04, -5.520347e-04, 1.782939e-04};
		Double_t BtoH[9] = {1.017306e+00, 7.403816e-03, 1.796262e-02, -1.432189e-02, -2.284996e-02, 6.013502e-03, 5.160505e-03, -6.922301e-04, -3.861758e-04};

		if( era == "BtoF" )
			for(Int_t i_order=0; i_order<=8; i_order++)
				Event_SF += BtoF[i_order] * TMath::Power(leadMu.eta, i_order);
		else if( era == "GtoH" )
			for(Int_t i_order=0; i_order<=8; i_order++)
				Event_SF += GtoH[i_order] * TMath::Power(leadMu.eta, i_order);
		else if( era == "BtoH" )
			for(Int_t i_order=0; i_order<=8; i_order++)
				Event_SF += BtoH[i_order] * TMath::Power(leadMu.eta, i_order);
		else
			cout << "[LeadEtaCorr] You chose wrong era!!" << endl;
	}

	//cout << "Event_SF : " << Event_SF << endl;
	return Event_SF;
}

void DYAnalyzer::SetupEfficiencyScaleFactor()
{
	TString Location_TnP = "/home/kplee/CommonCodes/DrellYanAnalysis/ROOTFile_TagProbeEfficiency_v20160329.root";
	cout << "[Tag&Probe efficiency is from " << Location_TnP << " (Default, 74X)]" << endl;
	
	TFile *f = new TFile( Location_TnP );
	TH2D *h_RecoID_data = (TH2D*)f->Get("h_2D_Eff_RecoID_Data");
	TH2D *h_RecoID_MC = (TH2D*)f->Get("h_2D_Eff_RecoID_MC");

	TH2D *h_Iso_data = (TH2D*)f->Get("h_2D_Eff_Iso_Data");
	TH2D *h_Iso_MC = (TH2D*)f->Get("h_2D_Eff_Iso_MC");

	TH2D *h_HLTv4p2_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_Data");
	TH2D *h_HLTv4p2_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_MC");

	TH2D *h_HLTv4p3_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_Data");
	TH2D *h_HLTv4p3_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_MC");


	Int_t nEtaBins = h_RecoID_data->GetNbinsX();
	Int_t nPtBins = h_RecoID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t RecoID_data = h_RecoID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t RecoID_MC = h_RecoID_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p2_data = h_HLTv4p2_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p2_MC = h_HLTv4p2_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p3_data = h_HLTv4p3_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p3_MC = h_HLTv4p3_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_RecoID_data[iter_x][iter_y] = RecoID_data;
			Eff_RecoID_MC[iter_x][iter_y] = RecoID_MC;

			Eff_Iso_data[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC[iter_x][iter_y] = Iso_MC;

			Eff_HLTv4p2_data[iter_x][iter_y] = HLTv4p2_data;
			Eff_HLTv4p2_MC[iter_x][iter_y] = HLTv4p2_MC;

			Eff_HLTv4p3_data[iter_x][iter_y] = HLTv4p3_data;
			Eff_HLTv4p3_MC[iter_x][iter_y] = HLTv4p3_MC;
		}
	}
	cout << "Setting for efficiency correction factors is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor(TString ROOTFileName)
{
	TString Location_TnP = "/home/kplee/CommonCodes/DrellYanAnalysis/"+ROOTFileName;
	cout << "[Tag&Probe efficiency is from " << Location_TnP << "]" << endl; 

	TFile *f = new TFile( Location_TnP );

	TH2D *h_RecoID_data = (TH2D*)f->Get("h_2D_Eff_RecoID_Data");
	TH2D *h_RecoID_MC = (TH2D*)f->Get("h_2D_Eff_RecoID_MC");

	TH2D *h_Iso_data = (TH2D*)f->Get("h_2D_Eff_Iso_Data");
	TH2D *h_Iso_MC = (TH2D*)f->Get("h_2D_Eff_Iso_MC");

	TH2D *h_HLTv4p2_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_Data");
	TH2D *h_HLTv4p2_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p2_MC");

	TH2D *h_HLTv4p3_data = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_Data");
	TH2D *h_HLTv4p3_MC = (TH2D*)f->Get("h_2D_Eff_HLTv4p3_MC");


	Int_t nEtaBins = h_RecoID_data->GetNbinsX();
	Int_t nPtBins = h_RecoID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t RecoID_data = h_RecoID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t RecoID_MC = h_RecoID_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p2_data = h_HLTv4p2_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p2_MC = h_HLTv4p2_MC->GetBinContent(i_etabin, i_ptbin);

			Double_t HLTv4p3_data = h_HLTv4p3_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLTv4p3_MC = h_HLTv4p3_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_RecoID_data[iter_x][iter_y] = RecoID_data;
			Eff_RecoID_MC[iter_x][iter_y] = RecoID_MC;

			Eff_Iso_data[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC[iter_x][iter_y] = Iso_MC;

			Eff_HLTv4p2_data[iter_x][iter_y] = HLTv4p2_data;
			Eff_HLTv4p2_MC[iter_x][iter_y] = HLTv4p2_MC;

			Eff_HLTv4p3_data[iter_x][iter_y] = HLTv4p3_data;
			Eff_HLTv4p3_MC[iter_x][iter_y] = HLTv4p3_MC;
		}
	}
	cout << "Setting for efficiency correction factors is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor_BtoF()
{
	//TString Location = "/home/dmpai/effSF_muon/";
	TString Location = "../../TOOL/effSF/effSF_muon/";
	cout << "[Tag&Probe efficiency is from " << Location+"*BtoF.root" << "]" << endl; 

	///////////////////
	// -- Reco SF -- //
	///////////////////
	TFile *f0 = new TFile( Location+"Tracking_SF_RunBtoF.root" );
	TGraphAsymmErrors *h_Reco_ratio = (TGraphAsymmErrors*)f0->Get("ratio_eff_eta3_dr030e030_corr");

	Int_t nEtaBins_reco = h_Reco_ratio->GetN();
	Int_t nPtBins_reco = 1;

	Double_t x_reco[15]; Double_t xx_reco[15];
	Double_t y_reco[15]; Double_t yy_reco[15];

	for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++)
	{
		for(Int_t i=0; i<nEtaBins_reco; i++)
		{
			h_Reco_ratio->GetPoint(i, x_reco[i], y_reco[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_reco; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_reco; k++)
			{
				if( xlow < x_reco[k] && x_reco[k] < xmin )
				{
					jj = k;
					xmin = x_reco[k];
				}
			}
			xx_reco[j] = x_reco[jj];
			yy_reco[j] = y_reco[jj];

			xlow = xmin;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
		{
			Eff_Reco_data_BtoF[iter_x][iter_y] = yy_reco[iter_x]; // actually, it is the scale factor.
			Eff_Reco_MC_BtoF[iter_x][iter_y] = 1;
		}
	}
	/////////////////
	// -- ID SF -- //
	/////////////////
	TFile *f1 = new TFile( Location+"ID_SF_RunBtoF.root" );
	//TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA"); //HighPt
	//TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesMC/abseta_pair_ne_MC");
	TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA"); //Tight
	TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

	Int_t nEtaBins1 = h_ID_data->GetNbinsX();
	Int_t nPtBins1 = h_ID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_ID_data_BtoF[iter_x][iter_y] = ID_data;
			Eff_ID_MC_BtoF[iter_x][iter_y] = ID_MC;
		}
	}
	//////////////////
	// -- Iso SF -- //
	//////////////////
	TFile *f2 = new TFile( Location+"ISO_SF_RunBtoF.root" );
	//TH2F *h_Iso_data = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA"); //trkiso
	//TH2F *h_Iso_MC = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesMC/abseta_pair_ne_MC"); //trkiso
	TH2F *h_Iso_data = (TH2F*)f2->Get("TightISO_TightID_pt_eta/efficienciesDATA/abseta_pt_DATA"); //pfiso
	TH2F *h_Iso_MC = (TH2F*)f2->Get("TightISO_TightID_pt_eta/efficienciesMC/abseta_pt_MC"); //pfiso

	Int_t nEtaBins3 = h_Iso_data->GetNbinsX(); //pfiso
	Int_t nPtBins3 = h_Iso_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins3; iter_x++) //pfiso
	{
		for(Int_t iter_y = 0; iter_y < nPtBins3; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_Iso_data_BtoF[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC_BtoF[iter_x][iter_y] = Iso_MC;
		}
	}
	//////////////////////
	// -- Trigger SF -- //
	//////////////////////
	TFile *f3 = new TFile( Location+"Trigger_SF_RunBtoF.root" );
	TH2F *h_HLT_data = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
	TH2F *h_HLT_MC = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/abseta_pt_MC");

	Int_t nEtaBins2 = h_HLT_data->GetNbinsX();
	Int_t nPtBins2 = h_HLT_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins2; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins2; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t HLT_data = h_HLT_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLT_MC = h_HLT_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_HLT_data_BtoF[iter_x][iter_y] = HLT_data;
			Eff_HLT_MC_BtoF[iter_x][iter_y] = HLT_MC;
		}
	}
	cout << "Setting for efficiency correction factors (BtoF) is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor_GtoH()
{
	//TString Location = "/home/dmpai/effSF_muon/";
	TString Location = "../../TOOL/effSF/effSF_muon/";
	cout << "[Tag&Probe efficiency is from " << Location+"*GtoH.root" << "]" << endl; 

	///////////////////
	// -- Reco SF -- //
	///////////////////
	TFile *f0 = new TFile( Location+"Tracking_SF_RunGtoH.root" );
	TGraphAsymmErrors *h_Reco_ratio = (TGraphAsymmErrors*)f0->Get("ratio_eff_eta3_dr030e030_corr");

	Int_t nEtaBins_reco = h_Reco_ratio->GetN();
	Int_t nPtBins_reco = 1;

	Double_t x_reco[15]; Double_t xx_reco[15];
	Double_t y_reco[15]; Double_t yy_reco[15];

	for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++)
	{
		for(Int_t i=0; i<nEtaBins_reco; i++)
		{
			h_Reco_ratio->GetPoint(i, x_reco[i], y_reco[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_reco; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_reco; k++)
			{
				if( xlow < x_reco[k] && x_reco[k] < xmin )
				{
					jj = k;
					xmin = x_reco[k];
				}
			}
			xx_reco[j] = x_reco[jj];
			yy_reco[j] = y_reco[jj];

			xlow = xmin;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
		{
			Eff_Reco_data_GtoH[iter_x][iter_y] = yy_reco[iter_x]; // actually, it is the scale factor.
			Eff_Reco_MC_GtoH[iter_x][iter_y] = 1;
		}
	}
	/////////////////
	// -- ID SF -- //
	/////////////////
	TFile *f1 = new TFile( Location+"ID_SF_RunGtoH.root" );
	//TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA"); //HighPt
	//TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/efficienciesMC/abseta_pair_ne_MC");
	TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA"); //Tight
	TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

	Int_t nEtaBins1 = h_ID_data->GetNbinsX();
	Int_t nPtBins1 = h_ID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_ID_data_GtoH[iter_x][iter_y] = ID_data;
			Eff_ID_MC_GtoH[iter_x][iter_y] = ID_MC;
		}
	}
	//////////////////
	// -- Iso SF -- //
	//////////////////
	TFile *f2 = new TFile( Location+"ISO_SF_RunGtoH.root" );
	//TH2F *h_Iso_data = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesDATA/abseta_pair_ne_DATA"); //trkiso
	//TH2F *h_Iso_MC = (TH2F*)f2->Get("tkLooseISO_highptID_newpt_eta/efficienciesMC/abseta_pair_ne_MC"); //trkiso
	TH2F *h_Iso_data = (TH2F*)f2->Get("TightISO_TightID_pt_eta/efficienciesDATA/abseta_pt_DATA"); //pfiso
	TH2F *h_Iso_MC = (TH2F*)f2->Get("TightISO_TightID_pt_eta/efficienciesMC/abseta_pt_MC"); //pfiso

	Int_t nEtaBins3 = h_Iso_data->GetNbinsX(); //pfiso
	Int_t nPtBins3 = h_Iso_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins3; iter_x++) //pfiso
	{
		for(Int_t iter_y = 0; iter_y < nPtBins3; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t Iso_data = h_Iso_data->GetBinContent(i_etabin, i_ptbin);
			Double_t Iso_MC = h_Iso_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_Iso_data_GtoH[iter_x][iter_y] = Iso_data;
			Eff_Iso_MC_GtoH[iter_x][iter_y] = Iso_MC;
		}
	}
	//////////////////////
	// -- Trigger SF -- //
	//////////////////////
	TFile *f3 = new TFile( Location+"Trigger_SF_RunGtoH.root" );
	TH2F *h_HLT_data = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
	TH2F *h_HLT_MC = (TH2F*)f3->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/abseta_pt_MC");

	Int_t nEtaBins2 = h_HLT_data->GetNbinsX();
	Int_t nPtBins2 = h_HLT_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins2; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins2; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t HLT_data = h_HLT_data->GetBinContent(i_etabin, i_ptbin);
			Double_t HLT_MC = h_HLT_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_HLT_data_GtoH[iter_x][iter_y] = HLT_data;
			Eff_HLT_MC_GtoH[iter_x][iter_y] = HLT_MC;
		}
	}
	cout << "Setting for efficiency correction factors (GtoH) is completed" << endl;
}


void DYAnalyzer::SetupEfficiencyScaleFactor_electron()
{
	//TString Location = "/home/dmpai/effSF_electron/";
	TString Location = "../../TOOL/effSF/effSF_electron/";
	cout << "[Tag&Probe efficiency is from " << Location+"*.root" << "]" << endl; 

	///////////////////
	// -- Reco SF -- //
	///////////////////
	TFile *f1 = new TFile( Location+"Reco_SF.root" );
	TGraphErrors *h_reco_sf = (TGraphErrors*)f1->Get("grSF1D_0");

	Int_t nEtaBins_reco = h_reco_sf->GetN();
	Int_t nPtBins_reco = 1;

	Double_t x_reco[30]; Double_t xx_reco[30];
	Double_t y_reco[30]; Double_t yy_reco[30];

	for(Int_t iter_y = 0; iter_y < nPtBins_reco; iter_y++)
	{
		for(Int_t i=0; i<nEtaBins_reco; i++)
		{
			h_reco_sf->GetPoint(i, x_reco[i], y_reco[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_reco; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_reco; k++)
			{
				if( xlow < x_reco[k] && x_reco[k] < xmin )
				{
					jj = k;
					xmin = x_reco[k];
				}
			}
			xx_reco[j] = x_reco[jj];
			yy_reco[j] = y_reco[jj];

			xlow = xmin;
//			cout << j << "  " << xx_reco[j] << "  " << yy_reco[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_reco; iter_x++)
		{
			Eff_Reco_data[iter_x][iter_y] = yy_reco[iter_x]; // actually, it is the scale factor.
			Eff_Reco_MC[iter_x][iter_y] = 1;
//			cout << "Reco: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << xx_reco[iter_x] << " sf = " << yy_reco[iter_x] << endl;
		}
	}
	/////////////////
	// -- ID SF -- //
	/////////////////
	TFile *f2 = new TFile( Location+"MediumID_SF.root" );
	TGraphErrors *h_id_sf_0 = (TGraphErrors*)f2->Get("grSF1D_0");
	TGraphErrors *h_id_sf_1 = (TGraphErrors*)f2->Get("grSF1D_1");
	TGraphErrors *h_id_sf_2 = (TGraphErrors*)f2->Get("grSF1D_2");
	TGraphErrors *h_id_sf_3 = (TGraphErrors*)f2->Get("grSF1D_3");
	TGraphErrors *h_id_sf_4 = (TGraphErrors*)f2->Get("grSF1D_4");

	Int_t nEtaBins_id = h_id_sf_0->GetN();
	Int_t nPtBins_id = 5;

	Double_t x_id[10]; Double_t xx_id[10];
	Double_t y_id[10]; Double_t yy_id[10];

	TGraphErrors *h_id_sf;
	for(Int_t iter_y = 0; iter_y < nPtBins_id; iter_y++)
	{
		if(iter_y == 0) h_id_sf = (TGraphErrors*)h_id_sf_0->Clone();
		else if(iter_y == 1) h_id_sf = (TGraphErrors*)h_id_sf_1->Clone();
		else if(iter_y == 2) h_id_sf = (TGraphErrors*)h_id_sf_2->Clone();
		else if(iter_y == 3) h_id_sf = (TGraphErrors*)h_id_sf_3->Clone();
		else if(iter_y == 4) h_id_sf = (TGraphErrors*)h_id_sf_4->Clone();

		for(Int_t i=0; i<nEtaBins_id; i++)
		{
			h_id_sf->GetPoint(i, x_id[i], y_id[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_id; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_id; k++)
			{
				if( xlow < x_id[k] && x_id[k] < xmin )
				{
					jj = k;
					xmin = x_id[k];
				}
			}
			xx_id[j] = x_id[jj];
			yy_id[j] = y_id[jj];

			xlow = xmin;
			//cout << j << "  " << xx_id[j] << "  " << yy_id[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_id; iter_x++)
		{
			Eff_ID_data[iter_x][iter_y] = yy_id[iter_x]; // actually, it is the scale factor.
			Eff_ID_MC[iter_x][iter_y] = 1;
			//cout << "ID: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << xx_id[iter_x] << " sf = " << yy_id[iter_x] << endl;
		}
	}
	/////////////////////////////
	// -- Trigger (leg1) SF -- //
	/////////////////////////////
	TFile *f3 = new TFile( Location+"Leg1_SF.root" );
	TGraphErrors *h_leg1_sf_0 = (TGraphErrors*)f3->Get("grSF1D_0");
	TGraphErrors *h_leg1_sf_1 = (TGraphErrors*)f3->Get("grSF1D_1");
	TGraphErrors *h_leg1_sf_2 = (TGraphErrors*)f3->Get("grSF1D_2");
	TGraphErrors *h_leg1_sf_3 = (TGraphErrors*)f3->Get("grSF1D_3");
	TGraphErrors *h_leg1_sf_4 = (TGraphErrors*)f3->Get("grSF1D_4");
	TGraphErrors *h_leg1_sf_5 = (TGraphErrors*)f3->Get("grSF1D_5");
	TGraphErrors *h_leg1_sf_6 = (TGraphErrors*)f3->Get("grSF1D_6");
	TGraphErrors *h_leg1_sf_7 = (TGraphErrors*)f3->Get("grSF1D_7");

	Int_t nEtaBins_leg1 = h_leg1_sf_0->GetN();
	Int_t nPtBins_leg1 = 8;

	Double_t x_leg1[10]; Double_t xx_leg1[10];
	Double_t y_leg1[10]; Double_t yy_leg1[10];

	TGraphErrors *h_leg1_sf;
	for(Int_t iter_y = 0; iter_y < nPtBins_leg1; iter_y++)
	{
		if(iter_y == 0) h_leg1_sf = (TGraphErrors*)h_leg1_sf_0->Clone();
		else if(iter_y == 1) h_leg1_sf = (TGraphErrors*)h_leg1_sf_1->Clone();
		else if(iter_y == 2) h_leg1_sf = (TGraphErrors*)h_leg1_sf_2->Clone();
		else if(iter_y == 3) h_leg1_sf = (TGraphErrors*)h_leg1_sf_3->Clone();
		else if(iter_y == 4) h_leg1_sf = (TGraphErrors*)h_leg1_sf_4->Clone();
		else if(iter_y == 5) h_leg1_sf = (TGraphErrors*)h_leg1_sf_5->Clone();
		else if(iter_y == 6) h_leg1_sf = (TGraphErrors*)h_leg1_sf_6->Clone();
		else if(iter_y == 7) h_leg1_sf = (TGraphErrors*)h_leg1_sf_7->Clone();

		for(Int_t i=0; i<nEtaBins_leg1; i++)
		{
			h_leg1_sf->GetPoint(i, x_leg1[i], y_leg1[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_leg1; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_leg1; k++)
			{
				if( xlow < x_leg1[k] && x_leg1[k] < xmin )
				{
					jj = k;
					xmin = x_leg1[k];
				}
			}
			xx_leg1[j] = x_leg1[jj];
			yy_leg1[j] = y_leg1[jj];

			xlow = xmin;
			//cout << j << "  " << xx_leg1[j] << "  " << yy_leg1[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_leg1; iter_x++)
		{
			Eff_HLT_Leg1_data[iter_x][iter_y] = yy_leg1[iter_x]; // actually, it is the scale factor.
			Eff_HLT_Leg1_MC[iter_x][iter_y] = 1;
			//cout << "ID: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << xx_leg1[iter_x] << " sf = " << yy_leg1[iter_x] << endl;
		}
	}
	/////////////////////////////
	// -- Trigger (leg2) SF -- //
	/////////////////////////////
	TFile *f4 = new TFile( Location+"Leg2_SF.root" );
	TGraphErrors *h_leg2_sf_0 = (TGraphErrors*)f4->Get("grSF1D_0");
	TGraphErrors *h_leg2_sf_1 = (TGraphErrors*)f4->Get("grSF1D_1");
	TGraphErrors *h_leg2_sf_2 = (TGraphErrors*)f4->Get("grSF1D_2");
	TGraphErrors *h_leg2_sf_3 = (TGraphErrors*)f4->Get("grSF1D_3");
	TGraphErrors *h_leg2_sf_4 = (TGraphErrors*)f4->Get("grSF1D_4");
	TGraphErrors *h_leg2_sf_5 = (TGraphErrors*)f4->Get("grSF1D_5");
	TGraphErrors *h_leg2_sf_6 = (TGraphErrors*)f4->Get("grSF1D_6");
	TGraphErrors *h_leg2_sf_7 = (TGraphErrors*)f4->Get("grSF1D_7");

	Int_t nEtaBins_leg2 = h_leg2_sf_0->GetN();
	Int_t nPtBins_leg2 = 8;

	Double_t x_leg2[10]; Double_t xx_leg2[10];
	Double_t y_leg2[10]; Double_t yy_leg2[10];

	TGraphErrors *h_leg2_sf;
	for(Int_t iter_y = 0; iter_y < nPtBins_leg2; iter_y++)
	{
		if(iter_y == 0) h_leg2_sf = (TGraphErrors*)h_leg2_sf_0->Clone();
		else if(iter_y == 1) h_leg2_sf = (TGraphErrors*)h_leg2_sf_1->Clone();
		else if(iter_y == 2) h_leg2_sf = (TGraphErrors*)h_leg2_sf_2->Clone();
		else if(iter_y == 3) h_leg2_sf = (TGraphErrors*)h_leg2_sf_3->Clone();
		else if(iter_y == 4) h_leg2_sf = (TGraphErrors*)h_leg2_sf_4->Clone();
		else if(iter_y == 5) h_leg2_sf = (TGraphErrors*)h_leg2_sf_5->Clone();
		else if(iter_y == 6) h_leg2_sf = (TGraphErrors*)h_leg2_sf_6->Clone();
		else if(iter_y == 7) h_leg2_sf = (TGraphErrors*)h_leg2_sf_7->Clone();

		for(Int_t i=0; i<nEtaBins_leg2; i++)
		{
			h_leg2_sf->GetPoint(i, x_leg2[i], y_leg2[i]);
		}

		// -- Rearrangement in order of x values -- //
		Double_t xmin;
		Double_t xlow = -9999;
		for(Int_t j=0; j<nEtaBins_leg2; j++)
		{
			Int_t jj = -9999;

			xmin = 9999;
			for(Int_t k=0; k<nEtaBins_leg2; k++)
			{
				if( xlow < x_leg2[k] && x_leg2[k] < xmin )
				{
					jj = k;
					xmin = x_leg2[k];
				}
			}
			xx_leg2[j] = x_leg2[jj];
			yy_leg2[j] = y_leg2[jj];

			xlow = xmin;
			//cout << j << "  " << xx_leg2[j] << "  " << yy_leg2[j] << endl;
		} // End of rearrangement

		for(Int_t iter_x = 0; iter_x < nEtaBins_leg2; iter_x++)
		{
			Eff_HLT_Leg2_data[iter_x][iter_y] = yy_leg2[iter_x]; // actually, it is the scale factor.
			Eff_HLT_Leg2_MC[iter_x][iter_y] = 1;
			//cout << "ID: bin# = [" << iter_x << ", " << iter_y << "]" << " eta = " << xx_leg2[iter_x] << " sf = " << yy_leg2[iter_x] << endl;
		}
	}

	cout << "Setting for efficiency correction factors is completed" << endl;
}
/*
Double_t DYAnalyzer::EfficiencySF_EventWeight_HLTv4p2(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];


	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLTv4p2_data[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_data = Eff_HLTv4p2_data[etabin2][ptbin2];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLTv4p2_MC[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_MC = Eff_HLTv4p2_MC[etabin2][ptbin2];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

	// cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin1][ptbin1], Eff_Iso_data[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin2][ptbin2], Eff_Iso_data[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin1][ptbin1], Eff_Iso_MC[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin2][ptbin2], Eff_Iso_MC[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLTv4p3(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];


	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLTv4p3_data[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_data = Eff_HLTv4p3_data[etabin2][ptbin2];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLTv4p3_MC[etabin1][ptbin1];
	Double_t Eff_Trig_muon2_MC = Eff_HLTv4p3_MC[etabin2][ptbin2];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

	// cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin1][ptbin1], Eff_Iso_data[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin2][ptbin2], Eff_Iso_data[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin1][ptbin1], Eff_Iso_MC[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin2][ptbin2], Eff_Iso_MC[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}
*/
Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 Reco)"); return -999; }
	if( ptbin1_ID == 9999 || etabin1_ID == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 ID)"); return -999; }
	if( ptbin1_Iso == 9999 || etabin1_Iso == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 Iso)"); return -999; }
	if( ptbin1_Trig == 9999 || etabin1_Trig == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 Trig)"); return -999; }
	if( ptbin2_Reco == 9999 || etabin2_Reco == 9999 ) { printf("ERROR! Wrong assigned bin number (mu2 Reco)"); return -999; }
	if( ptbin2_ID == 9999 || etabin2_ID == 9999) { printf("ERROR! Wrong assigned bin number (mu2 ID)"); return -999; }
	if( ptbin2_Iso == 9999 || etabin2_Iso == 9999 ) { printf("ERROR! Wrong assigned bin number (mu2 Iso)"); return -999; }
	if( ptbin2_Trig == 9999 || etabin2_Trig == 9999 ) { printf("ERROR! Wrong assigned bin number (mu2 Trig)"); return -999; }

	Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 Reco)"); return -999; }
	if( ptbin1_ID == 9999 || etabin1_ID == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 ID)"); return -999; }
	if( ptbin1_Iso == 9999 || etabin1_Iso == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 Iso)"); return -999; }
	if( ptbin1_Trig == 9999 || etabin1_Trig == 9999 ) { printf("ERROR! Wrong assigned bin number (mu1 Trig)"); return -999; }
	if( ptbin2_Reco == 9999 || etabin2_Reco == 9999 ) { printf("ERROR! Wrong assigned bin number (mu2 Reco)"); return -999; }
	if( ptbin2_ID == 9999 || etabin2_ID == 9999) { printf("ERROR! Wrong assigned bin number (mu2 ID)"); return -999; }
	if( ptbin2_Iso == 9999 || etabin2_Iso == 9999 ) { printf("ERROR! Wrong assigned bin number (mu2 Iso)"); return -999; }
	if( ptbin2_Trig == 9999 || etabin2_Trig == 9999 ) { printf("ERROR! Wrong assigned bin number (mu2 Trig)"); return -999; }

	Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

// -- Test SFs separately -- //
Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF_Reco(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	//Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	//Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco];
	Double_t Eff_MC_all = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco];

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF_RecoID(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	//Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	//Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID];
	Double_t Eff_MC_all = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID];

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF_RecoIDIso(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	Double_t Eff_muon2_data = Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_BtoF[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	//Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	//Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso] * Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_MC_all = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso] * Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso];

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_data_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_BtoF[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_BtoF[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_BtoF[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH_Reco(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	//Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	//Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco];
	Double_t Eff_MC_all = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco];

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH_RecoID(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	//Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	//Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID];
	Double_t Eff_MC_all = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID];

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH_RecoIDIso(Muon mu1, Muon mu2)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2_Reco = Find_muon_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_muon_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_muon_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_muon_EtaBin_ID( eta2 );

	Int_t ptbin2_Iso = Find_muon_PtBin_Iso( Pt2 );
	Int_t etabin2_Iso = Find_muon_EtaBin_Iso( eta2 );

	Int_t ptbin2_Trig = Find_muon_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_muon_EtaBin_Trig( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	Double_t Eff_muon2_data = Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_muon2_MC = Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];

	// -- Trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_data = Eff_HLT_data_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
	Double_t Eff_Trig_muon2_MC = Eff_HLT_MC_GtoH[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;

	//Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	//Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;
	Double_t Eff_data_all = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso] * Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso];
	Double_t Eff_MC_all = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso] * Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso];

	//cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);

		printf("[Data]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_data_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_data_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_data_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco], Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID],
										Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso]);
		printf("\t[Muon2] (Reco, ID, Iso): (%.3lf, %.3lf, %.3lf)\n", Eff_Reco_MC_GtoH[etabin2_Reco][ptbin2_Reco], Eff_ID_MC_GtoH[etabin2_ID][ptbin2_ID],
										Eff_Iso_MC_GtoH[etabin2_Iso][ptbin2_Iso]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_electron(Electron ele1, Electron ele2)
{
	Double_t weight = -999;

	// -- Electron1 -- //
	Double_t Pt1 = ele1.Pt;
	//Double_t eta1 = ele1.eta;
	Double_t eta1 = ele1.etaSC;

	Int_t ptbin1_Reco = Find_electron_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_electron_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_electron_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_electron_EtaBin_ID( eta1 );

	Int_t ptbin1_Trig = Find_electron_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_electron_EtaBin_Trig( eta1 );

	Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco] * Eff_ID_data[etabin1_ID][ptbin1_ID];
	Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC[etabin1_ID][ptbin1_ID];

	// -- Electron2 -- //
	Double_t Pt2 = ele2.Pt;
	//Double_t eta2 = ele2.eta;
	Double_t eta2 = ele2.etaSC;

	Int_t ptbin2_Reco = Find_electron_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_electron_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_electron_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_electron_EtaBin_ID( eta2 );

	Int_t ptbin2_Trig = Find_electron_PtBin_Trig( Pt2 );
	Int_t etabin2_Trig = Find_electron_EtaBin_Trig( eta2 );

	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

	// -- This is trigger part -- //
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Eff_EventTrig_data = Eff_HLT_Leg2_data[etabin1_Trig][ptbin1_Trig] * Eff_HLT_Leg2_data[etabin2_Trig][ptbin2_Trig];
	Eff_EventTrig_MC = Eff_HLT_Leg2_MC[etabin1_Trig][ptbin1_Trig] * Eff_HLT_Leg2_MC[etabin2_Trig][ptbin2_Trig];

	//cout << Eff_EventTrig_data << "\t" << Eff_EventTrig_MC << endl;

	Double_t Eff_data_all = Eff_ele1_data * Eff_ele2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_electron_Reco(Electron ele1, Electron ele2)
{
	Double_t weight = -999;

	// -- Electron1 -- //
	Double_t Pt1 = ele1.Pt;
	//Double_t eta1 = ele1.eta;
	Double_t eta1 = ele1.etaSC;

	Int_t ptbin1_Reco = Find_electron_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_electron_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_electron_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_electron_EtaBin_ID( eta1 );

	//Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco] * Eff_ID_data[etabin1_ID][ptbin1_ID];
	//Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC[etabin1_ID][ptbin1_ID];
	Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco]; // only reco SF
	Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1_Reco][ptbin1_Reco]; // only reco SF

	// -- Electron2 -- //
	Double_t Pt2 = ele2.Pt;
	//Double_t eta2 = ele2.eta;
	Double_t eta2 = ele2.etaSC;

	Int_t ptbin2_Reco = Find_electron_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_electron_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_electron_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_electron_EtaBin_ID( eta2 );

	//Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	//Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco]; // only reco SF
	Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco]; // only reco SF

	Double_t Eff_data_all = Eff_ele1_data * Eff_ele2_data;
	Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_electron_RecoID(Electron ele1, Electron ele2)
{
	Double_t weight = -999;

	// -- Electron1 -- //
	Double_t Pt1 = ele1.Pt;
	//Double_t eta1 = ele1.eta;
	Double_t eta1 = ele1.etaSC;

	Int_t ptbin1_Reco = Find_electron_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_electron_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_electron_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_electron_EtaBin_ID( eta1 );

	Double_t Eff_ele1_data = Eff_Reco_data[etabin1_Reco][ptbin1_Reco] * Eff_ID_data[etabin1_ID][ptbin1_ID];
	Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC[etabin1_ID][ptbin1_ID];

	// -- Electron2 -- //
	Double_t Pt2 = ele2.Pt;
	//Double_t eta2 = ele2.eta;
	Double_t eta2 = ele2.etaSC;

	Int_t ptbin2_Reco = Find_electron_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_electron_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_electron_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_electron_EtaBin_ID( eta2 );

	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

	Double_t Eff_data_all = Eff_ele1_data * Eff_ele2_data;
	Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("(pt1, eta1, pt2, eta2): (%.3lf, %.3lf, %.3lf, %.3lf)\n", Pt1, eta1, Pt2, eta2);
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_BtoF(Muon mu, Electron ele)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu.Pt;
	Double_t eta1 = mu.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	// -- Electron2 -- //
	Double_t Pt2 = ele.Pt;
	//Double_t eta2 = ele.eta;
	Double_t eta2 = ele.etaSC;

	Int_t ptbin2_Reco = Find_electron_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_electron_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_electron_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_electron_EtaBin_ID( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	//Muon1
	Double_t Eff_muon1_data = Eff_Reco_data_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_data_BtoF[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_BtoF[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_BtoF[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_BtoF[etabin1_Iso][ptbin1_Iso];

	//Electron2
	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

	//Trigger : We consider only SingleMuon trigger
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_BtoF[etabin1_Trig][ptbin1_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_BtoF[etabin1_Trig][ptbin1_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC;


	// -- emu event SF -- //
	Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 ) printf("[SF] Weight = %.3lf\n", weight);
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_emu_GtoH(Muon mu, Electron ele)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu.Pt;
	Double_t eta1 = mu.eta;

	Int_t ptbin1_Reco = Find_muon_PtBin_Reco( Pt1 );
	Int_t etabin1_Reco = Find_muon_EtaBin_Reco( eta1 );

	Int_t ptbin1_ID = Find_muon_PtBin_ID( Pt1 );
	Int_t etabin1_ID = Find_muon_EtaBin_ID( eta1 );

	Int_t ptbin1_Iso = Find_muon_PtBin_Iso( Pt1 );
	Int_t etabin1_Iso = Find_muon_EtaBin_Iso( eta1 );

	Int_t ptbin1_Trig = Find_muon_PtBin_Trig( Pt1 );
	Int_t etabin1_Trig = Find_muon_EtaBin_Trig( eta1 );

	// -- Electron2 -- //
	Double_t Pt2 = ele.Pt;
	//Double_t eta2 = ele.eta;
	Double_t eta2 = ele.etaSC;

	Int_t ptbin2_Reco = Find_electron_PtBin_Reco( Pt2 );
	Int_t etabin2_Reco = Find_electron_EtaBin_Reco( eta2 );

	Int_t ptbin2_ID = Find_electron_PtBin_ID( Pt2 );
	Int_t etabin2_ID = Find_electron_EtaBin_ID( eta2 );

	// -- Check about bin settings -- //
	if( ptbin1_Reco == 9999 || etabin1_Reco == 9999 || ptbin1_ID == 9999 || etabin1_ID == 9999 || ptbin1_Iso == 9999 || etabin1_Iso == 9999 || ptbin1_Trig == 9999 || etabin1_Trig == 9999 || ptbin2_Reco == 9999 || etabin2_Reco == 9999 || ptbin2_ID == 9999 || etabin2_ID == 9999)
	{
		printf("ERROR! Wrong assigned bin number ...");
		return -999;
	}

	//Muon1
	Double_t Eff_muon1_data = Eff_Reco_data_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_data_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_data_GtoH[etabin1_Iso][ptbin1_Iso];
	Double_t Eff_muon1_MC = Eff_Reco_MC_GtoH[etabin1_Reco][ptbin1_Reco] * Eff_ID_MC_GtoH[etabin1_ID][ptbin1_ID] * Eff_Iso_MC_GtoH[etabin1_Iso][ptbin1_Iso];

	//Electron2
	Double_t Eff_ele2_data = Eff_Reco_data[etabin2_Reco][ptbin2_Reco] * Eff_ID_data[etabin2_ID][ptbin2_ID];
	Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2_Reco][ptbin2_Reco] * Eff_ID_MC[etabin2_ID][ptbin2_ID];

	//Trigger : We consider only SingleMuon trigger
	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;

	Double_t Eff_Trig_muon1_data = Eff_HLT_data_GtoH[etabin1_Trig][ptbin1_Trig];
	Eff_EventTrig_data = Eff_Trig_muon1_data;

	Double_t Eff_Trig_muon1_MC = Eff_HLT_MC_GtoH[etabin1_Trig][ptbin1_Trig];
	Eff_EventTrig_MC = Eff_Trig_muon1_MC;


	// -- emu event SF -- //
	Double_t Eff_data_all = Eff_muon1_data * Eff_ele2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 ) printf("[SF] Weight = %.3lf\n", weight);
	return weight;
}

Int_t DYAnalyzer::Find_muon_PtBin_Reco(Double_t Pt)
{
	const Int_t nPtBins = 1;
	Double_t PtBinEdges[nPtBins+1] = {25, 500};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_muon_PtBin_ID(Double_t Pt)
{
	// -- HighPtID & TrkIso -- //
	//const Int_t nPtBins = 7;
	//Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 55, 60, 120};

	// -- TightID & PFIso -- //
	const Int_t nPtBins = 6;
	Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 60, 120};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_muon_PtBin_Iso(Double_t Pt)
{
	// -- HighPtID & TrkIso -- //
	//const Int_t nPtBins = 7;
	//Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 55, 60, 120};

	// -- TightID & PFIso -- //
	const Int_t nPtBins = 6;
	Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 60, 120};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_muon_PtBin_Trig(Double_t Pt)
{
	// -- IsoMu24_OR_IsoTkMu24 -- //
	const Int_t nPtBins = 7;
	Double_t PtBinEdges[nPtBins+1] = {26, 30, 40, 50, 60, 120, 200, 500};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_muon_EtaBin_Reco(Double_t eta)
{
	const Int_t nEtaBins = 15;
	Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		//if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::Find_muon_EtaBin_ID(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 4;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		//if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::Find_muon_EtaBin_Iso(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 4;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		//if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::Find_muon_EtaBin_Trig(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 4;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		//if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::Find_electron_PtBin_Reco(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 1;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {25, 500};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_electron_PtBin_ID(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 5;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {10, 20, 35, 50, 90, 150};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_electron_PtBin_Trig(Double_t Pt)
{
	//const Int_t nPtBins = 4;
	const Int_t nPtBins = 8;
	//Double_t PtBinEdges[nPtBins+1] = {10, 22, 40, 70, 250};
	Double_t PtBinEdges[nPtBins+1] = {10, 15, 20, 25, 30, 50, 90, 150, 500};

	Int_t ptbin = 9999;

	// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
	if( Pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	// -- if Pt is smaller than the smallest Pt bin edge, SF is same with the value for the first bin -- // updated at 14 Apr. 2017 by Dalmin Pai
	else if( Pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( Pt > PtBinEdges[i] && Pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_electron_EtaBin_Reco(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 30;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {-2.5,-2.45,-2.4,-2.3,-2.2,-2.0,-1.8,-1.63,-1.566,-1.4442,-1.2,-1.0,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,1.0,1.2,1.4442,1.566,1.63,1.8,2.0,2.2,2.3,2.4,2.45,2.5};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		//if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::Find_electron_EtaBin_ID(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 10;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		//if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

Int_t DYAnalyzer::Find_electron_EtaBin_Trig(Double_t eta)
{
	//const Int_t nEtaBins = 5;
	const Int_t nEtaBins = 10;
	//Double_t EtaBinEdges[nEtaBins+1] = {-2.4, -1.2, -0.3, 0.3, 1.2, 2.4};
	Double_t EtaBinEdges[nEtaBins+1] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( eta > EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
		//if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}

/*
Double_t DYAnalyzer::EfficiencySF_EventWeight(Muon mu1, Muon mu2, NtupleHandle *ntuple)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];

	Bool_t isHLTv4p2 = kFALSE;
	if( ntuple->runNum < 257932.5 )
		isHLTv4p2 = kTRUE;

	Double_t Eff_EventTrig_data = 0;
	Double_t Eff_EventTrig_MC = 0;
	if( isHLTv4p2 )
	{
		Double_t Eff_Trig_muon1_data = Eff_HLTv4p2_data[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_data = Eff_HLTv4p2_data[etabin2][ptbin2];
		Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

		Double_t Eff_Trig_muon1_MC = Eff_HLTv4p2_MC[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_MC = Eff_HLTv4p2_MC[etabin2][ptbin2];
		Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;
	}
	else
	{
		Double_t Eff_Trig_muon1_data = Eff_HLTv4p3_data[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_data = Eff_HLTv4p3_data[etabin2][ptbin2];
		Eff_EventTrig_data = Eff_Trig_muon1_data + Eff_Trig_muon2_data - Eff_Trig_muon1_data * Eff_Trig_muon2_data;

		Double_t Eff_Trig_muon1_MC = Eff_HLTv4p3_MC[etabin1][ptbin1];
		Double_t Eff_Trig_muon2_MC = Eff_HLTv4p3_MC[etabin2][ptbin2];
		Eff_EventTrig_MC = Eff_Trig_muon1_MC + Eff_Trig_muon2_MC - Eff_Trig_muon1_MC * Eff_Trig_muon2_MC;
	}

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data * Eff_EventTrig_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC * Eff_EventTrig_MC;

	// cout << "Eff_data_all: " << Eff_data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
	weight = Eff_data_all / Eff_MC_all;

	if( weight > 2 )
	{
		printf("[Data]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin1][ptbin1], Eff_Iso_data[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_data[etabin2][ptbin2], Eff_Iso_data[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_data, Eff_data_all);

		printf("[MC]\n");
		printf("\t[Muon1] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin1][ptbin1], Eff_Iso_MC[etabin1][ptbin1]);
		printf("\t[Muon2] (RecoID, Iso): (%.3lf, %.3lf)\n", Eff_RecoID_MC[etabin2][ptbin2], Eff_Iso_MC[etabin2][ptbin2]);
		printf("\t[Event] (TrigEvent, Total): (%.3lf, %.3lf)\n", Eff_EventTrig_MC, Eff_MC_all);

		printf("(ptbin1, etabin1, ptbin2, etabin2): (%d, %d, %d, %d)\n", ptbin1, etabin1, ptbin2, etabin2);
		
		printf("[SF] Weight = %.3lf\n", weight);
	}
	return weight;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_RecoIdIso(Muon mu1, Muon mu2, NtupleHandle *ntuple)
{
	Double_t weight = -999;

	// -- Muon1 -- //
	Double_t Pt1 = mu1.Pt;
	Double_t eta1 = mu1.eta;

	Int_t ptbin1 = FindPtBin( Pt1 );
	Int_t etabin1 = FindEtaBin( eta1 );

	if( ptbin1 == 9999 || etabin1 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin1, etabin1);
		return -999;
	}

	Double_t Eff_muon1_data = Eff_RecoID_data[etabin1][ptbin1] * Eff_Iso_data[etabin1][ptbin1];
	Double_t Eff_muon1_MC = Eff_RecoID_MC[etabin1][ptbin1] * Eff_Iso_MC[etabin1][ptbin1];

	// -- Muon2 -- //
	Double_t Pt2 = mu2.Pt;
	Double_t eta2 = mu2.eta;

	Int_t ptbin2 = FindPtBin( Pt2 );
	Int_t etabin2 = FindEtaBin( eta2 );

	if( ptbin2 == 9999 || etabin2 == 9999 )
	{
		printf("ERROR! Wrong assigned bin number ... (ptbin, etabin) = (%d, %d)\n", ptbin2, etabin2);
		return -999;
	}
	Double_t Eff_muon2_data = Eff_RecoID_data[etabin2][ptbin2] * Eff_Iso_data[etabin2][ptbin2];
	Double_t Eff_muon2_MC = Eff_RecoID_MC[etabin2][ptbin2] * Eff_Iso_MC[etabin2][ptbin2];

	Double_t Eff_data_all = Eff_muon1_data * Eff_muon2_data;
	Double_t Eff_MC_all = Eff_muon1_MC * Eff_muon2_MC;

	weight = Eff_data_all / Eff_MC_all;

	return weight;
}
*/
Bool_t DYAnalyzer::EventSelection(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15) //2018.03.02
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		//if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") ) //2018.03.02
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE ) //2018.03.02
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				//if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") ) //2018.03.02
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				//if( reco_M > 60 && reco_M < 120 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE ) //2018.03.02
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zpeak(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15) //2018.03.02
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		//if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") ) //2018.03.02
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			//if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE ) //2018.03.02
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				//if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") ) //2018.03.02
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

				//if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				if( reco_M > 60 && reco_M < 120 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE ) //2018.03.02
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}

// -- Test using the trigger without isolation condition: HLT_Mu50_v* -- //
Bool_t DYAnalyzer::EventSelection_Mu50(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, HLT) )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, HLT) )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}


Bool_t DYAnalyzer::EventSelection_minusDimuonVtxCut(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			// if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				// if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 )
				if( reco_M > 10 && Angle < TMath::Pi() - 0.005 )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}

// -- Event selection used for differential Z cross section measurement @ 13TeV -- // For HighPt muon id with 60-120GeV mass bin
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			// -- Opposite sign condition -- //
			Bool_t isOppositeSign = kFALSE;
			if( recolep1.charge != recolep2.charge )
				isOppositeSign = kTRUE;

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			// -- Dimuon Vtx cut -- //
			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

//			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
//			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//			if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && isOppositeSign == kTRUE )
			if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
			Double_t SumPt_BestPair = 0;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
				{
					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t SumPt_temp = 0;
								SumPt_temp = Mu.Pt + Mu_jth.Pt;

								// -- Find best pair by selecting largest Pt sum -- // 
								if( SumPt_BestPair < SumPt_temp )
								{
									SumPt_BestPair = SumPt_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( 0 < SumPt_BestPair ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				// -- Opposite sign condition -- //
				Bool_t isOppositeSign = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge )
					isOppositeSign = kTRUE;

				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- Dimuon Vtx cut -- //
				Double_t VtxProb = -999;
				Double_t VtxNormChi2 = 999;
				DimuonVertexProbNormChi2(ntuple, mu1_BestPair.Inner_pT, mu2_BestPair.Inner_pT, &VtxProb, &VtxNormChi2);

//				TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
//				TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

//				if( reco_M > 60 && reco_M < 120 && isOppositeSign == kTRUE )
				if( reco_M > 60 && reco_M < 120 && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}

// -- Event selection used for differential Z cross section measurement @ 13TeV -- // For all mass bin, and HighPt id
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		//if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
		//if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

// -- Event selection used for N-1 cuts -- // For all mass bin, and HighPt id
Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt1(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_isGLB() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_isGLB() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt2(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_muonHits() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_muonHits() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt3(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_nMatches() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_nMatches() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt4(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_dpT_over_pT() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_isPF() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt5(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_dxyVTX() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_dxyVTX() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt6(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_dzVTX() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt7(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_pixelHits() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_pixelHits() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt8(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon_minus_trackerLayers() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_trackerLayers() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt9(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt10(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt11(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon() )
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Zdiff_13TeV_HighPt12(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
							vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    //if( MuonCollection[j].isTightMuon() && MuonCollection[j].trkiso < 0.10)
	    if( MuonCollection[j].isTightMuon_minus_chi2dof() && MuonCollection[j].RelPFIso_dBeta < 0.15)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	if( nQMuons == 2)
	{
		Muon recolep1 = QMuonCollection[0];
		Muon recolep2 = QMuonCollection[1];

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( recolep1.charge != recolep2.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 );
			SelectedMuonCollection->push_back( recolep2 );
		}
	}
	else if( nQMuons > 2 )
	{
		// -- More then 2 Qualified Muon: Select the muons with highest pT -- // 
		Double_t Pt_leading = 0;
		Muon LeadingMuon;
		Double_t i_leading = 0;
		for(Int_t i_mu1 = 0; i_mu1 < nQMuons; i_mu1++)
		{
			Muon Mu = QMuonCollection[i_mu1];

			// printf("%dth Muon: Pt = %.3lf\n", i_mu1, Mu.Pt);

			if( Mu.Pt > Pt_leading )
			{
				Pt_leading = Mu.Pt;
				LeadingMuon	= Mu;
				i_leading = i_mu1;
			}
		}

		Double_t Pt_sub = 0;
		Muon SubMuon;
		for(Int_t i_mu2=0; i_mu2 < nQMuons; i_mu2++)
		{
			if( i_mu2 == i_leading ) continue;

			Muon Mu = QMuonCollection[i_mu2];

			if( Mu.Pt > Pt_sub )
			{
				Pt_sub = Mu.Pt;
				SubMuon	= Mu;
			}
		}

		// printf("\t(Pt_leading, Pt_sub) = (%.3lf, %.3lf)\n", LeadingMuon.Pt, SubMuon.Pt);

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(LeadingMuon, SubMuon);

		// -- Opposite sign condition -- //
		Bool_t isOppositeSign = kFALSE;
		if( LeadingMuon.charge != SubMuon.charge )
			isOppositeSign = kTRUE;

		Double_t reco_M = (LeadingMuon.Momentum + SubMuon.Momentum).M();

		// -- Dimuon Vtx cut -- //
		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, LeadingMuon.Inner_pT, SubMuon.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = LeadingMuon.Momentum_Inner;
		TLorentzVector inner_v2 = SubMuon.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = LeadingMuon.Momentum.Angle( SubMuon.Momentum.Vect() );

//		if( isPassAcc == kTRUE && isOppositeSign == kTRUE )
//		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOppositeSign == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( LeadingMuon );
			SelectedMuonCollection->push_back( SubMuon );
		}

	} // -- End of else if( nQMuons > 2 ) -- //

	return isPassEventSelection;
}



// test for emu event selection not using vtx cut
// Updated to use synchronized acceptance cut : 19 Jan. 2018
Bool_t DYAnalyzer::EventSelection_emu_method_test(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
						vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection)
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
		if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
			//&& MuonCollection[j].Pt > LeadPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut )
			&& MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut ) // pT>17 && |eta|<2.4
			QMuonCollection.push_back( MuonCollection[j] );
	}

	//Collect qualified electrons among electrons
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE
			//&& elec.Pt > SubPtCut && fabs(elec.etaSC) < 2.5 && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) ) // pT>17 && |eta|<2.4
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	Double_t Pt_mu = 0;
	Double_t Pt_el = 0;
	Muon mu_BestPair;
	Electron el_BestPair;

	// -- Select muon with highest pT -- //
	for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
	{
		Muon Mu = QMuonCollection[i_mu];

		// -- muon should be matched with HLT objects in emu best pair -- //
		if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
		{
			if( Mu.Pt > Pt_mu )
			{
				Pt_mu		= Mu.Pt;
				mu_BestPair	= Mu;
			}
		}
	}

	// -- Select electron with highest pT -- //
	for(Int_t j_el=0; j_el<nQElectrons; j_el++)
	{
		Electron El = QElectronCollection[j_el];

		if( El.Pt > Pt_el )
		{
			Pt_el		= El.Pt;
			el_BestPair	= El;
		}
	}

	//if( Pt_mu > 0 && Pt_el > 0 )
	if( Pt_mu > 0 && Pt_el > 0 && ( Pt_mu > LeadPtCut || Pt_el > LeadPtCut ) ) // At least one lepton has pT above 28 [GeV]
	{
		TLorentzVector reco_v1 = mu_BestPair.Momentum;
		TLorentzVector reco_v2 = el_BestPair.Momentum;
		Double_t reco_M = (reco_v1 + reco_v2).M();

		// -- 3D open angle -- //
		Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

		//if( reco_M > 10 && Angle < TMath::Pi() - 0.005 )
		if( reco_M > 10 )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( mu_BestPair );
			SelectedElectronCollection->push_back( el_BestPair );
		}
	}

	return isPassEventSelection;
}


Bool_t DYAnalyzer::EventSelection_emu_method(vector< Muon > MuonCollection, vector< Electron > ElectronCollection, NtupleHandle *ntuple,
						vector< Muon >* SelectedMuonCollection, vector< Electron >* SelectedElectronCollection)
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
		//if( MuonCollection[j].isHighPtMuon() && MuonCollection[j].trkiso < 0.10
		if( MuonCollection[j].isTightMuon() && MuonCollection[j].RelPFIso_dBeta < 0.15 //2018.03.02
			&& MuonCollection[j].Pt > SubPtCut && fabs(MuonCollection[j].eta) < LeadEtaCut ) // pT>17 && |eta|<2.4
			QMuonCollection.push_back( MuonCollection[j] );
	}

	//Collect qualified electrons among electrons
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < LeadEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) ) // pT>17 && |eta|<2.4
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQMuons = (Int_t)QMuonCollection.size();
	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	Double_t VtxProb_BestPair = -1;
	Double_t VtxNormChi2_BestPair = 999;
	Muon mu_BestPair;
	Electron el_BestPair;

	for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
	{
		Muon Mu = QMuonCollection[i_mu];

		// -- muon should be matched with HLT objects in emu best pair -- //
		if( Mu.isTrigMatched(ntuple, "HLT_IsoMu24_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu24_v*") )
		{
			// -- Start another loop for finding electron (for electron, we don't need to check about trigger) -- //
			for(Int_t j_el=0; j_el<nQElectrons; j_el++)
			{
				Electron El = QElectronCollection[j_el];

				Double_t VtxProb_temp = -999;
				Double_t VtxNormChi2_temp = 999;
				emuVertexProbNormChi2(ntuple, El.gsfpT, Mu.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

				// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
				if( VtxNormChi2_temp < VtxNormChi2_BestPair )
				{
					VtxNormChi2_BestPair = VtxNormChi2_temp;
					mu_BestPair = Mu;
					el_BestPair = El;
				}
			} // -- end of the loop for j_el (finding for electron)
		}
	} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

	if( VtxNormChi2_BestPair < 999 )
	{
		TLorentzVector reco_v1 = mu_BestPair.Momentum;
		TLorentzVector reco_v2 = el_BestPair.Momentum;
		Double_t reco_M = (reco_v1 + reco_v2).M();

		Bool_t isPassAcc = kFALSE;
		if( mu_BestPair.Pt > LeadPtCut || el_BestPair.Pt > LeadPtCut ) isPassAcc = kTRUE;

		// -- 3D open angle -- //
		Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

		//if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 )
		if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( mu_BestPair );
			SelectedElectronCollection->push_back( el_BestPair );
		}
	}

	return isPassEventSelection;
}

Bool_t DYAnalyzer::isPassAccCondition_Muon(Muon Mu1, Muon Mu2)
{
	Bool_t isPassAcc = kFALSE;
	Muon leadMu, subMu;
	CompareMuon(&Mu1, &Mu2, &leadMu, &subMu);
	if( leadMu.Pt > LeadPtCut && fabs(leadMu.eta) < LeadEtaCut && 
		subMu.Pt  > SubPtCut  && fabs(subMu.eta)  < SubEtaCut )
		isPassAcc = kTRUE;

	return isPassAcc;
}

Bool_t DYAnalyzer::isPassAccCondition_GenLepton(GenLepton genlep1, GenLepton genlep2)
{
	Bool_t isPassAcc = kFALSE;

	GenLepton leadGenLep, subGenLep;
	CompareGenLepton(&genlep1, &genlep2, &leadGenLep, &subGenLep);
	
	if( leadGenLep.Pt > LeadPtCut && fabs(leadGenLep.eta) < LeadEtaCut &&
		subGenLep.Pt  > SubPtCut  && fabs(subGenLep.eta) < SubEtaCut )
		isPassAcc = 1;

	return isPassAcc;
}

void DYAnalyzer::CompareMuon(Muon *Mu1, Muon *Mu2, Muon *leadMu, Muon *subMu)
{
    if( Mu1->Pt > Mu2->Pt )
    {
        *leadMu = *Mu1;
        *subMu = *Mu2;
    }
    else
    {
        *leadMu = *Mu2;
        *subMu = *Mu1;
    }
}

void DYAnalyzer::CompareGenLepton(GenLepton *genlep1, GenLepton *genlep2, GenLepton *leadgenlep, GenLepton *subgenlep)
{
	if( genlep1->Pt > genlep2->Pt )
	{
		*leadgenlep = *genlep1;
		*subgenlep = *genlep2;
	}
	else
	{
		*leadgenlep = *genlep2;
		*subgenlep = *genlep1;
	}
}

void DYAnalyzer::DimuonVertexProbNormChi2(NtupleHandle *ntuple, Double_t Pt1, Double_t Pt2, Double_t *VtxProb, Double_t *VtxNormChi2)
{
	vector<double> *PtCollection1 = ntuple->vtxTrkCkt1Pt;
	vector<double> *PtCollection2 = ntuple->vtxTrkCkt2Pt;
	vector<double> *VtxProbCollection = ntuple->vtxTrkProb;

	Int_t NPt1 = (Int_t)PtCollection1->size();
	Int_t NPt2 = (Int_t)PtCollection2->size();
	Int_t NProb = (Int_t)VtxProbCollection->size();

	if( NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb ) 
		cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

	// cout << "Pt1: " << Pt1 << " Pt2: " << Pt2 << endl;

	for(Int_t i=0; i<NProb; i++)
	{
		// cout << "\tPtCollection1->at("<< i << "): " << PtCollection1->at(i) << " PtCollection2->at("<< i << "): " << PtCollection2->at(i) << endl;
		if( ( PtCollection1->at(i) == Pt1 && PtCollection2->at(i) == Pt2 )  || ( PtCollection1->at(i) == Pt2 && PtCollection2->at(i) == Pt1 ) )
		{
			*VtxProb = VtxProbCollection->at(i);
			*VtxNormChi2 = ntuple->vtxTrkChi2->at(i) / ntuple->vtxTrkNdof->at(i);
			break;
		}
	}

	return;
}

void DYAnalyzer::emuVertexProbNormChi2(NtupleHandle *ntuple, Double_t ele_Pt, Double_t mu_Pt, Double_t *VtxProb, Double_t *VtxNormChi2)
{
	vector<double> *PtCollection1 = ntuple->vtxTrkEMu1Pt; //electron
	vector<double> *PtCollection2 = ntuple->vtxTrkEMu2Pt; //muon
	vector<double> *VtxProbCollection = ntuple->vtxTrkEMuProb;

	Int_t NPt1 = (Int_t)PtCollection1->size();
	Int_t NPt2 = (Int_t)PtCollection2->size();
	Int_t NProb = (Int_t)VtxProbCollection->size();

	if( NPt1 != NPt2 || NPt2 != NProb || NPt1 != NProb ) 
		cout << "NPt1: " << NPt1 << " NPt2: " << NPt2 << " Nprob: " << NProb << endl;

	for(Int_t i=0; i<NProb; i++)
	{
		if( PtCollection1->at(i) == ele_Pt && PtCollection2->at(i) == mu_Pt )
		{
			*VtxProb = VtxProbCollection->at(i);
			*VtxNormChi2 = ntuple->vtxTrkEMuChi2->at(i) / ntuple->vtxTrkEMuNdof->at(i);
			break;
		}
	}

	return;
}

// -- Event selecton for the electron channel (test) -- //
Bool_t DYAnalyzer::EventSelection_Electron(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		// cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
		if( elec.isMediumElectron_Spring25ns() && elec.ecalDriven == 1 && elec.Pt > 15 )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

// -- Event selecton for the electron channel (2016.02.11) -- // modified at 17 May 2017 by Dalmin Pai
Bool_t DYAnalyzer::EventSelection_ElectronChannel(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE // modified by Dalmin Pai
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		//if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel_Zpeak(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // input: All electrons in a event & NtupleHandle
						vector< Electron >* SelectedElectronCollection) // output: 2 electrons passing event selection conditions
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE // modified by Dalmin Pai
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

// -- Event selection for without energy scale correction -- //
Bool_t DYAnalyzer::EventSelection_ElectronChannel0(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.passMediumID == kTRUE // modified by Dalmin Pai
			&& elec.pTUnCorr > SubPtCut && fabs(elec.etaSCUnCorr) < SubEtaCut && !( fabs(elec.etaSCUnCorr) > 1.4442 && fabs(elec.etaSCUnCorr) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;

		Electron leadElec, subElec;
		if( recolep1.pTUnCorr > recolep2.pTUnCorr )
		{
			leadElec = recolep1;
			subElec = recolep2;
		}
		else
		{
			leadElec = recolep2;
			subElec = recolep1;
		}

		if( leadElec.pTUnCorr > LeadPtCut && fabs(leadElec.etaSCUnCorr) < LeadEtaCut &&
			!( fabs(leadElec.etaSCUnCorr) > 1.4442 && fabs(leadElec.etaSCUnCorr) < 1.566 ) &&
			subElec.pTUnCorr > SubPtCut && fabs(subElec.etaSCUnCorr) < SubEtaCut &&
			!( fabs(subElec.etaSCUnCorr) > 1.4442 && fabs(subElec.etaSCUnCorr) < 1.566 ) )
			isPassAcc = kTRUE;

		Double_t reco_M = (recolep1.Momentum_UnCorr + recolep2.Momentum_UnCorr).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

// -- N-1 cuts Event selecton for the electron channel -- // modified at 17 May 2017 by Dalmin Pai //
Bool_t DYAnalyzer::EventSelection_ElectronChannel1(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_Full5x5_SigmaIEtaIEta() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel2(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_dEtaInSeed() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel3(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_dPhiIn() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel4(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_HoverE() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel5(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_RelPFIso_Rho() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel6(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_InvEminusInvP() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel7(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_mHits() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::EventSelection_ElectronChannel8(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		if( elec.isMediumElectron_2016dataFor80X_minus_passConvVeto() == kTRUE
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		//if( reco_M > 10 && isPassAcc == kTRUE )
		if( reco_M > 60 && reco_M < 120 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}
// End of N-1 cuts Event Selection

// -- Event selecton for the electron channel (2016.02.11) -- //
Bool_t DYAnalyzer::EventSelection_ElectronChannel_NminusPFIso(vector< Electron > ElectronCollection, NtupleHandle *ntuple, // -- input: All electrons in a event & NtupleHandle -- //
						vector< Electron >* SelectedElectronCollection) // -- output: 2 electrons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	// -- Electron ID -- //
	vector< Electron > QElectronCollection;
	for(Int_t j=0; j<(int)ElectronCollection.size(); j++)
	{
		Electron elec = ElectronCollection[j];
		// cout << "elec.passConvVeto: " << elec.passConvVeto << endl;
		if( elec.isMediumElectron_Spring25ns_minus_PFIso() && elec.ecalDriven == 1 
			&& elec.Pt > SubPtCut && fabs(elec.etaSC) < SubEtaCut && !( fabs(elec.etaSC) > 1.4442 && fabs(elec.etaSC) < 1.566 ) )
			QElectronCollection.push_back( ElectronCollection[j] );
	}

	Int_t nQElectrons = (Int_t)QElectronCollection.size();

	// cout << "# qualified electrons: " << nQElectrons << endl;

	if( nQElectrons == 2 )
	{
		Electron recolep1 = QElectronCollection[0];
		Electron recolep2 = QElectronCollection[1];

		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Electron(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		if( reco_M > 10 && isPassAcc == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedElectronCollection->push_back( recolep1 );
			SelectedElectronCollection->push_back( recolep2 );
		}
	}
	return isPassEventSelection;

}

Bool_t DYAnalyzer::isPassAccCondition_Electron(Electron Elec1, Electron Elec2)
{
	Bool_t isPassAcc = kFALSE;
	Electron leadElec, subElec;
	CompareElectron(&Elec1, &Elec2, &leadElec, &subElec);
	if( leadElec.Pt > LeadPtCut && fabs(leadElec.etaSC) < LeadEtaCut && !( fabs(leadElec.etaSC) > 1.4442 && fabs(leadElec.etaSC) < 1.566 ) &&
		subElec.Pt  > SubPtCut  && fabs(subElec.etaSC)  < SubEtaCut && !( fabs(subElec.etaSC) > 1.4442 && fabs(subElec.etaSC) < 1.566 ) )
		isPassAcc = kTRUE;

	return isPassAcc;
}


Bool_t DYAnalyzer::isPassAccCondition_GenLepton_ECALGAP(GenLepton genlep1, GenLepton genlep2)
{
	Bool_t isPassAcc = kFALSE;

	GenLepton leadGenLep, subGenLep;
	CompareGenLepton(&genlep1, &genlep2, &leadGenLep, &subGenLep);
	
	if( leadGenLep.Pt > LeadPtCut && fabs(leadGenLep.eta) < LeadEtaCut && !( fabs(leadGenLep.eta) > 1.4442 && fabs(leadGenLep.eta) < 1.566 ) &&
		subGenLep.Pt  > SubPtCut  && fabs(subGenLep.eta) < SubEtaCut && !( fabs(subGenLep.eta) > 1.4442 && fabs(subGenLep.eta) < 1.566 ) )
		isPassAcc = 1;

	return isPassAcc;
}

void DYAnalyzer::CompareElectron(Electron *Elec1, Electron *Elec2, Electron *leadElec, Electron *subElec)
{
    if( Elec1->Pt > Elec2->Pt )
    {
        *leadElec = *Elec1;
        *subElec = *Elec2;
    }
    else
    {
        *leadElec = *Elec2;
        *subElec = *Elec1;
    }
}

void DYAnalyzer::PostToPreFSR_byDressedLepton(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection)
{
	TLorentzVector genlep_Mom_postFSR = genlep_postFSR->Momentum;

	TLorentzVector SumPhotonMom;
	SumPhotonMom.SetPxPyPzE(0,0,0,0);

	Int_t NGenOthers = ntuple->nGenOthers;
	for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
	{
		GenOthers genlep;
		genlep.FillFromNtuple(ntuple, i_gen);

		// -- Only for the photons whose mother is muon or anti-muon -- //
		if( fabs(genlep.ID) == 22 && fabs(genlep.Mother) == 13)
		{
			
			Double_t dR = Calc_dR_GenLepton_GenOthers(*genlep_postFSR, genlep);

			// -- Sum of all photon's momentum near the post-FSR muon -- //
			if( dR < dRCut )
			{
				SumPhotonMom  = SumPhotonMom + genlep.Momentum;
				GenPhotonCollection->push_back( genlep );
			}
		}
	}

	// -- Momentum(pre-FSR) = Momentum(post-FSR) + Sum of all Photon's momentum near the post-FSR muon -- //
	genlep_preFSR->Momentum = genlep_Mom_postFSR + SumPhotonMom;
	genlep_preFSR->Et = genlep_preFSR->Momentum.Et();
	genlep_preFSR->Pt = genlep_preFSR->Momentum.Pt();
	genlep_preFSR->eta = genlep_preFSR->Momentum.Eta();
	genlep_preFSR->phi = genlep_preFSR->Momentum.Phi();
	genlep_preFSR->Px = genlep_preFSR->Momentum.Px();
	genlep_preFSR->Py = genlep_preFSR->Momentum.Py();
	genlep_preFSR->Pz = genlep_preFSR->Momentum.Pz();
}

void DYAnalyzer::PostToPreFSR_byDressedLepton_AllPhotons(NtupleHandle *ntuple, GenLepton *genlep_postFSR, Double_t dRCut, GenLepton *genlep_preFSR, vector< GenOthers >* GenPhotonCollection)
{
	TLorentzVector genlep_Mom_postFSR = genlep_postFSR->Momentum;

	TLorentzVector SumPhotonMom;
	SumPhotonMom.SetPxPyPzE(0,0,0,0);

	Int_t NGenOthers = ntuple->nGenOthers;
	for(Int_t i_gen=0; i_gen<NGenOthers; i_gen++)
	{
		GenOthers genlep;
		genlep.FillFromNtuple(ntuple, i_gen);

		// -- all photons within dR < 0.1 -- //
		// if( fabs(genlep.ID) == 22 && fabs(genlep.Mother) == 13)
		if( fabs(genlep.ID) == 22 )
		{
			
			Double_t dR = Calc_dR_GenLepton_GenOthers(*genlep_postFSR, genlep);

			// -- Sum of all photon's momentum near the post-FSR muon -- //
			if( dR < dRCut )
			{
				SumPhotonMom  = SumPhotonMom + genlep.Momentum;
				GenPhotonCollection->push_back( genlep );
			}
		}
	}

	// -- Momentum(pre-FSR) = Momentum(post-FSR) + Sum of all Photon's momentum near the post-FSR muon -- //
	genlep_preFSR->Momentum = genlep_Mom_postFSR + SumPhotonMom;
	genlep_preFSR->Et = genlep_preFSR->Momentum.Et();
	genlep_preFSR->Pt = genlep_preFSR->Momentum.Pt();
	genlep_preFSR->eta = genlep_preFSR->Momentum.Eta();
	genlep_preFSR->phi = genlep_preFSR->Momentum.Phi();
	genlep_preFSR->Px = genlep_preFSR->Momentum.Px();
	genlep_preFSR->Py = genlep_preFSR->Momentum.Py();
	genlep_preFSR->Pz = genlep_preFSR->Momentum.Pz();
}

TString DYAnalyzer::DecideFSRType(GenLepton preFSR1, GenLepton preFSR2, GenLepton postFSR1, GenLepton postFSR2)
{
	TString FSRType = "";

	Bool_t isPassAcc_preFSREvent = kFALSE;
	isPassAcc_preFSREvent = isPassAccCondition_GenLepton(preFSR1, preFSR2);

	Bool_t isPassAcc_postFSREvent = kFALSE;
	isPassAcc_postFSREvent = isPassAccCondition_GenLepton(postFSR1, postFSR2);


	if( isPassAcc_preFSREvent == kFALSE && isPassAcc_postFSREvent == kTRUE )
		FSRType = "A";

	else if( isPassAcc_preFSREvent == kTRUE && isPassAcc_postFSREvent == kTRUE)
		FSRType = "B";
	
	else if( isPassAcc_preFSREvent == kTRUE && isPassAcc_postFSREvent == kFALSE)
		FSRType = "C";

	else if( isPassAcc_preFSREvent == kFALSE && isPassAcc_postFSREvent == kFALSE)
		FSRType = "D";
	else
	{
		cout << "ERROR: NO FSR TYPE CORRESPONDING TO THIS EVENT" << endl;
		FSRType = "NOTAssigned";
	}

	return FSRType;
}

Double_t DYAnalyzer::Calc_dR_GenLeptons( GenLepton genlep1, GenLepton genlep2 )
{
	Double_t eta1 = genlep1.eta;
	Double_t phi1 = genlep1.phi;

	Double_t eta2 = genlep2.eta;
	Double_t phi2 = genlep2.phi;

	Double_t diff_eta = eta1 - eta2;
	Double_t diff_phi = phi1 - phi2;

	Double_t dR = sqrt( diff_eta * diff_eta + diff_phi * diff_phi );
	return dR;
}

Double_t DYAnalyzer::Calc_dR_GenLepton_GenOthers( GenLepton genlep1, GenOthers genlep2 )
{
	Double_t eta1 = genlep1.eta;
	Double_t phi1 = genlep1.phi;

	Double_t eta2 = genlep2.eta;
	Double_t phi2 = genlep2.phi;

	Double_t diff_eta = eta1 - eta2;
	Double_t diff_phi = phi1 - phi2;

	Double_t dR = sqrt( diff_eta * diff_eta + diff_phi * diff_phi );
	return dR;
}

void DYAnalyzer::GenMatching(TString MuonType, NtupleHandle* ntuple, vector<Muon>* MuonCollection)
{
	vector<GenLepton> GenLeptonCollection;
	Int_t NGenLeptons = ntuple->gnpair;

	if( MuonType == "PromptFinalState" )
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isPromptFinalState )
				GenLeptonCollection.push_back( genlep );
		}
	}
	else if( MuonType == "fromTau")
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.isDirectPromptTauDecayProductFinalState )
				GenLeptonCollection.push_back( genlep );
		}

	}
	else if( MuonType == "fromHardProcess" )
	{
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( genlep.isMuon() && genlep.fromHardProcessFinalState )
				GenLeptonCollection.push_back( genlep );
		}
	}
	else
	{
		cout << "Incorrect MuonType!" << endl;
		return;
	}

	//Give Acceptance Cuts
	if( GenLeptonCollection.size() >= 2 )
	{
		GenLepton leadGenLep, subGenLep;
		CompareGenLepton(&GenLeptonCollection[0], &GenLeptonCollection[1], &leadGenLep, &subGenLep);
		if( !(leadGenLep.Pt > LeadPtCut && subGenLep.Pt > SubPtCut && abs(leadGenLep.eta) < LeadEtaCut && abs(subGenLep.eta) < SubEtaCut) )
			GenLeptonCollection.clear();
	}


	
	Int_t NMuons = (Int_t)MuonCollection->size();
	vector<Muon> RecoMuonCollection;
	//Copy all muons in MuonCollection into RecoMuonCollection
	for(Int_t i_mu=0; i_mu<NMuons; i_mu++)
		RecoMuonCollection.push_back( MuonCollection->at(i_mu) );

	MuonCollection->clear();

	Int_t NGen = (Int_t)GenLeptonCollection.size();
	for(Int_t i_gen=0; i_gen<NGen; i_gen++)
	{
		GenLepton genlep = GenLeptonCollection[i_gen];
		Double_t gen_Pt = genlep.Pt;
		Double_t gen_eta = genlep.eta;
		Double_t gen_phi = genlep.phi;

		Int_t i_matched = -1;
		Double_t dPtMin = 1e10;
		for(Int_t i_reco=0; i_reco<NMuons; i_reco++)
		{
			Muon mu = RecoMuonCollection[i_reco];
			Double_t reco_Pt = mu.Pt;
			Double_t reco_eta = mu.eta;
			Double_t reco_phi = mu.phi;

			Double_t dR = sqrt( (gen_eta-reco_eta)*(gen_eta-reco_eta) + (gen_phi-reco_phi)*(gen_phi-reco_phi) );
			Double_t dPt = fabs(gen_Pt - reco_Pt);
			if( dR < 0.3 )
			{
				if( dPt < dPtMin )
				{
					i_matched = i_reco;
					dPtMin = dPt;
				}
			}
		}

		if( i_matched != -1 )
			MuonCollection->push_back( RecoMuonCollection[i_matched] );
	}

	return;
}

void DYAnalyzer::ConvertToTunePInfo( Muon &mu )
{
	// -- Use TuneP information -- //
	mu.Pt = mu.TuneP_pT;
	mu.eta = mu.TuneP_eta;
	mu.phi = mu.TuneP_phi;

	/*Double_t Px = mu.TuneP_Px;
	Double_t Py = mu.TuneP_Py;
	Double_t Pz = mu.TuneP_Pz;
	Double_t E = sqrt( Px*Px + Py*Py + Pz*Pz + M_Mu*M_Mu );
	mu.Momentum.SetPxPyPzE( Px, Py, Pz, E );*/
	mu.Momentum.SetPtEtaPhiM( mu.Pt, mu.eta, mu.phi, M_Mu );
}

void DYAnalyzer::PrintOutDoubleMuInfo( Muon mu1, Muon mu2 )
{
	printf("\t[Muon1] (pT, eta, phi, charge) = (%10.5lf, %10.5lf, %10.5lf, %.1d)\n", mu1.Pt, mu1.eta, mu1.phi, mu1.charge);
	printf("\t[Muon2] (pT, eta, phi, charge) = (%10.5lf, %10.5lf, %10.5lf, %.1d)\n", mu2.Pt, mu2.eta, mu2.phi, mu2.charge);
	Double_t reco_M = ( mu1.Momentum + mu2.Momentum ).M();
	printf("\t\tDilepton Mass: %10.5lf\n", reco_M);

}

Bool_t DYAnalyzer::EventSelection_Dijet(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > PassingMuonCollection;
	vector< Muon > FailingMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() )
	    {
	    	if( MuonCollection[j].trkiso < 0.10 )
	    		PassingMuonCollection.push_back( MuonCollection[j] );
	    	else
	    		FailingMuonCollection.push_back( MuonCollection[j] );
	    }
	}

	Int_t nFailMuon = (Int_t)FailingMuonCollection.size();

	if( nFailMuon >= 2 ) // -- Dijet events: contains more than 2 failing muons regardless of # passing muons -- // 
	{
		if( nFailMuon == 2 )
		{
			Muon recolep1 = FailingMuonCollection[0];
			Muon recolep2 = FailingMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}

		} // -- end of if( nFailMuon == 2 ) -- //
		else // -- # failing muons > 2 -- // 
		{
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nFailMuon; i_mu++)
			{
				Muon Mu = FailingMuonCollection[i_mu];

				// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
				for(Int_t j_mu=0; j_mu<nFailMuon; j_mu++)
				{
					Muon Mu_jth = FailingMuonCollection[j_mu];

					if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
					{
						// -- Check that this pair is within acceptance -- //
						Bool_t isPassAcc = kFALSE;
						isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

						if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
						{
							Double_t VtxProb_temp = -999;
							Double_t VtxNormChi2_temp = 999;
							DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

							// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
							if( VtxNormChi2_temp < VtxNormChi2_BestPair )
							{
								VtxNormChi2_BestPair = VtxNormChi2_temp;
								mu1_BestPair = Mu;
								mu2_BestPair = Mu_jth;
							}
						}
					}
				} // -- end of the loop for j_mu (finding for second muon)
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- end of (# failing muons > 2) case -- //

	} // -- end of if( nFailMuon >= 2 ) -- //

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_Wjet(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > PassingMuonCollection;
	vector< Muon > FailingMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() )
	    {
	    	if( MuonCollection[j].trkiso < 0.10 )
	    		PassingMuonCollection.push_back( MuonCollection[j] );
	    	else
	    		FailingMuonCollection.push_back( MuonCollection[j] );
	    }
	}

	Int_t nFailMuon = (Int_t)FailingMuonCollection.size();
	Int_t nPassMuon = (Int_t)PassingMuonCollection.size();

	if( nFailMuon == 1 && nPassMuon == 1) // -- W+Jets events: exactly (# pass muon , # fail muon ) = (1, 1) -- //
	{
		Muon recolep1 = PassingMuonCollection[0]; // -- first one: passing muon -- //
		Muon recolep2 = FailingMuonCollection[1]; // -- second one: failing muon -- //

		// -- Check the Accpetance -- //
		Bool_t isPassAcc = kFALSE;
		isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

		Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

		Double_t VtxProb = -999;
		Double_t VtxNormChi2 = 999;
		DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

		TLorentzVector inner_v1 = recolep1.Momentum_Inner;
		TLorentzVector inner_v2 = recolep2.Momentum_Inner;

		// -- 3D open angle -- //
		Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

		Bool_t isOS = kFALSE;
		if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

		// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
		if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
		{
			isPassEventSelection = kTRUE;
			SelectedMuonCollection->push_back( recolep1 ); // -- first one: passing muon -- //
			SelectedMuonCollection->push_back( recolep2 ); // -- second one: failing muon -- //
		}
	}

	return isPassEventSelection;
}

Bool_t DYAnalyzer::EventSelection_CheckMoreThanOneDimuonCand(vector< Muon > MuonCollection, NtupleHandle *ntuple, // -- input: All muons in a event & NtupleHandle -- //
						vector< Muon >* SelectedMuonCollection, Bool_t& isMoreThanOneCand) // -- output: 2 muons passing event selection conditions -- //
{
	Bool_t isPassEventSelection = kFALSE;
	isMoreThanOneCand = kFALSE;

	//Collect qualified muons among muons
	vector< Muon > QMuonCollection;
	for(Int_t j=0; j<(int)MuonCollection.size(); j++)
	{
	    if( MuonCollection[j].isHighPtMuon_minus_dzVTX() && MuonCollection[j].trkiso < 0.10)
	        QMuonCollection.push_back( MuonCollection[j] );
	}

	// -- Check the existence of at least one muon matched with HLT-object -- //
	Bool_t isExistHLTMatchedMuon = kFALSE;
	for(Int_t i_mu=0; i_mu<(Int_t)QMuonCollection.size(); i_mu++)
	{
		Muon mu = QMuonCollection[i_mu];
		if( mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
		{
			isExistHLTMatchedMuon = kTRUE;
			break;
		}
	}

	if( isExistHLTMatchedMuon == kTRUE )
	{
		Int_t nQMuons = (Int_t)QMuonCollection.size();
		if( nQMuons == 2)
		{
			Muon recolep1 = QMuonCollection[0];
			Muon recolep2 = QMuonCollection[1];

			// -- Check the Accpetance -- //
			Bool_t isPassAcc = kFALSE;
			isPassAcc = isPassAccCondition_Muon(recolep1, recolep2);

			Double_t reco_M = (recolep1.Momentum + recolep2.Momentum).M();

			Double_t VtxProb = -999;
			Double_t VtxNormChi2 = 999;
			DimuonVertexProbNormChi2(ntuple, recolep1.Inner_pT, recolep2.Inner_pT, &VtxProb, &VtxNormChi2);

			TLorentzVector inner_v1 = recolep1.Momentum_Inner;
			TLorentzVector inner_v2 = recolep2.Momentum_Inner;

			// -- 3D open angle -- //
			Double_t Angle = recolep1.Momentum.Angle( recolep2.Momentum.Vect() );

			Bool_t isOS = kFALSE;
			if( recolep1.charge != recolep2.charge ) isOS = kTRUE;

			// if( reco_M > 10 && isPassAcc == kTRUE && Chi2/ndof(VTX) < 20 && Angle < TMath::Pi() - 0.005 )
			if( reco_M > 10 && isPassAcc == kTRUE && VtxNormChi2 < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
			{
				isPassEventSelection = kTRUE;
				SelectedMuonCollection->push_back( recolep1 );
				SelectedMuonCollection->push_back( recolep2 );
			}
		}
		else if( nQMuons > 2 )
		{
			isMoreThanOneCand = kTRUE;
			Double_t VtxProb_BestPair = -1;
			Double_t VtxNormChi2_BestPair = 999;
			Muon mu1_BestPair;
			Muon mu2_BestPair;

			for(Int_t i_mu=0; i_mu<nQMuons; i_mu++)
			{
				Muon Mu = QMuonCollection[i_mu];

				// -- at least 1 muon should be matched with HLT objects in best pair -- //
				if( Mu.isTrigMatched(ntuple, "HLT_IsoMu20_v*") || Mu.isTrigMatched(ntuple, "HLT_IsoTkMu20_v*") )
				{
					// -- Mu in this loop: QMuon Matched with HLT object -- //

					// -- Start another loop for finding second muon (for second muon, we don't need to check trigger matching) -- //
					for(Int_t j_mu=0; j_mu<nQMuons; j_mu++)
					{
						Muon Mu_jth = QMuonCollection[j_mu];

						if( j_mu != i_mu ) // -- do not calculate vertex variables(prob, chi2). with itself -- //
						{
							// -- Check that this pair is within acceptance -- //
							Bool_t isPassAcc = kFALSE;
							isPassAcc = isPassAccCondition_Muon(Mu, Mu_jth);

							if( isPassAcc == kTRUE ) // -- Find best pair ONLY for the pairs within acceptance -- //
							{
								Double_t VtxProb_temp = -999;
								Double_t VtxNormChi2_temp = 999;
								DimuonVertexProbNormChi2(ntuple, Mu.Inner_pT, Mu_jth.Inner_pT, &VtxProb_temp, &VtxNormChi2_temp);

								// -- Find best pair by selecting smallest Chi2/dnof(VTX) value -- // 
								if( VtxNormChi2_temp < VtxNormChi2_BestPair )
								{
									VtxNormChi2_BestPair = VtxNormChi2_temp;
									mu1_BestPair = Mu;
									mu2_BestPair = Mu_jth;
								}
							}
						}
					} // -- end of the loop for j_mu (finding for second muon)
				}
			} // -- end of the loop for i_mu (finding for the first muon matched with HLT matching)

			if( VtxNormChi2_BestPair < 999 ) // -- If at least one pair within acceptance & with at least one muon matched with HLT object exists -- //
			{
				TLorentzVector reco_v1 = mu1_BestPair.Momentum;
				TLorentzVector reco_v2 = mu2_BestPair.Momentum;
				Double_t reco_M = (reco_v1 + reco_v2).M();

				// -- 3D open angle is calculated using inner track information -- //
				// -- 3D open angle -- //
				Double_t Angle = reco_v1.Angle( reco_v2.Vect() );

				Bool_t isOS = kFALSE;
				if( mu1_BestPair.charge != mu2_BestPair.charge ) isOS = kTRUE;

				if( reco_M > 10 && VtxNormChi2_BestPair < 20 && Angle < TMath::Pi() - 0.005 && isOS == kTRUE )
				{
					isPassEventSelection = kTRUE;
					SelectedMuonCollection->push_back( mu1_BestPair );
					SelectedMuonCollection->push_back( mu2_BestPair );
				}
			}

		} // -- End of else if( nQMuons > 2 ) -- //

	} // -- End of if( isExistHLTMatchedMuon == kTRUE ) -- //

	return isPassEventSelection;
}





