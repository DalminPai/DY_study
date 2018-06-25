/////////////////////////////////////////////////////////////////////////////////////
// -- 2015.04.02: Remove duplicated declaration of gen-variables & Trigger
// -- 2016.01.02: Cleaning & Add Electron, Jet and MET informaiton
// -- 2017.04.27: Adding "Electron_passMediumID" (by Dalmin Pai)
// -- 2017.07.19: Updating Electron variables with EGM corrections (by Dalmin Pai)
// -- 2017.07.26: Varialbe type of "Electron_passConvVeto" and "Electron_passMediumID" were changed to "Bool_t" (by Dalmin Pai)
// -- 2017.10.31: Changing Electron variables (by Dalmin Pai)
// -- 2018.01.24: Add "Muon_PFSumPUIsoR04" to use "RelPFIso_dBeta" in muon
// -- 2018.02.20: Add "Mu50_OR_TkMu50" in isTriggered
/////////////////////////////////////////////////////////////////////////////////////
#pragma once

#define MaxN 50000
#include <TTree.h>
#include <TChain.h>
#include <vector>

class NtupleHandle
{
public:
	TChain *chain;

    //Event Informations
    Int_t nVertices;
    Int_t runNum;
    Int_t lumiBlock;
    //Int_t evtNum;
    unsigned long long evtNum;
    Int_t nPileUp;

    //Trigger variables
    Int_t HLT_ntrig;
    Int_t HLT_trigFired[MaxN];
    vector<string> *HLT_trigName;
    Double_t HLT_trigEta[MaxN];
    Double_t HLT_trigPhi[MaxN];

    //Generator level information
    Int_t gnpair;
    Double_t GenLepton_px[MaxN];
    Double_t GenLepton_py[MaxN];
    Double_t GenLepton_pz[MaxN];
    Double_t GenLepton_mother[MaxN];
    Double_t GenLepton_pT[MaxN];
    Double_t GenLepton_eta[MaxN];
    Double_t GenLepton_phi[MaxN];
    Int_t GenLepton_charge[MaxN];
    Int_t GenLepton_status[MaxN];
    Int_t GenLepton_ID[MaxN];
    Double_t GENEvt_weight;
    Int_t GenLepton_isPrompt[MaxN];
    Int_t GenLepton_isPromptFinalState[MaxN];
    Int_t GenLepton_isTauDecayProduct[MaxN];
    Int_t GenLepton_isPromptTauDecayProduct[MaxN];
    Int_t GenLepton_isDirectPromptTauDecayProductFinalState[MaxN];
    Int_t GenLepton_isHardProcess[MaxN];
    Int_t GenLepton_isLastCopy[MaxN];
    Int_t GenLepton_isLastCopyBeforeFSR[MaxN];
    Int_t GenLepton_isPromptDecayed[MaxN];
    Int_t GenLepton_isDecayedLeptonHadron[MaxN];
    Int_t GenLepton_fromHardProcessBeforeFSR[MaxN];
    Int_t GenLepton_fromHardProcessDecayed[MaxN];
    Int_t GenLepton_fromHardProcessFinalState[MaxN];

    Int_t nGenOthers;
    Double_t GenOthers_px[MaxN];
    Double_t GenOthers_py[MaxN];
    Double_t GenOthers_pz[MaxN];
    Double_t GenOthers_mother[MaxN];
    Double_t GenOthers_pT[MaxN];
    Double_t GenOthers_eta[MaxN];
    Double_t GenOthers_phi[MaxN];
    Int_t GenOthers_charge[MaxN];
    Int_t GenOthers_status[MaxN];
    Int_t GenOthers_ID[MaxN];
    Int_t GenOthers_isPrompt[MaxN];
    Int_t GenOthers_isPromptFinalState[MaxN];
    Int_t GenOthers_isTauDecayProduct[MaxN];
    Int_t GenOthers_isPromptTauDecayProduct[MaxN];
    Int_t GenOthers_isDirectPromptTauDecayProductFinalState[MaxN];
    Int_t GenOthers_isHardProcess[MaxN];
    Int_t GenOthers_isLastCopy[MaxN];
    Int_t GenOthers_isLastCopyBeforeFSR[MaxN];
    Int_t GenOthers_isPromptDecayed[MaxN];
    Int_t GenOthers_isDecayedLeptonHadron[MaxN];
    Int_t GenOthers_fromHardProcessBeforeFSR[MaxN];
    Int_t GenOthers_fromHardProcessDecayed[MaxN];
    Int_t GenOthers_fromHardProcessFinalState[MaxN];


    //////////////////////////////
    // -- Electron Variables -- //
    //////////////////////////////
    Int_t Nelectrons;
    Double_t Electron_Energy[MaxN];
    Double_t Electron_pT[MaxN];
    Double_t Electron_Px[MaxN];
    Double_t Electron_Py[MaxN];
    Double_t Electron_Pz[MaxN];
    Double_t Electron_eta[MaxN];
    Double_t Electron_phi[MaxN];
    Int_t Electron_charge[MaxN];
    Double_t Electron_gsfpT[MaxN];
    Double_t Electron_gsfPx[MaxN];
    Double_t Electron_gsfPy[MaxN];
    Double_t Electron_gsfPz[MaxN];
    Double_t Electron_gsfEta[MaxN];
    Double_t Electron_gsfPhi[MaxN];
    Double_t Electron_gsfCharge[MaxN];
    Double_t Electron_etaSC[MaxN];
    Double_t Electron_phiSC[MaxN];
    Double_t Electron_etaWidth[MaxN];
    Double_t Electron_phiWidth[MaxN];
    Double_t Electron_dEtaIn[MaxN];
    Double_t Electron_dEtaInSeed[MaxN]; // updated at 19 Jul. 2017 by Dalmin
    Double_t Electron_dPhiIn[MaxN];
    Double_t Electron_sigmaIEtaIEta[MaxN];
    Double_t Electron_Full5x5_SigmaIEtaIEta[MaxN]; // updated at 19 Jul. 2017 by Dalmin
    Double_t Electron_HoverE[MaxN];
    Double_t Electron_fbrem[MaxN];
    Double_t Electron_eOverP[MaxN];
    Double_t Electron_InvEminusInvP[MaxN];
    Double_t Electron_dxyVTX[MaxN];
    Double_t Electron_dzVTX[MaxN];
    Double_t Electron_dxy[MaxN];
    Double_t Electron_dz[MaxN];
    Double_t Electron_dxyBS[MaxN];
    Double_t Electron_dzBS[MaxN];
    Double_t Electron_chIso03[MaxN];
    Double_t Electron_nhIso03[MaxN];
    Double_t Electron_phIso03[MaxN];
    Double_t Electron_ChIso03FromPU[MaxN];
    Int_t Electron_mHits[MaxN];
    Double_t Electron_EnergySC[MaxN];
    Double_t Electron_preEnergySC[MaxN];
    Double_t Electron_rawEnergySC[MaxN];
    Double_t Electron_etSC[MaxN];
    Double_t Electron_E15[MaxN];
    Double_t Electron_E25[MaxN];
    Double_t Electron_E55[MaxN];
    Double_t Electron_RelPFIso_dBeta[MaxN];
    Double_t Electron_RelPFIso_Rho[MaxN];
    Double_t Electron_r9[MaxN];
    Int_t Electron_ecalDriven[MaxN];
//    Int_t Electron_passConvVeto[MaxN];
    Bool_t Electron_passConvVeto[MaxN]; // modified at 26 Jul. 2017 by Dalmin
//    Int_t Electron_passMediumID[MaxN]; // updated at 27 Apr. 2017 by Dalmin
    Bool_t Electron_passMediumID[MaxN]; // modified at 26 Jul. 2017 by Dalmin

    // -- Uncorrected Electrons -- //
    Int_t nUnCorrElectron;
    Double_t Electron_pTUnCorr[MaxN];
    Double_t Electron_etaUnCorr[MaxN];
    Double_t Electron_phiUnCorr[MaxN];
    Double_t Electron_PxUnCorr[MaxN];
    Double_t Electron_PyUnCorr[MaxN];
    Double_t Electron_PzUnCorr[MaxN];
    Double_t Electron_EnergyUnCorr[MaxN];
    Double_t Electron_EnergySCUnCorr[MaxN];
    Double_t Electron_etaSCUnCorr[MaxN];
    Double_t Electron_phiSCUnCorr[MaxN];
    Double_t Electron_etSCUnCorr[MaxN];

    //////////////////////////
    // -- Muon Variables -- //
    //////////////////////////
    // -- Physical Variables -- //
    Double_t Muon_pT[MaxN];
    Double_t Muon_eta[MaxN];
    Double_t Muon_phi[MaxN];
    
    // -- Cut variables -- //
    Int_t Muon_muonType[MaxN];
    Double_t Muon_chi2dof[MaxN];
    Int_t Muon_muonHits[MaxN];
    Int_t Muon_nSegments[MaxN];
    Int_t Muon_nMatches[MaxN];
    Int_t Muon_trackerLayers[MaxN];
    Int_t Muon_pixelHits[MaxN];
    Double_t Muon_dxyVTX[MaxN];
    Double_t Muon_dzVTX[MaxN];
    Double_t Muon_trkiso[MaxN];
    Int_t isGLBmuon[MaxN];
    Int_t isPFmuon[MaxN];
    Int_t isTRKmuon[MaxN];
    Int_t nMuon;
    
    // -- for invariant mass calculation -- //
    Double_t Muon_Px[MaxN];
    Double_t Muon_Py[MaxN];
    Double_t Muon_Pz[MaxN];
    
    // -- for the muon momentum corrections -- //
    Int_t Muon_charge[MaxN];
    
    //PF information
    Double_t Muon_PfChargedHadronIsoR04[MaxN];
    Double_t Muon_PfNeutralHadronIsoR04[MaxN];
    Double_t Muon_PfGammaIsoR04[MaxN];
    Double_t Muon_PFSumPUIsoR04[MaxN];

    //Dimuon variables
    std::vector<double> *CosAngle;
    std::vector<double> *vtxTrkChi2;
    std::vector<double> *vtxTrkProb;
    std::vector<double> *vtxTrkNdof;
    std::vector<double> *vtxTrkCkt1Pt;
    std::vector<double> *vtxTrkCkt2Pt;
    std::vector<double> *vtxTrkDiEChi2;
    std::vector<double> *vtxTrkDiEProb;
    std::vector<double> *vtxTrkDiENdof;
    std::vector<double> *vtxTrkDiE1Pt;
    std::vector<double> *vtxTrkDiE2Pt;
    std::vector<double> *vtxTrkEMuChi2;
    std::vector<double> *vtxTrkEMuProb;
    std::vector<double> *vtxTrkEMuNdof;
    std::vector<double> *vtxTrkEMu1Pt;
    std::vector<double> *vtxTrkEMu2Pt;

    // -- Various Track Information -- //
    Double_t Muon_Best_pT[MaxN];
    Double_t Muon_Best_pTError[MaxN];
    Double_t Muon_Best_Px[MaxN];
    Double_t Muon_Best_Py[MaxN];
    Double_t Muon_Best_Pz[MaxN];
    Double_t Muon_Best_eta[MaxN];
    Double_t Muon_Best_phi[MaxN];

    Double_t Muon_Inner_pT[MaxN];
    Double_t Muon_Inner_pTError[MaxN];
    Double_t Muon_Inner_Px[MaxN];
    Double_t Muon_Inner_Py[MaxN];
    Double_t Muon_Inner_Pz[MaxN];
    Double_t Muon_Inner_eta[MaxN];
    Double_t Muon_Inner_phi[MaxN];

    Double_t Muon_Outer_pT[MaxN];
    Double_t Muon_Outer_pTError[MaxN];
    Double_t Muon_Outer_Px[MaxN];
    Double_t Muon_Outer_Py[MaxN];
    Double_t Muon_Outer_Pz[MaxN];
    Double_t Muon_Outer_eta[MaxN];
    Double_t Muon_Outer_phi[MaxN];

    Double_t Muon_GLB_pT[MaxN];
    Double_t Muon_GLB_pTError[MaxN];
    Double_t Muon_GLB_Px[MaxN];
    Double_t Muon_GLB_Py[MaxN];
    Double_t Muon_GLB_Pz[MaxN];
    Double_t Muon_GLB_eta[MaxN];
    Double_t Muon_GLB_phi[MaxN];

    Double_t Muon_TuneP_pT[MaxN];
    Double_t Muon_TuneP_pTError[MaxN];
    Double_t Muon_TuneP_Px[MaxN];
    Double_t Muon_TuneP_Py[MaxN];
    Double_t Muon_TuneP_Pz[MaxN];
    Double_t Muon_TuneP_eta[MaxN];
    Double_t Muon_TuneP_phi[MaxN];

    // -- Photon Information -- //
    Int_t nPhotons;
    Int_t Photon_hasPixelSeed[MaxN];
    Double_t Photon_pT[MaxN];
    Double_t Photon_eta[MaxN];
    Double_t Photon_phi[MaxN];
    Double_t Photon_etaSC[MaxN];
    Double_t Photon_phiSC[MaxN];
    Double_t Photon_HoverE[MaxN];
    Double_t Photon_Full5x5_SigmaIEtaIEta[MaxN];
    Double_t Photon_ChIso[MaxN];
    Double_t Photon_NhIso[MaxN];
    Double_t Photon_PhIso[MaxN];
    Double_t Photon_ChIsoWithEA[MaxN];
    Double_t Photon_NhIsoWithEA[MaxN];
    Double_t Photon_PhIsoWithEA[MaxN];


    // -- Jet Information -- //
    Int_t Njets;
    Double_t Jet_pT[MaxN];
    Double_t Jet_eta[MaxN];
    Double_t Jet_phi[MaxN];
    Double_t Jet_Charge[MaxN];
    Int_t Jet_Flavor[MaxN];

    Double_t Jet_bTag[MaxN];
    Double_t Jet_CHfrac[MaxN];
    Double_t Jet_NHfrac[MaxN];
    Double_t Jet_NHEMfrac[MaxN];
    Double_t Jet_CHEMfrac[MaxN];
    Int_t Jet_CHmulti[MaxN];
    Int_t Jet_NHmulti[MaxN];

    // -- MET Information -- //
    Double_t pfMET_pT;
    Double_t pfMET_phi;
    Double_t pfMET_Px;
    Double_t pfMET_Py;
    Double_t pfMET_SumEt;

    Double_t pfMET_Type1_pT;
    Double_t pfMET_Type1_phi;
    Double_t pfMET_Type1_Px;
    Double_t pfMET_Type1_Py;
    Double_t pfMET_Type1_SumEt;



    // -- Constructor -- //
    NtupleHandle(TChain *chainptr)
    {
    	chain = chainptr;
        chain->SetBranchStatus("*", 0);

    	// -- Event Information -- //
    	chain->SetBranchStatus("nVertices", 1);
    	chain->SetBranchStatus("runNum", 1);
        chain->SetBranchStatus("lumiBlock", 1);
        chain->SetBranchStatus("evtNum", 1);
    	chain->SetBranchStatus("nPileUp", 1);

    	chain->SetBranchAddress("nVertices", &nVertices);
    	chain->SetBranchAddress("runNum", &runNum);
        chain->SetBranchAddress("lumiBlock", &lumiBlock);
        chain->SetBranchAddress("evtNum", &evtNum);
    	chain->SetBranchAddress("nPileUp", &nPileUp);

    	// -- Trigger Information -- //
    	chain->SetBranchStatus("HLT_trigName", 1);
    	chain->SetBranchStatus("HLT_ntrig", 1);
    	chain->SetBranchStatus("HLT_trigFired", 1);
    	chain->SetBranchStatus("HLT_trigEta", 1);
    	chain->SetBranchStatus("HLT_trigPhi", 1);

    	chain->SetBranchAddress("HLT_trigName", &HLT_trigName);
    	chain->SetBranchAddress("HLT_trigFired", HLT_trigFired);
    	chain->SetBranchAddress("HLT_ntrig", &HLT_ntrig);
    	chain->SetBranchAddress("HLT_trigEta", &HLT_trigEta);
    	chain->SetBranchAddress("HLT_trigPhi", &HLT_trigPhi);


    	// this->TurnOnBranches_GenLepton();
    	// this->TurnOnBranches_Muon();

    }

    void TurnOnBranches_GenLepton()
    {
    	//Generator level information
    	chain->SetBranchStatus("GENnPair", 1);
    	chain->SetBranchStatus("GENEvt_weight", 1);
        chain->SetBranchStatus("GENLepton_Px", 1);
        chain->SetBranchStatus("GENLepton_Py", 1);
        chain->SetBranchStatus("GENLepton_Pz", 1);
        chain->SetBranchStatus("GENLepton_mother", 1);
        chain->SetBranchStatus("GENLepton_pT", 1);
        chain->SetBranchStatus("GENLepton_eta", 1);
        chain->SetBranchStatus("GENLepton_phi", 1);
        chain->SetBranchStatus("GENLepton_charge", 1);
        chain->SetBranchStatus("GENLepton_status", 1);
        chain->SetBranchStatus("GENLepton_ID", 1);
        chain->SetBranchStatus("GENLepton_isPrompt", 1);
        chain->SetBranchStatus("GENLepton_isPromptFinalState", 1);
        chain->SetBranchStatus("GENLepton_isTauDecayProduct", 1);
        chain->SetBranchStatus("GENLepton_isPromptTauDecayProduct", 1);
        chain->SetBranchStatus("GENLepton_isDirectPromptTauDecayProductFinalState", 1);
        chain->SetBranchStatus("GENLepton_isHardProcess", 1);
        chain->SetBranchStatus("GENLepton_isLastCopy", 1);
        chain->SetBranchStatus("GENLepton_isLastCopyBeforeFSR", 1);
        chain->SetBranchStatus("GENLepton_isPromptDecayed", 1);
        chain->SetBranchStatus("GENLepton_isDecayedLeptonHadron", 1);
        chain->SetBranchStatus("GENLepton_fromHardProcessBeforeFSR", 1);
        chain->SetBranchStatus("GENLepton_fromHardProcessDecayed", 1);
        chain->SetBranchStatus("GENLepton_fromHardProcessFinalState", 1);

    	chain->SetBranchAddress("GENnPair", &gnpair);
    	chain->SetBranchAddress("GENLepton_Px", GenLepton_px);
    	chain->SetBranchAddress("GENLepton_Py", GenLepton_py);
    	chain->SetBranchAddress("GENLepton_Pz", GenLepton_pz);
    	chain->SetBranchAddress("GENLepton_mother", &GenLepton_mother);
    	chain->SetBranchAddress("GENLepton_pT", GenLepton_pT);
    	chain->SetBranchAddress("GENLepton_eta", GenLepton_eta);
    	chain->SetBranchAddress("GENLepton_phi", GenLepton_phi);
    	chain->SetBranchAddress("GENLepton_charge", &GenLepton_charge);
    	chain->SetBranchAddress("GENLepton_status", &GenLepton_status);
    	chain->SetBranchAddress("GENLepton_ID", &GenLepton_ID);
    	chain->SetBranchAddress("GENEvt_weight", &GENEvt_weight);
    	chain->SetBranchAddress("GENLepton_isPrompt", &GenLepton_isPrompt);
    	chain->SetBranchAddress("GENLepton_isPromptFinalState", &GenLepton_isPromptFinalState);
    	chain->SetBranchAddress("GENLepton_isTauDecayProduct", &GenLepton_isTauDecayProduct);
    	chain->SetBranchAddress("GENLepton_isPromptTauDecayProduct", &GenLepton_isPromptTauDecayProduct);
    	chain->SetBranchAddress("GENLepton_isDirectPromptTauDecayProductFinalState", &GenLepton_isDirectPromptTauDecayProductFinalState);
    	chain->SetBranchAddress("GENLepton_isHardProcess", &GenLepton_isHardProcess);
    	chain->SetBranchAddress("GENLepton_isLastCopy", &GenLepton_isLastCopy);
    	chain->SetBranchAddress("GENLepton_isLastCopyBeforeFSR", &GenLepton_isLastCopyBeforeFSR);
    	chain->SetBranchAddress("GENLepton_isPromptDecayed", &GenLepton_isPromptDecayed);
    	chain->SetBranchAddress("GENLepton_isDecayedLeptonHadron", &GenLepton_isDecayedLeptonHadron);
    	chain->SetBranchAddress("GENLepton_fromHardProcessBeforeFSR", &GenLepton_fromHardProcessBeforeFSR);
    	chain->SetBranchAddress("GENLepton_fromHardProcessDecayed", &GenLepton_fromHardProcessDecayed);
    	chain->SetBranchAddress("GENLepton_fromHardProcessFinalState", &GenLepton_fromHardProcessFinalState);
    }

    void TurnOnBranches_GenOthers()
    {
        chain->SetBranchStatus("nGenOthers", 1);
        chain->SetBranchStatus("GenOthers_Px", 1);
        chain->SetBranchStatus("GenOthers_Py", 1);
        chain->SetBranchStatus("GenOthers_Pz", 1);
        chain->SetBranchStatus("GenOthers_mother", 1);
        chain->SetBranchStatus("GenOthers_pT", 1);
        chain->SetBranchStatus("GenOthers_eta", 1);
        chain->SetBranchStatus("GenOthers_phi", 1);
        chain->SetBranchStatus("GenOthers_charge", 1);
        chain->SetBranchStatus("GenOthers_status", 1);
        chain->SetBranchStatus("GenOthers_ID", 1);
        chain->SetBranchStatus("GenOthers_isPrompt", 1);
        chain->SetBranchStatus("GenOthers_isPromptFinalState", 1);
        chain->SetBranchStatus("GenOthers_isTauDecayProduct", 1);
        chain->SetBranchStatus("GenOthers_isPromptTauDecayProduct", 1);
        chain->SetBranchStatus("GenOthers_isDirectPromptTauDecayProductFinalState", 1);
        chain->SetBranchStatus("GenOthers_isHardProcess", 1);
        chain->SetBranchStatus("GenOthers_isLastCopy", 1);
        chain->SetBranchStatus("GenOthers_isLastCopyBeforeFSR", 1);
        chain->SetBranchStatus("GenOthers_isPromptDecayed", 1);
        chain->SetBranchStatus("GenOthers_isDecayedLeptonHadron", 1);
        chain->SetBranchStatus("GenOthers_fromHardProcessBeforeFSR", 1);
        chain->SetBranchStatus("GenOthers_fromHardProcessDecayed", 1);
        chain->SetBranchStatus("GenOthers_fromHardProcessFinalState", 1);

        chain->SetBranchAddress("nGenOthers", &nGenOthers);
        chain->SetBranchAddress("GenOthers_Px", GenOthers_px);
        chain->SetBranchAddress("GenOthers_Py", GenOthers_py);
        chain->SetBranchAddress("GenOthers_Pz", GenOthers_pz);
        chain->SetBranchAddress("GenOthers_mother", &GenOthers_mother);
        chain->SetBranchAddress("GenOthers_pT", GenOthers_pT);
        chain->SetBranchAddress("GenOthers_eta", GenOthers_eta);
        chain->SetBranchAddress("GenOthers_phi", GenOthers_phi);
        chain->SetBranchAddress("GenOthers_charge", &GenOthers_charge);
        chain->SetBranchAddress("GenOthers_status", &GenOthers_status);
        chain->SetBranchAddress("GenOthers_ID", &GenOthers_ID);
        chain->SetBranchAddress("GenOthers_isPrompt", &GenOthers_isPrompt);
        chain->SetBranchAddress("GenOthers_isPromptFinalState", &GenOthers_isPromptFinalState);
        chain->SetBranchAddress("GenOthers_isTauDecayProduct", &GenOthers_isTauDecayProduct);
        chain->SetBranchAddress("GenOthers_isPromptTauDecayProduct", &GenOthers_isPromptTauDecayProduct);
        chain->SetBranchAddress("GenOthers_isDirectPromptTauDecayProductFinalState", &GenOthers_isDirectPromptTauDecayProductFinalState);
        chain->SetBranchAddress("GenOthers_isHardProcess", &GenOthers_isHardProcess);
        chain->SetBranchAddress("GenOthers_isLastCopy", &GenOthers_isLastCopy);
        chain->SetBranchAddress("GenOthers_isLastCopyBeforeFSR", &GenOthers_isLastCopyBeforeFSR);
        chain->SetBranchAddress("GenOthers_isPromptDecayed", &GenOthers_isPromptDecayed);
        chain->SetBranchAddress("GenOthers_isDecayedLeptonHadron", &GenOthers_isDecayedLeptonHadron);
        chain->SetBranchAddress("GenOthers_fromHardProcessBeforeFSR", &GenOthers_fromHardProcessBeforeFSR);
        chain->SetBranchAddress("GenOthers_fromHardProcessDecayed", &GenOthers_fromHardProcessDecayed);
        chain->SetBranchAddress("GenOthers_fromHardProcessFinalState", &GenOthers_fromHardProcessFinalState);
    }

    void TurnOnBranches_Muon()
    {
    	// chain->SetBranchStatus("Muon_*", 1);

        // -- # branches to be turned on is critical to the speed of the code: 
        // -- Usage of * is strongly not recommended! ... It would be better to exactly write down all branches to be used -- //
    	chain->SetBranchStatus("isPFmuon", 1);
    	chain->SetBranchStatus("isGLBmuon", 1);
        chain->SetBranchStatus("isTRKmuon", 1);
    	chain->SetBranchStatus("CosAngle", 1);
    	chain->SetBranchStatus("vtxTrkChi2", 1);
    	chain->SetBranchStatus("vtxTrkProb", 1);
    	chain->SetBranchStatus("vtxTrkNdof", 1);
    	chain->SetBranchStatus("vtxTrkCkt1Pt", 1);
    	chain->SetBranchStatus("vtxTrkCkt2Pt", 1);
    	chain->SetBranchStatus("vtxTrkDiEChi2", 1);
    	chain->SetBranchStatus("vtxTrkDiEProb", 1);
    	chain->SetBranchStatus("vtxTrkDiENdof", 1);
    	chain->SetBranchStatus("vtxTrkDiE1Pt", 1);
    	chain->SetBranchStatus("vtxTrkDiE2Pt", 1);
    	chain->SetBranchStatus("vtxTrkEMuChi2", 1);
    	chain->SetBranchStatus("vtxTrkEMuProb", 1);
    	chain->SetBranchStatus("vtxTrkEMuNdof", 1);
    	chain->SetBranchStatus("vtxTrkEMu1Pt", 1);
    	chain->SetBranchStatus("vtxTrkEMu2Pt", 1);

        chain->SetBranchStatus("nMuon", 1);
        chain->SetBranchStatus("Muon_pT", 1);
        chain->SetBranchStatus("Muon_eta", 1);
        chain->SetBranchStatus("Muon_phi", 1);
        chain->SetBranchStatus("Muon_muonType", 1);
        chain->SetBranchStatus("Muon_chi2dof", 1);
        chain->SetBranchStatus("Muon_muonHits", 1);
        chain->SetBranchStatus("Muon_nSegments", 1);
        chain->SetBranchStatus("Muon_nMatches", 1);
        chain->SetBranchStatus("Muon_trackerLayers", 1);
        chain->SetBranchStatus("Muon_pixelHits", 1);
        chain->SetBranchStatus("Muon_dxyVTX", 1);
        chain->SetBranchStatus("Muon_dzVTX", 1);
        chain->SetBranchStatus("Muon_trkiso", 1);
        
        chain->SetBranchStatus("Muon_Px", 1);
        chain->SetBranchStatus("Muon_Py", 1);
        chain->SetBranchStatus("Muon_Pz", 1);
        
        chain->SetBranchStatus("Muon_charge", 1);
        
        chain->SetBranchStatus("Muon_PfChargedHadronIsoR04", 1);
        chain->SetBranchStatus("Muon_PfNeutralHadronIsoR04" ,1);
        chain->SetBranchStatus("Muon_PfGammaIsoR04", 1);
	chain->SetBranchStatus("Muon_PFSumPUIsoR04", 1);

        chain->SetBranchStatus("Muon_Best_pT", 1);
        chain->SetBranchStatus("Muon_Best_pTError", 1);
        chain->SetBranchStatus("Muon_Best_Px", 1);
        chain->SetBranchStatus("Muon_Best_Py", 1);
        chain->SetBranchStatus("Muon_Best_Pz", 1);
        chain->SetBranchStatus("Muon_Best_eta", 1);
        chain->SetBranchStatus("Muon_Best_phi", 1);

        chain->SetBranchStatus("Muon_Inner_pT", 1);
        chain->SetBranchStatus("Muon_Inner_pTError", 1);
        chain->SetBranchStatus("Muon_Inner_eta", 1);
        chain->SetBranchStatus("Muon_Inner_phi", 1);
        chain->SetBranchStatus("Muon_Inner_Px", 1);
        chain->SetBranchStatus("Muon_Inner_Py", 1);
        chain->SetBranchStatus("Muon_Inner_Pz", 1);

        chain->SetBranchStatus("Muon_Outer_pT", 1);
        chain->SetBranchStatus("Muon_Outer_pTError", 1);
        chain->SetBranchStatus("Muon_Outer_Px", 1);
        chain->SetBranchStatus("Muon_Outer_Py", 1);
        chain->SetBranchStatus("Muon_Outer_Pz", 1);
        chain->SetBranchStatus("Muon_Outer_eta", 1);
        chain->SetBranchStatus("Muon_Outer_phi", 1);

        chain->SetBranchStatus("Muon_GLB_pT", 1);
        chain->SetBranchStatus("Muon_GLB_pTError", 1);
        chain->SetBranchStatus("Muon_GLB_Px", 1);
        chain->SetBranchStatus("Muon_GLB_Py", 1);
        chain->SetBranchStatus("Muon_GLB_Pz", 1);
        chain->SetBranchStatus("Muon_GLB_eta", 1);
        chain->SetBranchStatus("Muon_GLB_phi", 1);

        chain->SetBranchStatus("Muon_TuneP_pT", 1);
        chain->SetBranchStatus("Muon_TuneP_pTError", 1);
        chain->SetBranchStatus("Muon_TuneP_eta", 1);
        chain->SetBranchStatus("Muon_TuneP_phi", 1);
        chain->SetBranchStatus("Muon_TuneP_Px", 1);
        chain->SetBranchStatus("Muon_TuneP_Py", 1);
        chain->SetBranchStatus("Muon_TuneP_Pz", 1);


    	chain->SetBranchAddress("nMuon", &nMuon);
    	chain->SetBranchAddress("Muon_pT", Muon_pT);
    	chain->SetBranchAddress("Muon_eta", Muon_eta);
    	chain->SetBranchAddress("Muon_phi", Muon_phi);
    	chain->SetBranchAddress("Muon_muonType", Muon_muonType);
    	chain->SetBranchAddress("Muon_chi2dof", Muon_chi2dof);
    	chain->SetBranchAddress("Muon_muonHits", Muon_muonHits);
    	chain->SetBranchAddress("Muon_nSegments", Muon_nSegments);
    	chain->SetBranchAddress("Muon_nMatches", Muon_nMatches);
    	chain->SetBranchAddress("Muon_trackerLayers", Muon_trackerLayers);
    	chain->SetBranchAddress("Muon_pixelHits", Muon_pixelHits);
    	chain->SetBranchAddress("Muon_dxyVTX", Muon_dxyVTX);
    	chain->SetBranchAddress("Muon_dzVTX", Muon_dzVTX);
    	chain->SetBranchAddress("Muon_trkiso", Muon_trkiso);
    	chain->SetBranchAddress("isPFmuon", isPFmuon);
    	chain->SetBranchAddress("isGLBmuon", isGLBmuon);
        chain->SetBranchAddress("isTRKmuon", isTRKmuon);    	
    	chain->SetBranchAddress("Muon_Px", Muon_Px );
    	chain->SetBranchAddress("Muon_Py", Muon_Py );
    	chain->SetBranchAddress("Muon_Pz", Muon_Pz );
    	
    	chain->SetBranchAddress("Muon_charge", Muon_charge);
    	
    	chain->SetBranchAddress("Muon_PfChargedHadronIsoR04", Muon_PfChargedHadronIsoR04);
    	chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04", Muon_PfNeutralHadronIsoR04);
    	chain->SetBranchAddress("Muon_PfGammaIsoR04", Muon_PfGammaIsoR04);
	chain->SetBranchAddress("Muon_PFSumPUIsoR04", Muon_PFSumPUIsoR04);

    	chain->SetBranchAddress("CosAngle", &CosAngle);
    	chain->SetBranchAddress("vtxTrkChi2", &vtxTrkChi2);
    	chain->SetBranchAddress("vtxTrkProb", &vtxTrkProb);
    	chain->SetBranchAddress("vtxTrkNdof", &vtxTrkNdof);
    	chain->SetBranchAddress("vtxTrkCkt1Pt", &vtxTrkCkt1Pt);
    	chain->SetBranchAddress("vtxTrkCkt2Pt", &vtxTrkCkt2Pt);
    	chain->SetBranchAddress("vtxTrkDiEChi2", &vtxTrkDiEChi2);
    	chain->SetBranchAddress("vtxTrkDiEProb", &vtxTrkDiEProb);
    	chain->SetBranchAddress("vtxTrkDiENdof", &vtxTrkDiENdof);
    	chain->SetBranchAddress("vtxTrkDiE1Pt", &vtxTrkDiE1Pt);
    	chain->SetBranchAddress("vtxTrkDiE2Pt", &vtxTrkDiE2Pt);
    	chain->SetBranchAddress("vtxTrkEMuChi2", &vtxTrkEMuChi2);
    	chain->SetBranchAddress("vtxTrkEMuProb", &vtxTrkEMuProb);
    	chain->SetBranchAddress("vtxTrkEMuNdof", &vtxTrkEMuNdof);
    	chain->SetBranchAddress("vtxTrkEMu1Pt", &vtxTrkEMu1Pt);
    	chain->SetBranchAddress("vtxTrkEMu2Pt", &vtxTrkEMu2Pt);
    	
    	chain->SetBranchAddress("Muon_Best_pT", &Muon_Best_pT);
    	chain->SetBranchAddress("Muon_Best_pTError", &Muon_Best_pTError);
    	chain->SetBranchAddress("Muon_Best_Px", &Muon_Best_Px);
    	chain->SetBranchAddress("Muon_Best_Py", &Muon_Best_Py);
    	chain->SetBranchAddress("Muon_Best_Pz", &Muon_Best_Pz);
    	chain->SetBranchAddress("Muon_Best_eta", &Muon_Best_eta);
    	chain->SetBranchAddress("Muon_Best_phi", &Muon_Best_phi);

    	chain->SetBranchAddress("Muon_Inner_pT", &Muon_Inner_pT);
    	chain->SetBranchAddress("Muon_Inner_pTError", &Muon_Inner_pTError);
    	chain->SetBranchAddress("Muon_Inner_Px", &Muon_Inner_Px);
    	chain->SetBranchAddress("Muon_Inner_Py", &Muon_Inner_Py);
    	chain->SetBranchAddress("Muon_Inner_Pz", &Muon_Inner_Pz);
    	chain->SetBranchAddress("Muon_Inner_eta", &Muon_Inner_eta);
    	chain->SetBranchAddress("Muon_Inner_phi", &Muon_Inner_phi);

    	chain->SetBranchAddress("Muon_Outer_pT", &Muon_Outer_pT);
    	chain->SetBranchAddress("Muon_Outer_pTError", &Muon_Outer_pTError);
    	chain->SetBranchAddress("Muon_Outer_Px", &Muon_Outer_Px);
    	chain->SetBranchAddress("Muon_Outer_Py", &Muon_Outer_Py);
    	chain->SetBranchAddress("Muon_Outer_Pz", &Muon_Outer_Pz);
    	chain->SetBranchAddress("Muon_Outer_eta", &Muon_Outer_eta);
    	chain->SetBranchAddress("Muon_Outer_phi", &Muon_Outer_phi);

    	chain->SetBranchAddress("Muon_GLB_pT", &Muon_GLB_pT);
    	chain->SetBranchAddress("Muon_GLB_pTError", &Muon_GLB_pTError);
    	chain->SetBranchAddress("Muon_GLB_Px", &Muon_GLB_Px);
    	chain->SetBranchAddress("Muon_GLB_Py", &Muon_GLB_Py);
    	chain->SetBranchAddress("Muon_GLB_Pz", &Muon_GLB_Pz);
    	chain->SetBranchAddress("Muon_GLB_eta", &Muon_GLB_eta);
    	chain->SetBranchAddress("Muon_GLB_phi", &Muon_GLB_phi);

    	chain->SetBranchAddress("Muon_TuneP_pT", &Muon_TuneP_pT);
    	chain->SetBranchAddress("Muon_TuneP_pTError", &Muon_TuneP_pTError);
    	chain->SetBranchAddress("Muon_TuneP_Px", &Muon_TuneP_Px);
    	chain->SetBranchAddress("Muon_TuneP_Py", &Muon_TuneP_Py);
    	chain->SetBranchAddress("Muon_TuneP_Pz", &Muon_TuneP_Pz);
    	chain->SetBranchAddress("Muon_TuneP_eta", &Muon_TuneP_eta);
    	chain->SetBranchAddress("Muon_TuneP_phi", &Muon_TuneP_phi);
    }

    void TurnOnBranches_Electron()
    {
        chain->SetBranchStatus("Nelectrons", 1);
        chain->SetBranchStatus("Electron_Energy", 1);
        chain->SetBranchStatus("Electron_pT", 1);
        chain->SetBranchStatus("Electron_Px", 1);
        chain->SetBranchStatus("Electron_Py", 1);
        chain->SetBranchStatus("Electron_Pz", 1);
        chain->SetBranchStatus("Electron_eta", 1);
        chain->SetBranchStatus("Electron_phi", 1);
        chain->SetBranchStatus("Electron_charge", 1);
        chain->SetBranchStatus("Electron_gsfpT", 1);
        chain->SetBranchStatus("Electron_gsfPx", 1);
        chain->SetBranchStatus("Electron_gsfPy", 1);
        chain->SetBranchStatus("Electron_gsfPz", 1);
        chain->SetBranchStatus("Electron_gsfEta", 1);
        chain->SetBranchStatus("Electron_gsfPhi", 1);
        chain->SetBranchStatus("Electron_gsfCharge", 1);
        chain->SetBranchStatus("Electron_etaSC", 1);
        chain->SetBranchStatus("Electron_phiSC", 1);
        chain->SetBranchStatus("Electron_etaWidth", 1);
        chain->SetBranchStatus("Electron_phiWidth", 1);
        chain->SetBranchStatus("Electron_dEtaIn", 1);
        chain->SetBranchStatus("Electron_dEtaInSeed", 1); // updated at 19 Jul. 2017 by Dalmin
        chain->SetBranchStatus("Electron_dPhiIn", 1);
        chain->SetBranchStatus("Electron_sigmaIEtaIEta", 1);
        chain->SetBranchStatus("Electron_Full5x5_SigmaIEtaIEta", 1); // updated at 19 Jul. 2017 by Dalmin
        chain->SetBranchStatus("Electron_HoverE", 1);
        chain->SetBranchStatus("Electron_fbrem", 1);
        chain->SetBranchStatus("Electron_eOverP", 1);
        chain->SetBranchStatus("Electron_InvEminusInvP", 1);
        chain->SetBranchStatus("Electron_dxyVTX", 1);
        chain->SetBranchStatus("Electron_dzVTX", 1);
        chain->SetBranchStatus("Electron_dxy", 1);
        chain->SetBranchStatus("Electron_dz", 1);
        chain->SetBranchStatus("Electron_dxyBS", 1);
        chain->SetBranchStatus("Electron_dzBS", 1);
        chain->SetBranchStatus("Electron_chIso03", 1);
        chain->SetBranchStatus("Electron_nhIso03", 1);
        chain->SetBranchStatus("Electron_phIso03", 1);
        chain->SetBranchStatus("Electron_ChIso03FromPU", 1);

        chain->SetBranchStatus("Electron_mHits", 1);
        chain->SetBranchStatus("Electron_EnergySC", 1);
        chain->SetBranchStatus("Electron_preEnergySC", 1);
        chain->SetBranchStatus("Electron_rawEnergySC", 1);
        chain->SetBranchStatus("Electron_etSC", 1);
        chain->SetBranchStatus("Electron_E15", 1);
        chain->SetBranchStatus("Electron_E25", 1);
        chain->SetBranchStatus("Electron_E55", 1);
        chain->SetBranchStatus("Electron_RelPFIso_dBeta", 1);
        chain->SetBranchStatus("Electron_RelPFIso_Rho", 1);
        chain->SetBranchStatus("Electron_r9", 1);
        chain->SetBranchStatus("Electron_ecalDriven", 1);
        chain->SetBranchStatus("Electron_passConvVeto", 1);
        chain->SetBranchStatus("Electron_passMediumID", 1); // updated at 27 Apr. 2017 by Dalmin

    	// -- Uncorrected Electrons -- //
    	chain->SetBranchStatus("nUnCorrElectron", 1);
    	chain->SetBranchStatus("Electron_pTUnCorr", 1);
    	chain->SetBranchStatus("Electron_etaUnCorr", 1);
    	chain->SetBranchStatus("Electron_phiUnCorr", 1);
    	chain->SetBranchStatus("Electron_PxUnCorr", 1);
    	chain->SetBranchStatus("Electron_PyUnCorr", 1);
    	chain->SetBranchStatus("Electron_PzUnCorr", 1);
    	chain->SetBranchStatus("Electron_EnergyUnCorr", 1);
    	chain->SetBranchStatus("Electron_EnergySCUnCorr", 1);
    	chain->SetBranchStatus("Electron_etaSCUnCorr", 1);
    	chain->SetBranchStatus("Electron_phiSCUnCorr", 1);
    	chain->SetBranchStatus("Electron_etSCUnCorr", 1);

    	chain->SetBranchAddress("Nelectrons", &Nelectrons);
    	chain->SetBranchAddress("Electron_Energy", &Electron_Energy);
    	chain->SetBranchAddress("Electron_pT", &Electron_pT);
    	chain->SetBranchAddress("Electron_eta", &Electron_eta);
    	chain->SetBranchAddress("Electron_phi", &Electron_phi);
    	chain->SetBranchAddress("Electron_charge", &Electron_charge);
    	chain->SetBranchAddress("Electron_gsfpT", &Electron_gsfpT);
    	chain->SetBranchAddress("Electron_gsfPx", &Electron_gsfPx);
    	chain->SetBranchAddress("Electron_gsfPy", &Electron_gsfPy);
    	chain->SetBranchAddress("Electron_gsfPz", &Electron_gsfPz);
    	chain->SetBranchAddress("Electron_gsfEta", &Electron_gsfEta);
    	chain->SetBranchAddress("Electron_gsfPhi", &Electron_gsfPhi);
    	chain->SetBranchAddress("Electron_gsfCharge", &Electron_gsfCharge);
    	chain->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
    	chain->SetBranchAddress("Electron_phiSC", &Electron_phiSC);
    	chain->SetBranchAddress("Electron_etaWidth", &Electron_etaWidth);
    	chain->SetBranchAddress("Electron_phiWidth", &Electron_phiWidth);
    	chain->SetBranchAddress("Electron_dEtaIn", &Electron_dEtaIn);
    	chain->SetBranchAddress("Electron_dEtaInSeed", &Electron_dEtaInSeed); // updated at 19 Jul. 2017 by Dalmin
    	chain->SetBranchAddress("Electron_dPhiIn", &Electron_dPhiIn);
    	chain->SetBranchAddress("Electron_sigmaIEtaIEta", &Electron_sigmaIEtaIEta);
    	chain->SetBranchAddress("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta); // updated at 19 Jul. 2017 by Dalmin
    	chain->SetBranchAddress("Electron_HoverE", &Electron_HoverE);
    	chain->SetBranchAddress("Electron_fbrem", &Electron_fbrem);
    	chain->SetBranchAddress("Electron_eOverP", &Electron_eOverP);
    	chain->SetBranchAddress("Electron_InvEminusInvP", &Electron_InvEminusInvP);
    	chain->SetBranchAddress("Electron_dxyVTX", &Electron_dxyVTX);
    	chain->SetBranchAddress("Electron_dzVTX", &Electron_dzVTX);
    	chain->SetBranchAddress("Electron_dxy", &Electron_dxy);
    	chain->SetBranchAddress("Electron_dz", &Electron_dz);
    	chain->SetBranchAddress("Electron_dxyBS", &Electron_dxyBS);
    	chain->SetBranchAddress("Electron_dzBS", &Electron_dzBS);
    	chain->SetBranchAddress("Electron_chIso03", &Electron_chIso03);
    	chain->SetBranchAddress("Electron_nhIso03", &Electron_nhIso03);
    	chain->SetBranchAddress("Electron_phIso03", &Electron_phIso03);
    	chain->SetBranchAddress("Electron_ChIso03FromPU", &Electron_ChIso03FromPU);

    	chain->SetBranchAddress("Electron_mHits", &Electron_mHits);
    	chain->SetBranchAddress("Electron_EnergySC", &Electron_EnergySC);
    	chain->SetBranchAddress("Electron_preEnergySC", &Electron_preEnergySC);
    	chain->SetBranchAddress("Electron_rawEnergySC", &Electron_rawEnergySC);
    	chain->SetBranchAddress("Electron_etSC", &Electron_etSC);
    	chain->SetBranchAddress("Electron_E15", &Electron_E15);
    	chain->SetBranchAddress("Electron_E25", &Electron_E25);
    	chain->SetBranchAddress("Electron_E55", &Electron_E55);
    	chain->SetBranchAddress("Electron_RelPFIso_dBeta", &Electron_RelPFIso_dBeta);
    	chain->SetBranchAddress("Electron_RelPFIso_Rho", &Electron_RelPFIso_Rho);
    	chain->SetBranchAddress("Electron_r9", &Electron_r9);
    	chain->SetBranchAddress("Electron_ecalDriven", &Electron_ecalDriven);
        chain->SetBranchAddress("Electron_passConvVeto", &Electron_passConvVeto);
        chain->SetBranchAddress("Electron_passMediumID", &Electron_passMediumID); // updated at 27 Apr. 2017 by Dalmin

    	// -- Uncorrected Electrons -- //
    	chain->SetBranchAddress("nUnCorrElectron", &nUnCorrElectron);
    	chain->SetBranchAddress("Electron_pTUnCorr", &Electron_pTUnCorr);
    	chain->SetBranchAddress("Electron_etaUnCorr", &Electron_etaUnCorr);
    	chain->SetBranchAddress("Electron_phiUnCorr", &Electron_phiUnCorr);
    	chain->SetBranchAddress("Electron_PxUnCorr", &Electron_PxUnCorr);
    	chain->SetBranchAddress("Electron_PyUnCorr", &Electron_PyUnCorr);
    	chain->SetBranchAddress("Electron_PzUnCorr", &Electron_PzUnCorr);
    	chain->SetBranchAddress("Electron_EnergyUnCorr", &Electron_EnergyUnCorr);
    	chain->SetBranchAddress("Electron_EnergySCUnCorr", &Electron_EnergySCUnCorr);
    	chain->SetBranchAddress("Electron_etaSCUnCorr", &Electron_etaSCUnCorr);
    	chain->SetBranchAddress("Electron_phiSCUnCorr", &Electron_phiSCUnCorr);
    	chain->SetBranchAddress("Electron_etSCUnCorr", &Electron_etSCUnCorr);

    }

    void TurnOnBranches_Jet()
    {
    	chain->SetBranchStatus("Njets", 1);
        chain->SetBranchStatus("Jet_pT", 1);
        chain->SetBranchStatus("Jet_eta", 1);
        chain->SetBranchStatus("Jet_phi", 1);
        chain->SetBranchStatus("Jet_Charge", 1);
        chain->SetBranchStatus("Jet_Flavor", 1);

        chain->SetBranchStatus("Jet_bTag", 1);
        chain->SetBranchStatus("Jet_CHfrac", 1);
        chain->SetBranchStatus("Jet_NHfrac", 1);
        chain->SetBranchStatus("Jet_NHEMfrac", 1);
        chain->SetBranchStatus("Jet_CHEMfrac", 1);
        chain->SetBranchStatus("Jet_CHmulti", 1);
        chain->SetBranchStatus("Jet_NHmulti", 1);

    	chain->SetBranchAddress("Njets", &Njets);

    	chain->SetBranchAddress("Jet_pT", &Jet_pT);
    	chain->SetBranchAddress("Jet_eta", &Jet_eta);
    	chain->SetBranchAddress("Jet_phi", &Jet_phi);
    	chain->SetBranchAddress("Jet_Charge", &Jet_Charge);
    	chain->SetBranchAddress("Jet_Flavor", &Jet_Flavor);

    	chain->SetBranchAddress("Jet_bTag", &Jet_bTag);
    	chain->SetBranchAddress("Jet_CHfrac", &Jet_CHfrac);
    	chain->SetBranchAddress("Jet_NHfrac", &Jet_NHfrac);
    	chain->SetBranchAddress("Jet_NHEMfrac", &Jet_NHEMfrac);
    	chain->SetBranchAddress("Jet_CHEMfrac", &Jet_CHEMfrac);
    	chain->SetBranchAddress("Jet_CHmulti", &Jet_CHmulti);
    	chain->SetBranchAddress("Jet_NHmulti", &Jet_NHmulti);
    }

    void TurnOnBranches_Photon()
    {
    	chain->SetBranchStatus("nPhotons", 1);
        chain->SetBranchStatus("Photon_hasPixelSeed", 1);
        chain->SetBranchStatus("Photon_pT", 1);
        chain->SetBranchStatus("Photon_eta", 1);
        chain->SetBranchStatus("Photon_phi", 1);
        chain->SetBranchStatus("Photon_etaSC", 1);
        chain->SetBranchStatus("Photon_phiSC", 1);
        chain->SetBranchStatus("Photon_HoverE", 1);
        chain->SetBranchStatus("Photon_Full5x5_SigmaIEtaIEta", 1);
        chain->SetBranchStatus("Photon_ChIso", 1);
        chain->SetBranchStatus("Photon_NhIso", 1);
        chain->SetBranchStatus("Photon_PhIso", 1);
        chain->SetBranchStatus("Photon_ChIsoWithEA", 1);
        chain->SetBranchStatus("Photon_NhIsoWithEA", 1);
        chain->SetBranchStatus("Photon_PhIsoWithEA", 1);

    	chain->SetBranchAddress("nPhotons",&nPhotons);
    	chain->SetBranchAddress("Photon_hasPixelSeed",&Photon_hasPixelSeed);
    	chain->SetBranchAddress("Photon_pT",&Photon_pT);
    	chain->SetBranchAddress("Photon_eta",&Photon_eta);
    	chain->SetBranchAddress("Photon_phi",&Photon_phi);
    	chain->SetBranchAddress("Photon_etaSC",&Photon_etaSC);
    	chain->SetBranchAddress("Photon_phiSC",&Photon_phiSC);
    	chain->SetBranchAddress("Photon_HoverE",&Photon_HoverE);
    	chain->SetBranchAddress("Photon_Full5x5_SigmaIEtaIEta",&Photon_Full5x5_SigmaIEtaIEta);
    	chain->SetBranchAddress("Photon_ChIso",&Photon_ChIso);
    	chain->SetBranchAddress("Photon_NhIso",&Photon_NhIso);
    	chain->SetBranchAddress("Photon_PhIso",&Photon_PhIso);
    	chain->SetBranchAddress("Photon_ChIsoWithEA",&Photon_ChIsoWithEA);
    	chain->SetBranchAddress("Photon_NhIsoWithEA",&Photon_NhIsoWithEA);
    	chain->SetBranchAddress("Photon_PhIsoWithEA",&Photon_PhIsoWithEA);
    }

    void TurnOnBranches_MET()
    {
        chain->SetBranchStatus("pfMET_pT", 1);
        chain->SetBranchStatus("pfMET_phi", 1);
        chain->SetBranchStatus("pfMET_Px", 1);
        chain->SetBranchStatus("pfMET_Py", 1);
        chain->SetBranchStatus("pfMET_SumEt", 1);

        chain->SetBranchStatus("pfMET_Type1_pT", 1);
        chain->SetBranchStatus("pfMET_Type1_phi", 1);
        chain->SetBranchStatus("pfMET_Type1_Px", 1);
        chain->SetBranchStatus("pfMET_Type1_Py", 1);
        chain->SetBranchStatus("pfMET_Type1_SumEt", 1);
        
    	chain->SetBranchAddress("pfMET_pT", &pfMET_pT);
    	chain->SetBranchAddress("pfMET_phi", &pfMET_phi);
    	chain->SetBranchAddress("pfMET_Px", &pfMET_Px);
    	chain->SetBranchAddress("pfMET_Py", &pfMET_Py);
    	chain->SetBranchAddress("pfMET_SumEt", &pfMET_SumEt);

    	chain->SetBranchAddress("pfMET_Type1_pT", &pfMET_Type1_pT);
    	chain->SetBranchAddress("pfMET_Type1_phi", &pfMET_Type1_phi);
    	chain->SetBranchAddress("pfMET_Type1_Px", &pfMET_Type1_Px);
    	chain->SetBranchAddress("pfMET_Type1_Py", &pfMET_Type1_Py);
    	chain->SetBranchAddress("pfMET_Type1_SumEt", &pfMET_Type1_SumEt);
    }

    void GetEvent(Int_t i)
    {
        if(!chain) return;
        
        chain->GetEntry(i);
    }

    Bool_t isTriggered(TString HLT)
    {
        Bool_t isTrigger = false;
        if( HLT == "HLT_IsoMu20_v* || HLT_IsoTkMu20_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_IsoMu20_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_IsoTkMu20_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else if( HLT == "HLT_IsoMu24_v* || HLT_IsoTkMu24_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_IsoMu24_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_IsoTkMu24_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else if( HLT == "HLT_Mu50_v* || HLT_TkMu50_v*" )
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == "HLT_Mu50_v*" || (HLT_trigName->at((unsigned int)k)) == "HLT_TkMu50_v*" )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }
        else
        {
            for( Int_t k = 0; k < HLT_ntrig; k++ )
            {
                if( (HLT_trigName->at((unsigned int)k)) == HLT )
                {
                    if( HLT_trigFired[k] == 1 )
                    {
                        isTrigger = true;
                        break;
                    }
                }
            }
        }

        return isTrigger;
    }

};
