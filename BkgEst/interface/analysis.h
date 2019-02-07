#ifndef analysis_H
#define analysis_H

#include <TLorentzVector.h>
#include "NtupleMaker.h"

const double muon_mass     = 0.1056583715;
const double electron_mass = 0.000510998;

//const int binnum = 45;
//const double bins[46] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1200, 1500, 2000, 3000};

class PhysicsEvent : public NtupleEvent {
public:
	bool MuonSelection() {
		bool flag1 = false;
		bool flag2 = false;
		for( std::vector<NtupleGenParticle>::const_iterator par = genparticles.begin(); par != genparticles.end(); ++par ) {
			if( par->isHardProcess && par->id == 13 ) flag1 = true;
			else if( par->isHardProcess && par->id == -13 ) flag2 = true;
		}
		return flag1&&flag2;
	};
	bool TauSelection() {
		bool flag1 = false;
		bool flag2 = false;
		for( std::vector<NtupleGenParticle>::const_iterator par = genparticles.begin(); par != genparticles.end(); ++par ) {
			if( par->isHardProcess && par->id == 15 ) flag1 = true;
			else if( par->isHardProcess && par->id == -15 ) flag2 = true;
		}
		return flag1&&flag2;
	};
	bool TriggerSelection( std::string trigger_ = "HLT_IsoMu20_v") {
		bool isTriggered = false;
		for( std::vector<NtupleTrigger>::const_iterator trigger = triggers.begin(); trigger != triggers.end(); ++trigger ) {
			if( trigger->name.find(trigger_) == 0 ) {
				//cout<<trigger_<<", "<<trigger->name<<endl;
				//cout<<trigger->isFired<<endl;
				isTriggered = trigger->isFired;
				break;
			}
		}
		return isTriggered;
	}
    bool GenMass() {

        int check1 = 0;
        int check2 = 0;
        TLorentzVector mom1;
        TLorentzVector mom2;

        for(vector<NtupleGenParticle>::const_iterator gen = genparticles.begin(); gen != genparticles.end(); ++gen) {
            if( gen->id==6 && gen->isHardProcess ) {
                mom1.SetPtEtaPhiE(gen->pt,gen->eta,gen->phi,gen->energy);
            }
            else if( gen->id==-6 && gen->isHardProcess ) {
                mom2.SetPtEtaPhiE(gen->pt,gen->eta,gen->phi,gen->energy);
            }
        }

        double mass = (mom1+mom2).M();
        //cout<<"Mtt = "<<mass<<endl;
        if( mass > 700 ) return true;
        else return false;
    }
};

class PhysicsMuon : public NtupleMuon {
public:
	bool acceptance(double pt_, double eta_){
		if( pt > pt_ && fabs(eta) < eta_ ) return true;
		else return false;
	}
	TLorentzVector momentum(){
		TLorentzVector momentum_;
		momentum_.SetPtEtaPhiM(pt,eta,phi,muon_mass);
		return momentum_;
	}

  bool looseMuonID() {
    if( isPFMuon && (isGlobalMuon||isTrackerMuon) && fabs(dxyVTX) < 0.2 && fabs(dzVTX) < 0.5 )
      return true;
    else
      return false;
  }

  bool tightMuonID() {
    if( isPFMuon && isGlobalMuon && normalizedChi2 < 10 && nValidMuonHits>0 && nMatchedStations>1 && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX) < 0.2 && fabs(dzVTX) < 0.5 )
      return true;
    else
      return false;
  }

	bool highPtMuonID() {
		if( isGlobalMuon && nValidMuonHits>0 && nMatchedStations>1 && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX)<0.2 && muonBestTrack_ptError/muonBestTrack_pt<0.3 )
			return true;
		else
			return false;
	}
	bool isolation(double iso) {
		if( isolationR03_sumpt/pt < iso) return true;
		else return false;
	}
	bool miniIsolation(double iso) {
		if( miniIso < iso ) return true;
		else return false;
	}
};

class PhysicsElectron : public NtupleElectron {
public:
	bool acceptance(double pt_, double eta_){
		if( pt > pt_ && fabs(eta) < eta_  && (fabs(eta) < 1.4442 || fabs(eta) > 1.566) ) return true;
		else return false;
	}
	TLorentzVector momentum(){
		TLorentzVector momentum_;
		momentum_.SetPtEtaPhiM(pt,eta,phi,electron_mass);
		return momentum_;
	}	

  bool WPLoose() {
    bool flag = false;
    if(fabs(etaSC)<=1.479) { //Barrel
      if( sigmaIetaIeta < 0.0103
        && fabs(dEtaIn) < 0.0105
        && fabs(dPhiIn) < 0.115
        && hOverE < 0.104
        //&& isoRho < 0.0893
        && ooEmooP < 0.102
        && fabs(d0) < 0.0261
        && fabs(dz) < 0.41
        && expectedMissingInnerHits <= 2
        && passConversionVeto )
        flag = true;
    }
    else { //Endcap
      if( sigmaIetaIeta < 0.0301
        && fabs(dEtaIn) < 0.00814
        && fabs(dPhiIn) < 0.182
        && hOverE < 0.0897
        //&& isoRho < 0.121
        && ooEmooP < 0.126
        && fabs(d0) < 0.118
        && fabs(dz) < 0.822
        && expectedMissingInnerHits <= 1
        && passConversionVeto )
        flag = true;
    }
    return flag;
  }

	bool WPMedium() {
		bool flag = false;
		if(fabs(etaSC)<=1.479) { //Barrel
			if( sigmaIetaIeta < 0.0101 
				&& fabs(dEtaIn) < 0.0103 
				&& fabs(dPhiIn) < 0.0336 
				&& hOverE < 0.0876 
				&& isoRho < 0.0766 
				&& ooEmooP < 0.0174 
				&& fabs(d0) < 0.0118 
				&& fabs(dz) < 0.373 
				&& expectedMissingInnerHits <= 2 
				&& passConversionVeto ) 
				flag = true;
		}
		else { //Endcap
			if( sigmaIetaIeta < 0.0283 
				&& fabs(dEtaIn) < 0.00733 
				&& fabs(dPhiIn) < 0.114  
				&& hOverE < 0.0678 
				&& isoRho < 0.0678 
				&& ooEmooP < 0.0898 
				&& fabs(d0) < 0.0739 
				&& fabs(dz) < 0.602 
				&& expectedMissingInnerHits <= 1 
				&& passConversionVeto ) 
				flag = true;
		}
		return flag;
	} 

	bool WPTight() {
		bool flag = false;
		if(fabs(etaSC)<=1.479) { //Barrel
			if( sigmaIetaIeta < 0.0101 
				&& fabs(dEtaIn) < 0.00926 
				&& fabs(dPhiIn) < 0.0336 
				&& hOverE < 0.0597 
				&& isoRho < 0.0354 
				&& ooEmooP < 0.012 
				&& fabs(d0) < 0.0111 
				&& fabs(dz) < 0.0466 
				&& expectedMissingInnerHits <= 2 
				&& passConversionVeto ) 
				flag = true;
		}
		else { //Endcap
			if( sigmaIetaIeta < 0.0279 
				&& fabs(dEtaIn) < 0.00724 
				&& fabs(dPhiIn) < 0.0918 
				&& hOverE < 0.0615 
				&& isoRho < 0.0646 
				&& ooEmooP < 0.00999 
				&& fabs(d0) < 0.0351 
				&& fabs(dz) < 0.417 
				&& expectedMissingInnerHits <= 1 
				&& passConversionVeto ) 
				flag = true;
		}
		return flag;
	} 
	bool miniIsolation(double iso) {
		if( miniIso < iso ) return true;
		else return false;
	}
};

class PhysicsPhoton : public NtuplePhoton {
public:
	bool acceptance(double pt_, double eta_) {
		if( pt > pt_ && fabs(eta) < eta_  && (fabs(eta) < 1.4442 || fabs(eta) > 1.566) ) return true;
		else return false;
	}
	TLorentzVector momentum(){
		TLorentzVector momentum_;
		momentum_.SetPtEtaPhiM(pt,eta,phi,0);
		return momentum_;
	}
	
  bool loose() {
    if( fabs(eta) < 1.479
      && HoverE < 0.05
      && Full5x5_SigmaIEtaIEta < 0.0102
      && ChIso < 2.5
      && !hasPixelSeed )
      return true;
    else if( fabs(eta) > 1.479
      && HoverE < 0.05
      && Full5x5_SigmaIEtaIEta < 0.0274
      && ChIso < 2.5 
      && !hasPixelSeed )
      return true;
    else return false;
  }

	bool medium() {
		if( fabs(etaSC) < 1.479 
			&& HoverE < 0.05 
			&& Full5x5_SigmaIEtaIEta < 0.0102 
			&& ChIso < 1.37 
			&& NhIsoWithEA < 1.06 + 0.014 * pt + 0.000019 * pt * pt 
			&& PhIsoWithEA < 0.28 + 0.0053 * pt ) 
			return true;
		else if( fabs(etaSC) > 1.479 
			&& HoverE < 0.05 
			&& Full5x5_SigmaIEtaIEta < 0.0268 
			&& ChIso < 1.10 
			&& NhIsoWithEA < 2.69 + 0.0139 * pt + 0.000025 * pt * pt 
			&& PhIsoWithEA < 0.39 + 0.0034 * pt ) 
			return true;
		else return false;
	}
};

double deltaR(double eta1, double phi1, double eta2, double phi2) {
	double deta = fabs(eta1 - eta2);
	double dphi = fabs(phi1 - phi2);
	if(dphi>M_PI)
		dphi = 2*M_PI - dphi;
	return sqrt(deta*deta+dphi*dphi);
}

bool triggerMatch(std::vector<NtupleTriggerObject> triggerobjects, PhysicsMuon mu, std::string triggerName = "HLT_IsoMu20_v") {
	bool flag = false;
	for( std::vector<NtupleTriggerObject>::const_iterator trig = triggerobjects.begin(); trig != triggerobjects.end(); ++trig ) {
		if( trig->name.find(triggerName)!=std::string::npos) {
			double dR = deltaR(mu.eta, mu.phi, trig->eta, trig->phi);
			if(dR<0.2) {
				flag = true;
				break;
			}
		}
	}
	return flag;
}  

double openingAngle(std::vector<NtupleDimuon> dimuons, std::pair<PhysicsMuon,int> mu1, std::pair<PhysicsMuon,int> mu2) {
    bool fail = true;
	double angle = 9999;
	int index1, index2;
	if( mu1.second < mu2.second ) {
		index1 = mu1.second;
		index2 = mu2.second;
	}
	else {
		index1 = mu2.second;
		index2 = mu1.second;
	}
	for( unsigned i=0; i!=dimuons.size(); i++){
		if( dimuons.at(i).X==index1 && dimuons.at(i).Y==index2 ) {
            fail = false;
			angle = dimuons.at(i).openingAngle;
			break;
		}
	}
    if(fail) {
        cout<<"Wrong angle:"<<index1<<", "<<index2<<endl;
        for( unsigned i=0; i!=dimuons.size(); ++i) {
            cout<<dimuons.at(i).X<<", "<<dimuons.at(i).Y<<endl;
        }
        cout<<""<<endl;
    }
	return angle;
}

double openingAngle(vector<NtupleEmu> emus, std::pair<PhysicsElectron,int> el, std::pair<PhysicsMuon,int> mu) {
	double angle = 9999;
	int index1, index2;
	index1 = mu.second;
	index2 = el.second;

	for(unsigned i=0; i!=emus.size(); i++){
		if(emus.at(i).X==index1 && emus.at(i).Y==index2) {
			angle = emus.at(i).openingAngle;
			break;
		}
	}
	return angle;
}

double vertexFitChi2(std::vector<NtupleDimuon> dimuons, std::pair<PhysicsMuon,int> mu1, std::pair<PhysicsMuon,int> mu2) {
	double chi2 = 9999;
    bool fail = true;
	int index1, index2;
	if( mu1.second < mu2.second ) {
		index1 = mu1.second;
		index2 = mu2.second;
	}
	else {
		index1 = mu2.second;
		index2 = mu1.second;
	}
	for( unsigned i=0; i!=dimuons.size(); i++){
		if( dimuons.at(i).X==index1 && dimuons.at(i).Y==index2 ) {
			chi2 = dimuons.at(i).vertexFitChi2Ndof;
            fail = false;
			break;
		}
	}
    if(fail) cout<<"Wrong"<<endl;
	//cout<<"Dimuon chi2 = "<<chi2<<endl;
	return chi2;
}

double vertexFitChi2(std::vector<NtupleEmu> emus, std::pair<PhysicsElectron,int> el, std::pair<PhysicsMuon,int> mu) {
  double chi2 = 9999;
  int index1, index2;
  index1 = mu.second;
  index2 = el.second;

  for(unsigned i=0; i!=emus.size(); i++){
    if( emus.at(i).X==index1 && emus.at(i).Y==index2 ) {
      chi2 = emus.at(i).vertexFitChi2Ndof;
      break;
    }
  }
  //cout<<"Emu chi2 = "<<chi2<<endl;
  return chi2;
}

bool dimuonDY(std::vector<NtupleTriggerObject> triggerobjects, std::vector<NtupleDimuon> dimuons, vector<pair<PhysicsMuon,int>>* muons, pair<PhysicsMuon,PhysicsMuon>* dimuon) {
	bool flag = false;
	double min = 20;
	pair<PhysicsMuon,int> mu1;
	pair<PhysicsMuon,int> mu2;

	if(muons->size()>=2) {
		for(unsigned int i=0; i<muons->size(); i++) {
			for(unsigned int j=0; j<muons->size(); j++) {
				if(i>=j) continue;

				mu1 = muons->at(i);
				mu2 = muons->at(j);

				if( mu1.first.pt<=22 && mu2.first.pt<=22 ) continue;

				if(triggerMatch(triggerobjects, mu1.first,"HLT_IsoMu20_v")||triggerMatch(triggerobjects, mu1.first,"HLT_IsoTkMu20_v")||triggerMatch(triggerobjects, mu2.first,"HLT_IsoMu20_v")||triggerMatch(triggerobjects, mu2.first,"HLT_IsoTkMu20_v")) {
					if(openingAngle(dimuons,mu1,mu2)<M_PI-0.005) {
						if(vertexFitChi2(dimuons,mu1,mu2)<min) {
							flag = true;
							min = vertexFitChi2(dimuons,mu1,mu2);
							dimuon->first = mu1.first;
							dimuon->second = mu2.first;
						}
					}
				}
			}
		}
	}
	return flag;
}
bool emuDY(std::vector<NtupleTriggerObject> triggerobjects, std::vector<NtupleEmu> emus, vector<pair<PhysicsElectron,int>>* electrons, vector<pair<PhysicsMuon,int>>* muons, pair<PhysicsElectron,PhysicsMuon>* emu) {
	bool flag = false;
	double min = 20;
	pair<PhysicsMuon,int> mu;
	pair<PhysicsElectron,int> el;

	for(unsigned int i=0; i<muons->size(); i++) {
		for(unsigned int j=0; j<electrons->size(); j++) {
			mu = muons->at(i);
			el = electrons->at(j);

			if( mu.first.pt<=22 && el.first.etSC<=22 ) continue;

			if(!triggerMatch(triggerobjects, mu.first,"HLT_IsoMu20_v*")&&!triggerMatch(triggerobjects, mu.first,"HLT_IsoTkMu20_v*")) continue;

			if( openingAngle(emus,el,mu)<M_PI-0.005 ) {
				if( vertexFitChi2(emus,el,mu)<min ) {
					flag = true;
					min = vertexFitChi2(emus,el,mu);
					emu->first = el.first;
					emu->second = mu.first;
				}
			}
		}
	}
	return flag;
}

bool dijetDY(std::vector<NtupleDimuon> dimuons, vector<pair<PhysicsMuon,int>>* muons, pair<PhysicsMuon,PhysicsMuon>* dimuon) {
	bool flag = false;
	double min = 20;
	pair<PhysicsMuon,int> mu1;
	pair<PhysicsMuon,int> mu2;

	if(muons->size()>=2) {
		for(unsigned int i=0; i<muons->size(); i++) {
			for(unsigned int j=0; j<muons->size(); j++) {
				if(i>=j) continue;

				mu1 = muons->at(i);
				mu2 = muons->at(j);

				if( mu1.first.pt<=22 && mu2.first.pt<=22 ) continue;

				if(openingAngle(dimuons,mu1,mu2)<M_PI-0.005) {
					if(vertexFitChi2(dimuons,mu1,mu2)<min) {
						flag = true;
						min = vertexFitChi2(dimuons,mu1,mu2);
						dimuon->first = mu1.first;
						dimuon->second = mu2.first;
					}
				}
				
			}
		}
	}
	return flag;
}

bool wjetsDY(std::vector<NtupleDimuon> dimuons, pair<PhysicsMuon,int> mu1, pair<PhysicsMuon,int> mu2, pair<PhysicsMuon,PhysicsMuon>* dimuon) {
	bool flag = false;
	const double min = 20;


	if( mu1.first.pt>22 || mu2.first.pt>22 ) {
		if(openingAngle(dimuons,mu1,mu2)<M_PI-0.005) {
			if(vertexFitChi2(dimuons,mu1,mu2)<min) {
				flag = true;
				dimuon->first = mu1.first;
				dimuon->second = mu2.first;
			}
		}
	}

	return flag;
}

double FR_template(PhysicsMuon muon){
  double pT = muon.pt;
  double eta = muon.eta;
  double fakerate = -999;

  if(fabs(eta)<1.2) {
    double FR[] = {0.134681,0.13028,0.127178,0.127156,0.132708,0.133303,0.126991,0.142287,0.133668,0.151394,0.163657,0.187625,0.192752,0.242238,0.284034,0.307194,0.326136};
    //double FR[] = {0.130238,0.136145,0.132686,0.138025,0.121351,0.126273,0.142708,0.142825,0.149498,0.158084,0.146399,0.16534,0.244816,0.275748,0.24459,0.292898,0.639552};
    if(pT < 52) fakerate = FR[0];
    else if( pT > 52 && pT <60 ) fakerate = FR[1];
    else if( pT > 60 && pT <70 ) fakerate = FR[2];
    else if( pT > 70 && pT <80 ) fakerate = FR[3];
    else if( pT > 80 && pT <90 ) fakerate = FR[4];
    else if( pT > 90 && pT <100 ) fakerate = FR[5];
    else if( pT > 100 && pT <120 ) fakerate = FR[6];
    else if( pT > 120 && pT <140 ) fakerate = FR[7];
    else if( pT > 140 && pT <160 ) fakerate = FR[8];
    else if( pT > 160 && pT <180 ) fakerate = FR[9];
    else if( pT > 180 && pT <200 ) fakerate = FR[10];
    else if( pT > 200 && pT <250 ) fakerate = FR[11];
    else if( pT > 250 && pT <300 ) fakerate = FR[12];
    else if( pT > 300 && pT <350 ) fakerate = FR[13];
    else if( pT > 350 && pT <400 ) fakerate = FR[14];
    else if( pT > 400 && pT <450 ) fakerate = FR[15];
    else fakerate = FR[16];
  }
  else {
    double FR[] = {0.220145,0.229679,0.234881,0.233588,0.251451,0.2501,0.287229,0.291463,0.302949,0.264239,0.460734,0.359131,0.644995,0.929643,0.606793,0.606793,0.606793};
    //double FR[] = {0.207901,0.221861,0.210008,0.232894,0.221247,0.240958,0.249162,0.261152,0.324663,0.266702,0.327946,0.306969,0.322887,0.69845,0.658293,0.815957,0.686238};
    if(pT < 52) fakerate = FR[0];
    else if( pT > 52 && pT <60 ) fakerate = FR[1];
    else if( pT > 60 && pT <70 ) fakerate = FR[2];
    else if( pT > 70 && pT <80 ) fakerate = FR[3];
    else if( pT > 80 && pT <90 ) fakerate = FR[4];
    else if( pT > 90 && pT <100 ) fakerate = FR[5];
    else if( pT > 100 && pT <120 ) fakerate = FR[6];
    else if( pT > 120 && pT <140 ) fakerate = FR[7];
    else if( pT > 140 && pT <160 ) fakerate = FR[8];
    else if( pT > 160 && pT <180 ) fakerate = FR[9];
    else if( pT > 180 && pT <200 ) fakerate = FR[10];
    else if( pT > 200 && pT <250 ) fakerate = FR[11];
    else if( pT > 250 && pT <300 ) fakerate = FR[12];
    else if( pT > 300 && pT <350 ) fakerate = FR[13];
    else if( pT > 350 && pT <400 ) fakerate = FR[14];
    else if( pT > 400 && pT <450 ) fakerate = FR[15];
    else fakerate = FR[16];
  }
  return fakerate;
}

double FR_ratio(PhysicsMuon muon){
  double pT = muon.pt;
  double eta = muon.eta;
  double fakerate = -999;

  if(fabs(eta)<1.2) {
    double FR[] = {0.124991,0.119697,0.116249,0.11349,0.11714,0.115876,0.10989,0.120195,0.117077,0.136469,0.140139,0.159767,0.173791,0.2025,0.265724,0.289455,0.352381};
    //double FR[] = {0.117844,0.125069,0.12191,0.124546,0.105265,0.1117,0.124381,0.124537,0.12765,0.142944,0.124476,0.145098,0.197322,0.24474,0.283817,0.351978,0.322147};
    if(pT < 52) fakerate = FR[0];
    else if( pT > 52 && pT <60 ) fakerate = FR[1];
    else if( pT > 60 && pT <70 ) fakerate = FR[2];
    else if( pT > 70 && pT <80 ) fakerate = FR[3];
    else if( pT > 80 && pT <90 ) fakerate = FR[4];
    else if( pT > 90 && pT <100 ) fakerate = FR[5];
    else if( pT > 100 && pT <120 ) fakerate = FR[6];
    else if( pT > 120 && pT <140 ) fakerate = FR[7];
    else if( pT > 140 && pT <160 ) fakerate = FR[8];
    else if( pT > 160 && pT <180 ) fakerate = FR[9];
    else if( pT > 180 && pT <200 ) fakerate = FR[10];
    else if( pT > 200 && pT <250 ) fakerate = FR[11];
    else if( pT > 250 && pT <300 ) fakerate = FR[12];
    else if( pT > 300 && pT <350 ) fakerate = FR[13];
    else if( pT > 350 && pT <400 ) fakerate = FR[14];
    else if( pT > 400 && pT <450 ) fakerate = FR[15];
    else fakerate = FR[16];
  }
  else {
    double FR[] = {0.188425,0.200836,0.20233,0.196967,0.203623,0.198822,0.220961,0.223721,0.230471,0.199941,0.346051,0.260297,0.541437,0.885527,0.566429,0.389111,0.852681};
    //double FR[] = {0.173865,0.190552,0.179987,0.196405,0.18471,0.191826,0.193596,0.198557,0.242361,0.20703,0.256279,0.223002,0.250271,0.516487,0.456352,0.567045,0.683899};
    if(pT < 52) fakerate = FR[0];
    else if( pT > 52 && pT <60 ) fakerate = FR[1];
    else if( pT > 60 && pT <70 ) fakerate = FR[2];
    else if( pT > 70 && pT <80 ) fakerate = FR[3];
    else if( pT > 80 && pT <90 ) fakerate = FR[4];
    else if( pT > 90 && pT <100 ) fakerate = FR[5];
    else if( pT > 100 && pT <120 ) fakerate = FR[6];
    else if( pT > 120 && pT <140 ) fakerate = FR[7];
    else if( pT > 140 && pT <160 ) fakerate = FR[8];
    else if( pT > 160 && pT <180 ) fakerate = FR[9];
    else if( pT > 180 && pT <200 ) fakerate = FR[10];
    else if( pT > 200 && pT <250 ) fakerate = FR[11];
    else if( pT > 250 && pT <300 ) fakerate = FR[12];
    else if( pT > 300 && pT <350 ) fakerate = FR[13];
    else if( pT > 350 && pT <400 ) fakerate = FR[14];
    else if( pT > 400 && pT <450 ) fakerate = FR[15];
    else fakerate = FR[16];
  }
  return fakerate;
}

#endif
