TString variables(Double_t i_var)
{
	vector< pair<Double_t, TString> > var;

	var.push_back(make_pair(0,"mass"));
	var.push_back(make_pair(1,"mass_fine"));
	var.push_back(make_pair(2,"pT"));
	var.push_back(make_pair(2.1,"leadPt"));
	var.push_back(make_pair(2.2,"subPt"));
	var.push_back(make_pair(3,"eta"));
	var.push_back(make_pair(3.1,"leadEta"));
	var.push_back(make_pair(3.2,"subEta"));
	var.push_back(make_pair(3.3,"etaSC"));
	var.push_back(make_pair(4,"phi"));
	var.push_back(make_pair(4.1,"leadPhi"));
	var.push_back(make_pair(4.2,"subPhi"));
	var.push_back(make_pair(5,"diPt"));
	var.push_back(make_pair(6,"rapi"));
	var.push_back(make_pair(7,"nVTX_noPU"));
	var.push_back(make_pair(7.1,"nVTX_withPU"));
	var.push_back(make_pair(8,"isGLB"));
	var.push_back(make_pair(9,"muonHits"));
	var.push_back(make_pair(10,"nMatches"));
	var.push_back(make_pair(11,"dpT_over_pT"));
	var.push_back(make_pair(12,"dxyVTX"));
	var.push_back(make_pair(13,"dzVTX"));
	var.push_back(make_pair(14,"pixelHits"));
	var.push_back(make_pair(15,"trackerLayers"));
	var.push_back(make_pair(16,"VtxNormChi2"));
	var.push_back(make_pair(17,"3DAngle"));
	var.push_back(make_pair(18,"relTrkIso"));
	var.push_back(make_pair(19,"RelPFIso_dBeta"));
	var.push_back(make_pair(20,"Full5x5_SigmaIEtaIEta"));
	var.push_back(make_pair(21,"dEtaInSeed"));
	var.push_back(make_pair(22,"dPhiIn"));
	var.push_back(make_pair(23,"HoverE"));
	var.push_back(make_pair(24,"RelPFIso_Rho"));
	var.push_back(make_pair(25,"InvEminusInvP"));
	var.push_back(make_pair(26,"mHits"));
	var.push_back(make_pair(27,"passConvVeto"));

	TString variable = "null";
	for(Int_t i=0; i<var.size(); i++)
		if( var[i].first == i_var )
			variable = var[i].second;

	// -- Gen-level -- //
	//variable = "Gen_" + variable;

	// -- Barrel and End-cap -- //
	//variable = variable + "_BB";
	//variable = variable + "_BE";
	//variable = variable + "_EE";

	// -- Corrections -- //
	//variable = "before_PUCorr_" + variable;
	//variable = "before_RoccoR_" + variable;
	//variable = "before_EGMCorr_" + variable;
	//variable = "BEC_" + variable;
	//variable = "Reco_" + variable;
	//variable = "RecoID_" + variable;
	//variable = "RecoIDIso_" + variable;

	// -- Cut based -- //
	//variable = "PtCut_" + variable;
	//variable = variable + "_M20to30";
	//variable = variable + "_M30to45";
	//variable = variable + "_M60to120";
	//variable = variable + "_M120to200";
	//variable = variable + "_M200to1500";
	//variable = "M20to30_" + variable;
	//variable = "M30to45_" + variable;
	//variable = "M60to120_" + variable;
	//variable = "M120to200_" + variable;
	//variable = "M200to1500_" + variable;

	return 	variable;
}
