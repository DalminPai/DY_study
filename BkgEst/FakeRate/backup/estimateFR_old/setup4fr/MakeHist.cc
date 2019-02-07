#include <TFile.h>
#include <TH1D.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>
using namespace std;

void MakeHist(TString rootfile, TString var)
{
	//////////////////////////
	// -- Get histograms -- //
	//////////////////////////
	TFile f_input(rootfile, "read");
	TH1D *h[100]; TH1D* h_[21];

	//Data histogram
	h[0] = (TH1D*)f_input.Get(var+"_Data");
	h_[0] = (TH1D*)h[0]->Clone();
	
  //MC histograms
	h[1] = (TH1D*)f_input.Get(var+"_DY_M10to50_v1");
	h[2] = (TH1D*)f_input.Get(var+"_DY_M10to50_v2");
	h[3] = (TH1D*)f_input.Get(var+"_DY_M10to50_ext1v1");
	h[1]->Add(h[2]);
	h[1]->Add(h[3]);
	h_[1] = (TH1D*)h[1]->Clone(var+"_DY_M10to50");

	h[4] = (TH1D*)f_input.Get(var+"_DY_M50toInf");
  h_[2] = (TH1D*)h[4]->Clone();
	
  h[5] = (TH1D*)f_input.Get(var+"_ttbar");
	h[6] = (TH1D*)f_input.Get(var+"_ttbarBackup");
  h[5]->Add(h[6]);
  h_[3] = (TH1D*)h[5]->Clone();
	
  h[7] = (TH1D*)f_input.Get(var+"_WW");
  h_[4] = (TH1D*)h[7]->Clone();

	h[8] = (TH1D*)f_input.Get(var+"_WZ");
  h_[5] = (TH1D*)h[8]->Clone();

	h[9] = (TH1D*)f_input.Get(var+"_ZZ");
  h_[6] = (TH1D*)h[9]->Clone();
	
  h[10] = (TH1D*)f_input.Get(var+"_tW");
	h[11] = (TH1D*)f_input.Get(var+"_tbarW");
  h[10]->Add(h[11]);
  h_[7] = (TH1D*)h[10]->Clone();
	
  h[12] = (TH1D*)f_input.Get(var+"_WJetsToLNu");
	h[13] = (TH1D*)f_input.Get(var+"_WJetsToLNu_ext");
  h[12]->Add(h[13]);
  h_[8] = (TH1D*)h[12]->Clone(var+"_WJets");
	
  h[14] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt15to20");
  h_[9] = (TH1D*)h[14]->Clone();
	
  h[15] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt20to30");
	h_[10] = (TH1D*)h[15]->Clone();

  h[16] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt30to50");
  h_[11] = (TH1D*)h[16]->Clone();
	
  h[17] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt50to80");
	h_[12] = (TH1D*)h[17]->Clone();

  h[18] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt80to120");
	h[19] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt80to120_ext1");
  h[18]->Add(h[19]);
  h_[13] = (TH1D*)h[18]->Clone();
	
  h[20] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt120to170");
	h[21] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt120to170_backup");
  h[20]->Add(h[21]);
  h_[14] = (TH1D*)h[20]->Clone();
	
  h[22] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt170to300");
	h[23] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt170to300_ext1");
	h[24] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt170to300_backup");
  h[22]->Add(h[23]);
  h[22]->Add(h[24]);
  h_[15] = (TH1D*)h[22]->Clone();
	
  h[25] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt300to470");
	h[26] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt300to470_ext1");
	h[27] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt300to470_ext2");
  h[25]->Add(h[26]);
  h[25]->Add(h[27]);
  h_[16] = (TH1D*)h[25]->Clone();
	
  h[28] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt470to600");
	h[29] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt470to600_ext1");
  h[28]->Add(h[29]);
  h_[17] = (TH1D*)h[28]->Clone();
	
  h[30] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt600to800");
	h[31] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt600to800_ext1");
	h[32] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt600to800_backup");
  h[30]->Add(h[31]);
  h[30]->Add(h[32]);
  h_[18] = (TH1D*)h[30]->Clone();
	
  h[33] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt800to1000");
	h[34] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt800to1000_ext1");
	h[35] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt800to1000_ext2");
  h[33]->Add(h[34]);
  h[33]->Add(h[35]);
  h_[19] = (TH1D*)h[33]->Clone();
	
  h[36] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt1000toInf");
	h[37] = (TH1D*)f_input.Get(var+"_QCDMuEnriched_Pt1000toInf_ext1");
  h[36]->Add(h[37]);
  h_[20] = (TH1D*)h[36]->Clone();

	/////////////////////////////////////////////////
	// -- Normalization: Scale(lumi*xsec/nEvts) -- //
	/////////////////////////////////////////////////
	//Xsec and nEvent
	vector< pair< TString, pair<Double_t, Double_t> > > norm;
	norm.push_back( make_pair( "DY_M10to50", make_pair( 6016.88*3, 11372172.0+23921165.0+14460175.0+10929538.0+24025168.0+14926245.0 ) ) );
	norm.push_back( make_pair( "DY_M50toInf", make_pair( 1952.68432327*3, 81780984.0 ) ) );
	norm.push_back( make_pair( "ttbar", make_pair( 831.76, 77081149.0+77867729.0 ) ) );
	norm.push_back( make_pair( "WW", make_pair( 118.7, 6987123.0 ) ) );
	norm.push_back( make_pair( "WZ", make_pair( 47.13, 2995828.0 ) ) );
	norm.push_back( make_pair( "ZZ", make_pair( 16.523, 998034.0 ) ) );
	norm.push_back( make_pair( "tW", make_pair( 35.85, 6952830.0+6933093.0 ) ) );
	norm.push_back( make_pair( "WJets", make_pair( 61526.7, 15389794.0+28076023.0+14315916.0+28949965.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt15to20", make_pair( 720648000*0.00042, 2219296.0+1921955.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt20to30", make_pair( 1273190000*0.003, 15721703.0+15580377.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt30to50", make_pair( 139803000*0.01182, 14683358.0+15033813.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt50to80", make_pair( 19222500*0.02276, 9908530.0+9898384.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt80to120", make_pair( 2758420*0.03844, 6556825.0+4645907.0+6998498.0+5151336.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt120to170", make_pair( 469797*0.05362, 3993904.0+5991370.0+4048816.0+5946767.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt170to300", make_pair( 117989*0.07335, 3995209.0+4773157.0+9705028.0+3951949.0+4629913.0+9902747.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt300to470", make_pair( 7820.25*0.10196, 4041948.0+8199205.0+12237898.0+3895639.0+8253382.0+12367604.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt470to600", make_pair( 645.528*0.12242, 2003720.0+2909119.0+1847803.0+2754636.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt600to800", make_pair( 187.109*0.13412, 2032249.0+2858727.0+4874009.0+1977886.0+3112446.0+4882843.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt800to1000", make_pair( 32.3486*0.14552, 1927100.0+2752753.0+5048968.0+2035647.0+3085786.0+4917178.0 ) ) );
	norm.push_back( make_pair( "QCDMuEnriched_Pt1000toInf", make_pair( 10.4305*0.15544, 1896115.0+4641117.0+1965321.0+4968703.0 ) ) );

	Double_t lumi = 35867;
	for(int i=0;i<21;i++)
	{
		for(int j=0;j<norm.size();j++)
			if( h_[i]->GetName() == norm[j].first )
				h_[i]->Scale( lumi * norm[j].second.first / norm[j].second.second );
      else
        cout << "Something is wrong!"
	}

	///////////////////////////////////////////
	// -- Save histograms following types -- //
	///////////////////////////////////////////
	TH1D* h__[9];
  //Data
  h__[0] = (TH1D*)h_[0]->Clone(var+"_Data");
  //DY
  h_[1]->Add(h_[2]);
  h__[1] = (TH1D*)h_[1]->Clone(var+"_DY");
  //ttbar
  h__[2] = (TH1D*)h_[3]->Clone(var+"_ttbar");
  //WW
  h__[3] = (TH1D*)h_[4]->Clone(var+"_WW");
  //WZ
  h__[4] = (TH1D*)h_[5]->Clone(var+"_WZ");
  //ZZ
  h__[5] = (TH1D*)h_[6]->Clone(var+"_ZZ");
  //tW
  h__[6] = (TH1D*)h_[7]->Clone(var+"_tW");
  //WJets
  h__[7] = (TH1D*)h_[8]->Clone(var+"_WJets");
  //QCD
  for(int i=10;i<=20;i++)
    h_[9]->Add(h_[i]);
  h__[8] = (TH1D*)h_[9]->Clone(var+"_QCD");

	//Save with type names
	TFile* f_out = new TFile(var+"__"+rootfile,"recreate");
	for(int i=0; i<9; i++)
		h__[i]->Write();

}
