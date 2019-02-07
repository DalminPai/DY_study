#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "../../interface/MyCanvas.C"
#include "../../interface/variables.C"
using namespace std;

void removeNegativeBins( TH1D* hist );

void emuCheck( Double_t type = 1 )
{
	int c1 = 3;
	int c2 = 2;
	int c3 = 5;
	int c4 = 13;
	int c5 = 14;
	int c6 = 15;
	int c7 = 16;
	int c8 = 17;

	// -- Choose variable -- //
	TString var, varSS, objects;
    if( 1 <= type && type < 2 )
		objects = "emu";
	else
		objects = "mu";
		//objects = "ele";

	var = objects + "_" + variables(type);
	varSS = objects + "SS_" + variables(type);

	////////////////////////////
	// -- Set MC histogram -- //
	////////////////////////////
	TString inputname = "INPUT_ROOTFILE_MC.root"; //Choose your input root file for MC
	TFile f_input(inputname, "read");
	TH1D *emu_ttbar = (TH1D*)f_input.Get("h_"+var+"_ttbar");
	TH1D *emu_ttbarBackup = (TH1D*)f_input.Get("h_"+var+"_ttbarBackup");
	TH1D *emu_ttbar_M700to1000 = (TH1D*)f_input.Get("h_"+var+"_ttbar_M700to1000");
	TH1D *emu_ttbar_M1000toInf = (TH1D*)f_input.Get("h_"+var+"_ttbar_M1000toInf");
	TH1D *emu_DYtautau = (TH1D*)f_input.Get("h_"+var+"_DYTauTau");
	TH1D *emu_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_"+var+"_DYTauTau_M10to50_v1");
	TH1D *emu_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_"+var+"_DYTauTau_M10to50_v2");
	TH1D *emu_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_"+var+"_DYTauTau_M10to50_ext1v1");
	TH1D *emu_WW = (TH1D*)f_input.Get("h_"+var+"_WW");
	TH1D *emu_WZ = (TH1D*)f_input.Get("h_"+var+"_WZ");
	TH1D *emu_ZZ = (TH1D*)f_input.Get("h_"+var+"_ZZ");
	TH1D *emu_tW = (TH1D*)f_input.Get("h_"+var+"_tW");
	TH1D *emu_antitW = (TH1D*)f_input.Get("h_"+var+"_tbarW");

	TH1D *emuSS_ttbar = (TH1D*)f_input.Get("h_"+varSS+"_ttbar");
	TH1D *emuSS_ttbarBackup = (TH1D*)f_input.Get("h_"+varSS+"_ttbarBackup");
	TH1D *emuSS_ttbar_M700to1000 = (TH1D*)f_input.Get("h_"+varSS+"_ttbar_M700to1000");
	TH1D *emuSS_ttbar_M1000toInf = (TH1D*)f_input.Get("h_"+varSS+"_ttbar_M1000toInf");
	TH1D *emuSS_DYtautau = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau");
	TH1D *emuSS_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau_M10to50_v1");
	TH1D *emuSS_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau_M10to50_v2");
	TH1D *emuSS_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau_M10to50_ext1v1");
	TH1D *emuSS_WW = (TH1D*)f_input.Get("h_"+varSS+"_WW");
	TH1D *emuSS_WZ = (TH1D*)f_input.Get("h_"+varSS+"_WZ");
	TH1D *emuSS_ZZ = (TH1D*)f_input.Get("h_"+varSS+"_ZZ");
	TH1D *emuSS_tW = (TH1D*)f_input.Get("h_"+varSS+"_tW");
	TH1D *emuSS_antitW = (TH1D*)f_input.Get("h_"+varSS+"_tbarW");

	// -- Merge ttbar samples -- //
	emu_ttbar->Add(emu_ttbarBackup);
	emu_ttbar->Add(emu_ttbar_M700to1000);
	emu_ttbar->Add(emu_ttbar_M1000toInf);
	emuSS_ttbar->Add(emuSS_ttbarBackup);
	emuSS_ttbar->Add(emuSS_ttbar_M700to1000);
	emuSS_ttbar->Add(emuSS_ttbar_M1000toInf);

	// DY M50 + M10to50
	emu_DYtautau->Add(emu_DYtautau_M10to50_v1);
	emu_DYtautau->Add(emu_DYtautau_M10to50_v2);
	emu_DYtautau->Add(emu_DYtautau_M10to50_ext1v1);
	emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v1);
	emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v2);
	emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_ext1v1);

	// WW + WZ + ZZ
	TH1D* emu_diboson = emu_WW;
	emu_diboson->Add(emu_WZ);
	emu_diboson->Add(emu_ZZ);
	TH1D* emuSS_diboson = emuSS_WW;
	emuSS_diboson->Add(emuSS_WZ);
	emuSS_diboson->Add(emuSS_ZZ);

	// tW + antitW
	emu_tW->Add(emu_antitW);
	emuSS_tW->Add(emuSS_antitW);

	// -- Consider average weight of top pt re-weight -- //
	//emu_ttbar->Scale(1/AvgWeight);
	//emuSS_ttbar->Scale(1/AvgWeight);

	// -- Set Color -- //
	emu_ttbar->SetFillColor(c1);
	emu_DYtautau->SetFillColor(c2);
	emu_diboson->SetFillColor(c4);
	emu_tW->SetFillColor(c7);

	emuSS_ttbar->SetFillColor(c1);
	emuSS_DYtautau->SetFillColor(c2);
	emuSS_diboson->SetFillColor(c4);
	emuSS_tW->SetFillColor(c7);

	// -- No Stats -- //
	emu_ttbar->SetStats(kFALSE);
	emu_DYtautau->SetStats(kFALSE);
	emu_diboson->SetStats(kFALSE);
	emu_tW->SetStats(kFALSE);

	emuSS_ttbar->SetStats(kFALSE);
	emuSS_DYtautau->SetStats(kFALSE);
	emuSS_diboson->SetStats(kFALSE);
	emuSS_tW->SetStats(kFALSE);

	//////////////////////////////
	// -- Set Data histogram -- //
	//////////////////////////////
	TString inputname2 = "INPUT_ROOTFILE_DATA.root"; //Choose your input root file for data
	TFile f_input2(inputname2, "read");
	TH1D *emu_data = (TH1D*)f_input2.Get("h_"+var+"_Data");
	TH1D *emuSS_data = (TH1D*)f_input2.Get("h_"+varSS+"_Data");

	emu_data->SetMarkerStyle(20);
	emu_data->SetMarkerSize(0.8);
	emu_data->SetStats(kFALSE);

	emuSS_data->SetMarkerStyle(20);
	emuSS_data->SetMarkerSize(0.8);
	emuSS_data->SetStats(kFALSE);

	// -- Remove negative bins -- //
	removeNegativeBins( emu_DYtautau );
	removeNegativeBins( emuSS_DYtautau );

	// -- Rebin -- //
	if( 2 <= type && type < 7 )
	{
		Int_t nRebin = 10;
		//Int_t nRebin = 20;
		if( 2 <= type && type < 3 ) nRebin = 5;

		emu_diboson->Rebin(nRebin);
		emu_tW->Rebin(nRebin);
		emu_DYtautau->Rebin(nRebin);
		emu_ttbar->Rebin(nRebin);
		emu_data->Rebin(nRebin);

		emuSS_diboson->Rebin(nRebin);
		emuSS_tW->Rebin(nRebin);
		emuSS_DYtautau->Rebin(nRebin);
		emuSS_ttbar->Rebin(nRebin);
		emuSS_data->Rebin(nRebin);
	}

	// -- Stack histograms -- //
	THStack* emu_stackBkg = new THStack("emu_stackBkg","");
	emu_stackBkg->Add(emu_diboson);
	emu_stackBkg->Add(emu_tW);
	emu_stackBkg->Add(emu_DYtautau);
	emu_stackBkg->Add(emu_ttbar);

	THStack* emuSS_stackMC = new THStack("emuSS_stackMC","");
	emuSS_stackMC->Add(emuSS_diboson);
	emuSS_stackMC->Add(emuSS_tW);
	emuSS_stackMC->Add(emuSS_DYtautau);
	emuSS_stackMC->Add(emuSS_ttbar);

	// -- Legend -- //
	TLegend* legend = new TLegend(.75,.75,.95,.89);
	legend->AddEntry(emu_data,"Data");
	legend->AddEntry(emu_ttbar,"ttbar","F");
	legend->AddEntry(emu_DYtautau,"DY#tau#tau","F");
	legend->AddEntry(emu_tW,"tW+#bar{t}W","F");
	legend->AddEntry(emu_diboson,"VV","F");
	legend->SetBorderSize(0);  
	legend->SetFillStyle(0);  

	// -- Estimate emu QCD -- //
	TH1D* emu_QCD = (TH1D*)emuSS_data->Clone();
	emu_QCD->Add(emuSS_DYtautau,-1.0);
	emu_QCD->Add(emuSS_ttbar,-1.0);
	emu_QCD->Add(emuSS_diboson,-1.0);
	emu_QCD->Add(emuSS_tW,-1.0);
	emu_QCD->SetFillColor(7);

	const double RR = 0.57147108645;
	emu_QCD->Scale(1/RR);

	removeNegativeBins(emu_QCD);
	//emu_QCD->Smooth(Nsmooth,"R"); //Smoothing QCD
	//printf("Smoothing QCD : %d times\n", Nsmooth);

	// -- with QCD -- //	
	THStack* emu_stackBkg_QCD = new THStack("emu_stackBkg_QCD","");
	emu_stackBkg_QCD->Add(emu_QCD);
	emu_stackBkg_QCD->Add(emu_diboson);
	emu_stackBkg_QCD->Add(emu_tW);
	emu_stackBkg_QCD->Add(emu_DYtautau);
	emu_stackBkg_QCD->Add(emu_ttbar);

	TLegend* legend_QCD = (TLegend*)legend->Clone();
	legend_QCD->AddEntry(emu_QCD,"QCD","F");

    //////////////////
    // -- Result -- //
    //////////////////
	TString outputname = "./result/emu_check";

    cout << "==========================================================" << endl;
    cout << "MC input file   : " << inputname << endl;
    cout << "Data input file : " << inputname2 << endl;
    cout << "Output file     : " << outputname << endl;
    cout << endl << "Running for [" << var << "]..." << endl;
    cout << "==========================================================" << endl;

	// -- Check emu -- //
	cout<<"[Check emu events...!]"<<endl;
	cout<<endl;

	cout<<"data in emu: "<<emu_data->Integral()<<endl;
	cout<<"ttbar in emu: "<<emu_ttbar->Integral()<<endl;
	cout<<"DYtautau in emu: "<<emu_DYtautau->Integral()<<endl;
	cout<<"tW in emu: "<<emu_tW->Integral()<<endl;
	cout<<"diboson in emu: "<<emu_diboson->Integral()<<endl;
	cout<<"QCD in emu: "<<emu_QCD->Integral()<<endl;
	cout<<endl;

	double data = emu_data->Integral();
	cout<<"ttbar contribution in data: "<<(emu_ttbar->Integral())/data<<endl;
	cout<<"DYtautau contribution in data: "<<(emu_DYtautau->Integral())/data<<endl;
	cout<<"tW contribution in data: "<<(emu_tW->Integral())/data<<endl;
	cout<<"diboson contribution in data: "<<(emu_diboson->Integral())/data<<endl;
	cout<<"QCD contribution in data: "<<(emu_QCD->Integral())/data<<endl;
	cout<<endl;

	// -- Check emuSS -- //
	cout<<"[Check emuSS events...!]"<<endl;
	cout<<endl;

	cout<<"data in emuSS: "<<emuSS_data->Integral()<<endl;
	cout<<"ttbar in emuSS: "<<emuSS_ttbar->Integral()<<endl;
	cout<<"DYtautau in emuSS: "<<emuSS_DYtautau->Integral()<<endl;
	cout<<"tW in emuSS: "<<emuSS_tW->Integral()<<endl;
	cout<<"diboson in emuSS: "<<emuSS_diboson->Integral()<<endl;
	cout<<endl;

	double dataSS = emuSS_data->Integral();
	cout<<"ttbar contribution in data: "<<(emuSS_ttbar->Integral())/dataSS<<endl;
	cout<<"DYtautau contribution in data: "<<(emuSS_DYtautau->Integral())/dataSS<<endl;
	cout<<"tW contribution in data: "<<(emuSS_tW->Integral())/dataSS<<endl;
	cout<<"diboson contribution in data: "<<(emuSS_diboson->Integral())/dataSS<<endl;
	cout<<endl;


	// -- Make plot -- //
	TString obj = "#mu";
	TString obj2 = "e";
	if( objects == "mu" ) objects = obj;

	TString x_title = "";
	if( 1 <= type && type < 2 ) x_title = "M("+obj+obj2+") [GeV]";
	else if( 2 <= type && type < 3 ) x_title = "P_{T}("+objects+") [GeV]";
	else if( 3 <= type && type < 4 ) x_title = "#eta("+objects+")";
	else if( 4 <= type && type < 5 ) x_title = "#phi("+objects+")";

	MyCanvas *myc = new MyCanvas(outputname+"_"+var, x_title, "Number of events"); //emu without QCD
	//MyCanvas *myc = new MyCanvas(outputname+"_"+var+"_with_QCD", x_title, "Number of events"); //emu with QCD
	//MyCanvas *myc = new MyCanvas(outputname+"_"+varSS, x_title, "Number of events"); //emuSS
	if( 1 <= type && type < 2 )
	{
		myc->SetLogx();
		myc->SetLogy(0);
		myc->SetXRange(15, 3000);
		myc->SetYRange(0.1, 1e8);
		//myc->SetRatioRange(0.9, 1.1);
		myc->SetRatioRange(0.7, 1.3);
	}
	else if( 2 <= type && type < 3 )
	{ 
		myc->SetLogy(0);
		myc->SetXRange(0, 500);
		myc->SetYRange(0.1, 1e8);
		myc->SetRatioRange(0.7, 1.3);
	}
	else if( 3 <= type && type < 4 )
	{
		myc->SetLogy(0);
		myc->SetXRange(-2.5, 2.5);
		myc->SetYRange(0.1, 1e8);
		//myc->SetRatioRange(0.85, 1.15);
		//myc->SetRatioRange(0.9, 1.1);
		myc->SetRatioRange(0.7, 1.3);
	}
	else if( 4 <= type && type < 5 )
	{
		myc->SetLogy(0);
		myc->SetXRange(-3.5, 3.5);
		myc->SetYRange(0.1, 1e8);
		myc->SetRatioRange(0.85, 1.15);
	}
	myc->CanvasWithTHStackRatioPlot( emu_data, emu_stackBkg, legend, "Data/MC", 0); //emu without QCD
	//myc->CanvasWithTHStackRatioPlot( emu_data, emu_stackBkg_QCD, legend_QCD, "Data/MC", 0); //emu with QCD
	//myc->CanvasWithTHStackRatioPlot( emuSS_data, emuSS_stackBkg, legend, "Data/MC", 0); //emuSS
	myc->PrintCanvas();
}

void removeNegativeBins( TH1D* hist ) {

	for(int i=0; i<hist->GetNbinsX(); i++) {
		if(hist->GetBinContent(i+1)<0) {
			hist->SetBinContent(i+1,0);
			hist->SetBinError(i+1,0);
		}
	}   

}


