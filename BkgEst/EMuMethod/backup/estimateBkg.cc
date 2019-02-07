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

void fillSystematics( TH1D* data_driven, TH1D* stat, TH1D* systematic, TH1D* total );
void removeNegativeBins( TH1D* hist );

void estimateBkg( TString channel = "MuMu" )
{
    int c1 = 3;
    int c2 = 2;
    int c3 = 5;
    int c4 = 13;
    int c5 = 14;
    int c6 = 15;
    int c7 = 16;
    int c8 = 17;

    ////////////////////////////
    // -- Set MC histogram -- //
    ////////////////////////////
	TString inputname = "INPUT_ROOTFILE_MC.root"; //Choose your input root file for MC
    TFile f_input(inputname, "read");
    TH1D *emu_ttbar = (TH1D*)f_input.Get("h_emu_mass_ttbar");
    TH1D *emu_ttbarBackup = (TH1D*)f_input.Get("h_emu_mass_ttbarBackup");
    TH1D *emu_ttbar_M700to1000 = (TH1D*)f_input.Get("h_emu_mass_ttbar_M700to1000");
    TH1D *emu_ttbar_M1000toInf = (TH1D*)f_input.Get("h_emu_mass_ttbar_M1000toInf");
    TH1D *emu_DYtautau = (TH1D*)f_input.Get("h_emu_mass_DYTauTau");
    TH1D *emu_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_emu_mass_DYTauTau_M10to50_v1");
    TH1D *emu_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_emu_mass_DYTauTau_M10to50_v2");
    TH1D *emu_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_emu_mass_DYTauTau_M10to50_ext1v1");
    TH1D *emu_WW = (TH1D*)f_input.Get("h_emu_mass_WW");
    TH1D *emu_WZ = (TH1D*)f_input.Get("h_emu_mass_WZ");
    TH1D *emu_ZZ = (TH1D*)f_input.Get("h_emu_mass_ZZ");
    TH1D *emu_tW = (TH1D*)f_input.Get("h_emu_mass_tW");
    TH1D *emu_antitW = (TH1D*)f_input.Get("h_emu_mass_tbarW");

    TH1D *emuSS_ttbar = (TH1D*)f_input.Get("h_emuSS_mass_ttbar");
    TH1D *emuSS_ttbarBackup = (TH1D*)f_input.Get("h_emuSS_mass_ttbarBackup");
    TH1D *emuSS_ttbar_M700to1000 = (TH1D*)f_input.Get("h_emuSS_mass_ttbar_M700to1000");
    TH1D *emuSS_ttbar_M1000toInf = (TH1D*)f_input.Get("h_emuSS_mass_ttbar_M1000toInf");
    TH1D *emuSS_DYtautau = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau");
    TH1D *emuSS_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau_M10to50_v1");
    TH1D *emuSS_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau_M10to50_v2");
    TH1D *emuSS_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau_M10to50_ext1v1");
    TH1D *emuSS_WW = (TH1D*)f_input.Get("h_emuSS_mass_WW");
    TH1D *emuSS_WZ = (TH1D*)f_input.Get("h_emuSS_mass_WZ");
    TH1D *emuSS_ZZ = (TH1D*)f_input.Get("h_emuSS_mass_ZZ");
    TH1D *emuSS_tW = (TH1D*)f_input.Get("h_emuSS_mass_tW");
    TH1D *emuSS_antitW = (TH1D*)f_input.Get("h_emuSS_mass_tbarW");

	TString inputname3 = "INPUT_ROOTFILE_MC_MUON.root"; //Choose your muon channel input root file for MC
	if( channel == "EE" ) inputname3 = "INPUT_ROOTFILE_MC_ELECTRON.root"; //Choose your muon channel input root file for MC
    TFile f_input3(inputname3, "read");
    TH1D *dimu_ttbar = (TH1D*)f_input3.Get("h_mass_ttbar");
    TH1D *dimu_ttbarBackup = (TH1D*)f_input3.Get("h_mass_ttbarBackup");
    TH1D *dimu_ttbar_M700to1000 = (TH1D*)f_input3.Get("h_mass_ttbar_M700to1000");
    TH1D *dimu_ttbar_M1000toInf = (TH1D*)f_input3.Get("h_mass_ttbar_M1000toInf");
    TH1D *dimu_DYtautau = (TH1D*)f_input3.Get("h_mass_DYTauTau");
    TH1D *dimu_DYtautau_M10to50_v1 = (TH1D*)f_input3.Get("h_mass_DYTauTau_M10to50_v1");
    TH1D *dimu_DYtautau_M10to50_v2 = (TH1D*)f_input3.Get("h_mass_DYTauTau_M10to50_v2");
    TH1D *dimu_DYtautau_M10to50_ext1v1 = (TH1D*)f_input3.Get("h_mass_DYTauTau_M10to50_ext1v1");
    TH1D *dimu_WW = (TH1D*)f_input3.Get("h_mass_WW");
    TH1D *dimu_WZ = (TH1D*)f_input3.Get("h_mass_WZ");
    TH1D *dimu_ZZ = (TH1D*)f_input3.Get("h_mass_ZZ");
    TH1D *dimu_tW = (TH1D*)f_input3.Get("h_mass_tW");
    TH1D *dimu_antitW = (TH1D*)f_input3.Get("h_mass_tbarW");

    // -- Merge ttbar samples -- //
    emu_ttbar->Add(emu_ttbarBackup);
    emu_ttbar->Add(emu_ttbar_M700to1000);
    emu_ttbar->Add(emu_ttbar_M1000toInf);
    emuSS_ttbar->Add(emuSS_ttbarBackup);
    emuSS_ttbar->Add(emuSS_ttbar_M700to1000);
    emuSS_ttbar->Add(emuSS_ttbar_M1000toInf);
    dimu_ttbar->Add(dimu_ttbarBackup);
    dimu_ttbar->Add(dimu_ttbar_M700to1000);
    dimu_ttbar->Add(dimu_ttbar_M1000toInf);

    // DY M50 + M10to50
    emu_DYtautau->Add(emu_DYtautau_M10to50_v1);
    emu_DYtautau->Add(emu_DYtautau_M10to50_v2);
    emu_DYtautau->Add(emu_DYtautau_M10to50_ext1v1);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v1);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v2);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_ext1v1);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50_v1);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50_v2);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50_ext1v1);

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
    dimu_tW->Add(dimu_antitW);

    // -- Consider average weight of top pt re-weight -- //
    //emu_ttbar->Scale(1/AvgWeight);
    //emuSS_ttbar->Scale(1/AvgWeight);
    //dimu_ttbar->Scale(1/AvgWeight);

    // -- Set Color -- //
    emu_ttbar->SetFillColor(c1);
    emu_DYtautau->SetFillColor(c2);
    emu_diboson->SetFillColor(c4);
    emu_tW->SetFillColor(c7);

    emuSS_ttbar->SetFillColor(c1);
    emuSS_DYtautau->SetFillColor(c2);
    emuSS_diboson->SetFillColor(c4);
    emuSS_tW->SetFillColor(c7);

    dimu_ttbar->SetFillColor(c1);
    dimu_DYtautau->SetFillColor(c2);
    dimu_WW->SetFillColor(c4);
    dimu_tW->SetFillColor(c7);

    // -- No Stats -- //
    emu_ttbar->SetStats(kFALSE);
    emu_DYtautau->SetStats(kFALSE);
    emu_diboson->SetStats(kFALSE);
    emu_tW->SetStats(kFALSE);

    emuSS_ttbar->SetStats(kFALSE);
    emuSS_DYtautau->SetStats(kFALSE);
    emuSS_diboson->SetStats(kFALSE);
    emuSS_tW->SetStats(kFALSE);

    dimu_ttbar->SetStats(kFALSE);
    dimu_DYtautau->SetStats(kFALSE);
    dimu_WW->SetStats(kFALSE);
    dimu_tW->SetStats(kFALSE);

    // -- Set Data histogram -- //
	TString inputname2 = "INPUT_ROOTFILE_DATA.root"; //Choose your input root file for data
	TFile f_input2(inputname2, "read");
    TH1D *emu_data = (TH1D*)f_input2.Get("h_emu_mass_Data");
    TH1D *emuSS_data = (TH1D*)f_input2.Get("h_emuSS_mass_Data");

    emu_data->SetMarkerStyle(20);
    emu_data->SetMarkerSize(0.8);
    emu_data->SetStats(kFALSE);

    emuSS_data->SetMarkerStyle(20);
    emuSS_data->SetMarkerSize(0.8);
    emuSS_data->SetStats(kFALSE);

    // -- Remove negative bins -- //
    removeNegativeBins( emu_DYtautau );
    removeNegativeBins( emuSS_DYtautau );
    removeNegativeBins( dimu_DYtautau );

    // -- Stack histograms -- //
    TH1D* emu_sumBkg = new TH1D("emu_sumBkg","",binnum,bins);
    emu_sumBkg->Add(emu_DYtautau);	
    emu_sumBkg->Add(emu_ttbar);	
    emu_sumBkg->Add(emu_diboson);	
    emu_sumBkg->Add(emu_tW);	

    THStack* emu_stackBkg = new THStack("emu_stackBkg","");
    emu_stackBkg->Add(emu_diboson);
    emu_stackBkg->Add(emu_DYtautau);
    emu_stackBkg->Add(emu_tW);
    emu_stackBkg->Add(emu_ttbar);

    TH1D* emuSS_sumMC = new TH1D("emuSS_sumMC","",binnum,bins);
    emuSS_sumMC->Add(emuSS_DYtautau);  
    emuSS_sumMC->Add(emuSS_ttbar); 
    emuSS_sumMC->Add(emuSS_diboson); 
    emuSS_sumMC->Add(emuSS_tW);  

    THStack* emuSS_stackMC = new THStack("emuSS_stackMC","");
    emuSS_stackMC->Add(emuSS_diboson);
    emuSS_stackMC->Add(emuSS_DYtautau);
    emuSS_stackMC->Add(emuSS_tW);
    emuSS_stackMC->Add(emuSS_ttbar);

    // -- Legend -- //
    TLegend* legend = new TLegend(.75,.75,.95,.89);
    legend->AddEntry(emu_data,"Bkg(data-driven)");
    legend->AddEntry(dimu_ttbar,"ttbar","F");
    legend->AddEntry(dimu_DYtautau,"DY#tau#tau","F");
    legend->AddEntry(dimu_tW,"tW+#bar{t}W","F");
    legend->AddEntry(dimu_WW,"WW","F");
    legend->SetBorderSize(0);  
    legend->SetFillStyle(0);  

    TLegend* legend_ttbar = new TLegend(.75,.75,.95,.89);
    legend_ttbar->AddEntry(emu_data,"Bkg(data-driven)");
    legend_ttbar->AddEntry(dimu_ttbar,"ttbar","F");
    legend_ttbar->SetBorderSize(0);  
    legend_ttbar->SetFillStyle(0);  

    TLegend* legend_DYtautau = new TLegend(.75,.75,.95,.89);
    legend_DYtautau->AddEntry(emu_data,"Bkg(data-driven)");
    legend_DYtautau->AddEntry(dimu_DYtautau,"DY#tau#tau","F");
    legend_DYtautau->SetBorderSize(0);  
    legend_DYtautau->SetFillStyle(0);  

    TLegend* legend_tW = new TLegend(.75,.75,.95,.89);
    legend_tW->AddEntry(emu_data,"Bkg(data-driven)");
    legend_tW->AddEntry(dimu_tW,"tW+#bar{t}W","F");
    legend_tW->SetBorderSize(0);  
    legend_tW->SetFillStyle(0);  

    TLegend* legend_WW = new TLegend(.75,.75,.95,.89);
    legend_WW->AddEntry(emu_data,"Bkg(data-driven)");
    legend_WW->AddEntry(dimu_WW,"WW","F");
    legend_WW->SetBorderSize(0);  
    legend_WW->SetFillStyle(0);  

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

    //emu_QCD->Smooth(Nsmooth,"R"); //Smooting QCD
    //printf("Smoothing QCD : %d times\n", Nsmooth);

    // -- Ratio = emu_Data/emu_MC -- //
    TH1D* emu_data_nonQCD = (TH1D*)emu_data->Clone();
    emu_data_nonQCD->Add(emu_QCD,-1.0);

    TH1D* emu_ratio = (TH1D*)emu_data_nonQCD->Clone("emu_ratio");
    emu_ratio->Divide(emu_data_nonQCD,emu_sumBkg,1.0,1.0,"B");

    // -- Set data-driven histograms -- //
    removeNegativeBins(emu_ratio);

    TH1D* data_driven_ttbar = (TH1D*)emu_ratio->Clone("data_driven_ttbar");
    TH1D* data_driven_tW = (TH1D*)emu_ratio->Clone("data_driven_tW");
    TH1D* data_driven_WW = (TH1D*)emu_ratio->Clone("data_driven_WW");
    TH1D* data_driven_DYtautau = (TH1D*)emu_ratio->Clone("data_driven_DYtautau");

    data_driven_ttbar->Multiply(dimu_ttbar);
    data_driven_DYtautau->Multiply(dimu_DYtautau);
    data_driven_WW->Multiply(dimu_WW);
    data_driven_tW->Multiply(dimu_tW);

    // -- Replace the highest mass bins with MC (becasue data-driven entries are zero) -- //
    data_driven_ttbar->SetBinContent(binnum,dimu_ttbar->GetBinContent(binnum));
    data_driven_ttbar->SetBinError(binnum,dimu_ttbar->GetBinError(binnum)); 
    data_driven_tW->SetBinContent(binnum,dimu_tW->GetBinContent(binnum));
    data_driven_tW->SetBinError(binnum,dimu_tW->GetBinError(binnum)); 
    data_driven_WW->SetBinContent(binnum,dimu_WW->GetBinContent(binnum));
    data_driven_WW->SetBinError(binnum,dimu_WW->GetBinError(binnum)); 

    // -- Calculate systematic uncertainty -- //
    TH1D* ttbar_total = new TH1D("ttbar_total","",binnum,bins);
    TH1D* ttbar_systematic = new TH1D("ttbar_systematic","",binnum,bins);
    TH1D* ttbar_stat = new TH1D("ttbar_stat","",binnum,bins);

    TH1D* tW_total = new TH1D("tW_total","",binnum,bins);
    TH1D* tW_systematic = new TH1D("tW_systematic","",binnum,bins);
    TH1D* tW_stat = new TH1D("tW_stat","",binnum,bins);

    TH1D* WW_total = new TH1D("WW_total","",binnum,bins);
    TH1D* WW_systematic = new TH1D("WW_systematic","",binnum,bins);
    TH1D* WW_stat = new TH1D("WW_stat","",binnum,bins);

    TH1D* DYtautau_total = new TH1D("DYtautau_total","",binnum,bins);
    TH1D* DYtautau_systematic = new TH1D("DYtautau_systematic","",binnum,bins);
    TH1D* DYtautau_stat = new TH1D("DYtautau_stat","",binnum,bins);

    ttbar_systematic->Add(dimu_ttbar);
    ttbar_systematic->Add(data_driven_ttbar,-1.0);

    DYtautau_systematic->Add(dimu_DYtautau);
    DYtautau_systematic->Add(data_driven_DYtautau,-1.0);

    tW_systematic->Add(dimu_tW);
    tW_systematic->Add(data_driven_tW,-1.0);

    WW_systematic->Add(dimu_WW);
    WW_systematic->Add(data_driven_WW,-1.0);

    fillSystematics( data_driven_ttbar, ttbar_stat, ttbar_systematic, ttbar_total );
    fillSystematics( data_driven_DYtautau, DYtautau_stat, DYtautau_systematic, DYtautau_total );
    fillSystematics( data_driven_tW, tW_stat, tW_systematic, tW_total );
    fillSystematics( data_driven_WW, WW_stat, WW_systematic, WW_total );

    // -- Set Bkg histogram name -- //
    data_driven_ttbar->SetName("ttbar");  
    data_driven_DYtautau->SetName("DYtautau");
    data_driven_tW->SetName("tW");
    data_driven_WW->SetName("WW");

    dimu_ttbar->SetName("ttbar_MC");  
    dimu_DYtautau->SetName("DYtautau_MC");
    dimu_tW->SetName("tW_MC");
    dimu_WW->SetName("WW_MC");

    // -- Output ROOT file -- //
    TFile* g = new TFile("./result/estimatedBkg/emu.root","RECREATE");

    data_driven_ttbar->Write();
    data_driven_DYtautau->Write();
    data_driven_tW->Write();
    data_driven_WW->Write();

    dimu_ttbar->Write();
    dimu_DYtautau->Write();
    dimu_tW->Write();
    dimu_WW->Write();

    ttbar_systematic->Write();
    DYtautau_systematic->Write();
    tW_systematic->Write();
    WW_systematic->Write();

    ttbar_stat->Write();
    DYtautau_stat->Write();
    tW_stat->Write();
    WW_stat->Write();

    g->Close();

    // -- Stack MC-Bkg histograms -- //
    THStack* h_ttbar = new THStack("h_ttbar","");
    h_ttbar->Add(dimu_ttbar);

    THStack* h_DYtautau = new THStack("h_DYtautau","");
    h_DYtautau->Add(dimu_DYtautau);

    THStack* h_tW = new THStack("h_tW","");
    h_tW->Add(dimu_tW);

    THStack* h_WW = new THStack("h_WW","");
    h_WW->Add(dimu_WW);

    THStack* h_all = new THStack("h_all","");
    h_all->Add(dimu_WW);
    h_all->Add(dimu_tW);
    h_all->Add(dimu_DYtautau);
    h_all->Add(dimu_ttbar);

    // -- Merge Data -- //
    TList* datalist = new TList;
    datalist->Add(data_driven_ttbar);
    datalist->Add(data_driven_DYtautau);
    datalist->Add(data_driven_tW);
    datalist->Add(data_driven_WW);
    
    TH1D* data_driven_all = new TH1D("data_driven_all", "",binnum, bins);
    data_driven_all->Merge(datalist);

    // Latex style
    TString obj;
    if( channel == "MuMu" ) obj = "#mu";
    else if( channel == "EE" ) obj = "e";

    mkplot( h_all, data_driven_all, legend, "All" );

	// Make plot
	MyCanvas *myc_ttbar = new MyCanvas(outputname+"_"+var+"_ttbar", "M("+obj+obj+") [GeV]", "Number of events");
	myc_ttbar->SetLogx();
	myc_ttbar->SetLogy(0);
	myc_ttbar->SetXRange(15, 3000);
	myc_ttbar->SetYRange(0.1, 1e8);
	//myc_ttbar->SetRatioRange(0.9, 1.1);
	myc_ttbar->SetRatioRange(0.7, 1.3);
	myc_ttbar->CanvasWithTHStackRatioPlot( data_driven_ttbar, h_ttbar, legend_ttbar, "Data-driven/MC", 0);
	myc_ttbar->PrintCanvas();

	MyCanvas *myc_DYtautau = new MyCanvas(outputname+"_"+var+"_DYtautau", "M("+obj+obj+") [GeV]", "Number of events");
	myc_DYtautau->SetLogx();
	myc_DYtautau->SetLogy(0);
	myc_DYtautau->SetXRange(15, 3000);
	myc_DYtautau->SetYRange(0.1, 1e8);
	//myc_DYtautau->SetRatioRange(0.9, 1.1);
	myc_DYtautau->SetRatioRange(0.7, 1.3);
	myc_DYtautau->CanvasWithTHStackRatioPlot( data_driven_DYtautau, h_DYtautau, legend_DYtautau, "Data-driven/MC", 0);
	myc_DYtautau->PrintCanvas();

	MyCanvas *myc_tW = new MyCanvas(outputname+"_"+var+"_tW", "M("+obj+obj+") [GeV]", "Number of events");
	myc_tW->SetLogx();
	myc_tW->SetLogy(0);
	myc_tW->SetXRange(15, 3000);
	myc_tW->SetYRange(0.1, 1e8);
	//myc_tW->SetRatioRange(0.9, 1.1);
	myc_tW->SetRatioRange(0.7, 1.3);
	myc_tW->CanvasWithTHStackRatioPlot( data_driven_tW, h_tW, legend_tW, "Data-driven/MC", 0);
	myc_tW->PrintCanvas();

	MyCanvas *myc_WW = new MyCanvas(outputname+"_"+var+"_WW", "M("+obj+obj+") [GeV]", "Number of events");
	myc_WW->SetLogx();
	myc_WW->SetLogy(0);
	myc_WW->SetXRange(15, 3000);
	myc_WW->SetYRange(0.1, 1e8);
	//myc_WW->SetRatioRange(0.9, 1.1);
	myc_WW->SetRatioRange(0.7, 1.3);
	myc_WW->CanvasWithTHStackRatioPlot( data_driven_WW, h_WW, legend_WW, "Data-driven/MC", 0);
	myc_WW->PrintCanvas();

	MyCanvas *myc_all = new MyCanvas(outputname+"_"+var+"_all", "M("+obj+obj+") [GeV]", "Number of events");
	myc_all->SetLogx();
	myc_all->SetLogy(0);
	myc_all->SetXRange(15, 3000);
	myc_all->SetYRange(0.1, 1e8);
	//myc_all->SetRatioRange(0.9, 1.1);
	myc_all->SetRatioRange(0.7, 1.3);
	myc_all->CanvasWithTHStackRatioPlot( data_driven_all, h_all, legend, "Data-driven/MC", 0);
	myc_all->PrintCanvas();
}


void fillSystematics( TH1D* data_driven, TH1D* stat, TH1D* systematic, TH1D* total ) {

    double binSystematic = 0;
    double binStat = 0;
    double binTotal = 0;

    for(int i=0; i<data_driven->GetNbinsX(); i++) {

        if(data_driven->GetBinContent(i+1)!=0) {  
            binStat = data_driven->GetBinError(i+1);
            binSystematic = fabs(systematic->GetBinContent(i+1));
            binTotal = sqrt(binSystematic*binSystematic + binStat*binStat);
        }
        else{
            binSystematic = 0;
            binStat = 0;
            binTotal = 0;
        }

        systematic->SetBinContent(i+1,binSystematic);
        stat->SetBinContent(i+1,binStat);    
        total->SetBinContent(i+1,binTotal);

        data_driven->SetBinError(i+1,binTotal);

    }

}

void removeNegativeBins( TH1D* hist ) {

    for(int i=0; i<hist->GetNbinsX(); i++) {
        if(hist->GetBinContent(i+1)<0) {
            hist->SetBinContent(i+1,0);
            hist->SetBinError(i+1,0);
        }
    }   

}

