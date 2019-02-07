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
#include "../../interface/tdrstyle.C"
#include "../../interface/CMS_lumi.C"
using namespace std;

void setDataHist(TH1D* hist);
void setMCHist(TH1D* hist, const int& color);
//void setFRHist(TH1D* hist);
TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator);
TH1D* FRBytRatio(TH1D** numerator, TH1D** denominator);

// Binning for 2015 DY analysis
/*const int ptbinnum_endcap = 9;
double ptbin_endcap[ptbinnum_endcap+1] = {47,52,60,70,80,90,100,150,200,500};
const int ptbinnum = 17;
double ptbin[ptbinnum+1] = {47,52,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500};*/

// Binning for 2016 DY analysis
const int ptbinnum_endcap = 8;
double ptbin_endcap[ptbinnum_endcap+1] = {52,60,70,80,90,100,150,200,500};
const int ptbinnum = 16;
double ptbin[ptbinnum+1] = {52,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500};

const int etabinnum = 2;
double etabin[etabinnum+1] = {0,1.2,2.4};

void my_estimateFR()
{
	//Set types: start from 0
	vector< TString > types;
	types.push_back("Data");
	types.push_back("QCD");
	types.push_back("WJets");
	types.push_back("DY");
	types.push_back("ttbar");
	types.push_back("tW");
	types.push_back("WW");
	types.push_back("WZ");
	types.push_back("ZZ");

	//TFile* f = new TFile("../bin/histograms/Re_hist_test.root","READ");
	TFile* f = new TFile("../bin/histograms/Re_hist_TightID_PFIso_20181014.root","READ");

	TH1D* denominator_pt_fit_barrel[9];
	TH1D* denominator_pt_fit_endcap[9];
	TH1D* denominator_pt_xsec_barrel[9];
	TH1D* denominator_pt_xsec_endcap[9];

	TH1D* numerator_pt_barrel[9];
	TH1D* numerator_pt_endcap[9];

	TH1D* denominator_barrel[9];
	TH1D* denominator_endcap[9];

	Int_t icolor[9] = {0, 7, 4, 2, 3, 5, 20, 30, 40};

	for(int i=0;i<9;i++)
	{
		denominator_pt_fit_barrel[i] = (TH1D*)f->Get("denominator_pt_barrel_"+types[i])->Clone("denominator_pt_fit_barrel_"+types[i]);
		denominator_pt_xsec_barrel[i] = (TH1D*)f->Get("denominator_pt_barrel_"+types[i])->Clone("denominator_pt_xsec_barrel_"+types[i]);
		denominator_barrel[i] = (TH1D*)f->Get("denominator_barrel_"+types[i])->Clone("denominator_barrel_"+types[i]);
		numerator_pt_barrel[i] = (TH1D*)f->Get("numerator_pt_barrel_"+types[i])->Clone("numerator_pt_barrel_"+types[i]);
		denominator_pt_fit_endcap[i] = (TH1D*)f->Get("denominator_pt_endcap_"+types[i])->Clone("denominator_pt_fit_endcap_"+types[i]);
		denominator_pt_xsec_endcap[i] = (TH1D*)f->Get("denominator_pt_endcap_"+types[i])->Clone("denominator_pt_xsec_endcap_"+types[i]);
		denominator_endcap[i] = (TH1D*)f->Get("denominator_endcap_"+types[i])->Clone("denominator_endcap_"+types[i]);
		numerator_pt_endcap[i] = (TH1D*)f->Get("numerator_pt_endcap_"+types[i])->Clone("numerator_pt_endcap_"+types[i]);

		if( types[i] == "Data" ) //Data
		{
			setDataHist( denominator_pt_fit_barrel[i] );
			setDataHist( denominator_pt_xsec_barrel[i] );
			setDataHist( denominator_barrel[i] );
			setDataHist( numerator_pt_barrel[i] );
			setDataHist( denominator_pt_fit_endcap[i] );
			setDataHist( denominator_pt_xsec_endcap[i] );
			setDataHist( denominator_endcap[i] );
			setDataHist( numerator_pt_endcap[i] );
		}
		else
		{
			setMCHist( denominator_pt_fit_barrel[i], icolor[i] );
			setMCHist( denominator_pt_xsec_barrel[i], icolor[i] );
			setMCHist( denominator_barrel[i], icolor[i] );
			setMCHist( numerator_pt_barrel[i], icolor[i] );
			setMCHist( denominator_pt_fit_endcap[i], icolor[i] );
			setMCHist( denominator_pt_xsec_endcap[i], icolor[i] );
			setMCHist( denominator_endcap[i], icolor[i] );
			setMCHist( numerator_pt_endcap[i], icolor[i] );
		}

		TH1D* hist[8];
		hist[0] = denominator_pt_fit_barrel[i];
		hist[1] = denominator_pt_xsec_barrel[i];
		hist[2] = denominator_barrel[i];
		hist[3] = numerator_pt_barrel[i];
		hist[4] = denominator_pt_fit_endcap[i];
		hist[5] = denominator_pt_xsec_endcap[i];
		hist[6] = denominator_endcap[i];
		hist[7] = numerator_pt_endcap[i];

	// Rebin
	TString hName = hist->GetName();
	if( hName.Contains("pt") )
	{
		if( hName.Contains("endcap") )
			hist = (TH1D*)(hist->Rebin(ptbinnum_endcap, "", ptbin_endcap));
		else
			hist = (TH1D*)(hist->Rebin(ptbinnum, "", ptbin));
	}
	else hist->Rebin(2);






	}


	// Check of Rebin
	cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
	cout << numerator_pt_endcap[0]->GetNbinsX() << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++" << endl;

	double norm_fit_barrel[9];
	double norm_fit_endcap[9];
	//double fit_value_barrel[9] = {0, 2.1672e+07, 1.4593e+07, 4.6877e+06, 1.4304e+06, 6.7331e+04, 8.4353e+04, 3.4200e+04, 1.1022e+04}; //before Oct 2018
	//double fit_value_endcap[9] = {0, 1.0308e+07, 1.1846e+07, 3.3025e+06, 5.8866e+05, 2.6427e+04, 6.1413e+04, 2.3134e+04, 6.9664e+03}; //before Oct 2018
	double fit_value_barrel[9] = {0, 2.0850e+07, 1.4720e+07, 4.6443e+06, 1.4144e+06, 6.0201e+04, 8.3562e+04, 3.0658e+04, 9.8743e+03};
	double fit_value_endcap[9] = {0, 1.0130e+07, 1.1588e+07, 3.5864e+06, 5.8503e+05, 2.6256e+04, 5.9376e+04, 2.0971e+04, 6.9167e+03};

	for(int i=0; i<types.size(); i++)
	{
		if( types[i] == "Data" ) continue;

		norm_fit_barrel[i] = fit_value_barrel[i]/denominator_barrel[i]->Integral();
		norm_fit_endcap[i] = fit_value_endcap[i]/denominator_endcap[i]->Integral();

		denominator_barrel[i]->Scale(norm_fit_barrel[i]);
		denominator_pt_fit_barrel[i]->Scale(norm_fit_barrel[i]);

		denominator_endcap[i]->Scale(norm_fit_endcap[i]);
		denominator_pt_fit_endcap[i]->Scale(norm_fit_endcap[i]);
	}

	TH1D* FR_template_barrel = (TH1D*)FRByTemplate(numerator_pt_barrel, denominator_pt_fit_barrel);
	TH1D* FR_template_endcap = (TH1D*)FRByTemplate(numerator_pt_endcap, denominator_pt_fit_endcap);

	TH1D* FR_xsec_barrel = (TH1D*)FRBytRatio(numerator_pt_barrel, denominator_pt_xsec_barrel);
	TH1D* FR_xsec_endcap = (TH1D*)FRBytRatio(numerator_pt_endcap, denominator_pt_xsec_endcap);

	int W = 1200;
	int H = 1200;

	int H_ref = 1200;
	int W_ref = 1200;

	//references for T, B, L, R
	float T = 0.08*H_ref;
	float B = 0.12*H_ref;
	float L = 0.12*W_ref;
	float R = 0.04*W_ref;

	lumi_13TeV = "35.9 fb^{-1}";
	writeExtraText = true;
	extraText = "Preliminary";
	drawLogo = false;
	tdrGrid(true);
	lumiTextSize = 0.5;
	cmsTextSize = 0.75;

	//TH1D* ptFrame = new TH1D("ptFrame","",8,47,500);
	TH1D* ptFrame = new TH1D("ptFrame","",8,52,500);
	ptFrame->SetStats(kFALSE);
	ptFrame->GetXaxis()->SetTitle("p_{T}[GeV]");
	ptFrame->GetYaxis()->SetTitle("Fake Rate");

	ptFrame->SetMinimum(0);
	ptFrame->SetMaximum(1.0); 
	ptFrame->GetXaxis()->SetTitleOffset(1);
	ptFrame->GetYaxis()->SetTitleOffset(1.1);
	ptFrame->GetXaxis()->SetTitleSize(0.05);
	ptFrame->GetYaxis()->SetTitleSize(0.05);  
	ptFrame->GetXaxis()->SetLabelSize(0.035);
	ptFrame->GetYaxis()->SetLabelSize(0.035); 
	ptFrame->GetXaxis()->SetMoreLogLabels(); 

	setTDRStyle();
	tdrGrid(true);
	lumiTextSize = 0.5;
	cmsTextSize = 0.75;

	TCanvas* canv = new TCanvas("canv","",1200,1200);
	canv->Draw();
	canv->cd();

	canv->SetFillColor(0);
	canv->SetLeftMargin( L/W );
	canv->SetRightMargin( R/W );
	canv->SetTopMargin( T/H );
	canv->SetBottomMargin( B/H );

	//Int_t MS = 3;
	Int_t MS = 2;
	FR_template_barrel->SetMarkerSize(MS);
	FR_template_endcap->SetMarkerSize(MS);
	FR_xsec_barrel->SetMarkerSize(MS);
	FR_xsec_endcap->SetMarkerSize(MS);

	FR_template_barrel->SetLineColor(1);
	FR_template_barrel->SetLineWidth(2);
	FR_template_barrel->SetMarkerStyle(20);
	FR_template_barrel->SetMarkerColor(1);

	FR_xsec_barrel->SetLineColor(2);
	FR_xsec_barrel->SetLineWidth(2);
	FR_xsec_barrel->SetMarkerStyle(21);
	FR_xsec_barrel->SetMarkerColor(2);

	FR_template_endcap->SetLineColor(1);
	FR_template_endcap->SetLineWidth(2);
	FR_template_endcap->SetMarkerStyle(20);
	FR_template_endcap->SetMarkerColor(1);

	FR_xsec_endcap->SetLineColor(2);
	FR_xsec_endcap->SetLineWidth(2);
	FR_xsec_endcap->SetMarkerStyle(21);
	FR_xsec_endcap->SetMarkerColor(2);

	canv->cd();
	canv->SetLogx();
	
	TLegend* legend2 = new TLegend(.45,.65,.75,.89);
	legend2->AddEntry(FR_template_barrel,"Template fitting", "EP");
	legend2->AddEntry(FR_xsec_barrel,"Ratio method", "EP");
	legend2->SetBorderSize(0);
	legend2->SetFillStyle(0);

	ptFrame->Draw();
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->SetTicks();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	FR_template_barrel->Draw("EPSAME");
	FR_xsec_barrel->Draw("EPSAME");
	legend2->Draw("SAME");
	canv->Print("../bin/print/FR_Barrel.pdf");

	canv->Clear();
	ptFrame->Draw();
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->SetTicks();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	FR_template_endcap->Draw("EPSAME");
	FR_xsec_endcap->Draw("EPSAME");
	legend2->Draw("SAME");
	canv->Print("../bin/print/FR_Endcap.pdf");

	TFile* g = new TFile("../bin/result/fakerate.root","RECREATE");
	FR_template_barrel->Write();
	FR_template_endcap->Write();
	FR_xsec_barrel->Write();
	FR_xsec_endcap->Write();
	g->Close();

	cout<<"Template Barrel"<<endl;
	for(int i=1; i<ptbinnum+1; i++)
	{
		if(i!=1) cout<<",";
		cout<<FR_template_barrel->GetBinContent(i);
	}
	cout<<endl;
	cout<<"Template Endcap"<<endl;
	for(int i=1; i<ptbinnum_endcap+1; i++)
	{
		if(i!=1) cout<<",";
		cout<<FR_template_endcap->GetBinContent(i);
	}
	cout<<endl;
	cout<<"Ratio Barrel"<<endl;
	for(int i=1; i<ptbinnum+1; i++)
	{
		if(i!=1) cout<<",";
		cout<<FR_xsec_barrel->GetBinContent(i);
	}
	cout<<endl;
	cout<<"Ratio Endcap"<<endl;
	for(int i=1; i<ptbinnum_endcap+1; i++)
	{
		if(i!=1) cout<<",";
		cout<<FR_xsec_endcap->GetBinContent(i);
	}
	cout<<endl;
}

void setDataHist(TH1D* hist)
{
	// Rebin
	TString hName = hist->GetName();
	if( hName.Contains("pt") )
	{
		cout << "pass" << endl;
		if( hName.Contains("endcap") )
			hist = (TH1D*)(hist->Rebin(ptbinnum_endcap, "", ptbin_endcap));
		else
			hist = (TH1D*)(hist->Rebin(ptbinnum, "", ptbin));
	}
	else hist->Rebin(2);

	hist->SetLineWidth(2);
	hist->SetMarkerStyle(33);
	//hist->SetMarkerSize(3);
	hist->SetStats(kFALSE);
	hist->Sumw2();
}

void setMCHist(TH1D* hist, const int& color)
{
	// Rebin
	TString hName = hist->GetName();
	if( hName.Contains("pt") )
	{
		if( hName.Contains("endcap") )
			hist = (TH1D*)(hist->Rebin(ptbinnum_endcap, "", ptbin_endcap));
		else
			hist = (TH1D*)(hist->Rebin(ptbinnum, "", ptbin));
	}
	else hist->Rebin(2);

	hist->SetFillColor(color+2);
	hist->SetStats(kFALSE);
	hist->Sumw2();
}

TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator)
{
	Int_t i_data = 0, i_qcd = 1, i_wjet = 2, i_dy = 3, i_ttbar = 4, i_tw = 5, i_ww = 6, i_wz = 7, i_zz = 8;

	TString name = ( ((TString)(denominator[i_qcd]->GetName())).Contains("barrel") ) ? "FR_template_barrel" : "FR_template_endcap";

	TH1D* num = (TH1D*)numerator[i_qcd]->Clone(name);
	TH1D* den = (TH1D*)numerator[i_dy]->Clone(name+"_");

	num->Multiply(numerator[i_data]);

	den->Add(numerator[i_qcd]);
	den->Add(numerator[i_wjet]);
	den->Add(numerator[i_ttbar]);
	den->Add(numerator[i_tw]);
	den->Add(numerator[i_ww]);
	den->Add(numerator[i_wz]);
	den->Add(numerator[i_zz]);
	den->Multiply(denominator[i_qcd]);

	num->Divide(den);

	delete den;
	return num;
}

TH1D* FRBytRatio(TH1D** numerator, TH1D** denominator)
{
	Int_t i_data = 0, i_qcd = 1, i_wjet = 2, i_dy = 3, i_ttbar = 4, i_tw = 5, i_ww = 6, i_wz = 7, i_zz = 8;

	TString name = ( ((TString)(denominator[i_qcd]->GetName())).Contains("barrel") ) ? "FR_xsec_barrel" : "FR_xsec_endcap";

	TH1D* num = (TH1D*)denominator[i_dy]->Clone(name);
	TH1D* den = (TH1D*)numerator[i_dy]->Clone(name+"_");

	num->Add(denominator[i_qcd]);
	num->Add(denominator[i_wjet]);
	num->Add(denominator[i_ttbar]);
	num->Add(denominator[i_tw]);
	num->Add(denominator[i_ww]);
	num->Add(denominator[i_wz]);
	num->Add(denominator[i_zz]);
	num->Multiply(numerator[i_qcd]);
	num->Multiply(numerator[i_data]);

	den->Add(numerator[i_qcd]);
	den->Add(numerator[i_wjet]);
	den->Add(numerator[i_ttbar]);
	den->Add(numerator[i_tw]);
	den->Add(numerator[i_ww]);
	den->Add(numerator[i_wz]);
	den->Add(numerator[i_zz]);
	den->Multiply(denominator[i_qcd]);
	den->Multiply(denominator[i_data]);

	num->Divide(den);

	delete den;
	return num;
}
