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
# define MS 2

void my_estimateDijet()
{
	////////////////////////////
	// -- Setup for canvas -- //
	////////////////////////////
	int W = 1200;
	int H = 1200;

	int H_ref = 1200;
	int W_ref = 1200;

	// references for T, B, L, R
	float T = 0.08*H_ref;
	float B = 0.12*H_ref;
	float L = 0.12*W_ref;
	float R = 0.04*W_ref;

	lumi_13TeV = "35.9 fb^{-1}";
	lumiTextSize = 0.5;
	writeExtraText = true;
	extraText = "Preliminary";
	drawLogo = false;

	int binnum = 43;
	double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

	TH1D* massFrame = new TH1D("massFrame","",38,15,3000);
	massFrame->SetMinimum(0.001);
	massFrame->SetMaximum(1000000);
	massFrame->SetStats(kFALSE);
	massFrame->GetXaxis()->SetTitle("Mass[GeV]");
	//massFrame->GetXaxis()->CenterTitle(kTRUE);
	//massFrame->GetYaxis()->CenterTitle(kTRUE);
	massFrame->GetYaxis()->SetTitleOffset(1);
	massFrame->GetYaxis()->SetTitle("Number of events");
	massFrame->GetXaxis()->SetTitleSize(0);
	massFrame->GetYaxis()->SetTitleSize(0.05);
	massFrame->GetXaxis()->SetLabelSize(0);
	//massFrame->GetYaxis()->SetLabelSize(0.025); 
	massFrame->GetXaxis()->SetMoreLogLabels();

	//////////////////////////////////////////
	// -- Get histograms from input file -- //
	//////////////////////////////////////////
	TFile* f;
	//f = new TFile("../bin/histograms/Re_fake_test.root","READ");
	f = new TFile("../bin/histograms/Re_fake_test_20181010.root","READ");
	
	TH1D* dijet_template[14];
	TH1D* dijetSS_template[14];
	TH1D* dijet_ratio[14];
	TH1D* dijetSS_ratio[14];

	// opposite-sign dijet (template)
	dijet_template[5] = (TH1D*)f->Get("histDijet1_Data"); //from Data
	dijet_template[5]->SetFillColor(7);
	dijet_template[5]->SetLineColor(7);
	dijet_template[5]->SetStats(kFALSE);
	dijet_template[5]->Sumw2();

	dijet_template[0] = (TH1D*)f->Get("histDijet1_DY"); //from DY
	dijet_template[0]->SetFillColor(2);
	dijet_template[0]->SetLineColor(2);
	dijet_template[0]->SetStats(kFALSE);
	dijet_template[0]->Sumw2();

	dijet_template[1] = (TH1D*)f->Get("histDijet1_ttbar"); //from ttbar
	dijet_template[1]->SetFillColor(3);
	dijet_template[1]->SetLineColor(3);
	dijet_template[1]->SetStats(kFALSE);
	dijet_template[1]->Sumw2();

	// same-sign dijet (template)
	dijetSS_template[5] = (TH1D*)f->Get("histSameDijet1_Data"); //from Data
	dijetSS_template[5]->SetFillColor(7);
	dijetSS_template[5]->SetLineColor(7);
	dijetSS_template[5]->SetStats(kFALSE);
	dijetSS_template[5]->Sumw2();

	dijetSS_template[0] = (TH1D*)f->Get("histSameDijet1_DY"); //from DY
	dijetSS_template[1] = (TH1D*)f->Get("histSameDijet1_ttbar"); //from ttbar

	cout << "DY OS before = " << dijet_template[0]->Integral() << endl;
	cout << "DY SS before = " << dijetSS_template[0]->Integral() << endl;
	cout << "ttbar OS before = " << dijet_template[1]->Integral() << endl;
	cout << "ttbar SS before = " << dijetSS_template[1]->Integral() << endl;

	// opposite-sign dijet (ratio)
	dijet_ratio[5] = (TH1D*)f->Get("histDijet2_Data"); //from Data
	dijet_ratio[5]->SetFillColor(7);
	dijet_ratio[5]->SetLineColor(7);
	dijet_ratio[5]->SetStats(kFALSE);
	dijet_ratio[5]->Sumw2();

	dijet_ratio[0] = (TH1D*)f->Get("histDijet2_DY"); //from DY
	dijet_ratio[0]->SetFillColor(2);
	dijet_ratio[0]->SetLineColor(2);
	dijet_ratio[0]->SetStats(kFALSE);
	dijet_ratio[0]->Sumw2();

	dijet_ratio[1] = (TH1D*)f->Get("histDijet2_ttbar"); //from ttbar
	dijet_ratio[1]->SetFillColor(3);
	dijet_ratio[1]->SetLineColor(3);
	dijet_ratio[1]->SetStats(kFALSE);
	dijet_ratio[1]->Sumw2();

	// same-sign dijet (ratio)
	dijetSS_ratio[5] = (TH1D*)f->Get("histSameDijet2_Data"); //from Data
	dijetSS_ratio[5]->SetFillColor(7);
	dijetSS_ratio[5]->SetLineColor(7);
	dijetSS_ratio[5]->SetStats(kFALSE);
	dijetSS_ratio[5]->Sumw2();

	cout << "DY(template): " << dijet_template[0]->Integral() << endl;
	cout << "DY(ratio): " << dijet_ratio[0]->Integral() << endl;

	// -- Some reference histograms -- //
	// from template fitting
	TH1D* dijet_template_Data = (TH1D*)f->Get("histDijet1_Data");
	dijet_template_Data->SetFillColor(11);
	dijet_template_Data->SetLineColor(11);
	dijet_template_Data->SetStats(kFALSE);
	dijet_template_Data->Sumw2();

	// from ratio
	TH1D* dijet_ratio_Data = (TH1D*)f->Get("histDijet2_Data");
	dijet_ratio_Data->SetFillColor(11);
	dijet_ratio_Data->SetLineColor(11);
	dijet_ratio_Data->SetStats(kFALSE);
	dijet_ratio_Data->Sumw2();

	TH1D* dijet_Data_2f = (TH1D*)f->Get("histDijet_Data");
	dijet_Data_2f->SetFillColor(11);
	dijet_Data_2f->SetLineColor(11);
	dijet_Data_2f->SetStats(kFALSE);
	dijet_Data_2f->Sumw2();

	TH1D* dijet_DY_2f = (TH1D*)f->Get("histDijet_DY");
	dijet_DY_2f->SetFillColor(2);
	dijet_DY_2f->SetLineColor(2);
	dijet_DY_2f->SetStats(kFALSE);
	dijet_DY_2f->Sumw2();

	TH1D* dijet_ttbar_2f = (TH1D*)f->Get("histDijet_ttbar");
	dijet_ttbar_2f->SetFillColor(3);
	dijet_ttbar_2f->SetLineColor(3);
	dijet_ttbar_2f->SetStats(kFALSE);
	dijet_ttbar_2f->Sumw2();

	//////////////////////////////
	// -- Set histogram axis -- //
	//////////////////////////////
	// opposite-sign dijet (template)
	dijet_template[5]->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_template[5]->GetYaxis()->SetTitleOffset(1.5);
	dijet_template[5]->GetYaxis()->SetTitle("Number of events");
	dijet_template[5]->GetXaxis()->SetLabelSize(0.025);
	dijet_template[5]->GetYaxis()->SetLabelSize(0.025);
	dijet_template[5]->GetXaxis()->SetMoreLogLabels(); 

	dijet_template[0]->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_template[0]->GetYaxis()->SetTitleOffset(1.5);
	dijet_template[0]->GetYaxis()->SetTitle("Number of events");
	dijet_template[0]->GetXaxis()->SetLabelSize(0.025);
	dijet_template[0]->GetYaxis()->SetLabelSize(0.025);
	dijet_template[0]->GetXaxis()->SetMoreLogLabels();

	dijet_template[1]->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_template[1]->GetYaxis()->SetTitleOffset(1.5);
	dijet_template[1]->GetYaxis()->SetTitle("Number of events");
	dijet_template[1]->GetXaxis()->SetLabelSize(0.025);
	dijet_template[1]->GetYaxis()->SetLabelSize(0.025);
	dijet_template[1]->GetXaxis()->SetMoreLogLabels(); 

	// same-sign dijet (template)
	dijetSS_template[5]->GetXaxis()->SetTitle("Mass[GeV]");
	dijetSS_template[5]->GetYaxis()->SetTitleOffset(1.5);
	dijetSS_template[5]->GetYaxis()->SetTitle("Number of events");
	dijetSS_template[5]->GetXaxis()->SetLabelSize(0.025);
	dijetSS_template[5]->GetYaxis()->SetLabelSize(0.025);
	dijetSS_template[5]->GetXaxis()->SetMoreLogLabels();

	// opposite-sign dijet (ratio)
	dijet_ratio[5]->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_ratio[5]->GetYaxis()->SetTitleOffset(1.5);
	dijet_ratio[5]->GetYaxis()->SetTitle("Number of events");
	dijet_ratio[5]->GetXaxis()->SetLabelSize(0.025);
	dijet_ratio[5]->GetYaxis()->SetLabelSize(0.025);
	dijet_ratio[5]->GetXaxis()->SetMoreLogLabels(); 

	dijet_ratio[0]->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_ratio[0]->GetYaxis()->SetTitleOffset(1.5);
	dijet_ratio[0]->GetYaxis()->SetTitle("Number of events");
	dijet_ratio[0]->GetXaxis()->SetLabelSize(0.025);
	dijet_ratio[0]->GetYaxis()->SetLabelSize(0.025);
	dijet_ratio[0]->GetXaxis()->SetMoreLogLabels(); 

	dijet_ratio[1]->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_ratio[1]->GetYaxis()->SetTitleOffset(1.5);
	dijet_ratio[1]->GetYaxis()->SetTitle("Number of events");
	dijet_ratio[1]->GetXaxis()->SetLabelSize(0.025);
	dijet_ratio[1]->GetYaxis()->SetLabelSize(0.025);
	dijet_ratio[1]->GetXaxis()->SetMoreLogLabels(); 

	// same-sign dijet (ratio)
	dijetSS_ratio[5]->GetXaxis()->SetTitle("Mass[GeV]");
	dijetSS_ratio[5]->GetYaxis()->SetTitleOffset(1.5);
	dijetSS_ratio[5]->GetYaxis()->SetTitle("Number of events");
	dijetSS_ratio[5]->GetXaxis()->SetLabelSize(0.025);
	dijetSS_ratio[5]->GetYaxis()->SetLabelSize(0.025);
	dijetSS_ratio[5]->GetXaxis()->SetMoreLogLabels();

	// -- Some reference histograms -- //
	dijet_template_Data->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_template_Data->GetYaxis()->SetTitleOffset(1.5);
	dijet_template_Data->GetYaxis()->SetTitle("Number of events");
	dijet_template_Data->GetXaxis()->SetLabelSize(0.025);
	dijet_template_Data->GetYaxis()->SetLabelSize(0.025);
	dijet_template_Data->GetXaxis()->SetMoreLogLabels();

	dijet_ratio_Data->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_ratio_Data->GetYaxis()->SetTitleOffset(1.5);
	dijet_ratio_Data->GetYaxis()->SetTitle("Number of events");
	dijet_ratio_Data->GetXaxis()->SetLabelSize(0.025);
	dijet_ratio_Data->GetYaxis()->SetLabelSize(0.025);
	dijet_ratio_Data->GetXaxis()->SetMoreLogLabels();

	dijet_Data_2f->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_Data_2f->GetYaxis()->SetTitleOffset(1.5);
	dijet_Data_2f->GetYaxis()->SetTitle("Number of events");
	dijet_Data_2f->GetXaxis()->SetLabelSize(0.025);
	dijet_Data_2f->GetYaxis()->SetLabelSize(0.025);
	dijet_Data_2f->GetXaxis()->SetMoreLogLabels();

	dijet_DY_2f->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_DY_2f->GetYaxis()->SetTitleOffset(1.5);
	dijet_DY_2f->GetYaxis()->SetTitle("Number of events");
	dijet_DY_2f->GetXaxis()->SetLabelSize(0.025);
	dijet_DY_2f->GetYaxis()->SetLabelSize(0.025);
	dijet_DY_2f->GetXaxis()->SetMoreLogLabels();

	dijet_ttbar_2f->GetXaxis()->SetTitle("Mass[GeV]");
	dijet_ttbar_2f->GetYaxis()->SetTitleOffset(1.5);
	dijet_ttbar_2f->GetYaxis()->SetTitle("Number of events");
	dijet_ttbar_2f->GetXaxis()->SetLabelSize(0.025);
	dijet_ttbar_2f->GetYaxis()->SetLabelSize(0.025);
	dijet_ttbar_2f->GetXaxis()->SetMoreLogLabels();

	//////////////////////////////////
	// -- Background subtraction -- //
	//////////////////////////////////
	// Remove negative bins of DY for QCD dijet
	for(int i=1; i<46; i++)
	{
		// -- DY -- //
		if( dijet_template[0]->GetBinContent(i) < 0 )
		{
			dijet_template[0]->SetBinContent(i,0.0);
			dijet_template[0]->SetBinError(i,0.0);
		}

		if( dijet_ratio[0]->GetBinContent(i) < 0 )
		{
			dijet_ratio[0]->SetBinContent(i,0.0);
			dijet_ratio[0]->SetBinError(i,0.0);
		}

		if( dijet_DY_2f->GetBinContent(i) < 0 )
		{
			dijet_DY_2f->SetBinContent(i,0.0);
			dijet_DY_2f->SetBinError(i,0.0);
		}
	}

	// Data - ( DY + ttbar ) in QCD dijet
	dijet_template[5]->Add(dijet_template[0],-1.0);
	dijet_template[5]->Add(dijet_template[1],-1.0);

	dijet_ratio[5]->Add(dijet_ratio[0],-1.0);
	dijet_ratio[5]->Add(dijet_ratio[1],-1.0);

	// Remove negative bins in QCD dijet
	for(int i=1; i<46; i++)
	{
		if( dijet_template[5]->GetBinContent(i) < 0 )
		{
			dijet_template[5]->SetBinContent(i,0.0);
			dijet_template[5]->SetBinError(i,0.0);
		}

		if( dijet_ratio[5]->GetBinContent(i) < 0 )
		{
			dijet_ratio[5]->SetBinContent(i,0.0);
			dijet_ratio[5]->SetBinError(i,0.0);
		}
	}

	// Smooth function
	//dijet_template[5]->Smooth();

	//////////////////////
	// -- Make plots -- //
	//////////////////////
	setTDRStyle();
	tdrGrid(true);
	lumiTextSize = 0.5;
	cmsTextSize = 0.75;

	TCanvas* canv = new TCanvas("canv","",1200,1200);
	canv->SetFillColor(0);
	canv->SetLeftMargin( L/W );
	canv->SetRightMargin( R/W );
	canv->SetTopMargin( T/H );
	canv->SetBottomMargin( B/H );

	// -- Opposite sign QCD -- //
	TLegend *leg_QCD = new TLegend(0.6,0.8,.95,.85);
	leg_QCD->SetBorderSize(0);
	leg_QCD->SetFillStyle(0);
	leg_QCD->AddEntry(dijet_template[5],"QCD (Opposite sign)","F");

	dijet_template[5]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_QCD->Draw("SAME");
	canv->Print("../bin/print/dijet_template.pdf");
	canv->Clear();

	dijet_ratio[5]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_QCD->Draw("SAME");
	canv->Print("../bin/print/dijet_ratio.pdf");
	canv->Clear();

	// -- Same sign QCD -- //
	TLegend *leg_QCD_SS = new TLegend(0.6,0.8,.95,.85);
	leg_QCD_SS->SetBorderSize(0);
	leg_QCD_SS->SetFillStyle(0);
	leg_QCD_SS->AddEntry(dijetSS_template[5],"QCD (Same sign)","F");

	dijetSS_template[5]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_QCD_SS->Draw("SAME");
	canv->Print("../bin/print/dijetSS_template.pdf");
	canv->Clear();

	dijetSS_ratio[5]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_QCD_SS->Draw("SAME");
	canv->Print("../bin/print/dijetSS_ratio.pdf");
	canv->Clear();

	// -- Both QCD -- //
	dijetSS_template[5]->SetLineColor(1);
	dijetSS_template[5]->SetMarkerColor(1);
	//dijetSS_template[5]->SetMarkerStyle(22);
	dijetSS_template[5]->SetMarkerStyle(20);
	//dijetSS_template[5]->SetMarkerSize(3);
	dijetSS_template[5]->SetMarkerSize(MS);

	dijetSS_ratio[5]->SetLineColor(1);
	dijetSS_ratio[5]->SetMarkerColor(1);
	//dijetSS_ratio[5]->SetMarkerStyle(22);
	dijetSS_ratio[5]->SetMarkerStyle(20);
	//dijetSS_ratio[5]->SetMarkerSize(3);
	dijetSS_ratio[5]->SetMarkerSize(MS);

	TLegend* legg = new TLegend(.6,.8,.95,.85);
	legg->AddEntry(dijet_template[5],"Opposite sign", "F");
	legg->AddEntry(dijetSS_template[5],"Same sign", "EP");
	legg->SetBorderSize(0);
	legg->SetFillStyle(0);

	dijet_template[5]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	dijetSS_template[5]->Draw("EPSAME");
	legg->Draw("SAME");
	canv->Print("../bin/print/dijetBoth_template.pdf");
	canv->Clear();

	dijet_ratio[5]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	dijetSS_ratio[5]->Draw("EPSAME");
	legg->Draw("SAME");
	canv->Print("../bin/print/dijetBoth_ratio.pdf");
	canv->Clear();

	// -- Opposite sign DY -- //
	TLegend *leg_DY = new TLegend(0.6,0.8,.95,.85);
	leg_DY->SetBorderSize(0);
	leg_DY->SetFillStyle(0);
	leg_DY->AddEntry(dijet_template[0],"DY (Opposite sign)","F");

	dijet_template[0]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_DY->Draw("SAME");
	canv->Print("../bin/print/dijet_template_DY.pdf");
	canv->Clear();

	dijet_ratio[0]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_DY->Draw("SAME");
	canv->Print("../bin/print/dijet_ratio_DY.pdf");
	canv->Clear();

	// -- Opposite sign ttbar -- //
	TLegend *leg_ttbar = new TLegend(0.6,0.8,.95,.85);
	leg_ttbar->SetBorderSize(0);
	leg_ttbar->SetFillStyle(0);
	leg_ttbar->AddEntry(dijet_template[1],"ttbar (Opposite sign)","F");

	dijet_template[1]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_ttbar->Draw("SAME");
	canv->Print("../bin/print/dijet_template_ttbar.pdf");
	canv->Clear();

	dijet_ratio[1]->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_ttbar->Draw("SAME");
	canv->Print("../bin/print/dijet_ratio_ttbar.pdf");
	canv->Clear();

	// -- Some reference histograms -- //
	TLegend *leg_Data = new TLegend(0.6,0.8,.95,.85);
	leg_Data->SetBorderSize(0);
	leg_Data->SetFillStyle(0);
	leg_Data->AddEntry(dijet_template_Data,"Data (Opposite sign)","F");

	dijet_template_Data->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_Data->Draw("SAME");
	canv->Print("../bin/print/dijet_template_Data.pdf");
	canv->Clear();

	dijet_ratio_Data->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_Data->Draw("SAME");
	canv->Print("../bin/print/dijet_ratio_Data.pdf");
	canv->Clear();

	dijet_Data_2f->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_Data->Draw("SAME");
	canv->Print("../bin/print/dijet_Data_2f.pdf");
	canv->Clear();

	dijet_DY_2f->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_DY->Draw("SAME");
	canv->Print("../bin/print/dijet_DY_2f.pdf");
	canv->Clear();

	dijet_ttbar_2f->Draw("HIST");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	leg_ttbar->Draw("SAME");
	canv->Print("../bin/print/dijet_ttbar_2f.pdf");
	canv->Clear();

	//////////////////////////////////////////////////////////
	// -- Calculate the number of events and uncertainty -- //
	//////////////////////////////////////////////////////////
	// Check the number
	double error = 0;
	dijet_template[5]->IntegralAndError(1,45,error);
	cout<<"QCD(template) = "<<dijet_template[5]->Integral(1,45)<<"+-"<<error<<endl;
	error = 0;
	dijetSS_template[5]->IntegralAndError(1,45,error);
	cout<<"QCD(template) SS = "<<dijetSS_template[5]->Integral(1,45)<<"+-"<<error<<endl;
	error = 0;
	dijet_ratio[5]->IntegralAndError(1,45,error);
	cout<<"QCD(ratio) = "<<dijet_ratio[5]->Integral(1,45)<<"+-"<<error<<endl;
	error = 0;
	dijetSS_ratio[5]->IntegralAndError(1,45,error);
	cout<<"QCD(ratio) SS = "<<dijetSS_ratio[5]->Integral(1,45)<<"+-"<<error<<endl;

	TH1D* dijet = (TH1D*)dijet_template[5]->Clone();
	dijet->Sumw2();
	dijet->SetName("dijet");

	// -- Uncertainty calculation -- //
	TH1D* dijet_total      = new TH1D("dijet_total","",binnum,bins);
	TH1D* dijet_systematic = new TH1D("dijet_systematic","",binnum,bins);
	TH1D* dijet_stat       = new TH1D("dijet_stat","",binnum,bins);

	for(int i=1; i<binnum+1; i++)
	{
		double systematic = fabs( dijet->GetBinContent(i) - dijet_ratio[5]->GetBinContent(i) );
		double stat = dijet->GetBinError(i);
		double total = sqrt( systematic*systematic + stat*stat );
		
		if(dijet->GetBinContent(i)==0)
		{
		  systematic = 0;
		  stat = 0;
		  total = 0;
		}

		dijet_systematic->SetBinContent(i,systematic);
		dijet_stat->SetBinContent(i,stat);
		dijet_total->SetBinContent(i,total);

		dijet->SetBinError(i,total);
	}

	// Save the result
	TFile* gg = new TFile("../bin/result/dijet.root","RECREATE");
	dijet->Write();
	dijet_systematic->Write();
	dijet_stat->Write();
	gg->Close();

	// Uncertainty plot
	dijet_systematic->Divide(dijet);
	dijet_stat->Divide(dijet);
	dijet_total->Divide(dijet);

	dijet_total->SetMarkerStyle(20);
	//dijet_total->SetMarkerSize(3);
	dijet_total->SetMarkerSize(MS);
	dijet_total->SetMarkerColor(1);

	dijet_systematic->SetMarkerStyle(22);
	//dijet_systematic->SetMarkerSize(3);
	dijet_systematic->SetMarkerSize(MS);
	dijet_systematic->SetMarkerColor(2);

	dijet_stat->SetMarkerStyle(21);
	//dijet_stat->SetMarkerSize(3);
	dijet_stat->SetMarkerSize(MS);
	dijet_stat->SetMarkerColor(4);

	TLegend* leg = new TLegend(.8,.15,.95,.3);
	leg->AddEntry(dijet_total,"Total","P");
	leg->AddEntry(dijet_systematic,"Sys.","P");
	leg->AddEntry(dijet_stat,"Stat.","P");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	massFrame->GetYaxis()->SetTitle("Rel. Uncertainty");
	massFrame->SetMinimum(0);
	massFrame->SetMaximum(1);
	massFrame->GetXaxis()->SetLabelSize(0.025);
	massFrame->GetXaxis()->SetTitleSize(0.05);

	massFrame->Draw();
	dijet_total->Draw("HISTPSAME");
	dijet_systematic->Draw("HISTPSAME");
	dijet_stat->Draw("HISTPSAME");
	massFrame->Draw("AXISSAME");
	leg->Draw("SAME");
	CMS_lumi(canv,4,11);
	canv->Update();
	canv->RedrawAxis();
	canv->GetFrame()->Draw();
	canv->SetLogx();
	canv->Print("../bin/print/dijet_uncertainty.pdf");
}
