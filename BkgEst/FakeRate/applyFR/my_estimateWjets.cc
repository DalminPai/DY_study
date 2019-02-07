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

void my_estimateWjets()
{
    int W = 1200;
    int H = 1200;

    int H_ref = 1200;
    int W_ref = 1200;

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;

    // UPDATED IN 2017
    lumi_13TeV = "35.9 fb^{-1}";
    lumiTextSize = 0.5;
    writeExtraText = true;
    extraText = "Preliminary";
    drawLogo = false;

    const int binnum = 43;
    double bins[binnum+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000,1500,3000};

    TH1D* massFrame = new TH1D("massFrame","",38,15,3000);
    massFrame->SetMinimum(0.001);
    massFrame->SetMaximum(1000000);
    massFrame->SetStats(kFALSE);
    massFrame->GetXaxis()->SetTitle("Mass[GeV]");
    massFrame->GetYaxis()->SetTitleOffset(1);
    massFrame->GetYaxis()->SetTitle("Number of events");
    massFrame->GetXaxis()->SetTitleSize(0);
    massFrame->GetYaxis()->SetTitleSize(0.05);
    massFrame->GetXaxis()->SetLabelSize(0);
    massFrame->GetXaxis()->SetMoreLogLabels();

    TFile* f;
    //f = new TFile("../bin/histograms/Re_fake_test.root","READ");
    //f = new TFile("../bin/histograms/Re_fake_test_20181010.root","READ");
    f = new TFile("../bin/histograms/Re_fake_TightID_PFIso_EffSF_20181014.root","READ");

    TH1D* dijet_template[14];
    TH1D* dijetSS_template[14];
    TH1D* dijet_ratio[14];
    TH1D* dijetSS_ratio[14];
    TH1D* wjets_template[14];
    TH1D* wjetsSS_template[14];
    TH1D* wjets_ratio[14];
    TH1D* wjetsSS_ratio[14];

    // opposite sign QCD dijet (template)
    dijet_template[5] = (TH1D*)f->Get("histDijet1_Data")->Clone(); //from Data
    dijet_template[5]->SetFillColor(7);
    dijet_template[5]->SetStats(kFALSE);
    dijet_template[5]->Sumw2();

    dijet_template[0] = (TH1D*)f->Get("histDijet1_DY")->Clone(); //from DY
    dijet_template[0]->SetFillColor(2);
    dijet_template[0]->SetStats(kFALSE);
    dijet_template[0]->Sumw2();

    dijet_template[1] = (TH1D*)f->Get("histDijet1_ttbar")->Clone(); //from ttbar
    dijet_template[1]->SetFillColor(3);
    dijet_template[1]->SetStats(kFALSE);
    dijet_template[1]->Sumw2();

    // same sign QCD dijet (template)
    dijetSS_template[5] = (TH1D*)f->Get("histSameDijet1_Data")->Clone(); //from Data
    //dijetSS_template[5]->SetLineColor(2);
    //dijetSS_template[5]->SetMarkerColor(2);
    //dijetSS_template[5]->SetMarkerStyle(22);
    //dijetSS_template[5]->SetMarkerSize(3);
    dijetSS_template[5]->SetStats(kFALSE);
    dijetSS_template[5]->Sumw2();

    // opposite sign QCD dijet (ratio)
    dijet_ratio[5] = (TH1D*)f->Get("histDijet2_Data")->Clone(); //from Data
    dijet_ratio[5]->SetFillColor(7);
    dijet_ratio[5]->SetStats(kFALSE);
    dijet_ratio[5]->Sumw2();

    dijet_ratio[0] = (TH1D*)f->Get("histDijet2_DY")->Clone(); //from DY
    dijet_ratio[0]->SetFillColor(2);
    dijet_ratio[0]->SetStats(kFALSE);
    dijet_ratio[0]->Sumw2();

    dijet_ratio[1] = (TH1D*)f->Get("histDijet2_ttbar")->Clone(); //from ttbar
    dijet_ratio[1]->SetFillColor(3);
    dijet_ratio[1]->SetStats(kFALSE);
    dijet_ratio[1]->Sumw2();

    // same sign QCD dijet (ratio)
    dijetSS_ratio[5] = (TH1D*)f->Get("histSameDijet2_Data")->Clone(); //from Data
    //dijetSS_ratio[5]->SetLineColor(2);
    //dijetSS_ratio[5]->SetMarkerColor(2);
    //dijetSS_ratio[5]->SetMarkerStyle(22);
    //dijetSS_ratio[5]->SetMarkerSize(3);
    dijetSS_ratio[5]->SetStats(kFALSE);
    dijetSS_ratio[5]->Sumw2();

    // opposite sign Wjet (template)
    wjets_template[5] = (TH1D*)f->Get("histWJets1_Data")->Clone(); //from Data
    wjets_template[5]->SetFillColor(9);
    wjets_template[5]->SetStats(kFALSE);
    wjets_template[5]->Sumw2();

    wjets_template[0] = (TH1D*)f->Get("histWJets1_DY")->Clone(); //from DY
    wjets_template[0]->SetFillColor(2);
    wjets_template[0]->SetStats(kFALSE);
    wjets_template[0]->Sumw2();

    wjets_template[1] = (TH1D*)f->Get("histWJets1_ttbar")->Clone(); //from ttbar
    wjets_template[1]->SetFillColor(3);
    wjets_template[1]->SetStats(kFALSE);
    wjets_template[1]->Sumw2();

    // same sign Wjet (template)
    wjetsSS_template[5] = (TH1D*)f->Get("histSameWJets1_Data")->Clone(); //from Data, default
    //wjetsSS_template[5] = (TH1D*)f->Get("histSameWJets1_WJets")->Clone(); //by DM
    wjetsSS_template[5]->SetFillColor(9);
    wjetsSS_template[5]->SetLineColor(9);
    wjetsSS_template[5]->SetMarkerColor(2);
    wjetsSS_template[5]->SetMarkerStyle(22);
    wjetsSS_template[5]->SetMarkerSize(4);
    wjetsSS_template[5]->SetStats(kFALSE);
    wjetsSS_template[5]->Sumw2();

    wjetsSS_template[1] = (TH1D*)f->Get("histSameWJets1_ttbar")->Clone(); //from ttbar
	wjetsSS_template[1]->Scale(0.708540); //by DM

    // opposite sign Wjet (ratio)
    wjets_ratio[5] = (TH1D*)f->Get("histWJets2_Data")->Clone(); //from Data
    wjets_ratio[5]->SetFillColor(9);
    wjets_ratio[5]->SetStats(kFALSE);
    wjets_ratio[5]->Sumw2();

    wjets_ratio[0] = (TH1D*)f->Get("histWJets2_DY")->Clone(); //from DY
    wjets_ratio[0]->SetFillColor(2);
    wjets_ratio[0]->SetStats(kFALSE);
    wjets_ratio[0]->Sumw2();

    wjets_ratio[1] = (TH1D*)f->Get("histWJets2_ttbar")->Clone(); //from ttbar
    wjets_ratio[1]->SetFillColor(3);
    wjets_ratio[1]->SetStats(kFALSE);
    wjets_ratio[1]->Sumw2();

    // same sign Wjet (ratio)
    wjetsSS_ratio[5] = (TH1D*)f->Get("histSameWJets2_Data")->Clone(); //from Data, default
    //wjetsSS_ratio[5] = (TH1D*)f->Get("histSameWJets2_WJets")->Clone(); //by DM
    wjetsSS_ratio[5]->SetFillColor(9);
    wjetsSS_ratio[5]->SetLineColor(9);
    wjetsSS_ratio[5]->SetMarkerColor(2);
    wjetsSS_ratio[5]->SetMarkerStyle(22);
    wjetsSS_ratio[5]->SetMarkerSize(4);
    wjetsSS_ratio[5]->SetStats(kFALSE);
    wjetsSS_ratio[5]->Sumw2();

    wjetsSS_ratio[1] = (TH1D*)f->Get("histSameWJets2_ttbar")->Clone(); //from ttbar
	wjetsSS_ratio[1]->Scale(0.871850); //by DM

    double norm1[14];
    double norm2[14];

    // -- Wjet fitting results -- //
    // by template:
    //double n_DYJets = 3.4288e+05;
    //double n_QCD    = 1.0054e+04;
    //double n_WJets  = 5.3206e+02;
    //double n_ttbar  = 1.0091e+04;
	// ++ using WJets MC ++
    double n_DYJets = 3.4297e+05;
    double n_QCD    = 6.8905e+03;
    double n_WJets  = 2.9273e+03;
    double n_ttbar  = 1.0815e+04;

    // by ratio:
    //double nn_DYJets = 2.8378e+05;
    //double nn_QCD    = 7.2674e+03;
    //double nn_WJets  = 1.3191e+03;
    //double nn_ttbar  = 8.3471e+03;
	// ++ using WJets MC ++
    double nn_DYJets = 2.8371e+05;
    double nn_QCD    = 4.9082e+03;
    double nn_WJets  = 2.4212e+03;
    double nn_ttbar  = 9.6425e+03;

    // -- Smooth function -- //
    /*wjets_template[0]->Smooth();
    wjets_ratio[0]->Smooth();
    wjets_template[1]->Smooth();
    wjets_ratio[1]->Smooth();
    wjets_template[5]->Smooth();
    wjets_ratio[5]->Smooth();
    dijet_template[5]->Smooth();
    dijet_ratio[5]->Smooth();
    wjetsSS_template[5]->Smooth();
    wjetsSS_ratio[5]->Smooth();*/

    // Data - ( DY + ttbar ) in QCD dijet
    dijet_template[5]->Add(dijet_template[0],-1.0);
    dijet_template[5]->Add(dijet_template[1],-1.0);

    dijet_ratio[5]->Add(dijet_ratio[0],-1.0);
    dijet_ratio[5]->Add(dijet_ratio[1],-1.0);

    // Data - ( 2 * QCD + ttbar ) in same sign Wjet
    wjetsSS_template[5]->Add(dijetSS_template[5],-2.0); //default
    wjetsSS_template[5]->Add(wjetsSS_template[1],-1.0); //default
    cout<<"Same_template = "<<wjetsSS_template[5]->Integral(1,45)<<endl;

    wjetsSS_ratio[5]->Add(dijetSS_ratio[5],-2.0); //default
    wjetsSS_ratio[5]->Add(wjetsSS_ratio[1],-1.0); //default
    cout<<"Same_ratio = "<<wjetsSS_ratio[5]->Integral(1,45)<<endl;


    /*cout<<"ttbar(template): "<<wjets_template[1]->Integral()<<endl;
    cout<<"ttbar(ratio): "<<wjets_ratio[1]->Integral()<<endl;
    cout<<"DY(template): "<<wjets_template[0]->Integral()<<endl;
    cout<<"DY(ratio): "<<wjets_ratio[0]->Integral()<<endl;*/

    // -- Edge? -- //
    //cout<<"Edge="<<wjets_template[0]->GetBinLowEdge(31)<<endl;

    // by template
    norm1[0] = n_DYJets/wjets_template[0]->Integral(1,30);
    norm1[1] = n_ttbar/wjets_template[1]->Integral(1,30);
    norm1[2] = n_WJets/wjetsSS_template[5]->Integral(1,30);
    norm1[5] = n_QCD/dijet_template[5]->Integral(1,30);
    
    // by ratio
    norm2[0] = nn_DYJets/wjets_ratio[0]->Integral(1,30);
    norm2[1] = nn_ttbar/wjets_ratio[1]->Integral(1,30);
    norm2[2] = nn_WJets/wjetsSS_ratio[5]->Integral(1,30);
    norm2[5] = nn_QCD/dijet_ratio[5]->Integral(1,30);

    /*wjets_template[6] = (TH1D*)wjets_template[5]->Clone(); //from Data
    wjets_ratio[6] = (TH1D*)wjets_ratio[5]->Clone(); //from Data

    wjets_template[6]->SetLineColor(1);
    wjets_template[6]->SetMarkerColor(1);
    wjets_template[6]->SetLineWidth(3);
    wjets_template[6]->SetMarkerSize(2);
    wjets_template[6]->SetMarkerStyle(20);

    wjets_ratio[6]->SetLineColor(1);
    wjets_ratio[6]->SetMarkerColor(1);
    wjets_ratio[6]->SetLineWidth(3);
    wjets_ratio[6]->SetMarkerSize(2);
    wjets_ratio[6]->SetMarkerStyle(20);*/

    // -- OS Wjet from SS Wjet template -- //
    wjets_template[5] = (TH1D*)wjetsSS_template[5]->Clone();
    wjets_ratio[5] = (TH1D*)wjetsSS_ratio[5]->Clone();

    wjets_template[0]->Scale(norm1[0]); //DY
    wjets_template[1]->Scale(norm1[1]); //ttbar
    dijet_template[5]->Scale(norm1[5]); //estimated opposite sign QCD dijet
    wjets_template[5]->Scale(norm1[2]); //estimated opposite sign Wjet

    wjets_ratio[0]->Scale(norm2[0]); //DY
    wjets_ratio[1]->Scale(norm2[1]); //ttbar
    dijet_ratio[5]->Scale(norm2[5]); //estimated opposite sign QCD dijet
    wjets_ratio[5]->Scale(norm2[2]); //estimated opposite sign Wjet


    setTDRStyle();
    tdrGrid(true);

    lumiTextSize = 0.6;
    cmsTextSize = 1.0;

    wjets_template[5]->GetXaxis()->SetTitle("Mass[GeV]");
    wjets_template[5]->GetYaxis()->SetTitleOffset(1.5);
    wjets_template[5]->GetYaxis()->SetTitle("Number of events");
    wjets_template[5]->GetXaxis()->SetLabelSize(0.025);
    wjets_template[5]->GetYaxis()->SetLabelSize(0.025);
    wjets_template[5]->GetXaxis()->SetMoreLogLabels();

    wjets_ratio[5]->GetXaxis()->SetTitle("Mass[GeV]");
    wjets_ratio[5]->GetYaxis()->SetTitleOffset(1.5);
    wjets_ratio[5]->GetYaxis()->SetTitle("Number of events");
    wjets_ratio[5]->GetXaxis()->SetLabelSize(0.025);
    wjets_ratio[5]->GetYaxis()->SetLabelSize(0.025);
    wjets_ratio[5]->GetXaxis()->SetMoreLogLabels();

    wjetsSS_template[5]->GetXaxis()->SetTitle("Mass[GeV]");
    wjetsSS_template[5]->GetYaxis()->SetTitleOffset(1.5);
    wjetsSS_template[5]->GetYaxis()->SetTitle("Number of events");
    wjetsSS_template[5]->GetXaxis()->SetLabelSize(0.025);
    wjetsSS_template[5]->GetYaxis()->SetLabelSize(0.025);
    wjetsSS_template[5]->GetXaxis()->SetMoreLogLabels();

    wjetsSS_ratio[5]->GetXaxis()->SetTitle("Mass[GeV]");
    wjetsSS_ratio[5]->GetYaxis()->SetTitleOffset(1.5);
    wjetsSS_ratio[5]->GetYaxis()->SetTitle("Number of events");
    wjetsSS_ratio[5]->GetXaxis()->SetLabelSize(0.025);
    wjetsSS_ratio[5]->GetYaxis()->SetLabelSize(0.025);
    wjetsSS_ratio[5]->GetXaxis()->SetMoreLogLabels();

    TLegend* legg = new TLegend(.6,.65,.95,.89);
    legg->AddEntry(wjets_template[5],"Opposite sign", "F");
    legg->AddEntry(wjetsSS_template[5],"Same sign", "P");

    // Remove negative bins in Wjet
    for(int i=1; i<46; i++)
    {
        if( wjets_template[5]->GetBinContent(i) < 0 )
        {
            wjets_template[5]->SetBinContent(i,0.0);
            wjets_template[5]->SetBinError(i,0.0);
        }

        if( wjets_ratio[5]->GetBinContent(i) < 0 )
        {
            wjets_ratio[5]->SetBinContent(i,0.0);
            wjets_ratio[5]->SetBinError(i,0.0);
        }

        if( wjetsSS_template[5]->GetBinContent(i) < 0 )
        {
            wjetsSS_template[5]->SetBinContent(i,0.0);
            wjetsSS_template[5]->SetBinError(i,0.0);
        }

        if( wjetsSS_ratio[5]->GetBinContent(i) < 0 ) //added by DM
        {
            wjetsSS_ratio[5]->SetBinContent(i,0.0);
            wjetsSS_ratio[5]->SetBinError(i,0.0);
        }
    }

    // Smooth function
    //wjets_template[5]->Smooth();
    //wjets_ratio[5]->Smooth();
    //wjetsSS_template[5]->Smooth();

	//////////////////////
	// -- Make plots -- //
	//////////////////////
    TCanvas* canv = new TCanvas("canv","",1200,1200);
    canv->SetFillColor(0);
    canv->SetLeftMargin( L/W );
    canv->SetRightMargin( R/W );
    canv->SetTopMargin( T/H );
    canv->SetBottomMargin( B/H );

	// -- Opposite sign WJets -- //
	TLegend *leg_wjets = new TLegend(0.6,0.8,.95,.85);
	leg_wjets->SetBorderSize(0);
	leg_wjets->SetFillStyle(0);
	leg_wjets->AddEntry(wjetsSS_template[5],"W+jets (Opposite sign)","F");

	// -- Same sign WJets -- //
	TLegend *leg_wjets_SS = new TLegend(0.6,0.8,.95,.85);
	leg_wjets_SS->SetBorderSize(0);
	leg_wjets_SS->SetFillStyle(0);
	leg_wjets_SS->AddEntry(wjetsSS_template[5],"W+jets (Same sign)","F");

    wjets_template[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    leg_wjets->Draw("SAME");
    canv->Print("../bin/print/wjets_template.pdf");
    wjetsSS_template[5]->Draw("HISTPSAME");
    legg->Draw("SAME");
    canv->Print("../bin/print/wjetsBoth_template.pdf");
    canv->Clear();

    wjetsSS_template[5]->SetFillColor(9);
    wjetsSS_template[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
	leg_wjets_SS->Draw("same");
    canv->Print("../bin/print/wjetsSS_template.pdf");
    canv->Clear();

    wjets_ratio[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    leg_wjets->Draw("SAME");
    canv->Print("../bin/print/wjets_ratio.pdf");
    legg->Draw("SAME");
    wjetsSS_ratio[5]->Draw("HISTPSAME");
    canv->Print("../bin/print/wjetsBoth_ratio.pdf");
    canv->Clear();

    wjetsSS_ratio[5]->SetFillColor(9);
    wjetsSS_ratio[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
	leg_wjets_SS->Draw("same");
    canv->Print("../bin/print/wjetsSS_ratio.pdf");
    canv->Clear();


    double error = 0;
    wjets_template[5]->IntegralAndError(1,45,error);
    cout<<"QCD(template) = "<<wjets_template[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;
    wjetsSS_template[5]->IntegralAndError(1,45,error);
    cout<<"QCD(template) SS = "<<wjetsSS_template[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;
    wjets_ratio[5]->IntegralAndError(1,45,error);
    cout<<"QCD(ratio) = "<<wjets_ratio[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;
    wjetsSS_ratio[5]->IntegralAndError(1,45,error);
    cout<<"QCD(ratio) SS = "<<wjetsSS_ratio[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;

    TH1D* wjets = (TH1D*)wjets_template[5]->Clone();
    TH1D* wjets_control = (TH1D*)wjets_ratio[5]->Clone();
    wjets->SetName("wjets");

    // -- Uncertainty calculation -- //
    TH1D* wjets_total      = new TH1D("wjets_total","",binnum,bins);
    TH1D* wjets_systematic = new TH1D("wjets_systematic","",binnum,bins);
    TH1D* wjets_stat       = new TH1D("wjets_stat","",binnum,bins);

    for(int i=1; i<binnum+1; i++)
    {
        double systematic = fabs( wjets->GetBinContent(i) - wjets_control->GetBinContent(i) );
        double stat = wjets->GetBinError(i);
        double total = sqrt( systematic*systematic + stat*stat );
        
        if(wjets->GetBinContent(i)==0)
        {
            systematic = 0;
            stat = 0;
            total = 0;
        }

        wjets_systematic->SetBinContent(i,systematic);
        wjets_stat->SetBinContent(i,stat);
        wjets_total->SetBinContent(i,total);

        wjets->SetBinError(i,total);
    }

    // Save the result
    TFile* gg = new TFile("../bin/result/wjets.root","RECREATE");
    wjets->Write();
    wjets_systematic->Write();
    wjets_stat->Write();
    gg->Close();

    // Uncertainty plot
    wjets_systematic->Divide(wjets);
    wjets_stat->Divide(wjets);
    wjets_total->Divide(wjets);

    wjets_total->SetMarkerStyle(20);
    wjets_total->SetMarkerSize(3);
    wjets_total->SetMarkerColor(1);

    wjets_systematic->SetMarkerStyle(22);
    wjets_systematic->SetMarkerSize(3);
    wjets_systematic->SetMarkerColor(2);

    wjets_stat->SetMarkerStyle(21);
    wjets_stat->SetMarkerSize(3);
    wjets_stat->SetMarkerColor(4);

    TLegend* leg = new TLegend(.8,.15,.95,.3);
    leg->AddEntry(wjets_total,"Total","P");
    leg->AddEntry(wjets_systematic,"Sys.","P");
    leg->AddEntry(wjets_stat,"Stat.","P");

    massFrame->GetYaxis()->SetTitle("Unceratinty");
    massFrame->SetMinimum(0);
    massFrame->SetMaximum(1);
    massFrame->GetXaxis()->SetLabelSize(0.025);
    massFrame->GetXaxis()->SetTitleSize(0.05);

    massFrame->Draw();
    wjets_total->Draw("HISTPSAME");
    wjets_systematic->Draw("HISTPSAME");
    wjets_stat->Draw("HISTPSAME");
    massFrame->Draw("AXISSAME");
    leg->Draw("SAME");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("../bin/print/wjets_uncertainty.pdf");
}
