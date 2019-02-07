#include "../../interface/CMS_lumi.C"
#include "../../interface/MyCanvas.C"
using namespace RooFit;

void my_fit_wjets(TString type = "template")
{
	TString n_type;
	if( type == "template" )
		n_type = "1";
	else if( type == "ratio" )
		n_type = "2";
	else
		cout << "ERROR: Wron type!!" << endl;

	// Get ROOT Files
	//TFile *f = new TFile("../bin/histograms/Re_fake_test.root", "read");
	//TFile *f = new TFile("../bin/histograms/Re_fake_test_20181010.root", "read");
	TFile *f = new TFile("../bin/histograms/Re_fake_TightID_PFIso_EffSF_20181014.root","READ");

	// Get Histograms
	TH1D *h_data = (TH1D*)f->Get( "fitWJets"+n_type+"_Data" );
	TH1D *h_ttbar = (TH1D*)f->Get( "fitWJets"+n_type+"_ttbar" );
	TH1D *h_DYJets = (TH1D*)f->Get( "fitWJets"+n_type+"_DY" );

	// default
	/*TH1D *h_WJets = (TH1D*)f->Get( "fitSameWJets"+n_type+"_Data" );
	TH1D *h_ttbar_SS = (TH1D*)f->Get( "fitSameWJets"+n_type+"_ttbar" );
	TH1D *h_QCD_SS = (TH1D*)f->Get( "fitSameDijet"+n_type+"_Data" );
	h_WJets->Add(h_QCD_SS,-2.0);
	h_WJets->Add(h_ttbar_SS,-1.0);*/

	// by DM
	TH1D* h_WJets = (TH1D*)f->Get( "fitSameWJets"+n_type+"_WJets" );



	TH1D *h_QCD = (TH1D*)f->Get( "fitDijet"+n_type+"_Data" );
	TH1D *h_DYJets_Dijet = (TH1D*)f->Get( "fitDijet"+n_type+"_DY" );
	TH1D *h_ttbar_Dijet = (TH1D*)f->Get( "fitDijet"+n_type+"_ttbar" );
	h_QCD->Add(h_DYJets_Dijet,-1.0);
	h_QCD->Add(h_ttbar_Dijet,-1.0);

	// Remove negative bins: added by DM (only for default method)
    /*for(int i=1; i<50; i++)
    {
        if( h_WJets->GetBinContent(i) < 0 )
        {
            h_WJets->SetBinContent(i,0.0);
            h_WJets->SetBinError(i,0.0);
        }

		//if( h_QCD->GetBinContent(i) < 0 )
		//{
        //    h_QCD->SetBinContent(i,0.0);
        //    h_QCD->SetBinError(i,0.0);
		//}
	}*/

	// Convert TH1D to RooDataHist
	RooRealVar mass("mass", "Dimuon mass [GeV]", 15,200);

	RooDataHist *RooHist_ttbar = new RooDataHist("RooHist_ttbar", "RooHistogram_ttbar", mass, h_ttbar);
	RooDataHist *RooHist_DYJets = new RooDataHist("RooHist_DYJets", "RooHistogram_DYJets", mass, h_DYJets);
	RooDataHist *RooHist_WJets = new RooDataHist("RooHist_WJets", "RooHistogram_WJets", mass, h_WJets);
	RooDataHist *RooHist_QCD = new RooDataHist("RooHist_QCD", "RooHistogram_QCD", mass, h_QCD);
	RooDataHist *RooHist_data = new RooDataHist("RooHist_data", "RooHistogram_data", mass, h_data);

	// Convert RooDataHist to RooHistPdf
	RooHistPdf *pdf_ttbar = new RooHistPdf("pdf_ttbar", "Template from ttbar MC", mass, *RooHist_ttbar, 0);
	RooHistPdf *pdf_DYJets = new RooHistPdf("pdf_DYJets", "Template from DYJets MC", mass, *RooHist_DYJets, 0);
	RooHistPdf *pdf_WJets = new RooHistPdf("pdf_WJets", "Template from same-sign WJets", mass, *RooHist_WJets, 0);
	RooHistPdf *pdf_QCD = new RooHistPdf("pdf_QCD", "Template from data-driven dijet", mass, *RooHist_QCD, 0);

	// Construct model = n_ttbar * ttbar + n_WJets * WJets
	double N_ttbar = h_ttbar->Integral();
	double N_DYJets = h_DYJets->Integral();
	double N_WJets = h_WJets->Integral()*3;
	double N_QCD = h_QCD->Integral()*2;
	double ratio = h_data->Integral()/(N_ttbar + N_DYJets + N_WJets + N_QCD);
	N_ttbar *= ratio;
	N_DYJets *= ratio;
	N_WJets *= ratio;
	N_QCD *= ratio;

	RooRealVar n_ttbar("n_ttbar", "n_ttbar", N_ttbar, N_ttbar*0.9, N_ttbar*1.1);
	RooRealVar n_DYJets("n_DYJets", "n_DYJets", N_DYJets, N_DYJets*0.9, N_DYJets*1.1);
	RooRealVar n_WJets("n_WJets", "n_WJets", N_WJets, N_WJets*0.7, N_WJets*1.3);
	RooRealVar n_QCD("n_QCD", "n_QCD", N_QCD, N_QCD*0.7, N_QCD*1.3);
	RooAddPdf model( "model", "model", RooArgList(*pdf_WJets, *pdf_QCD, *pdf_DYJets, *pdf_ttbar), RooArgList(n_WJets, n_QCD, n_DYJets, n_ttbar) );

	// Fit to data
	RooFitResult* r = model.fitTo( *RooHist_data, Save() );
	r->Print();
	RooAbsReal *chi2 = model.createChi2(*RooHist_data);
	cout << "chi2ndof: " << chi2->getVal() / ((Double_t)h_data->GetNbinsX()) << endl;

	/////////////////////
	// -- Make plot -- //
	/////////////////////
	int W = 1200;
	int H = 1200;

	int W_ref = 1200;
	int H_ref = 1200;

	// references for T, B, L, R
	float T = 0.08*H_ref;
	//float T = 0.12*H_ref;
	float B = 0.12*H_ref;
	float L = 0.12*W_ref;
	float R = 0.04*W_ref;

	lumi_13TeV = "35.9 fb^{-1}";
	writeExtraText = true;
	extraText = "Preliminary";
	drawLogo = false;

	TCanvas *c_fit = new TCanvas("c_fit", "", 1200, 1200);
	c_fit->cd();
	c_fit->SetFillColor(0);
	c_fit->SetLeftMargin( L/W );
	c_fit->SetRightMargin( R/W );
	c_fit->SetTopMargin( T/H );
	c_fit->SetBottomMargin( B/H );

	//Top Pad
	TPad *c1_1 = new TPad("padc1_1","padc1_1",0.01,0.01,0.99,0.99);
	c1_1->Draw();
	c1_1->cd();
	//c1_1->SetTopMargin(0.01);
	c1_1->SetTopMargin(0.07);
	c1_1->SetBottomMargin(0.25);
	c1_1->SetRightMargin(0.03);
	c1_1->SetLeftMargin(0.11);
	//c1_1->SetFillStyle(0.01);
	//c1_1->SetLogx();
	c1_1->SetLogy();
	c1_1->SetTicks();

	RooPlot* frame1 = mass.frame( Title(" ") ) ;
	RooHist_data->plotOn(frame1, DataError(RooAbsData::SumW2));
	model.plotOn(frame1, Components("pdf_WJets,pdf_ttbar,pdf_QCD,pdf_DYJets"), LineColor(0), FillColor(2), DrawOption("F") );
	//model.plotOn(frame1, Components("pdf_WJets,pdf_ttbar,pdf_QCD"), LineColor(0), FillColor(7), DrawOption("F") );
	//model.plotOn(frame1, Components("pdf_WJets,pdf_ttbar"), LineColor(0), FillColor(3), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_WJets,pdf_QCD,pdf_ttbar"), LineColor(0), FillColor(3), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_WJets,pdf_QCD"), LineColor(0), FillColor(7), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_WJets"), LineColor(0), FillColor(9), DrawOption("F") );
	RooHist_data->plotOn(frame1, DataError(RooAbsData::SumW2));
	frame1->Draw();
	r->Print();

	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(2);
	//h_data->SetMarkerSize(1.5);
	h_data->Draw("SAMEEP");

	TLegend *leg1 = new TLegend(0.7,0.7,.95,.9);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg1->AddEntry(frame1->nameOf(0),"Data", "EP");
	leg1->AddEntry(frame1->nameOf(1),"DY","F");
	//leg1->AddEntry(frame1->nameOf(2),"QCD","F");
	//leg1->AddEntry(frame1->nameOf(3),"ttbar","F");
	leg1->AddEntry(frame1->nameOf(2),"ttbar","F");
	leg1->AddEntry(frame1->nameOf(3),"QCD","F");
	leg1->AddEntry(frame1->nameOf(4),"W+jets","F");
	leg1->Draw();

	frame1->GetXaxis()->SetLabelSize(0);
	frame1->GetXaxis()->SetTitle("");
	frame1->GetXaxis()->SetTitleSize(0);
	frame1->GetYaxis()->SetTitle("Number of events");
	frame1->GetXaxis()->SetRangeUser(0, 0.5);
	//frame1->GetYaxis()->SetRangeUser(0.1, 1e8);
	frame1->GetYaxis()->SetRangeUser(10, 1e7);
	frame1->GetYaxis()->SetTitleOffset(1.2);
	frame1->GetYaxis()->SetTitleSize(0.04);
	frame1->GetYaxis()->SetLabelSize(0.025);

	TH1D *h_MC = (TH1D*)model.createHistogram("h_MC", mass);
	// h_MC->Sumw2();
	// TCanvas *c_MC = new TCanvas("c_MC", "", 700, 700);
	
	Double_t Ndata = h_data->Integral();
	Double_t NMC = h_MC->Integral();
	h_MC->Scale( Ndata / NMC );
	// cout << "# data: " << Ndata << endl;
	// h_MC->Draw(); 

	//Bottom Pad
	TPad *c1_2 = new TPad("padc1_2","padc1_2",0.01,0.01,0.99,0.25);
	c1_2->Draw();
	c1_2->cd();
	c1_2->SetTopMargin(0.1);
	c1_2->SetBottomMargin(0.30);
	c1_2->SetRightMargin(0.02);
	//c1_2->SetLeftMargin(0.08);
	c1_2->SetLeftMargin(0.102);
	c1_2->SetFillStyle(0);
	c1_2->SetGrid();
	//c1_2->SetLogx();

	//Make ratio plot
	TH1D *h_ratio = (TH1D*)h_data->Clone();
	h_data->Sumw2(); h_MC->Sumw2();
	h_ratio->Divide(h_data, h_MC);
	h_ratio->SetTitle("");
	h_ratio->GetXaxis()->SetMoreLogLabels();
	h_ratio->GetXaxis()->SetNoExponent();
	h_ratio->GetXaxis()->SetTitle( "Mass [GeV]" );
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetXaxis()->SetTitleSize(0.13);
	h_ratio->GetYaxis()->SetTitleSize(0.09);
	h_ratio->GetYaxis()->SetTitleOffset(0.4);
	h_ratio->GetXaxis()->SetLabelSize(0.11);
	h_ratio->GetYaxis()->SetLabelSize(0.07);
	h_ratio->GetYaxis()->SetTickLength(0.015);
	h_ratio->GetXaxis()->SetRangeUser(0, 0.5);
	h_ratio->SetMaximum( 1.3 );
	h_ratio->SetMinimum( 0.7 );
	h_ratio->SetMarkerStyle(20);
	//h_ratio->SetMarkerSize(2);
	h_ratio->SetMarkerSize(1.5);
	h_ratio->SetStats(kFALSE);

	//h_ratio->Draw("e1p");
	h_ratio->Draw("EP");
	
	TH1D *h_line = (TH1D*)h_data->Clone();
	h_line->Reset("ICES");
	Int_t Nbins = h_line->GetNbinsX();
	for(Int_t i_bin=0; i_bin< Nbins; i_bin++)
		h_line->SetBinContent(i_bin+1, 1);

	h_line->SetLineColor(kRed);
	h_line->Draw("LSAME");

	CMS_lumi(c_fit,4,11);
	c_fit->Update();
	c_fit->RedrawAxis();
	//c_fit->GetFrame()->Draw();
	c_fit->Print("../bin/print/fit_wjets_"+type+".pdf");
}
