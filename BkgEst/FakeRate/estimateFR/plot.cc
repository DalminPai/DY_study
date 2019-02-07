//#include "../../interface/tdrstyle.C"
#include "../../interface/CMS_lumi.C"
#include "../../interface/MyCanvas.C"

// -- Plots for ratio method -- //
// ex category: "denominator_pt_barrel"
//              "denominator_pt_endcap"
//              "numerator_pt_barrel"
//              "numerator_pt_endcap"
void plot(const TString& category)
{
	Int_t icolor[9] = {0, 7, 4, 2, 3, 5, 20, 30, 40};

	//Get ROOT Files
	//TFile *f = new TFile("../bin/histograms/Re_hist_test.root");
	//TFile *f = new TFile("../bin/histograms/Re_hist_test_20181010.root");
	TFile *f = new TFile("../bin/histograms/Re_hist_TightID_PFIso_20181014.root");

	//Get Histograms
	TH1D *h_data = (TH1D*)f->Get( category + "_Data" );
	TH1D *h_QCD = (TH1D*)f->Get( category + "_QCD" );
	TH1D *h_WJets = (TH1D*)f->Get( category + "_WJets" );
	TH1D *h_DYJets = (TH1D*)f->Get( category + "_DY" );
	TH1D *h_ttbar = (TH1D*)f->Get( category + "_ttbar" );
	TH1D *h_tW = (TH1D*)f->Get( category + "_tW" );
	TH1D *h_WW = (TH1D*)f->Get( category + "_WW" );
	TH1D *h_WZ = (TH1D*)f->Get( category + "_WZ" );
	TH1D *h_ZZ = (TH1D*)f->Get( category + "_ZZ" );

	// Rebin
	if( category.Contains("pt") )
	{
		Int_t nRebin = 6;
		h_data->Rebin(nRebin);
		h_QCD->Rebin(nRebin);
		h_WJets->Rebin(nRebin);
		h_DYJets->Rebin(nRebin);
		h_ttbar->Rebin(nRebin);
		h_tW->Rebin(nRebin);
		h_WW->Rebin(nRebin);
		h_WZ->Rebin(nRebin);
		h_ZZ->Rebin(nRebin);
	}
	
	// Set Marker
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.8);

	// Set Fill Color
	h_QCD->SetFillColor(icolor[1]);
	h_WJets->SetFillColor(icolor[2]);
	h_DYJets->SetFillColor(icolor[3]);
	h_ttbar->SetFillColor(icolor[4]);
	h_tW->SetFillColor(icolor[5]);
	h_WW->SetFillColor(icolor[6]);
	h_WZ->SetFillColor(icolor[7]);
	h_ZZ->SetFillColor(icolor[8]);

	// Set Line Color
	h_QCD->SetLineColor(icolor[1]);
	h_WJets->SetLineColor(icolor[2]);
	h_DYJets->SetLineColor(icolor[3]);
	h_ttbar->SetLineColor(icolor[4]);
	h_tW->SetLineColor(icolor[5]);
	h_WW->SetLineColor(icolor[6]);
	h_WZ->SetLineColor(icolor[7]);
	h_ZZ->SetLineColor(icolor[8]);

	// No Stats
	h_data->SetStats(kFALSE);
	h_QCD->SetStats(kFALSE);
	h_WJets->SetStats(kFALSE);
	h_DYJets->SetStats(kFALSE);
	h_ttbar->SetStats(kFALSE);
	h_tW->SetStats(kFALSE);
	h_WW->SetStats(kFALSE);
	h_WZ->SetStats(kFALSE);
	h_ZZ->SetStats(kFALSE);

	// Stack histograms
	THStack* h_stack = new THStack("h_stack","");
	h_stack->Add(h_ZZ);
	h_stack->Add(h_WZ);
	h_stack->Add(h_WW);
	h_stack->Add(h_tW);
	h_stack->Add(h_ttbar);
	h_stack->Add(h_DYJets);
	h_stack->Add(h_WJets);
	h_stack->Add(h_QCD);

	TLegend *leg1 = new TLegend(0.7,0.65,.95,.9);
	//leg1->SetFillColor(kWhite);
	//leg1->SetLineColor(kWhite);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg1->AddEntry(h_data,"Data", "EP");
	leg1->AddEntry(h_QCD,"QCD","F");
	leg1->AddEntry(h_WJets,"WJets","F");
	leg1->AddEntry(h_DYJets,"DYJets","F");
	leg1->AddEntry(h_ttbar,"ttbar","F");
	leg1->AddEntry(h_tW,"tW","F");
	leg1->AddEntry(h_WW,"WW","F");
	leg1->AddEntry(h_WZ,"WZ","F");
	leg1->AddEntry(h_ZZ,"ZZ","F");

	TH1D* h_top = (TH1D*)h_ttbar->Clone();
	h_top->Add(h_tW);
	TH1D* h_VV = (TH1D*)h_WW->Clone();
	h_VV->Add(h_WZ);
	h_VV->Add(h_ZZ);

	THStack* h_stack2 = new THStack("h_stack2","");
	h_stack2->Add(h_VV);
	h_stack2->Add(h_top);
	h_stack2->Add(h_DYJets);
	h_stack2->Add(h_WJets);
	h_stack2->Add(h_QCD);

	TLegend *leg2 = new TLegend(0.7,0.67,.95,.92);
	//leg2->SetFillColor(kWhite);
	//leg2->SetLineColor(kWhite);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->AddEntry(h_data,"Data", "EP");
	leg2->AddEntry(h_QCD,"QCD","F");
	leg2->AddEntry(h_WJets,"WJets","F");
	leg2->AddEntry(h_DYJets,"DYJets","F");
	leg2->AddEntry(h_ttbar,"ttbar+tW","F");
	leg2->AddEntry(h_WW,"VV","F");

	MyCanvas *myc = new MyCanvas("../bin/print/"+category, "P_{T}(#mu)", "Number of events");
    myc->SetLogy(0);
    //myc->SetYRange(0.1, 1e8);
    myc->SetYRange(10, 5e8);
    myc->SetXRange(50, 500);
    myc->SetRatioRange(0.7, 1.3);
    myc->CanvasWithTHStackRatioPlot( h_data, h_stack2, leg2, "Data/MC", 0);
    myc->PrintCanvas();

}

