#include "../interface/tdrstyle.C"
#include "../interface/CMS_lumi.C"
using namespace RooFit;

void my_fitTemplates(TString category)
{
	/////////////////////
	// -- Set types -- //
	/////////////////////
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

	Int_t i_data, i_dy, i_ttbar, i_ww, i_wz, i_zz, i_wjet, i_qcd, i_tw;
	for(int i=0; i<types.size(); i++)
	{
		if( types[i] == "Data" ) i_data = i;
		else if( types[i] == "DY" ) i_dy = i;
		else if( types[i] == "ttbar" ) i_ttbar = i;
		else if( types[i] == "WW" ) i_ww = i;
		else if( types[i] == "WZ" ) i_wz = i;
		else if( types[i] == "ZZ" ) i_zz = i;
		else if( types[i] == "WJets" ) i_wjet = i;
		else if( types[i] == "QCD" ) i_qcd = i;
		else if( types[i] == "tW" ) i_tw = i;
	}

	//////////////////////////
	// -- Get histograms -- //
	//////////////////////////
	TFile f_input("bin/histograms/Re_hist_test.root", "read");
	vector<TH1D*> hist_;
	for(int i=0; i<types.size(); i++)
	{
		TH1D *h_med;
		h_med = (TH1D*)f_input.Get("denominator_"+category+"_"+types[i]);
		hist_->push_back(h_med);
	}

	//Check
	//for(int i=0;i<hist_.size();i++) cout << hist_[i]->GetName() << endl;

	////////////////////////////
	// -- Template fitting -- //
	////////////////////////////
	//Initializing
	RooRealVar obs("obs", "TrkIso/p_{T}", 0, 5);
	Double_t N_hist_[types.size()], N_tot = 0, N_data = 0;

	//Convert TH1D to RooDataHist
	RooDataHist *RooHist_Data = new RooDataHist("RooHist_" + types[i_data], "RooHistogram_" + types[i_data], obs, hist_[i_data]);
	RooDataHist *RooHist_QCD = new RooDataHist("RooHist_" + types[i_qcd], "RooHistogram_" + types[i_qcd], obs, hist_[i_qcd]);
	RooDataHist *RooHist_WJets = new RooDataHist("RooHist_" + types[i_wjet], "RooHistogram_" + types[i_wjet], obs, hist_[i_wjet]);
	RooDataHist *RooHist_DY = new RooDataHist("RooHist_" + types[i_dy], "RooHistogram_" + types[i_dy], obs, hist_[i_dy]);
	RooDataHist *RooHist_ttbar = new RooDataHist("RooHist_" + types[i_ttbar], "RooHistogram_" + types[i_ttbar], obs, hist_[i_ttbar]);
	RooDataHist *RooHist_tW = new RooDataHist("RooHist_" + types[i_tw], "RooHistogram_" + types[i_tw], obs, hist_[i_tw]);
	RooDataHist *RooHist_WW = new RooDataHist("RooHist_" + types[i_ww], "RooHistogram_" + types[i_ww], obs, hist_[i_ww]);
	RooDataHist *RooHist_WZ = new RooDataHist("RooHist_" + types[i_wz], "RooHistogram_" + types[i_wz], obs, hist_[i_wz]);
	RooDataHist *RooHist_ZZ = new RooDataHist("RooHist_" + types[i_zz], "RooHistogram_" + types[i_zz], obs, hist_[i_zz]);

	for(int i=0; i<types.size(); i++)
	{
		N_hist_[i] = hist_[i]->Integral(); //cout << "N_" + types[i] + ": " << N_hist_[i] << endl;

		if( types[i] == "Data" )
		{
			N_data = N_hist_[i];
			continue;
		}

		N_tot += N_hist_[i];
	}

	//Convert RooDataHist to RooHistPdf
	RooHistPdf *pdf_QCD = new RooHistPdf("pdf_QCD", "Template from QCD MC", obs, *RooHist_QCD, 0);
	RooHistPdf *pdf_WJets = new RooHistPdf("pdf_WJets", "Template from WJets MC", obs, *RooHist_WJets, 0);
	RooHistPdf *pdf_DY = new RooHistPdf("pdf_DY", "Template from DY MC", obs, *RooHist_DY, 0);
	RooHistPdf *pdf_ttbar = new RooHistPdf("pdf_ttbar", "Template from ttbar MC", obs, *RooHist_ttbar, 0);
	RooHistPdf *pdf_tW = new RooHistPdf("pdf_tW", "Template from tW MC", obs, *RooHist_tW, 0);
	RooHistPdf *pdf_WW = new RooHistPdf("pdf_WW", "Template from WW MC", obs, *RooHist_WW, 0);
	RooHistPdf *pdf_WZ = new RooHistPdf("pdf_WZ", "Template from WZ MC", obs, *RooHist_WZ, 0);
	RooHistPdf *pdf_ZZ = new RooHistPdf("pdf_ZZ", "Template from ZZ MC", obs, *RooHist_ZZ, 0);

	//Construct model
	for(int i=0; i<types.size(); i++)
	{
		if( types[i] == "Data" ) continue;

		N_hist_[i] = N_data * N_hist_[i] / N_tot;
	}

	RooRealVar n_QCD("n_" + types[i_qcd], "n_" + types[i_qcd], N_hist_[i_qcd], N_hist_[i_qcd]*0.95, N_hist_[i_qcd]*1.05);
	RooRealVar n_WJets("n_" + types[i_wjet], "n_" + types[i_wjet], N_hist_[i_wjet], N_hist_[i_wjet]*0.95, N_hist_[i_wjet]*1.05);
	RooRealVar n_DY("n_" + types[i_dy], "n_" + types[i_dy], N_hist_[i_dy], N_hist_[i_dy]*0.95, N_hist_[i_dy]*1.05);
	RooRealVar n_ttbar("n_" + types[i_ttbar], "n_" + types[i_ttbar], N_hist_[i_ttbar], N_hist_[i_ttbar]*0.95, N_hist_[i_ttbar]*1.05);
	RooRealVar n_tW("n_" + types[i_tw], "n_" + types[i_tw], N_hist_[i_tw], N_hist_[i_tw]*0.95, N_hist_[i_tw]*1.05);
	RooRealVar n_WW("n_" + types[i_ww], "n_" + types[i_ww], N_hist_[i_ww], N_hist_[i_ww]*0.95, N_hist_[i_ww]*1.05);
	RooRealVar n_WZ("n_" + types[i_wz], "n_" + types[i_wz], N_hist_[i_wz], N_hist_[i_wz]*0.95, N_hist_[i_wz]*1.05);
	RooRealVar n_ZZ("n_" + types[i_zz], "n_" + types[i_zz], N_hist_[i_zz], N_hist_[i_zz]*0.95, N_hist_[i_zz]*1.05);

  RooAddPdf model( "model", "model", RooArgList(*pdf_QCD, *pdf_WJets, *pdf_DY, *pdf_ttbar, *pdf_tW, *pdf_WW, *pdf_WZ, *pdf_ZZ),
										RooArgList(n_QCD, n_WJets, n_DY, n_ttbar, n_tW, n_WW, n_WZ, n_ZZ) );

	//fit to Data
  RooFitResult* r = model.fitTo( *RooHist_Data, Save() );

  TCanvas *c_fit = new TCanvas("c_fit", "", 800, 800);
  c_fit->cd();

  //Top Pad
	TPad *c1_1 = new TPad("padc1_1","padc1_1",0.01,0.01,0.99,0.99);
	c1_1->Draw();
	c1_1->cd();
	c1_1->SetTopMargin(0.01);
	c1_1->SetBottomMargin(0.25);
	c1_1->SetRightMargin(0.03);
	c1_1->SetLeftMargin(0.09);
	//c1_1->SetFillStyle(0.01);
	c1_1->SetLogy();

	RooPlot* frame1 = obs.frame( Title(" ") ) ;
	// pdf_ttbar->plotOn( frame1, LineColor(kOrange) );
	// pdf_WJets->plotOn(frame1, LineColor(kGreen) );
	RooHist_Data->plotOn(frame1, DataError(RooAbsData::SumW2));

	vector< TString > plots;
	for(int i=0; i<types.size(); i++)
	{
		if( types[i] == "Data" ) continue;

		TString plot_list = "";
		for(int j=i; j<types.size(); j++)
		{
			if( j!=i ) plot_list += ",";
			plot_list += "pdf_"+types[j];
		}
		plots.push_back(plot_list);
	}

	//Color
	vector< Int_t > i_color;
	Int_t array_color[types.size()] = {7,4,2,3,13,14,15};
	for(int i=0; i<types.size(); i++)
		i_color.push_back(array_color[i]);

	for(int i=0; i<plots.size(); i++)
		model.plotOn(frame1, Components(plots[i]), LineColor(0), FillColor(i_color[i]), DrawOption("F") );

	RooHist_Data->plotOn(frame1, DataError(RooAbsData::SumW2));
	frame1->Draw();
	r->Print();

	TLegend *leg1 = new TLegend(0.65,0.7,.95,.97);
	leg1->SetFillColor(kWhite);
	leg1->SetLineColor(kWhite);
	for(int i=0; i<types.size(); i++)
	{
		TString draw_opt = "EP";
		if( types[i] == "Data" ) draw_opt = "F";
		leg1->AddEntry(frame1->nameOf(i),types[i], draw_opt);
	}
	leg1->Draw();

	frame1->GetYaxis()->SetTitle("Entry");
	frame1->GetXaxis()->SetLabelSize(0);

	TH1D *h_MC = (TH1D*)model.createHistogram("h_MC", obs);
	// h_MC->Sumw2();
	// TCanvas *c_MC = new TCanvas("c_MC", "", 700, 700);
	
	Double_t N_MC = h_MC->Integral();
	h_MC->Scale(N_data / N_MC);
	// cout << "# data: " << Ndata << endl;
	// h_MC->Draw(); 
	// h_data->Draw("SAMEEP");

	//Bottom Pad
	TPad *c1_2 = new TPad("padc1_2","padc1_2",0.01,0.01,0.99,0.25);
	c1_2->Draw();
	c1_2->cd();
	c1_2->SetTopMargin(0.1);
	c1_2->SetBottomMargin(0.30);
	c1_2->SetRightMargin(0.02);
	c1_2->SetLeftMargin(0.08);
	c1_2->SetFillStyle(0);
	c1_2->SetGrid();

	//Make ratio plot
	TH1D *h_ratio = (TH1D*)hist_[i_data]->Clone();
	hist_[i_data]->Sumw2(); h_MC->Sumw2();
	h_ratio->Divide(hist_[i_data], h_MC);
	h_ratio->SetTitle("");
	h_ratio->GetXaxis()->SetMoreLogLabels();
	h_ratio->GetXaxis()->SetNoExponent();
	h_ratio->GetXaxis()->SetTitle( "TrkIso/p_{T}" );
	h_ratio->GetYaxis()->SetTitle("data/MC");
	h_ratio->GetXaxis()->SetTitleSize(0.13);
	h_ratio->GetYaxis()->SetTitleSize(0.09);
	h_ratio->GetYaxis()->SetTitleOffset(0.4);
	h_ratio->GetXaxis()->SetLabelSize(0.11);
	h_ratio->GetYaxis()->SetLabelSize(0.07);
	h_ratio->GetYaxis()->SetTickLength(0.015);
	h_ratio->SetMaximum( 1.3 );
	h_ratio->SetMinimum( 0.7 );
	h_ratio->SetMarkerSize(0.3);
	h_ratio->SetStats(kFALSE);

	h_ratio->Draw("e1p");
	
	TH1D *h_line = (TH1D*)hist_[i_data]->Clone();
	h_line->Reset("ICES");
	Int_t Nbins = h_line->GetNbinsX();
	for(Int_t i_bin=0; i_bin< Nbins; i_bin++)
		h_line->SetBinContent(i_bin+1, 1);

	h_line->SetLineColor(kRed);
	h_line->Draw("LSAME");

	//leg1->Draw();

	RooAbsReal *chi2 = model.createChi2(*RooHist_Data);
	cout << "chi2: " << chi2->getVal() << endl;
	cout << "Normalized chi2: " << chi2->getVal() / ((Double_t)hist_[i_data]->GetNbinsX()) << endl;

	c_fit->Print("bin/print/fit_"+category+".pdf");

}
