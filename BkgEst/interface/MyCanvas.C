#pragma once

#include <TLegend.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPaletteAxis.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TPad.h>
#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <vector>

#include "tdrstyle.C"

const int binnum = 43;
const Double_t bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

using namespace std;


TH1D* DYMassBin( TH1D* h_1D );
TH1D* AbsBin( TH1D* h_1D );
TH2D* doTranspose( TH2D* h_2D );
void ErrorCalc( TH1D* h1, TH1D* h2, TH1D* h3, Int_t nBins );

class MyCanvas
{
public:
	TCanvas *c;
	TLegend *legend;
	TPad *TopPad;
	TPad *BottomPad;

	TString TitleX;
	TString TitleY;
	Bool_t isLogX;
	Bool_t isLogY;
	Bool_t isLogZ;
	Bool_t isSetNoExpo_MoreLogLabels_X;
	Bool_t isSetNoExpo_MoreLogLabels_Y;
	
	Double_t LowerEdge_X;
	Double_t UpperEdge_X;
	Double_t LowerEdge_Y;
	Double_t UpperEdge_Y;
	Double_t LowerEdge_Z; // -- for 2D histogram -- //
	Double_t UpperEdge_Z;

	Double_t Legend_x1;
	Double_t Legend_y1;
	Double_t Legend_x2;
	Double_t Legend_y2;

	// -- for multiple histograms & graphs -- //
	vector< Int_t > Colors;

	// -- For Ratio plot -- //
	TH1D *h_ratio;
	TH1D *h_ratio1; TH1D *h_ratio2; TH1D *h_ratio3; TH1D *h_ratio4; // -- for the canvas with multiple ratio plot -- //
	TGraphAsymmErrors *g_ratio;
	Double_t LowerEdge_Ratio;
	Double_t UpperEdge_Ratio;

	Bool_t DoTranspose; // -- for 2D histogram -- //

	MyCanvas(TString CanvasName, TString XTitle, TString YTitle, Int_t length_x = 800, Int_t length_y = 800)
	{
		setTDRStyle();
		// tdrStyle->SetErrorX(0.5);
		c = new TCanvas(CanvasName, "", length_x, length_y);

		TitleX = XTitle;
		TitleY = YTitle;
		isLogX = kFALSE;
		isLogY = kFALSE;
		isLogZ = kFALSE;
		LowerEdge_X = 0;
		UpperEdge_X = 0;
		LowerEdge_Y = 0;
		UpperEdge_Y	= 0;
		LowerEdge_Z = 0;
		UpperEdge_Z	= 0;
		Legend_x1 = 0.65;
		Legend_y1 = 0.85;
		Legend_x2 = 0.95;
		Legend_y2 = 0.95;

		LowerEdge_Ratio = 0.5;
		UpperEdge_Ratio = 1.5;
		isSetNoExpo_MoreLogLabels_X = kTRUE;
		isSetNoExpo_MoreLogLabels_Y = kFALSE;

		for(Int_t i=1; i<=9; i++)
			//Colors.push_back( i );
			Colors.push_back( i+1 );
		Colors.push_back( kOrange );
		Colors.push_back( kViolet );

		DoTranspose = kFALSE;
	}

	~MyCanvas()
	{
		cout << "Call destructor.." << endl;
	}

	void SetLegendPosition(Double_t _x1, Double_t _y1, Double_t _x2, Double_t _y2)
	{
		Legend_x1 = _x1;
		Legend_y1 = _y1;
		Legend_x2 = _x2;
		Legend_y2 = _y2;
	}

	void SetTranspose( Bool_t _DoTranspose )
	{
		DoTranspose = _DoTranspose;
		cout << "DoTranspose: " << DoTranspose << endl;
	}

	void PrintVariables()
	{
		cout << "TitleX: " << TitleX << endl;
		cout << "TitleY: " << TitleY << endl;
		cout << "isLogX: " << isLogX << endl;
		cout << "isLogY: " << isLogY << endl;
		cout << "LowerEdge_X: " << LowerEdge_X << endl;
		cout << "UpperEdge_X: " << UpperEdge_X << endl;
		cout << "LowerEdge_Y: " << LowerEdge_Y << endl;
		cout << "UpperEdge_Y: " << UpperEdge_Y << endl;
		cout << "LowerEdge_Ratio: " << LowerEdge_Ratio << endl;
		cout << "UpperEdge_Ratio: " << UpperEdge_Ratio << endl;
	}

	void SetLogx(Bool_t isNoExpo_MoeLogLabels = kTRUE)
	{
		isLogX = kTRUE;
		if( isNoExpo_MoeLogLabels == kTRUE ) isSetNoExpo_MoreLogLabels_X = kTRUE;
	}

	void SetLogy(Bool_t isNoExpo_MoeLogLabels = kTRUE)
	{
		isLogY = kTRUE;
		if( isNoExpo_MoeLogLabels == kTRUE ) isSetNoExpo_MoreLogLabels_Y = kTRUE;
	}

	void SetLogz()
	{
		isLogZ = kTRUE;
	}

	void SetXRange(Double_t min, Double_t max)
	{
		LowerEdge_X = min;
		UpperEdge_X = max;
	}

	void SetYRange(Double_t min, Double_t max)
	{
		LowerEdge_Y = min;
		UpperEdge_Y = max;
	}

	void SetZRange(Double_t min, Double_t max)
	{
		LowerEdge_Z = min;
		UpperEdge_Z = max;
	}

	void SetRatioRange(Double_t min, Double_t max)
	{
		LowerEdge_Ratio = min;
		UpperEdge_Ratio = max;
	}

	void CanvasWithHistogramsRatioPlot(TH1D *h1, TH1D *h2,
										TString Name1, TString Name2, TString Name_Ratio,
										Int_t color1 = kBlack, Int_t color2 = kRed,
										Bool_t isFill1 = kFALSE, Bool_t isFill2 = kFALSE,
										TString DrawOp1 = "EP", TString DrawOp2 = "EPSAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		gPad->SetGridx(0);
		gPad->SetGridy(0);

		// -- ensure additional space at the bottom side for ratio plot -- //
		TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		h1->Draw(DrawOp1);
		h2->Draw(DrawOp2);
		h1->Draw("AXISSAME");
		h1->Draw(DrawOp1+"SAME");

		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(color1);
		h1->SetLineWidth(1);
		h1->SetMarkerColor(color1);
		h1->SetMarkerSize(1);
		h1->SetMarkerStyle(20);
		if( isFill1 == kTRUE ) h1->SetFillColor(color1);
		else h1->SetFillColorAlpha(kWhite, 0);

		// -- Setting for h2 -- //
		h2->SetLineColor(color2);
		h2->SetLineWidth(1);
		h2->SetMarkerColor(color2);
		h2->SetMarkerSize(1);
		h2->SetMarkerStyle(20);
		if( isFill2 == kTRUE ) h2->SetFillColor(color2);
		else h2->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetLabelSize(0);
		h1->GetXaxis()->SetTitleSize(0);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(h1, Name1);
		legend->AddEntry(h2, Name2);
		legend->Draw();


		// -- Bottom Pad -- //
		BottomPad = new TPad("BottomPad","BottomPad",0.01,0.01,0.99,0.3);
		BottomPad->Draw();
		BottomPad->cd();

		BottomPad->SetBottomMargin(0.4);
		BottomPad->SetRightMargin(0.04);
		BottomPad->SetLeftMargin(0.15);	

		if( isLogX == kTRUE ) gPad->SetLogx();

		//gPad->SetGridx(0);
		//gPad->SetGridy(0);

		// -- Make Ratio plot & Draw it -- //
		h_ratio = (TH1D*)h1->Clone();
		h_ratio->SetName("h_ratio");
		h1->Sumw2(); h2->Sumw2();
		h_ratio->Divide(h1, h2);
		h_ratio->Draw("EP");


		// -- General Setting -- //
		//h_ratio->SetLineColor(kBlack);
		h_ratio->SetLineColor(kBlue);
		h_ratio->SetMarkerStyle(20);
		h_ratio->SetMarkerSize(1);
		//h_ratio->SetMarkerColor(kBlack);
		h_ratio->SetMarkerColor(kBlue);
		h_ratio->SetStats(kFALSE);

		// -- X-axis Setting -- //		
		h_ratio->GetXaxis()->SetTitle( TitleX );
		h_ratio->GetXaxis()->SetTitleOffset( 0.9 );
		h_ratio->GetXaxis()->SetTitleSize( 0.2 );
		h_ratio->GetXaxis()->SetLabelColor(1);
		h_ratio->GetXaxis()->SetLabelFont(42);
		h_ratio->GetXaxis()->SetLabelOffset(0.007);
		h_ratio->GetXaxis()->SetLabelSize(0.15);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h_ratio->GetXaxis()->SetMoreLogLabels(); h_ratio->GetXaxis()->SetNoExponent(); }

		// -- Y-axis Setting -- //
		h_ratio->GetYaxis()->SetTitle( Name_Ratio );
		h_ratio->GetYaxis()->SetTitleOffset( 0.4 );
		h_ratio->GetYaxis()->SetTitleSize( 0.1);
		h_ratio->GetYaxis()->SetLabelSize( 0.07 );
		h_ratio->GetYaxis()->SetRangeUser( LowerEdge_Ratio, UpperEdge_Ratio );

		// -- flat line = 1.00 -- //
		TF1 *f_line = new TF1("f_line", "1", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("SAME");

		h_ratio->Draw("EPSAME");
	}

	void CanvasWithHistogramsRatioPlot_DYstyle(TH1D *h1, TH1D *h2,
												TString Name1, TString Name2, TString Name_Ratio,
												Int_t color1 = kBlack, Int_t color2 = kRed,
												Bool_t isFill1 = kFALSE, Bool_t isFill2 = kTRUE,
												TString DrawOp1 = "EP", TString DrawOp2 = "HIST SAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		// -- ensure additional space at the bottom side for ratio plot -- //
		TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		h1->Draw(DrawOp1);
		h2->Draw(DrawOp2);
		h1->Draw("AXISSAME");
		h1->Draw(DrawOp1+"SAME");

		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(color1);
		h1->SetLineWidth(1);
		h1->SetMarkerColor(color1);
		h1->SetMarkerSize(1);
		h1->SetMarkerStyle(20);
		if( isFill1 == kTRUE ) h1->SetFillColor(color1);
		else h1->SetFillColorAlpha(kWhite, 0);

		// -- Setting for h2 -- //
		h2->SetLineColor(kOrange+3);
		h2->SetLineWidth(1);
		h2->SetMarkerColor(color2);
		h2->SetMarkerSize(1);
		h2->SetMarkerStyle(20);
		if( isFill2 == kTRUE ) h2->SetFillColor(kOrange-2);
		else h2->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetLabelSize(0);
		h1->GetXaxis()->SetTitleSize(0);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(h1, Name1);
		legend->AddEntry(h2, Name2, "F");
		legend->Draw();

		// -- Draw Grid -- //
		this->DrawGrid(h2);


		// -- Bottom Pad -- //
		BottomPad = new TPad("BottomPad","BottomPad",0.01,0.01,0.99,0.3);
		BottomPad->Draw();
		BottomPad->cd();

		BottomPad->SetBottomMargin(0.4);
		BottomPad->SetRightMargin(0.04);
		BottomPad->SetLeftMargin(0.15);	

		if( isLogX == kTRUE ) gPad->SetLogx();

		// -- Make Ratio plot & Draw it -- //
		h_ratio = (TH1D*)h1->Clone();
		h_ratio->SetName("h_ratio");
		h1->Sumw2(); h2->Sumw2();
		h_ratio->Divide(h1, h2);
		h_ratio->Draw("EP");

		// -- General Setting -- //
		h_ratio->SetLineColor(kBlack);
		h_ratio->SetMarkerStyle(20);
		h_ratio->SetMarkerSize(1);
		h_ratio->SetMarkerColor(kBlack);
		h_ratio->SetStats(kFALSE);

		// -- X-axis Setting -- //		
		h_ratio->GetXaxis()->SetTitle( TitleX );
		h_ratio->GetXaxis()->SetTitleOffset( 0.9 );
		h_ratio->GetXaxis()->SetTitleSize( 0.2 );
		h_ratio->GetXaxis()->SetLabelColor(1);
		h_ratio->GetXaxis()->SetLabelFont(42);
		h_ratio->GetXaxis()->SetLabelOffset(0.007);
		h_ratio->GetXaxis()->SetLabelSize(0.15);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h_ratio->GetXaxis()->SetMoreLogLabels(); h_ratio->GetXaxis()->SetNoExponent(); }

		// -- Y-axis Setting -- //
		h_ratio->GetYaxis()->SetTitle( Name_Ratio );
		h_ratio->GetYaxis()->SetTitleOffset( 0.4 );
		h_ratio->GetYaxis()->SetTitleSize( 0.1);
		h_ratio->GetYaxis()->SetLabelSize( 0.07 );
		h_ratio->GetYaxis()->SetRangeUser( LowerEdge_Ratio, UpperEdge_Ratio );

		// -- flat line = 1.00 -- //
		TF1 *f_line = new TF1("f_line", "1", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("SAME");

		h_ratio->Draw("EPSAME");
	}

	void CanvasWithTHStackRatioPlot(TH1D *h1, THStack *h_stack, TLegend* legend, TString Name_Ratio, Int_t isGrid = 0,
									Int_t color1 = kBlack, Bool_t isFill1 = kFALSE,
									TString DrawOp1 = "EP", TString DrawOp2 = "HIST SAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		//gPad->SetGridx(0);
		//gPad->SetGridy(0);

		// -- ensure additional space at the bottom side for ratio plot -- //
		TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		// -- Set THStack and h2 -- //
		TH1D *h2 = (TH1D*)h1->Clone();
		h2->Add(h1,-1.0);
		TList *list = (TList*)h_stack->GetHists();
		h2->Merge(list);
		h_stack->SetHistogram(h2);

		h1->Draw(DrawOp1);
		h_stack->Draw(DrawOp2);
		h1->Draw("AXISSAME");
		h1->Draw(DrawOp1+"SAME");

		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(color1);
		h1->SetLineWidth(1);
		h1->SetMarkerColor(color1);
		h1->SetMarkerSize(1);
		h1->SetMarkerStyle(20);
		if( isFill1 == kTRUE ) h1->SetFillColor(color1);
		else h1->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetLabelSize(0);
		h1->GetXaxis()->SetTitleSize(0);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend->Draw();

		// -- Draw Grid -- //
		if( isGrid == 1 ) this->DrawGrid(h2);

		// -- Draw latex -- //
		this->DrawLatex();


		// -- Bottom Pad -- //
		BottomPad = new TPad("BottomPad","BottomPad",0.01,0.01,0.99,0.3);
		BottomPad->Draw();
		BottomPad->cd();

		BottomPad->SetBottomMargin(0.4);
		BottomPad->SetRightMargin(0.04);
		BottomPad->SetLeftMargin(0.153);	

		if( isLogX == kTRUE ) gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();

		// -- Make Ratio plot & Draw it -- //
		h_ratio = (TH1D*)h1->Clone();
		h_ratio->SetName("h_ratio");
		h1->Sumw2(); h2->Sumw2();
		h_ratio->Divide(h1, h2);
		h_ratio->Draw("EP");

		// -- General Setting -- //
		h_ratio->SetLineColor(kBlack);
		h_ratio->SetMarkerStyle(20);
		h_ratio->SetMarkerSize(1);
		h_ratio->SetMarkerColor(kBlack);
		h_ratio->SetStats(kFALSE);

		// -- X-axis Setting -- //		
		h_ratio->GetXaxis()->SetTitle( TitleX );
		h_ratio->GetXaxis()->SetTitleOffset( 0.9 );
		h_ratio->GetXaxis()->SetTitleSize( 0.2 );
		h_ratio->GetXaxis()->SetLabelColor(1);
		h_ratio->GetXaxis()->SetLabelFont(42);
		h_ratio->GetXaxis()->SetLabelOffset(0.007);
		h_ratio->GetXaxis()->SetLabelSize(0.15);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h_ratio->GetXaxis()->SetMoreLogLabels(); h_ratio->GetXaxis()->SetNoExponent(); }

		// -- Y-axis Setting -- //
		h_ratio->GetYaxis()->SetTitle( Name_Ratio );
		h_ratio->GetYaxis()->SetTitleOffset( 0.4 );
		h_ratio->GetYaxis()->SetTitleSize( 0.1);
		h_ratio->GetYaxis()->SetLabelSize( 0.07 );
		h_ratio->GetYaxis()->SetRangeUser( LowerEdge_Ratio, UpperEdge_Ratio );

		// -- flat line = 1.00 -- //
		TF1 *f_line = new TF1("f_line", "1", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("SAME");

		h_ratio->Draw("EPSAME");
	}

	void CanvasWithTHStackRatioPlot_jets(TH1D *h1, THStack *h_stack, vector< TH1D* > h_j,
									TLegend* legend, TString Name_Ratio, Int_t isGrid = 0, Bool_t isLogY_Ratio = kFALSE,
									Int_t color1 = kBlack, Bool_t isFill1 = kFALSE,
									TString DrawOp1 = "EP", TString DrawOp2 = "HIST SAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		gPad->SetGridx(0);
		gPad->SetGridy(0);

		// -- ensure additional space at the bottom side for ratio plot -- //
		TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		// -- Set THStack and h2 -- //
		TH1D *h2 = (TH1D*)h1->Clone();
		h2->Add(h1,-1.0);
		TList *list = (TList*)h_stack->GetHists();
		h2->Merge(list);
		h_stack->SetHistogram(h2);

		h1->Draw(DrawOp1);
		h_stack->Draw(DrawOp2);
		h1->Draw("AXISSAME");
		h1->Draw(DrawOp1+"SAME");

		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(color1);
		h1->SetLineWidth(1);
		h1->SetMarkerColor(color1);
		h1->SetMarkerSize(1);
		h1->SetMarkerStyle(20);
		if( isFill1 == kTRUE ) h1->SetFillColor(color1);
		else h1->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetLabelSize(0);
		h1->GetXaxis()->SetTitleSize(0);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend->Draw();

		// -- Draw Grid -- //
		if( isGrid == 1 ) this->DrawGrid(h2);

		// -- Draw latex -- //
		//this->DrawLatex();
		this->DrawLatex_simulation();


		// -- Bottom Pad -- //
		BottomPad = new TPad("BottomPad","BottomPad",0.01,0.01,0.99,0.3);
		BottomPad->Draw();
		BottomPad->cd();

		BottomPad->SetBottomMargin(0.4);
		BottomPad->SetRightMargin(0.04);
		BottomPad->SetLeftMargin(0.15);	

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY_Ratio == kTRUE && LowerEdge_Ratio > 0 ) gPad->SetLogy();

		// -- Make Ratio plot & Draw it -- //
		vector<TH1D*> h_r;
		Int_t nHist = (Int_t)h_j.size();
		for(Int_t i_hist=0; i_hist<nHist; i_hist++)
		{
			TH1D *h = h_j[i_hist];
			//h->Sumw2(); h1->Sumw2();

			h_ratio = (TH1D*)h->Clone();
			h_ratio->Divide(h, h1);

			if( i_hist == 0 ) h_ratio->Draw("EP");
			else if( i_hist != 0 ) h_ratio->Draw("EPSAME");

			// -- General Setting -- //
			h_ratio->SetMarkerStyle(20);
			h_ratio->SetMarkerSize(1);
			h_ratio->SetStats(kFALSE);
			h_ratio->SetFillColorAlpha(kWhite, 0);

			h_r.push_back(h_ratio);
		}

		// -- X-axis Setting -- //		
		h_r[0]->GetXaxis()->SetTitle( TitleX );
		h_r[0]->GetXaxis()->SetTitleOffset( 0.9 );
		h_r[0]->GetXaxis()->SetTitleSize( 0.2 );
		h_r[0]->GetXaxis()->SetLabelColor(1);
		h_r[0]->GetXaxis()->SetLabelFont(42);
		h_r[0]->GetXaxis()->SetLabelOffset(0.007);
		h_r[0]->GetXaxis()->SetLabelSize(0.15);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h_r[0]->GetXaxis()->SetMoreLogLabels(); h_r[0]->GetXaxis()->SetNoExponent(); }
		h_r[0]->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h_r[0]->GetYaxis()->SetTitle( Name_Ratio );
		h_r[0]->GetYaxis()->SetTitleOffset( 0.4 );
		h_r[0]->GetYaxis()->SetTitleSize( 0.1);
		h_r[0]->GetYaxis()->SetLabelSize( 0.07 );
		if( isLogY_Ratio == kTRUE && LowerEdge_Ratio > 0 ) { h_r[0]->GetYaxis()->SetMoreLogLabels(); h_r[0]->GetYaxis()->SetNoExponent(); }
		h_r[0]->GetYaxis()->SetRangeUser( LowerEdge_Ratio, UpperEdge_Ratio );

		// -- flat line = 1.00 -- //
		TF1 *f_line = new TF1("f_line", "1", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("SAME");

		for(Int_t i_hist=0; i_hist<nHist; i_hist++)
			h_r[i_hist]->Draw("PSAME LHIST");
	}

	void CanvasTHStackPlot(TH1D *h1, THStack *h_stack, TLegend* legend, Int_t isGrid = 0, TString DrawOp1 = "EP", TString DrawOp2 = "HIST SAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		gPad->SetGridx(0);
		gPad->SetGridy(0);

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		// -- ensure additional space at the bottom side for ratio plot -- //
		//TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		// -- Set THStack and h2 -- //
		TH1D *h2 = (TH1D*)h1->Clone();
		h2->Add(h1,-1.0);
		TList *list = (TList*)h_stack->GetHists();
		h2->Merge(list);
		h_stack->SetHistogram(h2);

		h1->Draw(DrawOp1);
		h_stack->Draw(DrawOp2);
		h1->Draw("AXISSAME");

		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(0);
		h1->SetLineWidth(0);
		h1->SetMarkerColor(0);
		h1->SetMarkerSize(0);
		h1->SetMarkerStyle(20);
		h1->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetTitleOffset( 0.9 );
		h1->GetXaxis()->SetTitleSize( 0.05 );
		h1->GetXaxis()->SetLabelColor(1);
		h1->GetXaxis()->SetLabelFont(42);
		h1->GetXaxis()->SetLabelOffset(0.007);
		h1->GetXaxis()->SetLabelSize(0.035);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend->Draw();

		// -- Draw Grid -- //
		if( isGrid == 1 ) this->DrawGrid(h2);

		// -- Draw latex -- //
		this->DrawLatex();
	}

	void CanvasWithTHStackPlot(TH1D *h1, THStack *h_stack, TLegend* legend, Int_t isGrid = 0,
									Int_t color1 = kBlack, Bool_t isFill1 = kFALSE,
									TString DrawOp1 = "EP", TString DrawOp2 = "HIST SAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		gPad->SetGridx(0);
		gPad->SetGridy(0);

		// -- ensure additional space at the bottom side for ratio plot -- //
		//TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		// -- Set THStack and h2 -- //
		TH1D *h2 = (TH1D*)h1->Clone();
		h2->Add(h1,-1.0);
		TList *list = (TList*)h_stack->GetHists();
		h2->Merge(list);
		h_stack->SetHistogram(h2);

		h1->Draw(DrawOp1);
		h_stack->Draw(DrawOp2);
		h1->Draw("AXISSAME");
		h1->Draw(DrawOp1+"SAME");


		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(color1);
		h1->SetLineWidth(1);
		h1->SetMarkerColor(color1);
		h1->SetMarkerSize(1);
		h1->SetMarkerStyle(20);
		if( isFill1 == kTRUE ) h1->SetFillColor(color1);
		else h1->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetTitleOffset( 0.9 );
		h1->GetXaxis()->SetTitleSize( 0.05 );
		h1->GetXaxis()->SetLabelColor(1);
		h1->GetXaxis()->SetLabelFont(42);
		h1->GetXaxis()->SetLabelOffset(0.007);
		h1->GetXaxis()->SetLabelSize(0.035);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend->Draw();

		// -- Draw Grid -- //
		if( isGrid == 1 ) this->DrawGrid(h2);

		// -- Draw latex -- //
		this->DrawLatex();
	}


	void CanvasWithThreeHistogramsRatioPlot(TH1D *h1, TH1D *h2, TH1D* h_ref, 
											TString Name1, TString Name2, TString Name_ref, TString Name_Ratio,
											Int_t color1 = kBlue, Int_t color2 = kGreen+1, Int_t color_ref = kRed,
											Bool_t isFill1 = kFALSE, Bool_t isFill2 = kFALSE, Bool_t isFill_ref = kFALSE,
											TString DrawOp1 = "EP", TString DrawOp2 = "EPSAME", TString DrawOp_ref = "EPSAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		gPad->SetGridx(0);
		gPad->SetGridy(0);

		// -- ensure additional space at the bottom side for ratio plot -- //
		TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		h1->Draw(DrawOp1);
		h2->Draw(DrawOp2);
		h_ref->Draw(DrawOp_ref);

		// -- General Setting for both plots -- //
		h1->SetXTitle( TitleX );
		h1->SetYTitle( TitleY );
		h1->SetStats(kFALSE);
		h2->SetStats(kFALSE);
		h_ref->SetStats(kFALSE);

		// -- Setting for h1 -- //
		h1->SetLineColor(color1);
		h1->SetLineWidth(1);
		h1->SetMarkerColor(color1);
		h1->SetMarkerSize(1);
		h1->SetMarkerStyle(20);
		if( isFill1 == kTRUE ) h1->SetFillColor(color1);
		else h1->SetFillColorAlpha(kWhite, 0);

		// -- Setting for h2 -- //
		h2->SetLineColor(color2);
		h2->SetLineWidth(1);
		h2->SetMarkerColor(color2);
		h2->SetMarkerSize(1);
		h2->SetMarkerStyle(20);
		if( isFill2 == kTRUE ) h2->SetFillColor(color2);
		else h2->SetFillColorAlpha(kWhite, 0);

		// -- Setting for h_ref -- //
		h_ref->SetLineColor(color_ref);
		h_ref->SetLineWidth(1);
		h_ref->SetMarkerColor(color_ref);
		h_ref->SetMarkerSize(1);
		h_ref->SetMarkerStyle(20);
		if( isFill_ref == kTRUE ) h_ref->SetFillColor(color_ref);
		else h_ref->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		h1->GetXaxis()->SetLabelSize(0);
		h1->GetXaxis()->SetTitleSize(0);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h1->GetXaxis()->SetNoExponent(); h1->GetXaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h1->GetYaxis()->SetTitleSize(0.06);
		h1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h1->GetYaxis()->SetNoExponent(); h1->GetYaxis()->SetMoreLogLabels(); }
		if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(h_ref, Name_ref);
		legend->AddEntry(h1, Name1);
		legend->AddEntry(h2, Name2);
		legend->Draw();

		// -- Bottom Pad -- //
		BottomPad = new TPad("BottomPad","BottomPad",0.01,0.01,0.99,0.3);
		BottomPad->Draw();
		BottomPad->cd();

		BottomPad->SetBottomMargin(0.4);
		BottomPad->SetRightMargin(0.04);
		BottomPad->SetLeftMargin(0.15);	

		if( isLogX == kTRUE ) gPad->SetLogx();

		// -- Make Ratio plot & Draw it -- //
		//h1->Sumw2(); h2->Sumw2(); h_ref->Sumw2();
		h_ratio1 = (TH1D*)h1->Clone();
		h_ratio1->SetName("h_ratio1");
		h_ratio1->Divide(h1, h_ref);

		h_ratio2 = (TH1D*)h1->Clone();
		h_ratio2->SetName("h_ratio2");
		h_ratio2->Divide(h2, h_ref);

		h_ratio1->Draw("EP");
		//h_ratio1->Draw("EPSAME");

		// -- General Setting for 1st ratio plot -- //
		h_ratio1->SetLineColor(color1);
		h_ratio1->SetMarkerStyle(20);
		h_ratio1->SetMarkerSize(1);
		h_ratio1->SetMarkerColor(color1);
		h_ratio1->SetStats(kFALSE);

		// -- General Setting for 2nd ratio plot -- //
		h_ratio2->SetLineColor(color2);
		h_ratio2->SetMarkerStyle(20);
		h_ratio2->SetMarkerSize(1);
		h_ratio2->SetMarkerColor(color2);
		h_ratio2->SetStats(kFALSE);

		// -- X-axis Setting -- //		
		h_ratio1->GetXaxis()->SetTitle( TitleX );
		h_ratio1->GetXaxis()->SetTitleOffset( 0.9 );
		h_ratio1->GetXaxis()->SetTitleSize( 0.2 );
		h_ratio1->GetXaxis()->SetLabelColor(1);
		h_ratio1->GetXaxis()->SetLabelFont(42);
		h_ratio1->GetXaxis()->SetLabelOffset(0.007);
		h_ratio1->GetXaxis()->SetLabelSize(0.15);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h_ratio1->GetXaxis()->SetMoreLogLabels(); h_ratio1->GetXaxis()->SetNoExponent(); }

		// -- Y-axis Setting -- //
		h_ratio1->GetYaxis()->SetTitle( Name_Ratio );
		h_ratio1->GetYaxis()->SetTitleOffset( 0.4 );
		h_ratio1->GetYaxis()->SetTitleSize( 0.1);
		h_ratio1->GetYaxis()->SetLabelSize( 0.07 );
		h_ratio1->GetYaxis()->SetRangeUser( LowerEdge_Ratio, UpperEdge_Ratio );

		// -- flat line = 1.00 -- //
		TF1 *f_line = new TF1("f_line", "1", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("SAME");

		h_ratio1->Draw("EPSAME");
		h_ratio2->Draw("EPSAME");
	}

	void CanvasWithMultipleHistograms(vector< TH1D* > Histos, vector< TString > Names, TString DrawOp = "EP")
	{
		Int_t nHist = (Int_t)Histos.size();
		Int_t nColors = (Int_t)Colors.size();
		if( nHist > nColors )
		{
			cout << "# histograms > # colors! ... add more colors: MyCanvas->Colors.push_back( 'Color' )" << endl;
			return;
		}

		c->cd();

		gPad->SetGridx(0);
		gPad->SetGridy(0);

		gPad->SetRightMargin(0.05);
		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);

		for(Int_t i_hist=0; i_hist<nHist; i_hist++)
		{
			TH1D *h = Histos[i_hist];
			Int_t color = Colors[i_hist];

			h->Draw(DrawOp+"SAME");

			// -- Axis Setting: It is enough to do on the first histogram -- //
			if( i_hist == 0 )
			{
				h->SetXTitle( TitleX );
				h->SetYTitle( TitleY );

				// -- X-axis Setting -- //
				h->GetXaxis()->SetTitleSize(0.06);
				h->GetXaxis()->SetTitleOffset(0.9);
				h->GetXaxis()->SetLabelSize(0.04);
				if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h->GetXaxis()->SetNoExponent(); h->GetXaxis()->SetMoreLogLabels(); }
				if( !(LowerEdge_X == 0 && UpperEdge_X == 0) ) h->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

				// -- Y-axis Setting -- //
				h->GetYaxis()->SetTitleSize(0.06);
				h->GetYaxis()->SetTitleOffset(1.25);
				h->GetYaxis()->SetLabelSize(0.04);
				if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h->GetYaxis()->SetNoExponent(); h->GetYaxis()->SetMoreLogLabels(); }
				if( !(LowerEdge_Y == 0 && UpperEdge_Y == 0) ) h->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );
			}

			h->SetStats(kFALSE);
			h->SetLineColor(color);
			h->SetLineWidth(1);
			h->SetMarkerColor(color);
			h->SetMarkerSize(1);
			h->SetMarkerStyle(20);
			h->SetFillColorAlpha(kWhite, 0);

			legend->AddEntry( h, Names[i_hist] );

		}
		legend->Draw();


	}

	void CanvasWithMultipleGraphs(vector< TGraphAsymmErrors* > Graphs, vector< TString > Names, TString DrawOp = "EP")
	{
		Int_t nGraphs = (Int_t)Graphs.size();
		Int_t nColors = (Int_t)Colors.size();
		if( nGraphs > nColors )
		{
			cout << "# Graphs > # colors! ... add more colors: MyCanvas->Colors.push_back( 'Color' )" << endl;
			return;
		}

		c->cd();

		gPad->SetRightMargin(0.05);
		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);

		for(Int_t i_graph=0; i_graph<nGraphs; i_graph++)
		{
			TGraphAsymmErrors *g = Graphs[i_graph];
			Int_t color = Colors[i_graph];

			if( i_graph == 0 )
				g->Draw(DrawOp+"A");
			else
				g->Draw(DrawOp+"SAME");

			// -- Axis Setting: It is enough to do on the first graph -- //
			if( i_graph == 0 )
			{
				g->GetXaxis()->SetTitle( TitleX );
				g->GetYaxis()->SetTitle( TitleY );

				// -- X-axis Setting -- //
				if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { g->GetXaxis()->SetNoExponent(); g->GetXaxis()->SetMoreLogLabels(); }
				if( LowerEdge_X != 0 && UpperEdge_X != 0 ) g->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

				// -- Y-axis Setting -- //
				g->GetYaxis()->SetTitleSize(0.06);
				g->GetYaxis()->SetTitleOffset(1.25);
				if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { g->GetYaxis()->SetNoExponent(); g->GetYaxis()->SetMoreLogLabels(); }
				if( LowerEdge_Y != 0 && UpperEdge_Y != 0 ) g->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );
			}

			g->SetLineColor(color);
			g->SetLineWidth(1);
			g->SetMarkerColor(color);
			g->SetMarkerSize(1);
			g->SetMarkerStyle(20);
			g->SetFillColorAlpha(kWhite, 0);

			legend->AddEntry( g, Names[i_graph] );

		}
		legend->Draw();

		
	}

	void CanvasWithGraphRatioPlot(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2,
									TString Name1, TString Name2, TString Name_Ratio,
									Int_t color1 = kBlack, Int_t color2 = kRed,
									TString DrawOp1 = "EP", TString DrawOp2 = "EPSAME")
	{
		c->cd();
		// -- Top Pad -- //
		TopPad = new TPad("TopPad", "TopPad", 0.01, 0.01, 0.99, 0.99);
		TopPad->Draw();
		TopPad->cd();

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();


		// -- ensure additional space at the bottom side for ratio plot -- //
		TopPad->SetBottomMargin(0.32);
		TopPad->SetRightMargin(0.05);

		g1->Draw(DrawOp1+"A");
		g2->Draw(DrawOp2);

		// -- Setting for g1 -- //
		g1->SetLineColor(color1);
		g1->SetLineWidth(1);
		g1->SetMarkerColor(color1);
		g1->SetMarkerSize(1);
		g1->SetMarkerStyle(20);
		g1->SetFillColorAlpha(kWhite, 0);

		// -- Setting for g2 -- //
		g2->SetLineColor(color2);
		g2->SetLineWidth(1);
		g2->SetMarkerColor(color2);
		g2->SetMarkerSize(1);
		g2->SetMarkerStyle(20);
		g2->SetFillColorAlpha(kWhite, 0);

		// -- X-axis Setting -- //
		g1->GetXaxis()->SetTitle( TitleX );
		g1->GetXaxis()->SetLabelSize(0);
		g1->GetXaxis()->SetTitleSize(0);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { g1->GetXaxis()->SetNoExponent(); g1->GetXaxis()->SetMoreLogLabels(); }
		// if( LowerEdge_X != 0 && UpperEdge_X != 0 ) g1->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );
		if( LowerEdge_X != 0 && UpperEdge_X != 0 ) g1->GetXaxis()->SetLimits( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		g1->GetYaxis()->SetTitle( TitleY );
		g1->GetYaxis()->SetTitleSize(0.06);
		g1->GetYaxis()->SetTitleOffset(1.25);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { g1->GetYaxis()->SetNoExponent(); g1->GetYaxis()->SetMoreLogLabels(); }
		if( LowerEdge_Y != 0 && UpperEdge_Y != 0 ) g1->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		// -- Add Legend -- //
		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(g1, Name1);
		legend->AddEntry(g2, Name2);
		legend->Draw();


		// -- Bottom Pad -- //
		BottomPad = new TPad("BottomPad","BottomPad",0.01,0.01,0.99,0.3);
		BottomPad->Draw();
		BottomPad->cd();

		BottomPad->SetBottomMargin(0.4);
		BottomPad->SetRightMargin(0.04);
		BottomPad->SetLeftMargin(0.15);	

		if( isLogX == kTRUE ) gPad->SetLogx();

		// -- Make Ratio plot & Draw it -- //
		g_ratio = (TGraphAsymmErrors*)g1->Clone();
		MakeRatioGraph(g_ratio, g1, g2);
		g_ratio->Draw("AEP");

		// -- General Setting -- //
		g_ratio->SetLineColor(kBlack);
		g_ratio->SetMarkerStyle(20);
		g_ratio->SetMarkerSize(1);
		g_ratio->SetMarkerColor(kBlack);

		// -- X-axis Setting -- //		
		g_ratio->GetXaxis()->SetTitle( TitleX );
		g_ratio->GetXaxis()->SetTitleOffset( 0.9 );
		g_ratio->GetXaxis()->SetTitleSize( 0.2 );
		g_ratio->GetXaxis()->SetLabelColor(1);
		g_ratio->GetXaxis()->SetLabelFont(42);
		g_ratio->GetXaxis()->SetLabelOffset(0.007);
		g_ratio->GetXaxis()->SetLabelSize(0.15);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { g_ratio->GetXaxis()->SetMoreLogLabels(); g_ratio->GetXaxis()->SetNoExponent(); }
		// if( LowerEdge_X != 0 && UpperEdge_X != 0 ) g_ratio->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );
		if( LowerEdge_X != 0 && UpperEdge_X != 0 ) g_ratio->GetXaxis()->SetLimits( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		g_ratio->GetYaxis()->SetTitle( Name_Ratio );
		g_ratio->GetYaxis()->SetTitleOffset( 0.4 );
		g_ratio->GetYaxis()->SetTitleSize( 0.1);
		g_ratio->GetYaxis()->SetLabelSize( 0.07 );
		g_ratio->GetYaxis()->SetRangeUser( LowerEdge_Ratio, UpperEdge_Ratio );

		// -- flat line = 1.00 -- //
		TF1 *f_line = new TF1("f_line", "1", -10000, 10000);
		f_line->SetLineColor(kAzure+1);
		f_line->SetLineWidth(1);
		f_line->Draw("SAME");

		g_ratio->Draw("PSAME");		
	}

	void MakeRatioGraph(TGraphAsymmErrors *g_ratio, TGraphAsymmErrors *g1, TGraphAsymmErrors *g2)
	{
		g_ratio->Set(0); // -- Remove all points (reset) -- //

		Int_t NPoints = g1->GetN();
		for(Int_t i_p=0; i_p<NPoints; i_p++)
		{
			// -- Get g1 point -- //
			Double_t x1, y1;
			g1->GetPoint(i_p, x1, y1);
			Double_t error1 = ReturnLargerValue( g1->GetErrorYhigh(i_p), g1->GetErrorYlow(i_p) );

			// -- Get g2 point -- //
			Double_t x2, y2;
			g2->GetPoint(i_p, x2, y2);
			Double_t error2 = ReturnLargerValue( g2->GetErrorYhigh(i_p), g2->GetErrorYlow(i_p) );

			Double_t ratio;
			Double_t ratio_error;
			if(y1 != 0 && error1 != 0 && y2 != 0 && error2 != 0)
			{
				// -- calculate ratio & error -- //
				ratio = y1 / y2;
				ratio_error = Error_PropagatedAoverB(y1, error1, y2, error2);
			}
			else if( y1 != 0 && y2 != 0 && (error1 == 0 || error2 == 0) )
			{
				ratio = y1 / y2;
				ratio_error = 0;
			}
			else
			{
				ratio = 0;
				ratio_error = 0;
			}

			// -- Set Central value -- //
			g_ratio->SetPoint(i_p, x1, ratio);

			// -- Set the error -- //
			Double_t error_XLow = g1->GetErrorXlow(i_p);
			Double_t error_Xhigh = g1->GetErrorXhigh(i_p);
			g_ratio->SetPointError(i_p, error_XLow, error_Xhigh, ratio_error, ratio_error);

		}
	}

	void CanvasWitOneHistogram(TH1D *h, TString Name, TString DrawOp = "EP", Int_t color = kRed)
	{
		c->cd();

		gPad->SetRightMargin(0.05);
		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		//gPad->SetGridx(0);
		//gPad->SetGridy(0);
		gPad->SetGridx(1);
		gPad->SetGridy(1);

		h->Draw(DrawOp+"SAME");

		h->SetXTitle( TitleX );
		h->SetYTitle( TitleY );

		// -- X-axis Setting -- //
		h->GetXaxis()->SetTitleSize(0.06);
		h->GetXaxis()->SetTitleOffset(0.9);
		//h->GetXaxis()->SetLabelSize(0.06);
		h->GetXaxis()->SetLabelSize(0.04);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h->GetXaxis()->SetNoExponent(); h->GetXaxis()->SetMoreLogLabels(); }
		if( LowerEdge_X != 0 && UpperEdge_X != 0 ) h->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h->GetYaxis()->SetTitleSize(0.06);
		h->GetYaxis()->SetTitleOffset(1.25);
		//h->GetYaxis()->SetLabelSize(0.06);
		h->GetYaxis()->SetLabelSize(0.04);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h->GetYaxis()->SetNoExponent(); h->GetYaxis()->SetMoreLogLabels(); }
		if( LowerEdge_Y != 0 && UpperEdge_Y != 0 ) h->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		h->SetStats(kFALSE);
		h->SetLineColor(color);
		h->SetLineWidth(1);
		h->SetMarkerColor(color);
		h->SetMarkerSize(1);
		h->SetMarkerStyle(20);
		//h->SetFillColorAlpha(kWhite, 0);
		h->SetFillColor(color);

		// -- Add Legend -- //
		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(h, Name);
		legend->Draw();

		// -- Draw latex -- //
		this->DrawLatex();
	}

	void CanvasWitOneHistogramAndFunc(TH1D *h, TF1 *f, TString Name, TString DrawOp = "EP", Int_t color = kRed)
	{
		c->cd();

		gPad->SetRightMargin(0.05);
		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();

		h->Draw(DrawOp+"SAME");

		h->SetXTitle( TitleX );
		h->SetYTitle( TitleY );

		// -- X-axis Setting -- //
		h->GetXaxis()->SetTitleSize(0.06);
		h->GetXaxis()->SetTitleOffset(0.9);
		//h->GetXaxis()->SetLabelSize(0.06);
		h->GetXaxis()->SetLabelSize(0.04);
		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h->GetXaxis()->SetNoExponent(); h->GetXaxis()->SetMoreLogLabels(); }
		if( LowerEdge_X != 0 && UpperEdge_X != 0 ) h->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		// -- Y-axis Setting -- //
		h->GetYaxis()->SetTitleSize(0.06);
		h->GetYaxis()->SetTitleOffset(1.25);
		//h->GetYaxis()->SetLabelSize(0.06);
		h->GetYaxis()->SetLabelSize(0.04);
		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h->GetYaxis()->SetNoExponent(); h->GetYaxis()->SetMoreLogLabels(); }
		if( LowerEdge_Y != 0 && UpperEdge_Y != 0 ) h->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );

		//h->SetStats(kFALSE);
		h->SetLineColor(color);
		h->SetLineWidth(1);
		h->SetMarkerColor(color);
		h->SetMarkerSize(1);
		h->SetMarkerStyle(20);
		h->SetFillColorAlpha(kWhite, 0);

		// -- Add Legend -- //
		legend = new TLegend(Legend_x1, Legend_y1, Legend_x2, Legend_y2);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(h, Name);
		legend->Draw();

		// -- Draw function -- //
		gStyle->SetOptFit();
		f->SetLineWidth(1);
		f->SetLineColor(kBlue);
		f->Draw("SAME");
	}


	void CanvasWith2DHistogram( TH2D *_h_2D, TString DrawOp="COLZTEXT")
	{
		// gStyle->Reset();
		// gStyle->SetPalette(1);
		TH2D* h_2D = (TH2D*)_h_2D->Clone();


		if( DoTranspose )
			h_2D = this->Transpose(h_2D);

		
		c->cd();
		gPad->SetLeftMargin(0.14);
		gPad->SetRightMargin(0.12);
		gPad->SetTopMargin(0.03);
		gPad->SetBottomMargin(0.10);
		// gPad->SetLeftMargin(0.12);
		// gPad->SetRightMargin(0.15);

		if( isLogX == kTRUE ) gPad->SetLogx();
		if( isLogY == kTRUE ) gPad->SetLogy();
		if( isLogZ == kTRUE ) gPad->SetLogz();


		h_2D->Draw(DrawOp);


		h_2D->SetStats(kFALSE);
		h_2D->SetXTitle(TitleX);
		h_2D->SetYTitle(TitleY);


		h_2D->GetYaxis()->SetNoExponent();
		h_2D->GetYaxis()->SetMoreLogLabels();
		h_2D->GetYaxis()->SetTitleOffset(2);
		h_2D->GetZaxis()->SetRangeUser(LowerEdge_Z, UpperEdge_Z);
		// h_2D->SetMinimum(LowerEdge_Z);
		// h_2D->SetMaximum(UpperEdge_Z);

		// TPaletteAxis *palette = (TPaletteAxis*)h_2D->GetListOfFunctions()->FindObject("palette");
		// if( palette != NULL )
		// {
		// 	// the following lines moe the paletter. Choose the values you need for the position.
		// 	palette->SetX1NDC(0.885);

		// 	palette->SetX2NDC(0.87);
		// 	palette->SetY1NDC(0.10);
		// 	palette->SetY2NDC(0.97);
		// 	gPad->Modified();
		// 	gPad->Update();
		// }

		if( isSetNoExpo_MoreLogLabels_X == kTRUE ) { h_2D->GetXaxis()->SetNoExponent(); h_2D->GetXaxis()->SetMoreLogLabels(); }
		if( LowerEdge_X != 0 || UpperEdge_X != 0 ) h_2D->GetXaxis()->SetRangeUser( LowerEdge_X, UpperEdge_X );

		if( isSetNoExpo_MoreLogLabels_Y == kTRUE ) { h_2D->GetYaxis()->SetNoExponent(); h_2D->GetYaxis()->SetMoreLogLabels(); }
		if( LowerEdge_Y != 0 || UpperEdge_Y != 0 ) h_2D->GetYaxis()->SetRangeUser( LowerEdge_Y, UpperEdge_Y );
	}

	TH2D* Transpose( TH2D* h_2D )
	{
		TH2D* h_2D_Transposed = (TH2D*)h_2D->Clone();
		h_2D_Transposed->SetDirectory(0);

		Int_t nBin_X = h_2D->GetNbinsX();
		Int_t nBin_Y = h_2D->GetNbinsY();

		if(nBin_X != nBin_Y )
		{
			printf("(nBin_X, nBin_Y) = (%d, %d): this is not squared matrix! ... transpose for this matrix is not provided by this function\n", nBin_X, nBin_Y);
			return NULL;
		}

		for(Int_t iter_x=0; iter_x<nBin_X; iter_x++)
		{
			Int_t i_binX = iter_x+1;
			for(Int_t iter_y=0; iter_y<nBin_Y; iter_y++)
			{
				Int_t i_binY = iter_y+1;

				Double_t value_before = h_2D->GetBinContent(i_binX, i_binY);

				h_2D_Transposed->SetBinContent(i_binY, i_binX, value_before);
			}
		}

		return h_2D_Transposed;
	}

	//void PrintCanvas()
	//{
	//	TString CanvasName = c->GetName();
	//	c->SaveAs(CanvasName+".pdf");
	//}
	void PrintCanvas( TString name = "" )
	{
		TString CanvasName;
		if( name == "" )
			CanvasName = c->GetName();
		else
			CanvasName = name;
		c->SaveAs(CanvasName+".pdf");
	}

	void PrintCanvas_ROOT()
	{
		TString CanvasName = c->GetName();
		c->SaveAs(CanvasName+".root");
	}

	void PrintCanvas_C()
	{
		TString CanvasName = c->GetName();
		c->SaveAs(CanvasName+".C");
	}

	void DrawLatex()
	{
		TLatex *tex = new TLatex(0.185,0.9,"#font[62]{#scale[0.75]{CMS}}");
		tex->SetNDC();
		tex->SetLineWidth(2);
		tex->Draw();

		tex = new TLatex(0.185,0.865,"#font[42]{#it{#scale[0.7]{Preliminary}}}");
		tex->SetNDC();
		tex->SetLineWidth(2);
		tex->Draw();

		tex = new TLatex(0.75,0.96,"#font[42]{#scale[0.6]{35.9 fb^{-1} (13TeV)}}");
		tex->SetNDC();
		tex->SetLineWidth(2);
		tex->Draw();
	}

	void DrawLatex_simulation()
	{
		TLatex *tex = new TLatex(0.185,0.9,"#font[62]{#scale[0.75]{CMS}}");
		tex->SetNDC();
		tex->SetLineWidth(2);
		tex->Draw();

		tex = new TLatex(0.185,0.865,"#font[42]{#it{#scale[0.7]{Simulation}}}");
		tex->SetNDC();
		tex->SetLineWidth(2);
		tex->Draw();

		tex = new TLatex(0.85,0.96,"#font[42]{#scale[0.6]{13TeV}}");
		tex->SetNDC();
		tex->SetLineWidth(2);
		tex->Draw();
	}

	void DrawGrid(TH1D* h)
	{
		// -- Grid -- //
		TLine grid_;
		grid_.SetLineColor(kGray+2);
		grid_.SetLineStyle(kSolid);
		for( size_t ii=0; ii<44; ii++ )
			grid_.DrawLine(bins[ii],0.1,bins[ii],h->GetBinContent(ii+1)); // using Drell-Yan mass bins
	}

protected:
	Double_t Error_PropagatedAoverB(Double_t A, Double_t sigma_A, Double_t B, Double_t sigma_B)
	{
		Double_t ratio_A = (sigma_A) / A;
		Double_t ratio_B = (sigma_B) / B;

		Double_t errorSquare = ratio_A * ratio_A + ratio_B * ratio_B;

		return (A/B) * sqrt(errorSquare);
	}

	Double_t ReturnLargerValue(Double_t a, Double_t b)
	{
		if( a > b )
			return a;
		else
			return b;
	}
};

TH1D* DYMassBin( TH1D* h_1D ) // Please use Rebin() !!
{
	TH1D *h_1D_DYMassBin = new TH1D("h_1D_DYMassBin", "", 43, bins);
	h_1D_DYMassBin->SetDirectory(0);

	Int_t nBin_X = h_1D->GetNbinsX();

	for(Int_t iter_x=1; iter_x<=nBin_X; iter_x++)
	{
		Double_t weight = h_1D->GetBinContent(iter_x);
		Double_t xvalue = h_1D->GetBinCenter(iter_x);
		h_1D_DYMassBin->Fill(xvalue, weight);
	}

	for(Int_t i_x=1; i_x<=43; i_x++)
	{
		Double_t error = sqrt(h_1D_DYMassBin->GetBinContent(i_x));
		h_1D_DYMassBin->SetBinError(i_x, error);
	}

	return h_1D_DYMassBin;
}

TH1D* AbsBin( TH1D* h_1D ) // Please use Rebin() !!
{
	TH1D *h_1D_AbsBin = new TH1D("h_1D_AbsBin", "", 500, 0, 5);
	h_1D_AbsBin->SetDirectory(0);

	Int_t nBin_X = h_1D->GetNbinsX();

	for(Int_t iter_x=1; iter_x<=nBin_X; iter_x++)
	{
		Double_t weight = h_1D->GetBinContent(iter_x);
		Double_t xvalue = h_1D->GetBinCenter(iter_x);
		h_1D_AbsBin->Fill(fabs(xvalue), weight);
	}

	for(Int_t i_x=1; i_x<=500; i_x++)
	{
		Double_t error = sqrt(h_1D_AbsBin->GetBinContent(i_x));
		h_1D_AbsBin->SetBinError(i_x, error);
	}

	return h_1D_AbsBin;
}

void ErrorCalc( TH1D* h1, TH1D* h2, TH1D* h3, Int_t nBins )
{
	for(Int_t i_x=1; i_x<=nBins; i_x++)
	{
		if( h1->GetBinContent(i_x) != 0 && h2->GetBinContent(i_x) != 0 )
		{
			Double_t RelErr1 = (h1->GetBinError(i_x))/(h1->GetBinContent(i_x));
			Double_t RelErr2 = (h2->GetBinError(i_x))/(h2->GetBinContent(i_x));
			Double_t RelErr = sqrt(RelErr1*RelErr1 + RelErr2*RelErr2);
			h3->SetBinError( i_x, RelErr * ( h3->GetBinContent(i_x) ) );
		}
	}
}

TH2D* doTranspose( TH2D* h_2D )
{
	TH2D* h_2D_Transposed = (TH2D*)h_2D->Clone();
	h_2D_Transposed->SetDirectory(0);

	Int_t nBin_X = h_2D->GetNbinsX();
	Int_t nBin_Y = h_2D->GetNbinsY();

	if(nBin_X != nBin_Y )
	{
		printf("(nBin_X, nBin_Y) = (%d, %d): this is not squared matrix! ... transpose for this matrix is not provided by this function\n", nBin_X, nBin_Y);
		return NULL;
	}

	for(Int_t iter_x=0; iter_x<nBin_X; iter_x++)
	{
		Int_t i_binX = iter_x+1;
		for(Int_t iter_y=0; iter_y<nBin_Y; iter_y++)
		{
			Int_t i_binY = iter_y+1;

			Double_t value_before = h_2D->GetBinContent(i_binX, i_binY);

			h_2D_Transposed->SetBinContent(i_binY, i_binX, value_before);
		}
	}

	return h_2D_Transposed;
}
