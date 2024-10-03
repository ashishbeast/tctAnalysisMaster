//This macro calculates gain for multiples files
//Gain is defined as the ratio of charge collected by LGAD to charge collected by PIN sensor

#include "AnalyzeTCTData.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPaveStats.h"
#include "TF1.h"

using namespace std;
void SetStyle(Bool_t threeD = false);
Float_t CalcCFD(AnalyzeTCTData *lgad, TH1F *his, Float_t thr);

int main()
{ 
  //Set Canvas Settings
  SetStyle(false);

  //Define Canvas and Legend
  TCanvas *canvas = new TCanvas("canvas","",1800,1800);
  canvas->Divide(2,2);
  TLegend *leg1 = new TLegend(0.2078464,0.625,0.4482471,0.914959,NULL,"brNDC");
  leg1->SetHeader("Test Timing Setup","C");
  TLegend *leg2 = new TLegend(0.5078464,0.725,0.8482471,0.914959,NULL,"brNDC");
  leg2->SetHeader("Test Timing Setup","C");

  Int_t col[12] = {kRed,kBlue,kGreen+2,kBlack,13,6,7,8,9,28,34,49};
  Int_t ms[12] = {20,21,22,23,29,33,34,39,41,43,45,47};

  //Define the filename
  
  TString file = "../timingLaser/amplitudeTest_XY_IR_c.tct";
  //TString file = "../timingLaser/amplitudeTest_VS_IR_c_Avg.tct";
  //TString file = "../timingLaser/amplitudeTest_VS_IR_c_NoAvg.tct";

  TObjArray *info;
  info = file.Tokenize("/");
  TString fileName = ((TObjString*)(info->At(2)))->String();
  const char* fileOutName = fileName.ReplaceAll(".tct","");
  
  //Read data
  AnalyzeTCTData *lgad = new AnalyzeTCTData(file, 7.2);
  lgad->CorrectBaseline();
  lgad->CalcNoise();
  lgad->CalculateSignalProperties();
  lgad->SaveSignalShape();

  Int_t nEvents = lgad->_events; //Total number of Events
  Int_t nV = lgad->_nV1;

  //cout<<nEvents<< '\t' << nV <<endl;
  
  Float_t amp1, amp2, cfd1, cfd2;
  Float_t t0 = lgad->_tct->t0; //nanoseconds
  Float_t tA = lgad->_tA; //in nanoseconds
  Float_t nPoints = lgad->_tct->nPoints;
  Float_t dt = lgad->_tct->dt; //picoseconds

  TH1F *deltaT = new TH1F("deltaT", "", 100, 26, 27);
  TH1F *deltaA = new TH1F("deltaA", "", 50, 0, 0.2);

  if(nV>1)
    {
      Float_t voltage[nV], ampP1[nV], ampP2[nV], ratio[nV];
      Float_t step = TMath::Abs(lgad->_tct->V1[1]-lgad->_tct->V1[0]);

      TH1F *hSignal, *hSignal1, *hSignal2; 
      for(Int_t i=0; i<nEvents; ++i)
	{
	  hSignal1 = new TH1F("hSignal1", "", nPoints/2, 0, (t0+(nPoints/2)*dt)*1e9);
	  hSignal2 = new TH1F("hSignal2", "", nPoints/2, 0, (t0+(nPoints/2)*dt)*1e9);
	  
	  hSignal = (TH1F*) lgad->_histo[0][i]->Clone();
	  hSignal->Scale(-1);
	  
	  //Split main signal into two individual signals
	  for(Int_t j=0; j<nPoints/2; ++j)
	    {
	      hSignal1->SetBinContent(j+1, hSignal->GetBinContent(j+1));
	      hSignal2->SetBinContent(j+1, hSignal->GetBinContent(nPoints/2+j+1));
	    }
	  
	  amp1 = hSignal1->GetMaximum();
	  amp2 = hSignal2->GetMaximum();
	  
	  cfd1 = CalcCFD(lgad, hSignal1, 0.5);
	  cfd2 = CalcCFD(lgad, hSignal2, 0.5);

	  //Fill arrays for the ratio 
	  voltage[i] = i*step;
	  ampP1[i] = amp1;
	  ampP2[i] = amp2;
	  ratio[i] = abs(amp1-amp2)/amp1;
	  //ratio[i] = amp1/amp2;

	  //cout<<hSignal->GetMaximum()<<'\t'<<amp1<<'\t'<<amp2<<'\t'<<ratio[i]<<endl;

	  deltaT->Fill(cfd1-cfd2);
	}
      TGraphErrors *gr[3];
      gr[0] = new TGraphErrors(nEvents,voltage,ampP1,0,0);
      gr[1] = new TGraphErrors(nEvents,voltage,ampP2,0,0);
      gr[2] = new TGraphErrors(nEvents,voltage,ratio,0,0);
      
      //Save full signal shape in Canvas 1 and splitted signal shape in Canvas 2
      canvas->cd(1);	  
      hSignal->Draw("HIST");
      hSignal->GetXaxis()->SetTitle("Time (nS)");
      hSignal->GetYaxis()->SetTitle("Amplitude (mV)");
      hSignal->SetLineColor(kRed);
      hSignal->SetLineWidth(3);
      hSignal->GetYaxis()->SetRangeUser(-5, 100);
      //hSignal->GetXaxis()->SetRangeUser(0, 5);

      canvas->cd(2);
      hSignal1->Draw("");
      hSignal1->GetXaxis()->SetTitle("Time (nS)");
      hSignal1->GetYaxis()->SetTitle("Amplitude (mV)");
      hSignal1->SetLineColor(kBlue);
      hSignal1->SetLineWidth(3);
      hSignal1->GetYaxis()->SetRangeUser(-5, 100);      

      hSignal2->Draw("lsame");
      hSignal2->GetXaxis()->SetTitle("Time (nS)");
      hSignal2->GetYaxis()->SetTitle("Amplitude (mV)");
      hSignal2->SetLineColor(kRed);
      hSignal2->SetLineWidth(3);
      hSignal2->GetYaxis()->SetRangeUser(-5, 100);

      for(Int_t j=0; j<3;++j)
	{
	  gr[j]->SetMarkerStyle(20);
	  gr[j]->SetMarkerColor(col[j]);
	  gr[j]->SetMarkerSize(1.5);
	  gr[j]->SetLineColor(col[j]);
	  switch(j)
	    {
	    case 0:
	      canvas->cd(3);
	      gr[j]->Draw("apl");
	      gr[j]->GetXaxis()->SetTitle("V_{bias} (V)");
	      gr[j]->GetYaxis()->SetTitle("Amplitude (mV)");
	      gr[j]->GetYaxis()->SetRangeUser(0,100);
	      gr[j]->GetXaxis()->SetRangeUser(0,360);
	      break;
	    case 1:
	      canvas->cd(3);
	      gr[j]->Draw("plsame");
	      gr[j]->GetYaxis()->SetRangeUser(0,100);
	      gr[j]->GetXaxis()->SetRangeUser(0,360);
	      break;
	    case 2:
	      canvas->cd(4);
	      gr[j]->Draw("apl");
	      gr[j]->GetXaxis()->SetTitle("V_{bias} (V)");
	      gr[j]->GetYaxis()->SetTitle("Ratio");
	      gr[j]->GetYaxis()->SetRangeUser(-0.05,0.15);
	      gr[j]->GetXaxis()->SetRangeUser(100,360);
	      break;
	    }
	}
      leg1->AddEntry(gr[0], "Pulse_{1}", "epl");
      leg1->AddEntry(gr[1], "Pulse_{2}", "epl");
      leg2->AddEntry(gr[2], "#frac{Pulse_{1}-Pulse_{2}}{Pulse_{1}}", "epl");
      canvas->cd(3);
      leg1->Draw();
      canvas->cd(4);
      leg2->Draw();
    }
  else
    {
      for(Int_t i=0; i<nEvents; ++i)
	{
	  TH1F *hSignal;
	  TH1F *hSignal1 = new TH1F("hSignal1", "", nPoints/2, 0, (t0+(nPoints/2)*dt)*1e9);
	  TH1F *hSignal2 = new TH1F("hSignal2", "", nPoints/2, 0, (t0+(nPoints/2)*dt)*1e9);
	  
	  hSignal = (TH1F*) lgad->_histo[0][i]->Clone();
	  hSignal->Scale(-1);
	  
	  //Split main signal into two individual signals
	  for(Int_t j=0; j<nPoints/2; ++j)
	    {
	      hSignal1->SetBinContent(j+1, hSignal->GetBinContent(j+1));
	      hSignal2->SetBinContent(j+1, hSignal->GetBinContent(nPoints/2+j+1));
	    }
	  
	  amp1 = hSignal1->GetMaximum();
	  amp2 = hSignal2->GetMaximum();
	  
	  cfd1 = CalcCFD(lgad, hSignal1, 0.5);
	  cfd2 = CalcCFD(lgad, hSignal2, 0.5);

	  //Fill Delta_amplitude and Delta_time  histogram
	  if(hSignal->GetMaximum() > 5)
	    {
	      deltaA->Fill(abs(amp1-amp2)/amp2);
	      deltaT->Fill(cfd1-cfd2);
	    }

	  if(i==0)
	    {
	      canvas->cd(1);	  
	      hSignal->Draw("");
	      hSignal->GetXaxis()->SetTitle("Time (nS)");
	      hSignal->GetYaxis()->SetTitle("Amplitude (mV)");
	      
	      canvas->cd(2);
	      hSignal1->Draw("");
	      hSignal1->GetXaxis()->SetTitle("Time (nS)");
	      hSignal1->GetYaxis()->SetTitle("Amplitude (mV)");
	      hSignal2->Draw("lsame");
	      hSignal2->GetXaxis()->SetTitle("Time (nS)");
	      hSignal2->GetYaxis()->SetTitle("Amplitude (mV)");
	    }
	  else
	    {
	      canvas->cd(1);	  
	      hSignal->Draw("lsame");
	      
	      canvas->cd(2);
	      hSignal1->Draw("lsame");
	      hSignal2->Draw("lsame"); 
	    }
	  hSignal->SetLineColor(kRed);
	  hSignal->SetLineWidth(3);
	  hSignal->GetYaxis()->SetRangeUser(-5, 100);
	  hSignal1->SetLineColor(kBlue);
	  hSignal1->SetLineWidth(3);
	  hSignal2->SetLineColor(kRed);
	  hSignal2->SetLineWidth(3);
	  hSignal1->GetYaxis()->SetRangeUser(-5, 100);
	  hSignal2->GetYaxis()->SetRangeUser(-5, 100);
	}

      Double_t height, mean, sigma, sigmaError, fitLow, fitHigh;
      
      canvas->cd(3);
      //gPad->SetLogy();
      height =  deltaA->GetBinContent(deltaA->GetMaximumBin());
      mean = deltaA->GetMean();
      sigma = deltaA->GetRMS();
      fitLow = mean - 0.1;
      fitHigh = mean + 0.1;
      TF1 *GaussA = new TF1("GaussA", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", fitLow, fitHigh);  
      GaussA->SetParameters(height, mean, sigma);
      deltaA->Draw("");
      deltaA->Fit("GaussA", "RMS+");
      deltaA->SetLineColor(kRed);
      deltaA->SetLineWidth(3);
      deltaA->GetXaxis()->SetTitle("(Amp_{2} -Amp_{1})/Amp_{2}");
      deltaA->GetYaxis()->SetTitle("Counts");
      deltaA->GetXaxis()->SetRangeUser(fitLow-1, fitHigh+1);
      TPaveStats *ptstatsA = new TPaveStats(0.6597008,0.5910006,0.9694723,0.9511307,"brNDC");
      ptstatsA->SetBorderSize(1);
      ptstatsA->SetFillColor(0);
      TText *ptstatsA_LaTex = ptstatsA->AddText("#sigma_{amplitude}");
      ptstatsA_LaTex = ptstatsA->AddText(Form("Events: %d",nEvents));
      ptstatsA_LaTex = ptstatsA->AddText(Form("#chi^{2}/ndf: %2.0lf / %d", GaussA->GetChisquare(), GaussA->GetNDF()));
      ptstatsA_LaTex = ptstatsA->AddText(Form("Mean : %2.2e #pm %2.2e", GaussA->GetParameter(1),GaussA->GetParError(1)));
      ptstatsA_LaTex = ptstatsA->AddText(Form("Sigma: %2.2e #pm %2.2e", TMath::Abs(GaussA->GetParameter(2)),GaussA->GetParError(2)));
      ptstatsA->Draw();
    
      canvas->cd(4);
      //gPad->SetLogy();
      height =  deltaT->GetBinContent(deltaT->GetMaximumBin());
      mean = deltaT->GetMean();
      sigma = deltaT->GetRMS();
      fitLow = mean - 0.2;
      fitHigh = mean + 0.2;
      TF1 *GaussT = new TF1("GaussT", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", fitLow, fitHigh);  
      GaussT->SetParameters(height, mean, sigma);
      deltaT->Draw("");
      deltaT->Fit("GaussT", "RMS+");
      deltaT->SetLineColor(kRed);
      deltaT->SetLineWidth(3);
      deltaT->GetXaxis()->SetTitle("#Delta_{time} (ns)");
      deltaT->GetYaxis()->SetTitle("Counts");
      deltaT->GetXaxis()->SetRangeUser(fitLow, fitHigh);
      TPaveStats *ptstatsT = new TPaveStats(0.6597008,0.5910006,0.9694723,0.9511307,"brNDC");
      ptstatsT->SetBorderSize(1);
      ptstatsT->SetFillColor(0);
      TText *ptstatsT_LaTex = ptstatsT->AddText("#sigma_{Jitter}");
      ptstatsT_LaTex = ptstatsT->AddText(Form("Events: %d",nEvents));
      ptstatsT_LaTex = ptstatsT->AddText(Form("#chi^{2}/ndf: %2.0lf / %d", GaussT->GetChisquare(), GaussT->GetNDF()));
      ptstatsT_LaTex = ptstatsT->AddText(Form("Mean : %2.2e #pm %2.2e ns", GaussT->GetParameter(1),GaussT->GetParError(1)));
      ptstatsT_LaTex = ptstatsT->AddText(Form("Sigma: %2.2e #pm %2.2e ns", TMath::Abs(GaussT->GetParameter(2)),GaussT->GetParError(2)));
      ptstatsT->Draw();
    }
  canvas->SaveAs(Form("../figures/ampVariationTimingSetup_%s.png", fileOutName),"png");
  return 0;
}

Float_t CalcCFD(AnalyzeTCTData *lgad, TH1F *his, Float_t thr)
{
  //Check Threshold
  if(thr > 1)
    thr /= 100;
    
  Float_t maxAmp = his->GetMaximum();
  Float_t t;
    
  vector<Float_t> xData;
  vector<Float_t> yData;
    
  for(Int_t i=0; i<his->GetNbinsX(); ++i)
    {
      xData.push_back(his->GetXaxis()->GetBinCenter(i+1));
      yData.push_back(his->GetBinContent(i+1));
    }
    
  t = lgad->LinearInterpolation(xData, yData, thr*maxAmp);
    
  xData.clear();
  yData.clear();
  return t;
}


void SetStyle(Bool_t threeD)
{
  gErrorIgnoreLevel=kError; //Removes annoying Potential memory leak warnings
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTopMargin(0.04854282);
  gStyle->SetPadBottomMargin(0.1353861);
  gStyle->SetPadLeftMargin(0.1418293);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kBlue);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.005,"y");
  gStyle->SetLabelOffset(0.010,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  if(threeD)
    {
      gStyle->SetPadRightMargin(0.18);
      gStyle->SetTitleOffset(1.15,"y");
    }
  else
    gStyle->SetTitleOffset(1.10,"y");
  gStyle->SetTitleOffset(0.95,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.03,"X");
  gStyle->SetTickLength(0.03,"Y");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.04);
  gStyle->SetNdivisions(505,"xy");
}
