#include "AnalyzeTCTData.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TLine.h"
#include "TF1.h"
#include "TPaveStats.h"
#include <map>
#include "TMap.h"

using namespace std;

void SetStyle(Bool_t threeD = false);

void analyzeChargeErrors();
void analyzeChargePINs();

// Define markers and colors
Int_t mCol[12] = {kRed, kRed, kRed, kBlue, kBlue, kBlue, kGreen+2, kGreen+2, kGreen+2, kMagenta, kMagenta, kMagenta};
Int_t mStyle[12] = {71, 20, 20, 72, 21, 21, 74, 33, 33, 75, 34, 34};
Int_t lStyle[12] = {7, 1, 7, 7, 1, 7, 7, 1, 7, 7, 1, 7};
Float_t mSize[12] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};

int main()
{
  //analyzeChargeErrors();
  analyzeChargePINs();
  return 0;
}  

void analyzeChargeErrors()
{
  SetStyle(false);
  Int_t nTotFiles = 9;
  TString tctFiles[nTotFiles] = {
    "../chargeData/lgad_VS_IR_b.tct",
    "../chargeData/lgad_VS_IR_c.tct",
    "../chargeData/lgad_VS_IR_d.tct",
    "../chargeData/lgad_VS_IR_e.tct",
    "../chargeData/lgad_VS_IR_f.tct",
    "../chargeData/lgad_VS_IR_g.tct",
    "../chargeData/lgad_VS_IR_h.tct",
    "../chargeData/lgad_VS_IR_i.tct",
    "../chargeData/lgad_VS_IR_j.tct"};
    
  TH1F *histC = new TH1F("histC", "", 200, 0, 2);
  TH1F *histN = new TH1F("histN", "", 200, 0, 1.5);
  TH1F *histRef = new TH1F("histRef", "", 10, 0.1, 0.2);
  
  AnalyzeTCTData *lgadRef = new AnalyzeTCTData("../chargeData/lgad_VS_IR_a.tct");
  lgadRef->CorrectBaseline();
  lgadRef->CalcNoise();
  lgadRef->SetIntegralLimits(0, 50.);
  lgadRef->CalculateWaveformProperties();

  Int_t nevRef=0, nev=0; 
  for(Int_t k=0;k<nTotFiles;k++)
    {
      if(gSystem->AccessPathName(tctFiles[k]))
	{
	  cout<<"\033[1;33m File not Found: "<<tctFiles[k]<<"\033[0m"<<endl;
	  cout<<"\033[1;32m SKIPPING THE ANALYSIS FOR THIS FILE!!!!\033[0m"<<endl;
	  continue;
	}
      //Read data  

      AnalyzeTCTData *lgad = new AnalyzeTCTData(tctFiles[k]);
      lgad->CorrectBaseline();
      lgad->CalcNoise();
      lgad->SetIntegralLimits(0, 50.);
      lgad->CalculateWaveformProperties();

      for(Int_t i=0; i<lgad->_events; ++i)
	{
	  // Store Laser reference amplitude
	  if(lgad->_bmON)
	    {
	      histRef->Fill(lgad->_bmValue[i]);
	      nevRef++;
	    }

	  // Fill noise Histogram Normalized to the Reference file
	  histN->Fill(lgad->_sigNoise[0][i]/lgadRef->_sigNoise[0][i]);
	  
	  // Fill charge Histogram Normalized to the Reference file
	  histC->Fill(lgad->_sigNormCharge[0][i]/lgadRef->_sigNormCharge[0][i]);
	 
	  nev++;
	}
      	    
    }

  //histC->Scale(1/nev);
  Float_t mean, height, stdDev;
  TF1 *Gaus;
  
  TCanvas *canN = new TCanvas("canN","",1200,1000);
  histN->GetXaxis()->SetTitle("Noise (mV)");
  histN->GetYaxis()->SetTitle("Counts");
  histN->SetLineColor(kBlue);
  histN->SetLineWidth(2);
  histN->Draw();
  mean = histN->GetMean();
  height= histN->GetBinContent(histN->GetMaximumBin());
  stdDev = histN->GetStdDev();
  Gaus = new TF1("Gaus", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", -2, 2);
  Gaus->SetParameters(height, mean, stdDev);
  histN->Fit("Gaus", "RMS+");
  Gaus->SetLineColor(kRed);
  Gaus->SetLineWidth(2);
  Gaus->Draw("LSAME");
  TPaveStats *ptstats1 = new TPaveStats(0.6,0.75,0.9,0.951,"brNDC");
  ptstats1->SetBorderSize(1);
  ptstats1->SetFillColor(0);
  TText *ptstats1_LaTex = ptstats1->AddText("NOISE");
  ptstats1_LaTex = ptstats1->AddText(Form("Events: %d",nev));
  ptstats1_LaTex = ptstats1->AddText(Form("#chi^{2}/ndf: %2.0lf / %d", Gaus->GetChisquare(), Gaus->GetNDF()));
  ptstats1_LaTex = ptstats1->AddText(Form("Mean : %2.2e #pm %2.2e mV", Gaus->GetParameter(1),Gaus->GetParError(1)));
  ptstats1_LaTex = ptstats1->AddText(Form("Sigma: %2.2e #pm %2.2e mV", TMath::Abs(Gaus->GetParameter(2)),Gaus->GetParError(2)));
  ptstats1->Draw();
  canN->SaveAs("../figures/spaceLgadsNoiseHistogram.png","png");
  
  TCanvas *canC = new TCanvas("canC","",1200,1000);
  histC->GetXaxis()->SetTitle("Norm. charge (arb.)");
  histC->GetYaxis()->SetTitle("Counts");  
  histC->SetLineColor(kBlue);
  histC->SetLineWidth(2);
  histC->Draw();
  histC->GetXaxis()->SetRangeUser(0.5, 1.5);
  mean = histC->GetMean();
  height= histC->GetBinContent(histC->GetMaximumBin());
  stdDev = histC->GetStdDev();
  Gaus = new TF1("Gaus", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", 0, 2);
  Gaus->SetParameters(height, mean, stdDev);
  histC->Fit("Gaus", "RMS+");
  Gaus->SetLineColor(kRed);
  Gaus->SetLineWidth(2);
  Gaus->Draw("LSAME");
  TPaveStats *ptstats2 = new TPaveStats(0.6,0.75,0.9,0.951,"brNDC");
  ptstats2->SetBorderSize(1);
  ptstats2->SetFillColor(0);
  TText *ptstats2_LaTex = ptstats2->AddText("CHARGE");
  ptstats2_LaTex = ptstats2->AddText(Form("Events: %d",nev));
  ptstats2_LaTex = ptstats2->AddText(Form("#chi^{2}/ndf: %2.0lf / %d", Gaus->GetChisquare(), Gaus->GetNDF()));
  ptstats2_LaTex = ptstats2->AddText(Form("Mean : %2.2e #pm %2.2e arb.", Gaus->GetParameter(1),Gaus->GetParError(1)));
  ptstats2_LaTex = ptstats2->AddText(Form("Sigma: %2.2e #pm %2.2e arb.", TMath::Abs(Gaus->GetParameter(2)),Gaus->GetParError(2)));
  ptstats2->Draw();
  canC->SaveAs("../figures/spaceLgadsChargeHistogram.png","png");

  //Laser stability check
  TCanvas *canLas = new TCanvas("canLas","",1200,1000);
  TLegend *legLas = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
  histRef->GetXaxis()->SetTitle("Reference amplitude (arb.)");
  histRef->GetYaxis()->SetTitle("Counts");
  histRef->SetLineColor(kBlue);
  histRef->SetLineWidth(2);
  histRef->Draw();
  mean = histRef->GetMean();
  height = histRef->GetBinContent(histRef->GetMaximumBin());
  stdDev = histRef->GetStdDev();
  Gaus = new TF1("Gaus", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", -2, 2);
  Gaus->SetParameters(height, mean, stdDev);
  histRef->Fit("Gaus", "RMS+");
  Gaus->SetLineColor(kRed);
  Gaus->SetLineWidth(2);
  Gaus->Draw("LSAME");
  TPaveStats *ptstats3 = new TPaveStats(0.6,0.75,0.9,0.951,"brNDC");
  ptstats3->SetBorderSize(1);
  ptstats3->SetFillColor(0);
  TText *ptstats3_LaTex = ptstats3->AddText("Laser stability");
  ptstats3_LaTex = ptstats3->AddText(Form("Events: %d",nevRef));
  ptstats3_LaTex = ptstats3->AddText(Form("#chi^{2}/ndf: %2.0lf / %d", Gaus->GetChisquare(), Gaus->GetNDF()));
  ptstats3_LaTex = ptstats3->AddText(Form("Mean : %2.2e #pm %2.2e arb.", Gaus->GetParameter(1),Gaus->GetParError(1)));
  ptstats3_LaTex = ptstats3->AddText(Form("Sigma: %2.2e #pm %2.2e arb.", TMath::Abs(Gaus->GetParameter(2)),Gaus->GetParError(2)));
  ptstats3->Draw();
  canLas->SaveAs("../figures/spaceLgadsLaserStabilityHistogram.png","png");
}

void analyzeChargePINs()
{
  SetStyle(false);
  Int_t nPinFiles = 6;
  TString tctFiles[nPinFiles] = {
    "../pinData/pin_W_3_type_1_pad_C_a.tct", "../pinData/pin_W_3_type_1_pad_C_b.tct",
    "../pinData/pin_W_9_type_1_pad_C_a.tct", "../pinData/pin_W_9_type_1_pad_C_b.tct",
    "../pinData/pin_W_12_type_1_pad_C_a.tct", "../pinData/pin_W_12_type_1_pad_C_b.tct"};
 
  Float_t thickness[3] = {50.0, 100.0, 150.0};
  Float_t charge[3];
  Float_t chargeErr[3];
  
  for(Int_t k=0;k<nPinFiles;k+=2)
    {
      if(gSystem->AccessPathName(tctFiles[k]))
	{
	  cout<<"\033[1;33m File not Found: "<<tctFiles[k]<<"\033[0m"<<endl;
	  cout<<"\033[1;32m SKIPPING THE ANALYSIS FOR THIS FILE!!!!\033[0m"<<endl;
	  continue;
	}

      //Read data  
      AnalyzeTCTData *pin1 = new AnalyzeTCTData(tctFiles[k], 0.5);
      pin1->SetAveragesInOscilloscope(512);
      pin1->CorrectBaseline();
      pin1->CalcNoise();
      pin1->SetIntegralLimits(0, 30.);
      pin1->CalculateWaveformProperties();

      AnalyzeTCTData *pin2 = new AnalyzeTCTData(tctFiles[k+1], 0.5);
      pin2->SetAveragesInOscilloscope(512);
      pin2->CorrectBaseline();
      pin2->CalcNoise();
      pin2->SetIntegralLimits(0, 30.);
      pin2->CalculateWaveformProperties();

      //Number of Voltage points for PIN
      Int_t nVP1 = pin1->_nV1;
      Int_t nVP2 = pin2->_nV1;
      Int_t stepP = TMath::Abs(pin1->_tct->V1[1]-pin1->_tct->V1[0]);

      // vector<Float_t> chargeP(nVP1+nVP2-20);
      // vector<Float_t> errP(nVP1+nVP2-20);
      
      vector<Float_t> chargeP(nVP1-10);
      vector<Float_t> errP(nVP1-10);
      
      for(Int_t i=10; i< nVP1; ++i)
	{
	  chargeP.push_back(pin1->_sigNormCharge[0][i]);
	  //chargeP.push_back(pin2->_sigNormCharge[0][i]);
	  errP.push_back(pin1->_sigChargeError[0][i]);
	  //errP.push_back(pin2->_sigChargeError[0][i]);
	}

      //Average Charge in Pin
      Float_t avgQ = TMath::Mean(chargeP.begin(),chargeP.end());

      //Error in Average Charge
      //Float_t errAvg = TMath::StdDev(chargeP.begin(),chargeP.end())/TMath::Sqrt(nVP1+nVP2-20);
      Float_t errAvg = TMath::StdDev(chargeP.begin(),chargeP.end())/TMath::Sqrt(nVP1-10);

      charge[k/2] = avgQ;
      chargeErr[k/2] = errAvg;
    }

  TCanvas *canvasC = new TCanvas("canvasC","",1200,1000);
  TGraphErrors *gr = new TGraphErrors(3, thickness, charge, 0, chargeErr);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerColor(kRed);
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->Draw("apl");
  gr->GetXaxis()->SetTitle("Thickness [#mum]");  
  gr->GetYaxis()->SetTitle("Norm. charge (arb.)");
  canvasC->SaveAs("../figures/spaceLgadsChargePin_PadC.png","png");
}

void SetStyle(Bool_t threeD)
{
  gErrorIgnoreLevel=kError; //Removes annoying Potential memory leak warnings
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
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

/*
  TGraphErrors *grRTPin = new TGraphErrors(sizeof(rtPin[]),thickness,rtPin,0,rtPinError);
  grRTPin->SetMarkerStyle(20);
  grRTPin->SetMarkerColor(kBlue);
  grRTPin->SetLineColor(kBlue);
  grRTPin->SetMarkerSize(1.5);
  grRTPin->Draw("apl");
  grRTPin->GetXaxis()->SetTitle("Thickness (#mum)");
  grRTPin->GetYaxis()->SetTitle("Jitter|_{V_{bias}=240 V} (ps)");
  grRTPin->GetYaxis()->SetRangeUser(0, 3000);
  legJPin->AddEntry(grRTPin, "Rise-Time", "epl");
  legJPin->Draw();

  TGraphErrors *grJPin = new TGraphErrors(3,thickness,JitterPin,0,JitterPinError);
  grJPin->SetMarkerStyle(21);
  grJPin->SetMarkerColor(kRed);
  grJPin->SetLineColor(kRed);
  grJPin->SetMarkerSize(1.5);
  grJPin->Draw("EPLSAME");
      
  legJPin->AddEntry(grJPin, "Jitter", "epl");
*/
  
/*
//================================================Gain=======================================================
if(processGain)
{
TCanvas *canG = new TCanvas("canG","",1200,1000);
TLegend *legG = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
legG->SetHeader("Pad area - 1 cm^{2}", "C");
      
TGraphErrors *gr;
TH1F *histG = new TH1F("histG", "", 700, 0, 700);
      
TString gainFiles[nTotFiles] = {
"../data/slapp_1MIP_W_1_pad_C_type_0_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_1_pad_C_type_2_gain_0_VS_IR_a.tct",
"../data/slapp_1MIP_W_1_pad_C_type_1_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_1_pad_C_type_2_gain_0_VS_IR_a.tct",
"../data/slapp_1MIP_W_1_pad_C_type_2_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_1_pad_C_type_2_gain_0_VS_IR_a.tct",

"../data/slapp_1MIP_W_9_pad_C_type_0_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_9_pad_C_type_2_gain_0_VS_IR_a.tct",
"../data/slapp_1MIP_W_9_pad_C_type_1_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_9_pad_C_type_2_gain_0_VS_IR_a.tct",
"../data/slapp_1MIP_W_9_pad_C_type_2_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_9_pad_C_type_2_gain_0_VS_IR_a.tct",

"../data/slapp_1MIP_W_12_pad_C_type_0_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_12_pad_C_type_2_gain_0_VS_IR_a.tct",
"../data/slapp_1MIP_W_12_pad_C_type_1_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_12_pad_C_type_2_gain_0_VS_IR_a.tct",
"../data/slapp_1MIP_W_12_pad_C_type_2_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_12_pad_C_type_2_gain_0_VS_IR_a.tct"};
	
for(Int_t k=0;k<nTotFiles;k+=2)
{
info = gainFiles[k].Tokenize("_");
wafer = ((TObjString*)(info->At(3)))->String();
pad   = ((TObjString*)(info->At(5)))->String();
type  = ((TObjString*)(info->At(7)))->String();
sensor= ((TObjString*)(info->At(9)))->String();
legend = "Wafer-" + wafer;

cout<<"\033[1;32mProcessing File: "<<gainFiles[k]<<"\033[0m"<<endl;

//Read LGAD data
AnalyzeTCTData lgad(gainFiles[k]);
lgad.CorrectBaseline();
lgad.CalcNoise();
lgad.SetIntegralLimits(0, 50.);
lgad.CalculateWaveformProperties();
        
cout<<"\033[1;32mProcessing File: "<<gainFiles[k+1]<<"\033[0m"<<endl;
//Read PIN data
AnalyzeTCTData pin(gainFiles[k+1], 0.9);
pin.CorrectBaseline();
pin.CalcNoise();
pin.SetIntegralLimits(0, 50.);
pin.CalculateWaveformProperties();
        
Int_t nVL = lgad._nV1; //Number of Voltage points for LGAD
Int_t nVP = pin._nV1; //Number of Voltage points for PIN

vector<Float_t> chargeL(nVL);
vector<Float_t> errL(nVL);
vector<Float_t> chargeP(nVP);
vector<Float_t> errP(nVP);

//Int_t arrS = (nVL-10)/5 + 1;
Int_t arrS = nVL;
Float_t voltage[arrS];
Float_t gain[arrS];
Float_t errG[arrS];
        
for(Int_t i=0; i< nVP; ++i)
{
chargeP[i] = pin._sigNormCharge[0][i];
errP[i] = pin._sigChargeError[0][i];
histRef->Fill(pin._bmValue[i]);
nEv++;
}

//Average Charge in Pin
Float_t avgQ = TMath::Mean(chargeP.begin()+10,chargeP.end());
//Error in Average Charge
Float_t errAvg = TMath::StdDev(chargeP.begin()+10,chargeP.end())/TMath::Sqrt(nVP-10);
        
Float_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);


// for(Int_t i=10; i<nVL; i+=5)
// {
// voltage[(i-10)/5] = i*step;
// chargeL[i] = lgad._sigNormCharge[0][i];
// errL[i] = lgad._sigChargeError[0][i];
// gain[(i-10)/5] = lgad._sigNormCharge[0][i]/avgQ;
// errG[(i-10)/5] = gain[(i-10)/5]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
// }

//Store Gain for each file to plot Jitter v/s Gain
Gain[k/2] = new Float_t[nVL];
GainError[k/2] = new Float_t[nVL];
      
for(Int_t i=0; i<arrS; ++i)
{
voltage[i] = i*step;
chargeL[i] = lgad._sigNormCharge[0][i];
errL[i] = lgad._sigChargeError[0][i];
gain[i] = lgad._sigNormCharge[0][i]/avgQ;
errG[i] = gain[i]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
Gain[k/2][i] = gain[i];
GainError[k/2][i] = errG[i];
}

gr = new TGraphErrors(arrS,voltage,gain,0,errG);
gr->SetMarkerStyle(mStyle[k/2]);
gr->SetMarkerColor(mCol[k/2]);
gr->SetLineStyle(lStyle[k/2]);
gr->SetLineColor(mCol[k/2]);
gr->SetMarkerSize(mSize[k/2]);
      
if(k==0)
{
TGraphErrors *grT = new TGraphErrors(arrS,voltage,gain,0,errG);
grT->SetMarkerColor(kBlack);
grT->SetLineColor(kBlack);
grT->SetLineStyle(1);
grT->SetLineWidth(4);

TGraphErrors *hrT = new TGraphErrors(arrS,voltage,gain,0,errG);
hrT->SetMarkerColor(kBlack);
hrT->SetLineColor(kBlack);
hrT->SetLineStyle(7);
hrT->SetLineWidth(4);

histG->Draw();
histG->GetXaxis()->SetTitle("V_{bias} (V)");
histG->GetYaxis()->SetTitle("Gain");
gr->Draw("eplsame");
legT->Draw();
}
else
gr->Draw("eplsame");
if(type=="1")
legG->AddEntry(gr, legend, "epl");
}
histG->GetYaxis()->SetRangeUser(0,120);
histG->GetXaxis()->SetRangeUser(0,650);
legG->Draw();
canG->SaveAs("../figures/spaceLgadsGain_PadC.png","png");
}

*/
