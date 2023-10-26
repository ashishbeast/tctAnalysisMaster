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

void analyzeSlappPads(const char *pad, const char *area);

// Define markers and colors
Int_t mCol[12] = {kRed, kRed, kRed, kBlue, kBlue, kBlue, kGreen+2, kGreen+2, kGreen+2, kMagenta, kMagenta, kMagenta};
Int_t mStyle[12] = {71, 20, 20, 72, 21, 21, 74, 33, 33, 75, 34, 34};
Int_t lStyle[12] = {7, 1, 7, 7, 1, 7, 7, 1, 7, 7, 1, 7};
Float_t mSize[12] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};

// Save info of each file while processing the data
TObjArray *info;
TString wafer, type, padType, legend, sensor;

// Define no of wafers to be processed
const Int_t nWaf = 3;
const Int_t nTypes = 3;
const Int_t nFiles = nWaf * nTypes;
const Int_t nTotFiles = nWaf * nTypes * 3;

Int_t wafers[nWaf] = {1, 9, 12};
  
// Save Jitter and Gain from the files
Int_t arrSize[nFiles];
Float_t *Jitter[nFiles];
Float_t *JitterError[nFiles];
Float_t *Gain[nFiles];
Float_t *GainError[nFiles];
Float_t thickness[nWaf] = {50.0, 100.0, 150.0};
Float_t ChargePin[nWaf];
Float_t ChargePinError[nWaf];
Int_t arrSizePin[nWaf];
Float_t *vBiasPin[nWaf];
Float_t *noisePin[nWaf];
Float_t *JitterPin[nWaf];
Float_t *JitterPinError[nWaf];
Float_t *rtPin[nWaf];
Float_t *rtPinError[nWaf];

// Switches for plotting and DEBUGGING
Bool_t processJitter          = kTRUE;
Bool_t plotLaserStability     = kTRUE;
Bool_t plotSignals            = kFALSE;
Bool_t plotJitterVsGain       = kFALSE;
Bool_t plotPinChargeThickness = kTRUE;

// Save PulseShape at 200 V for each file
TH1F *pulseShape[nFiles];

// Save the reference amplitude of beam monitor into an histogram
TH1F *histRef = new TH1F("histRef","",400, 1, 1.5);
Int_t nEv = 0;

// Global Legend for the layout types
TLegend *legT, *legTJ;

int main()
{ 
  // Three different pad sizes in space LGADS
  // A: pad with active area of 6.25 mm^2
  // B: pad with active area of 25 mm^2
  // C: pad with active area of 100 mm^2
  TMap *padArea = new TMap();
  padArea->Add(new TObjString("A"), new TObjString("Pad area - 6.25 mm^{2}"));
  padArea->Add(new TObjString("B"), new TObjString("Pad area - 25 mm^{2}"));
  padArea->Add(new TObjString("C"), new TObjString("Pad area - 100 mm^{2}"));
  
  TIter iter(padArea->MakeIterator());
  for(auto dict = iter.Next(); dict != nullptr; dict = iter.Next())
    {
      TObjString *padKey = dynamic_cast<TObjString*>(dict);
      const char *pad = padKey->GetString().Data();
      const char *area = static_cast<const char *>(padArea->GetValue(pad)->GetName());

      SetStyle(false);
      legT = new TLegend(0.15,0.8,0.4,0.95,NULL,"brNDC");
      legTJ = new TLegend(0.40,0.8,0.55,0.95,NULL,"brNDC");
      legT->SetHeader("Layout type", "L");
      legTJ->SetHeader("Layout type", "L");
      
      analyzeSlappPads(pad, area);
    }
  return 0;
}

void analyzeSlappPads(const char *pad, const char *area)
{
  if(pad == nullptr || area == nullptr)
    return;
  
  SetStyle(false);
  if(processJitter)
    {
      TString tctFiles[nTotFiles] = {
	Form("../data/slapp_1MIP_W_1_pad_%s_type_0_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_1_pad_%s_type_0_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_1_pad_%s_type_0_gain_0_VS_IR_a.tct", pad),
	Form("../data/slapp_1MIP_W_1_pad_%s_type_1_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_1_pad_%s_type_1_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_1_pad_%s_type_1_gain_0_VS_IR_a.tct", pad),
	Form("../data/slapp_1MIP_W_1_pad_%s_type_2_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_1_pad_%s_type_2_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_1_pad_%s_type_2_gain_0_VS_IR_a.tct", pad),
	
	Form("../data/slapp_1MIP_W_9_pad_%s_type_0_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_9_pad_%s_type_0_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_9_pad_%s_type_0_gain_0_VS_IR_a.tct", pad),
	Form("../data/slapp_1MIP_W_9_pad_%s_type_1_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_9_pad_%s_type_1_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_9_pad_%s_type_1_gain_0_VS_IR_a.tct", pad),
	Form("../data/slapp_1MIP_W_9_pad_%s_type_2_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_9_pad_%s_type_2_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_9_pad_%s_type_0_gain_0_VS_IR_a.tct", pad),

	Form("../data/slapp_1MIP_W_12_pad_%s_type_0_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_12_pad_%s_type_0_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_12_pad_%s_type_0_gain_0_VS_IR_a.tct", pad),
	Form("../data/slapp_1MIP_W_12_pad_%s_type_1_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_12_pad_%s_type_1_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_12_pad_%s_type_1_gain_0_VS_IR_a.tct", pad),
	Form("../data/slapp_1MIP_W_12_pad_%s_type_2_gain_1_VS_IR_a.tct", pad), Form("../data/slapp_1MIP_W_12_pad_%s_type_2_gain_1_VS_IR_a_Noise.tct", pad), Form("../data/slapp_1MIP_W_12_pad_%s_type_2_gain_0_VS_IR_a.tct", pad)};

      TH1F *histJ[7];
      for(Int_t i=0; i<7;++i)
	histJ[i] = new TH1F("histJ", "", 700, 0, 700);
      
      TGraphErrors *grN, *grA, *grC, *grSNR, *grRT, *grJ, *grG;

      TCanvas *canN = new TCanvas("canN","",1200,1000);
      TLegend *legN = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legN->SetHeader(area, "C");

      TCanvas *canA = new TCanvas("canA","",1200,1000);
      TLegend *legA = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legA->SetHeader(area, "C");

      TCanvas *canC = new TCanvas("canC","",1200,1000);
      TLegend *legC = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legC->SetHeader(area, "C");

      TCanvas *canSNR = new TCanvas("canSNR","",1200,1000);
      TLegend *legSNR = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legSNR->SetHeader(area, "C");

      TCanvas *canRT = new TCanvas("canRT","",1200,1000);
      TLegend *legRT = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legRT->SetHeader(area, "C");

      TCanvas *canJ = new TCanvas("canJ","",1200,1000);
      TLegend *legJ = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legJ->SetHeader(area, "C");

      TCanvas *canG = new TCanvas("canG","",1200,1000);
      TLegend *legG = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legG->SetHeader(area, "C");
      
      for(Int_t k=0;k<nTotFiles;k+=3)
	{
	  if(gSystem->AccessPathName(tctFiles[k]))
	    {
	      cout<<"\033[1;33m File not Found: "<<tctFiles[k]<<"\033[0m"<<endl;
	      cout<<"\033[1;32m SKIPPING THE ANALYSIS FOR THIS FILE!!!!\033[0m"<<endl;
	      continue;
	    }
	  info = tctFiles[k].Tokenize("_");
	  wafer = ((TObjString*)(info->At(3)))->String();
	  //pad   = ((TObjString*)(info->At(5)))->String();
	  type  = ((TObjString*)(info->At(7)))->String();
	  sensor= ((TObjString*)(info->At(9)))->String();
	  legend = "Wafer-" + wafer;

	  //Read data  
	  // lgad File -> Amplitude, Charge, SlewRate, Jitter
	  cout<<"\033[1;32m Processing File (LGAD): "<<tctFiles[k]<<"\033[0m"<<endl;
	  AnalyzeTCTData lgad(tctFiles[k], 5);
	  lgad.CorrectBaseline();
	  lgad.CalcNoise();
	  lgad.SetIntegralLimits(0, 50.);
	  lgad.CalculateWaveformProperties();
        
	  // noise File -> Noise for Jitter
	  cout<<"\033[1;32m Processing File (NOISE): "<<tctFiles[k+1]<<"\033[0m"<<endl;
	  AnalyzeTCTData noise(tctFiles[k+1], 7.2);
	  noise.CorrectBaseline();
	  noise.CalcNoise();
	  noise.CalculateWaveformProperties();

	  // pin File -> Charge
	  cout<<"\033[1;32m Processing File (PIN): "<<tctFiles[k+2]<<"\033[0m"<<endl;
	  AnalyzeTCTData pin(tctFiles[k+2], 0.9);
	  pin.CorrectBaseline();
	  pin.CalcNoise();
	  pin.SetIntegralLimits(0, 50.);
	  pin.CalculateWaveformProperties();

	  //Number of Voltage points for LGAD
	  Int_t nVL = lgad._nV1; 
	  Int_t stepL = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);

	  //Number of Voltage points for PIN
	  Int_t nVP = pin._nV1;
	  Int_t stepP = TMath::Abs(pin._tct->V1[1]-pin._tct->V1[0]);

	  // Array of applied voltage
	  Float_t vBias[nVL];
	  Float_t noiseSig[nVL];
	  Float_t amp[nVL];
	  Float_t snr[nVL];
	  Float_t charge[nVL];
	  Float_t chargeErr[nVL];
	  Float_t gain[nVL];
	  Float_t gainErr[nVL];
	  
	  // Reduced number of points
	  Int_t arrS = (nVL-20)/5 + 1;
	  Float_t vBiasReduced[arrS];
	  Float_t rt[arrS];
	  Float_t rtErr[arrS];
	  Float_t jitter[arrS];
	  Float_t jitterErr[arrS];
	 
	  vector<Float_t> chargeP(nVP);
	  vector<Float_t> errP(nVP);
	  
	  for(Int_t i=0; i< nVP; ++i)
	    {
	      //chargeP[i] = pin._sigNormCharge[0][i];
	      chargeP[i] = pin._sigCharge[0][i];
	      errP[i] = pin._sigChargeError[0][i];
	      histRef->Fill(pin._bmValue[i]);
	      nEv++;
	    }
	  //Average Charge in Pin
	  Float_t avgQ = TMath::Mean(chargeP.begin()+10,chargeP.end());
	  //Error in Average Charge
	  Float_t errAvg = TMath::StdDev(chargeP.begin()+10,chargeP.end())/TMath::Sqrt(nVP-10);

	  for(Int_t i=0; i<nVL; ++i)
	    {
	      // Store Laser reference amplitude
	      if(lgad._bmON)
		{
		  histRef->Fill(lgad._bmValue[i]);
		  nEv++;
		}
	      if(noise._bmON)
		{
		  histRef->Fill(noise._bmValue[i]);
		  nEv++;
		}

	      // Store bias voltages
	      vBias[i] = i*stepL;

	      // Store noise in the pre signal region at the given bias voltage
	      noiseSig[i] = noise._sigNoise[0][i];

	      // Store the signal peak value at the given bias voltage
	      amp[i] = lgad._sigNormAmplitude[0][i];

	      // Store signal charge at the given bias voltage
	      charge[i] = lgad._sigNormCharge[0][i];
	      chargeErr[i] = lgad._sigChargeError[0][i];

	      // Calculate and store the signal to noise ratio
	      snr[i] = amp[i]/noiseSig[i];

	      // Evaluate gain of the sensor
	      gain[i] = lgad._sigNormCharge[0][i]/avgQ;
	      gainErr[i] = gain[i]*TMath::Sqrt(TMath::Power(lgad._sigChargeError[0][i]/lgad._sigNormCharge[0][i],2)+TMath::Power(errAvg/avgQ,2));
	      
	      if( i>=20 && i%5==0)
		{
		  // Store bias voltages
		  vBiasReduced[(i-20)/5] = i*stepL;

		  // Store the rise time of the signal and the error in RT
		  rt[(i-20)/5] = lgad._sigRiseTime[0][i];
		  rtErr[(i-20)/5] = lgad._sigRiseTimeError[0][i];

		  // Evaluate jitter at the given bias voltage
		  if(lgad._sigSlewRate[0][i] == 0)
		    {
		      jitter[(i-20)/5] = 0;
		      jitterErr[(i-20)/5] = 0;
		    }
		  else
		    {
		      jitter[(i-20)/5] = 1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
		      jitterErr[(i-20)/5] = (jitter[(i-20)/5] / lgad._sigSlewRate[0][i]) * lgad._sigSlewRateError[0][i];
		    }
		}
	    }

	  /*
	    for(Int_t i=20; i<nVL; i+=5)
	    {
	    // Store bias voltages
	    vBiasReduced[(i-20)/5] = i*stepL;

	    // Store noise in the pre signal region at the given bias voltage
	    noiseSig[(i-20)/5] = noise._sigNoise[0][i];

	    // Store the signal peak value at the given bias voltage
	    amp[(i-20)/5] = lgad._sigNormAmplitude[0][i];

	    // Store signal charge at the given bias voltage
	    charge[(i-20)/5] = lgad._sigNormCharge[0][i];
	    chargeErr[(i-20)/5] = lgad._sigChargeError[0][i];
	  
	    // Store the rise time of the signal and the error in RT
	    rt[(i-20)/5] = lgad._sigRiseTime[0][i];
	    rtErr[(i-20)/5] = lgad._sigRiseTimeError[0][i];

	    // Calculate and store the signal to noise ratio
	    snr[(i-20)/5] = amp[(i-20)/5]/noiseSig[(i-20)/5];

	    // Evaluate jitter at the given bias voltage
	    if(lgad._sigSlewRate[0][i] == 0)
	    {
	    jitter[(i-20)/5] = 0;
	    jitterErr[(i-20)/5] = 0;
	    }
	    else
	    {
	    jitter[(i-20)/5] = 1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
	    jitterErr[(i-20)/5] = (jitter[(i-20)/5] / lgad._sigSlewRate[0][i]) * lgad._sigSlewRateError[0][i];
	    }

	    // Evaluate gain of the sensor
	    gain[(i-20)/5] = lgad._sigNormCharge[0][i]/avgQ;
	    gainErr[(i-20)/5] = gain[(i-20)/5]*TMath::Sqrt(TMath::Power(lgad._sigChargeError[0][i]/lgad._sigNormCharge[0][i],2)+TMath::Power(errAvg/avgQ,2));
	    }
	  */

	  // Store Jitter for each file to plot Jitter v/s Gain
	  arrSize[k/3] = nVL;
	  Jitter[k/3] = new Float_t[nVL];
	  JitterError[k/3] = new Float_t[nVL];

	  // Store Gain for each file to plot Jitter v/s Gain
	  Gain[k/3] = new Float_t[nVL];
	  GainError[k/3] = new Float_t[nVL];
	  
	  if(type=="2")
	    {
	      arrSizePin[k/9] = nVP;
	      ChargePin[k/9] = avgQ;
	      ChargePinError[k/9] = errAvg;
	      vBiasPin[k/9] = new Float_t[nVP];
	      noisePin[k/9] = new Float_t[nVP];
	      JitterPin[k/9] = new Float_t[nVP];
	      JitterPinError[k/9] = new Float_t[nVP];
	      rtPin[k/9] = new Float_t[nVP];
	      rtPinError[k/9] = new Float_t[nVP];
	    }

	  for(Int_t i=0; i<nVL; ++i)
	    {
	      Jitter[k/3][i] = 1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
	      JitterError[k/3][i] = (Jitter[k/3][i]/ lgad._sigSlewRate[0][i]) * lgad._sigSlewRateError[0][i];
	      Gain[k/3][i] = lgad._sigNormCharge[0][i]/avgQ;
	      GainError[k/3][i] = Gain[k/3][i]*TMath::Sqrt(TMath::Power(lgad._sigChargeError[0][i]/lgad._sigNormCharge[0][i],2)+TMath::Power(errAvg/avgQ,2));
	      //Store Signal at 240 V
	      if(i*stepL == 240)
		{
		  pulseShape[k/3] = (TH1F*) lgad._histo[0][i]->Clone();
		  pulseShape[k/3]->Scale(lgad._polarity);
		  //pulseShape[k/3]->Scale(1/lgad._bmValue[i]);
		}
	    }
	  
	  // Estimate the jitter of PIN sensors
	  if(type=="2")
	    for(Int_t i=0; i<nVP; ++i)
	      {
		vBiasPin[k/9][i] = i*stepP;
		noisePin[k/9][i] = pin._sigNoise[0][i];
		JitterPin[k/9][i] = (1000 * pin._sigNoise[0][i] * pin._sigRiseTime[0][i]) / (0.8*pin._sigAmplitude[0][i]);
		JitterPinError[k/9][i] = 0;
		rtPin[k/9][i] = 1000 * pin._sigRiseTime[0][i];
		rtPinError[k/9][i] = TMath::Sqrt(2) * pin._deltaT;
	      }

	  /*
	    if(type=="2")
	    {
	    ChargePin[k/9] = avgQ;
	    ChargePinError[k/9] = errAvg;
	    if(i*stepP == 240)
	    {
	    //cout<<Form("\033[1;31mPin Noise: %.2lf mV \t Pin RT: %.2lf ps \033[0m\n", pin._sigNoise[0][i], 1000 * pin._sigRiseTime[0][i]);
	    JitterPin[k/9] = (1000 * pin._sigNoise[0][i] * pin._sigRiseTime[0][i]) / (0.8*pin._sigAmplitude[0][i]);
	    JitterPinError[k/9] = 0;
	    rtPin[k/9] = 1000 * pin._sigRiseTime[0][i];
	    rtPinError[k/9] = TMath::Sqrt(2) * pin._deltaT;
	    }
	    }
	  */
	    
	
	  grN = new TGraphErrors(nVL,vBias,noiseSig,0,0);
	  grN->SetMarkerStyle(mStyle[k/3]);
	  grN->SetMarkerColor(mCol[k/3]);
	  grN->SetLineColor(mCol[k/3]);
	  grN->SetLineStyle(lStyle[k/3]);
	  grN->SetLineWidth(2);
	  grN->SetMarkerSize(mSize[k/3]);

	  grA = new TGraphErrors(nVL,vBias,amp,0,noiseSig);
	  grA->SetMarkerStyle(mStyle[k/3]);
	  grA->SetMarkerColor(mCol[k/3]);
	  grA->SetLineColor(mCol[k/3]);
	  grA->SetLineStyle(lStyle[k/3]);
	  grA->SetLineWidth(2);
	  grA->SetMarkerSize(mSize[k/3]);

	  grC = new TGraphErrors(nVL,vBias,charge,0,chargeErr);
	  grC->SetMarkerStyle(mStyle[k/3]);
	  grC->SetMarkerColor(mCol[k/3]);
	  grC->SetLineColor(mCol[k/3]);
	  grC->SetLineStyle(lStyle[k/3]);
	  grC->SetLineWidth(2);
	  grC->SetMarkerSize(mSize[k/3]);

	  grSNR = new TGraphErrors(nVL,vBias,snr,0,0);
	  grSNR->SetMarkerStyle(mStyle[k/3]);
	  grSNR->SetMarkerColor(mCol[k/3]);
	  grSNR->SetLineColor(mCol[k/3]);
	  grSNR->SetLineStyle(lStyle[k/3]);
	  grSNR->SetLineWidth(2);
	  grSNR->SetMarkerSize(mSize[k/3]);
     
	  grRT = new TGraphErrors(arrS,vBiasReduced,rt,0,rtErr);
	  grRT->SetMarkerStyle(mStyle[k/3]);
	  grRT->SetMarkerColor(mCol[k/3]);
	  grRT->SetLineColor(mCol[k/3]);
	  grRT->SetLineStyle(lStyle[k/3]);
	  grRT->SetLineWidth(2);
	  grRT->SetMarkerSize(mSize[k/3]);

	  grJ = new TGraphErrors(arrS,vBiasReduced,jitter,0,0);
	  grJ->SetMarkerStyle(mStyle[k/3]);
	  grJ->SetMarkerColor(mCol[k/3]);
	  grJ->SetLineColor(mCol[k/3]);
	  grJ->SetLineStyle(lStyle[k/3]);
	  grJ->SetLineWidth(2);
	  grJ->SetMarkerSize(mSize[k/3]);

	  grG = new TGraphErrors(nVL,vBias,gain,0,gainErr);
	  grG->SetMarkerStyle(mStyle[k/3]);
	  grG->SetMarkerColor(mCol[k/3]);
	  grG->SetLineColor(mCol[k/3]);
	  grG->SetLineStyle(lStyle[k/3]);
	  grG->SetLineWidth(2);
	  grG->SetMarkerSize(mSize[k/3]);

	  if(k==0)
	    {
	      TGraphErrors *grT = new TGraphErrors(nVL,vBias,noiseSig,0,0);
	      grT->SetMarkerStyle(71);
	      grT->SetMarkerSize(2);
	      grT->SetMarkerColor(kBlack);
	      grT->SetLineColor(kBlack);
	      grT->SetLineStyle(7);
	      grT->SetLineWidth(4);

	      TGraphErrors *hrT = new TGraphErrors(nVL,vBias,noiseSig,0,0);
	      hrT->SetMarkerColor(kBlack);
	      hrT->SetLineColor(kBlack);
	      hrT->SetLineStyle(1);
	      hrT->SetLineWidth(4);

	      TGraphErrors *irT = new TGraphErrors(nVL,vBias,noiseSig,0,0);
	      irT->SetMarkerColor(kBlack);
	      irT->SetLineColor(kBlack);
	      irT->SetLineStyle(7);
	      irT->SetLineWidth(4);

	      legT->AddEntry(grT,"OM+CR", "pl");
	      legT->AddEntry(hrT,"FM+CR", "l");
	      legT->AddEntry(irT,"FM+CR+CG", "l");

	      legTJ->AddEntry(grT,"OM+CR", "pl");
	      legTJ->AddEntry(hrT,"FM+CR", "l");
	      legTJ->AddEntry(irT,"FM+CR+CG", "l");

	      canN->cd();
	      histJ[0]->Draw();
	      histJ[0]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[0]->GetYaxis()->SetTitle("Noise (mV)");
	      grN->Draw("eplsame");
	      legT->Draw();
	      TLine *lineN = new TLine(50, 2, 610, 2);
	      lineN->SetLineColor(kBlack);
	      lineN->SetLineWidth(4);
	      lineN->SetLineStyle(9);
	      lineN->Draw("lsame");
	  
	      canA->cd();
	      histJ[1]->Draw();
	      histJ[1]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[1]->GetYaxis()->SetTitle("Norm. amplitude (arb.)");
	      grA->Draw("eplsame");
	      legT->Draw();

	      canC->cd();
	      histJ[2]->Draw();
	      histJ[2]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[2]->GetYaxis()->SetTitle("Norm. charge (arb.)");
	      grC->Draw("eplsame");
	      legT->Draw();

	      canSNR->cd();
	      histJ[3]->Draw();
	      histJ[3]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[3]->GetYaxis()->SetTitle("SNR");
	      grSNR->Draw("eplsame");
	      legT->Draw();
	  
	      canRT->cd();
	      histJ[4]->Draw();
	      histJ[4]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[4]->GetYaxis()->SetTitle("Rise-Time (ns)");
	      grRT->Draw("eplsame");
	      legT->Draw();
	      TLine *lineRT[2];
	      for(Int_t i=0;i<3;++i)
		{
		  lineRT[i] = new TLine(50, 0.5*(i+1), 610, 0.5*(i+1));
		  lineRT[i]->SetLineColor(kBlack);
		  lineRT[i]->SetLineWidth(4);
		  lineRT[i]->SetLineStyle(9);
		  lineRT[i]->Draw("lsame");
		}
	  
	      canJ->cd();
	      histJ[5]->Draw();
	      histJ[5]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[5]->GetYaxis()->SetTitle("Jitter (ps)");
	      grJ->Draw("eplsame");
	      legTJ->Draw();
	      TLine *lineJ[2];
	      for(Int_t i=0;i<2;++i)
		{
		  lineJ[i] = new TLine(50, 50*(i+1), 610, 50*(i+1));
		  lineJ[i]->SetLineColor(kBlack);
		  lineJ[i]->SetLineWidth(4);
		  lineJ[i]->SetLineStyle(9);
		  lineJ[i]->Draw("lsame");
		}

	      canG->cd();
	      histJ[6]->Draw();
	      histJ[6]->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJ[6]->GetYaxis()->SetTitle("Gain");
	      grSNR->Draw("eplsame");
	      legT->Draw();
	    }
	  else
	    {
	      canN->cd();
	      grN->Draw("eplsame");
	      canA->cd();
	      grA->Draw("eplsame");
	      canC->cd();
	      grC->Draw("eplsame");
	      canSNR->cd();
	      grSNR->Draw("eplsame");
	      canRT->cd();
	      grRT->Draw("eplsame");
	      canJ->cd();
	      grJ->Draw("eplsame");
	      canG->cd();
	      grG->Draw("eplsame");

	    }
	  canN->cd();
	  histJ[0]->GetYaxis()->SetRangeUser(0,10);
	  histJ[0]->GetXaxis()->SetRangeUser(0,610);
	  if(type=="1")
	    legN->AddEntry(grN, legend, "epl");

	  canA->cd();
	  histJ[1]->GetYaxis()->SetRangeUser(0,1600);
	  histJ[1]->GetXaxis()->SetRangeUser(0,610);
	  if(type=="1")
	    legA->AddEntry(grA, legend, "epl");

	  canC->cd();
	  histJ[2]->GetYaxis()->SetRangeUser(0,300);
	  histJ[2]->GetXaxis()->SetRangeUser(0,610);
	  if(type=="1")
	    legC->AddEntry(grC, legend, "epl");

	  canSNR->cd();
	  histJ[3]->GetYaxis()->SetRangeUser(0,1000);
	  histJ[3]->GetXaxis()->SetRangeUser(0,610);
	  if(type=="1")
	    legSNR->AddEntry(grSNR, legend, "epl");

	  canRT->cd();
	  histJ[4]->GetYaxis()->SetRangeUser(0,4);
	  histJ[4]->GetXaxis()->SetRangeUser(50,610);
	  if(type=="1")
	    legRT->AddEntry(grRT, legend, "epl");

	  canJ->cd();
	  histJ[5]->GetYaxis()->SetRangeUser(0,500);
	  histJ[5]->GetXaxis()->SetRangeUser(50,610);
	  if(type=="1")
	    legJ->AddEntry(grJ, legend, "epl");

	  canG->cd();
	  histJ[6]->GetYaxis()->SetRangeUser(0,120);
	  histJ[6]->GetXaxis()->SetRangeUser(0,610);
	  if(type=="1")
	    legG->AddEntry(grJ, legend, "epl");

	}

      canN->cd();
      legN->Draw();
      canN->SaveAs(Form("../figures/spaceLgadsNoise_Pad%s.png", pad),"png");
  
      canA->cd();
      legA->Draw();
      canA->SaveAs(Form("../figures/spaceLgadsAmplitude_Pad%s.png", pad),"png");

      canC->cd();
      legC->Draw();
      canC->SaveAs(Form("../figures/spaceLgadsCharge_Pad%s.png", pad),"png");

      canSNR->cd();
      legSNR->Draw();
      canSNR->SaveAs(Form("../figures/spaceLgadsSNR_Pad%s.png", pad),"png");

      canRT->cd();
      legRT->Draw();
      canRT->SaveAs(Form("../figures/spaceLgadsRiseTime_Pad%s.png", pad),"png");
    
      canJ->cd();
      legJ->Draw();
      canJ->SaveAs(Form("../figures/spaceLgadsJitter_Pad%s.png", pad),"png");

      canG->cd();
      legG->Draw();
      canG->SaveAs(Form("../figures/spaceLgadsGain_Pad%s.png", pad),"png");
    }
      
  //==============================================PulseShape=======================================================
  // Plot PulseShape shape at 200 V for all the data
  if(processJitter && plotSignals)
    {
      TCanvas *canSignal = new TCanvas("canSignal","",1200,1000);
      TLegend *legSignal = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legSignal->SetHeader(area, "C");
      for(Int_t k=0; k<nFiles; ++k)
	{
	  //pulseShape[k]->SetMarkerStyle(mStyle[k]);
	  //pulseShape[k]->SetMarkerColor(mCol[k]);
	  //pulseShape[k]->SetMarkerSize(mSize[k]);
	  pulseShape[k]->SetLineStyle(lStyle[k]);
	  pulseShape[k]->SetLineColor(mCol[k]);
	  pulseShape[k]->SetLineWidth(2);
	  if(k==0)
	    {
	      pulseShape[k]->Draw("HIST PL");
	      pulseShape[k]->GetXaxis()->SetTitle("Time (ns)");
	      pulseShape[k]->GetYaxis()->SetTitle("Voltage (mV)");
	    }
	  else
	    pulseShape[k]->Draw("HIST PLSAME");
	  if(k%nWaf==0)
	    legSignal->AddEntry(pulseShape[k], Form("Wafer-%d", wafers[k/nWaf]), "l");
	  pulseShape[k]->GetYaxis()->SetRangeUser(0, 100);
	  pulseShape[k]->GetXaxis()->SetRangeUser(0, 10);
	}
      legSignal->Draw();
      canSignal->SaveAs(Form("../figures/spaceLgadsSignal_Pad%s.png", pad),"png");
    }

  //===========================================Jitter VS Gain=================================================
  // Plot Jitter as a function of Gain
  // Charge of PIN sensors as a function of thickness
  if(processJitter && plotJitterVsGain)
    {
      TCanvas *canJG = new TCanvas("canJG","",1200,1000);
      TLegend *legJG = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legJG->SetHeader(area, "C");

      TGraphErrors *grJG[nFiles];
      TH1F *histJG = new TH1F("histJG","",200, 0, 200);
      for(Int_t k=0; k<nFiles; ++k)
	{
	  grJG[k] = new TGraphErrors(arrSize[k], &Gain[k][0], &Jitter[k][0], &GainError[k][0], 0);
	  grJG[k]->SetMarkerStyle(mStyle[k]);
	  grJG[k]->SetMarkerColor(mCol[k]);
	  grJG[k]->SetLineStyle(lStyle[k]);
	  grJG[k]->SetLineColor(mCol[k]);
	  grJG[k]->SetMarkerSize(mSize[k]);
	  if(k==0)
	    {
	      histJG->Draw();
	      histJG->GetXaxis()->SetTitle("Gain");
	      histJG->GetYaxis()->SetTitle("Jitter (ps)");
	      grJG[k]->Draw("EPLSAME");

	      legTJ->Draw();
	    }
	  else
	    grJG[k]->Draw("EPLSAME");
	  if(k%nWaf==0)
	    legJG->AddEntry(grJG[k], Form("Wafer-%d", wafers[k/nWaf]), "epl");
	}
      histJG->GetYaxis()->SetRangeUser(0, 1000);
      histJG->GetXaxis()->SetRangeUser(0, 100);
      legJG->Draw();
      canJG->SaveAs(Form("../figures/spaceLgadsJitterVsGain_Pad%s.png", pad),"png");
    }
  //============================================Laser Stability================================================
  //Laser stability check
  //Draw a histogram of reference sensor amplitude
  if(processJitter && plotLaserStability)
    {
      TCanvas *canLas = new TCanvas("canLas","",1200,1000);
      TLegend *legLas = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legLas->SetHeader("Laser Intensity", "C");

      histRef->GetXaxis()->SetTitle("Reference amplitude (mV)");
      histRef->GetYaxis()->SetTitle("Counts");
      histRef->SetLineColor(kRed);
      histRef->SetLineWidth(2);
      histRef->Draw();
      histRef->GetXaxis()->SetRangeUser(1.1,1.25);
      Float_t mean = histRef->GetMean();
      Float_t height = histRef->GetBinContent(histRef->GetMaximumBin());
      Float_t stdDev = histRef->GetStdDev();
      TF1 *Gaus = new TF1("Gaus", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", 1.13, 1.2);
      Gaus->SetParameters(height, mean, stdDev);
      histRef->Fit("Gaus", "RMS+");
      Gaus->SetLineColor(kBlue);
      Gaus->SetLineWidth(2);
      Gaus->Draw("LSAME");
      TPaveStats *ptstats4 = new TPaveStats(0.6,0.75,0.9,0.951,"brNDC");
      ptstats4->SetBorderSize(1);
      ptstats4->SetFillColor(0);
      TText *ptstats4_LaTex = ptstats4->AddText("Laser stability");
      ptstats4_LaTex = ptstats4->AddText(Form("Events: %d",nEv));
      ptstats4_LaTex = ptstats4->AddText(Form("#chi^{2}/ndf: %2.0lf / %d", Gaus->GetChisquare(), Gaus->GetNDF()));
      ptstats4_LaTex = ptstats4->AddText(Form("Mean : %2.2e #pm %2.2e mV", Gaus->GetParameter(1),Gaus->GetParError(1)));
      ptstats4_LaTex = ptstats4->AddText(Form("Sigma: %2.2e #pm %2.2e mV", TMath::Abs(Gaus->GetParameter(2)),Gaus->GetParError(2)));
      ptstats4->Draw();
      canLas->SaveAs("../figures/spaceLgadsLaserStability.png","png");
    }
  //=============================================PIN Charge============================================================
  // Charge of PIN sensors as a function of thickness
  if(plotPinChargeThickness)
    {
      TCanvas *canQPin = new TCanvas("canQPin","",2400,2000);
      canQPin->Divide(2,2);

      // Noise Canvas
      canQPin->cd(1);
      TLegend *legNPin = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legNPin->SetHeader(Form("PIN (%s)", area), "C");
      TH1F *histNPin = new TH1F("histJPin", "", 700, 0, 700);
      TGraphErrors *grNPin[nWaf];
      for(Int_t k=0; k<nWaf; ++k)
	{
	  grNPin[k] = new TGraphErrors(arrSizePin[k], &vBiasPin[k][0], &noisePin[k][0], 0);
	  grNPin[k]->SetMarkerStyle(mStyle[3*k]);
	  grNPin[k]->SetMarkerColor(mCol[3*k]);
	  grNPin[k]->SetLineStyle(lStyle[3*k]);
	  grNPin[k]->SetLineColor(mCol[3*k]);
	  grNPin[k]->SetMarkerSize(mSize[3*k]);
	  if(k==0)
	    {
	      histNPin->Draw();
	      histNPin->GetXaxis()->SetTitle("V_{bias} (V)");
	      histNPin->GetYaxis()->SetTitle("Noise (mV)");
	      histNPin->GetYaxis()->SetRangeUser(0, 5);
	      histNPin->GetXaxis()->SetRangeUser(0, 400);
	      grNPin[k]->Draw("EPLSAME");
	    }
	  else
	    grNPin[k]->Draw("EPLSAME");
	  legNPin->AddEntry(grNPin[k], Form("Wafer-%d", wafers[k]), "epl");
	}
      legNPin->Draw();

      // Charge canvas
      canQPin->cd(2);
      TLegend *legQPin = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legQPin->SetHeader(Form("PIN (%s)", area), "C");
      TGraphErrors *grQPin = new TGraphErrors(3,thickness,ChargePin,0,ChargePinError);
      grQPin->SetMarkerStyle(21);
      grQPin->SetMarkerColor(kRed);
      grQPin->SetLineColor(kRed);
      grQPin->SetMarkerSize(1.5);
      grQPin->Draw("apl");
      grQPin->GetXaxis()->SetTitle("Thickness (#mum)");
      //grQPin->GetYaxis()->SetTitle("Norm. charge (arb.)");
      grQPin->GetYaxis()->SetTitle("Charge (fC)");
      //legQPin->AddEntry(grQPin, "Type 2", "epl");
      //legQPin->Draw();

      // Jitter Canvas
      canQPin->cd(3);
      TLegend *legJPin = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legJPin->SetHeader(Form("PIN (%s)", area), "C");
      TH1F *histJPin = new TH1F("histJPin", "", 700, 0, 700);
      TGraphErrors *grJPin[nWaf];
      for(Int_t k=0; k<nWaf; ++k)
	{
	  grJPin[k] = new TGraphErrors(arrSizePin[k], &vBiasPin[k][0], &JitterPin[k][0], &JitterPinError[k][0]);
	  grJPin[k]->SetMarkerStyle(mStyle[3*k]);
	  grJPin[k]->SetMarkerColor(mCol[3*k]);
	  grJPin[k]->SetLineStyle(lStyle[3*k]);
	  grJPin[k]->SetLineColor(mCol[3*k]);
	  grJPin[k]->SetMarkerSize(mSize[3*k]);
	  if(k==0)
	    {
	      histJPin->Draw();
	      histJPin->GetXaxis()->SetTitle("V_{bias} (V)");
	      histJPin->GetYaxis()->SetTitle("Jitter (ps)");
	      histJPin->GetYaxis()->SetRangeUser(0, 4000);
	      histJPin->GetXaxis()->SetRangeUser(0, 400);
	      grJPin[k]->Draw("EPLSAME");
	    }
	  else
	    grJPin[k]->Draw("EPLSAME");
	  legJPin->AddEntry(grJPin[k], Form("Wafer-%d", wafers[k]), "epl");
	}
      legJPin->Draw();

      // Rise Time Canvas
      canQPin->cd(4);
      TLegend *legRTPin = new TLegend(0.85,0.95,0.65,0.65,NULL,"brNDC");
      legRTPin->SetHeader(Form("PIN (%s)", area), "C");
      TH1F *histRTPin = new TH1F("histRTPin", "", 700, 0, 700);
      TGraphErrors *grRTPin[nWaf];
      for(Int_t k=0; k<nWaf; ++k)
	{
	  grRTPin[k] = new TGraphErrors(arrSizePin[k], &vBiasPin[k][0], &rtPin[k][0], &rtPinError[k][0]);
	  grRTPin[k]->SetMarkerStyle(mStyle[3*k]);
	  grRTPin[k]->SetMarkerColor(mCol[3*k]);
	  grRTPin[k]->SetLineStyle(lStyle[3*k]);
	  grRTPin[k]->SetLineColor(mCol[3*k]);
	  grRTPin[k]->SetMarkerSize(mSize[3*k]);
	  if(k==0)
	    {
	      histRTPin->Draw();
	      histRTPin->GetXaxis()->SetTitle("V_{bias} (V)");
	      histRTPin->GetYaxis()->SetTitle("Rise-Time (ps)");
	      histRTPin->GetYaxis()->SetRangeUser(0, 4000);
	      histRTPin->GetXaxis()->SetRangeUser(0, 400);
	      grRTPin[k]->Draw("EPLSAME");
	    }
	  else
	    grRTPin[k]->Draw("EPLSAME");
	  legRTPin->AddEntry(grRTPin[k], Form("Wafer-%d", wafers[k]), "epl");
	}
      legRTPin->Draw();
      canQPin->SaveAs(Form("../figures/spaceLgadsChargeVsThickness_Pad%s_PinSensors.png", pad),"png");
    }
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
