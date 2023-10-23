#include "AnalyzeTCTData.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TObjString.h"

using namespace std;
void SetStyle(Bool_t threeD = false);

int main()
{
  SetStyle(false);

  TCanvas *can = new TCanvas("can","",1800,1200);
  can->Divide(3,2);  
  TLegend *leg11 = new TLegend(0.45,0.9,0.15,0.714959,NULL,"brNDC");
  leg11->SetHeader("Wafer 11","C");
  TLegend *leg12 = new TLegend(0.65,0.9,0.45,0.714959,NULL,"brNDC");
  leg12->SetHeader("Wafer 12","C");
  TLegend *leg14 = new TLegend(0.85,0.9,0.75,0.714959,NULL,"brNDC");
  leg14->SetHeader("Wafer 14","C");
  TLegend *leg9 = new TLegend(0.85,0.5,0.45,0.314959,NULL,"brNDC");
  leg9->SetHeader("Wafer 9","C");
  TLegend *leg1 = new TLegend(0.45,0.5,0.15,0.314959,NULL,"brNDC");
  leg1->SetHeader("Wafer 1","C");

  TLegend *legG11 = new TLegend(0.6078464,0.625,0.2482471,0.914959,NULL,"brNDC");
  legG11->SetHeader("Wafer 11","C");

  TLegend *legG12 = new TLegend(0.6078464,0.625,0.2482471,0.914959,NULL,"brNDC");
  legG12->SetHeader("Wafer 12","C");

  TLegend *legG9 = new TLegend(0.6078464,0.625,0.2482471,0.914959,NULL,"brNDC");
  legG9->SetHeader("Wafer 9","C");

  TLegend *legG1 = new TLegend(0.6078464,0.625,0.2482471,0.914959,NULL,"brNDC");
  legG1->SetHeader("Wafer 1","C");

  Int_t col[12] = {kRed,kBlue,kGreen,kBlack,13,6,7,8,9,28,34,49};
  Int_t ms[12] = {20,21,22,23,29,33,34,39,41,43,45,47};
    
  //gPad->SetLogy();
  TObjArray *info;
  TString wafer, pad, type, padType, legend, sensor;
 
  //Wafer 12
  TString filesW12[2]={"../w12Data/w_12_type_2_area_25_av_256", "../w12Data/w_12_type_2_area_25"};
  
  for(Int_t k=0;k<2;k+=2)
    {
      cout<<"Processing File: "<<filesW12[k]<<endl;
      info = filesW12[k].Tokenize("_");
      wafer = ((TObjString*)(info->At(1)))->String();
      type  = ((TObjString*)(info->At(3)))->String();
      pad   = ((TObjString*)(info->At(5)))->String();

      if(type == "1")
	type = "0";
      else if(type =="2")
	type = "1";
      else
	type = "2";
      
      if(pad == "6")
	pad = "A";
      else if(pad =="25")
	pad = "B";
      else
	pad = "C";
      
      legend = pad+", Type-"+ type;
        
      //slewRate File
      AnalyzeTCTData lgad(filesW12[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      //noise File
      AnalyzeTCTData noise(filesW12[k+1]);
      noise.CorrectBaseline();
      noise.CalcNoise();
      noise.CalculateWaveformProperties();
        
      Int_t nV = lgad._nV1; //Number of Voltage points for LGAD
      Int_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
      //Int_t arrS = nV;
      Int_t arrS = (nV-20)/10 + 1;
        
      Float_t jitter[arrS];
      Float_t vBias[arrS];
      Float_t noiseSig[arrS];
      Float_t amplitude[arrS];
      Float_t rt[arrS];
    
      for(Int_t i=20; i<nV; i+=10)
        {
	  vBias[(i-20)/10] = i*step;
	  noiseSig[(i-20)/10] = noise._sigNoise[0][i];
	  amplitude[(i-20)/10] = lgad._sigAmplitude[0][i];
	  rt[(i-20)/10] = lgad._sigRiseTime[0][i];
	  
	  if(lgad._sigSlewRate[0][i] == 0)
	    jitter[(i-20)/10] = 0;
	  else
	    jitter[(i-20)/10] =  1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
        }

      TGraphErrors *grN = new TGraphErrors(arrS,vBias,noiseSig,0,0);
      TGraphErrors *grA = new TGraphErrors(arrS,vBias,amplitude,0,0);
      TGraphErrors *grRT = new TGraphErrors(arrS,vBias,rt,0,0);
      TGraphErrors *grJ = new TGraphErrors(arrS,vBias,jitter,0,0);

      Int_t mStyle, mCol=kGreen+2;
      
      if(type == "0")
	mStyle = 25;
      else if(type =="1")
	mStyle = 21;
      else
	mStyle = 49;
            	  
      grN->SetMarkerStyle(mStyle);
      grN->SetMarkerColor(mCol);
      grN->SetLineColor(mCol);
      grN->SetLineStyle(4);
      grN->SetMarkerSize(1.5);
      
      grA->SetMarkerStyle(mStyle);
      grA->SetMarkerColor(mCol);
      grA->SetLineColor(mCol);
      grA->SetLineStyle(4);
      grA->SetMarkerSize(1.5);
      
      grRT->SetMarkerStyle(mStyle);
      grRT->SetMarkerColor(mCol);
      grRT->SetLineColor(mCol);
      grRT->SetLineStyle(4);
      grRT->SetMarkerSize(1.5);
      
      grJ->SetMarkerStyle(mStyle);
      grJ->SetMarkerColor(mCol);
      grJ->SetLineColor(mCol);
      grJ->SetLineStyle(4);
      grJ->SetMarkerSize(1.5);

      if(k==0)
        {
	  can->cd(1);
	  grN->Draw("apl");
	  grN->GetXaxis()->SetTitle("V_{bias} (V)");
	  grN->GetYaxis()->SetTitle("Noise (mV)");

	  can->cd(2);
	  grA->Draw("apl");
	  grA->GetXaxis()->SetTitle("V_{bias} (V)");
	  grA->GetYaxis()->SetTitle("Amplitude (mV)");

	  can->cd(3);
	  grRT->Draw("apl");
	  grRT->GetXaxis()->SetTitle("V_{bias} (V)");
	  grRT->GetYaxis()->SetTitle("Rise-Time (ns)");

	  can->cd(4);
	  grJ->Draw("apl");
	  grJ->GetXaxis()->SetTitle("V_{bias} (V)");
	  grJ->GetYaxis()->SetTitle("Jitter (ps)");
	  
        }
      else
	{
	  can->cd(1);
	  grN->Draw("plsame");
	  can->cd(2);
	  grA->Draw("plsame");
	  can->cd(3);
	  grRT->Draw("plsame");
	  can->cd(4);
	  grJ->Draw("plsame");
	  
	}

      grN->GetYaxis()->SetRangeUser(0,5);
      grN->GetXaxis()->SetRangeUser(0,610);

      grA->GetYaxis()->SetRangeUser(0,80);
      grA->GetXaxis()->SetRangeUser(0,610);

      grRT->GetYaxis()->SetRangeUser(1,3);
      grRT->GetXaxis()->SetRangeUser(0,610);
      
      grJ->GetYaxis()->SetRangeUser(0,500);
      grJ->GetXaxis()->SetRangeUser(0,610);

      leg12->AddEntry(grJ, legend, "p");
    }

  //Wafer 9
  TString filesW9[2]={"../w9Data/w_9_type_2_area_25_av_256", "../w9Data/w_9_type_2_area_25"};
  
  for(Int_t k=0;k<2;k+=2)
    {
      cout<<"Processing File: "<<filesW9[k]<<endl;
      info = filesW12[k].Tokenize("_");
      wafer = ((TObjString*)(info->At(1)))->String();
      type  = ((TObjString*)(info->At(3)))->String();
      pad   = ((TObjString*)(info->At(5)))->String();

      if(type == "1")
	type = "0";
      else if(type =="2")
	type = "1";
      else
	type = "2";
      
      if(pad == "6")
	pad = "A";
      else if(pad =="25")
	pad = "B";
      else
	pad = "C";
      
      legend = pad+", Type-"+ type;
        
      //slewRate File
      AnalyzeTCTData lgad(filesW9[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      //noise File
      AnalyzeTCTData noise(filesW9[k+1]);
      noise.CorrectBaseline();
      noise.CalcNoise();
      noise.CalculateWaveformProperties();
        
      Int_t nV = lgad._nV1; //Number of Voltage points for LGAD
      Int_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
      //Int_t arrS = nV;
      Int_t arrS = (nV-20)/10 + 1;

      Float_t jitter[arrS];
      Float_t vBias[arrS];
      Float_t noiseSig[arrS];
      Float_t amplitude[arrS];
      Float_t rt[arrS];
    
      for(Int_t i=20; i<nV; i+=10)
        {
	  vBias[(i-20)/10] = i*step;
	  noiseSig[(i-20)/10] = noise._sigNoise[0][i];
	  amplitude[(i-20)/10] = lgad._sigAmplitude[0][i];
	  rt[(i-20)/10] = lgad._sigRiseTime[0][i];
	  
	  if(lgad._sigSlewRate[0][i] == 0)
	    jitter[(i-20)/10] = 0;
	  else
	    jitter[(i-20)/10] =  1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
        }

      TGraphErrors *grN = new TGraphErrors(arrS,vBias,noiseSig,0,0);
      TGraphErrors *grA = new TGraphErrors(arrS,vBias,amplitude,0,0);
      TGraphErrors *grRT = new TGraphErrors(arrS,vBias,rt,0,0);
      TGraphErrors *grJ = new TGraphErrors(arrS,vBias,jitter,0,0);

      Int_t mStyle, mCol=kMagenta;
      
      if(type == "0")
	mStyle = 25;
      else if(type =="1")
	mStyle = 21;
      else
	mStyle = 49;
            	  
      grN->SetMarkerStyle(mStyle);
      grN->SetMarkerColor(mCol);
      grN->SetLineColor(mCol);
      grN->SetLineStyle(4);
      grN->SetMarkerSize(1.5);
      
      grA->SetMarkerStyle(mStyle);
      grA->SetMarkerColor(mCol);
      grA->SetLineColor(mCol);
      grA->SetLineStyle(4);
      grA->SetMarkerSize(1.5);
      
      grRT->SetMarkerStyle(mStyle);
      grRT->SetMarkerColor(mCol);
      grRT->SetLineColor(mCol);
      grRT->SetLineStyle(4);
      grRT->SetMarkerSize(1.5);
      
      grJ->SetMarkerStyle(mStyle);
      grJ->SetMarkerColor(mCol);
      grJ->SetLineColor(mCol);
      grJ->SetLineStyle(4);
      grJ->SetMarkerSize(1.5);
      
      can->cd(1);
      grN->Draw("plsame");
      can->cd(2);
      grA->Draw("plsame");
      can->cd(3);
      grRT->Draw("plsame");
      can->cd(4);
      grJ->Draw("plsame");

      leg9->AddEntry(grJ, legend, "p");
    }

  //Wafer 1
  TString filesW1[2]={"../w1Data/w_1_type_2_area_25_gain_vsweep_avg", "../w1Data/w_1_type_2_area_25_gain_vsweep"};
    
  for(Int_t k=0;k<2;k+=2)
    {
      cout<<"Processing File: "<<filesW1[k]<<endl;
      info = filesW12[k].Tokenize("_");
      wafer = ((TObjString*)(info->At(1)))->String();
      type  = ((TObjString*)(info->At(3)))->String();
      pad   = ((TObjString*)(info->At(5)))->String();

      if(type == "1")
	type = "0";
      else if(type =="2")
	type = "1";
      else
	type = "2";
      
      if(pad == "6")
	pad = "A";
      else if(pad =="25")
	pad = "B";
      else
	pad = "C";
      
      legend = pad+", Type-"+ type;
        
      //slewRate File
      AnalyzeTCTData lgad(filesW1[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      //noise File
      AnalyzeTCTData noise(filesW1[k+1]);
      noise.CorrectBaseline();
      noise.CalcNoise();
      noise.CalculateWaveformProperties();
        
      Int_t nV = lgad._nV1; //Number of Voltage points for LGAD
      Int_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
      //Int_t arrS = nV;
      Int_t arrS = (nV-20)/10 + 1;

      Float_t jitter[arrS];
      Float_t vBias[arrS];
      Float_t noiseSig[arrS];
      Float_t amplitude[arrS];
      Float_t rt[arrS];
    
      for(Int_t i=20; i<nV; i+=10)
        {
	  vBias[(i-20)/10] = i*step;
	  noiseSig[(i-20)/10] = noise._sigNoise[0][i];
	  amplitude[(i-20)/10] = lgad._sigAmplitude[0][i];
	  rt[(i-20)/10] = lgad._sigRiseTime[0][i];
	  
	  if(lgad._sigSlewRate[0][i] == 0)
	    jitter[(i-20)/10] = 0;
	  else
	    jitter[(i-20)/10] =  1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
        }

      TGraphErrors *grN = new TGraphErrors(arrS,vBias,noiseSig,0,0);
      TGraphErrors *grA = new TGraphErrors(arrS,vBias,amplitude,0,0);
      TGraphErrors *grRT = new TGraphErrors(arrS,vBias,rt,0,0);
      TGraphErrors *grJ = new TGraphErrors(arrS,vBias,jitter,0,0);

      Int_t mStyle, mCol=kBlack;
      
      if(type == "0")
	mStyle = 25;
      else if(type =="1")
	mStyle = 21;
      else
	mStyle = 49;
            	  
      grN->SetMarkerStyle(mStyle);
      grN->SetMarkerColor(mCol);
      grN->SetLineColor(mCol);
      grN->SetLineStyle(4);
      grN->SetMarkerSize(1.5);
      
      grA->SetMarkerStyle(mStyle);
      grA->SetMarkerColor(mCol);
      grA->SetLineColor(mCol);
      grA->SetLineStyle(4);
      grA->SetMarkerSize(1.5);
      
      grRT->SetMarkerStyle(mStyle);
      grRT->SetMarkerColor(mCol);
      grRT->SetLineColor(mCol);
      grRT->SetLineStyle(4);
      grRT->SetMarkerSize(1.5);
      
      grJ->SetMarkerStyle(mStyle);
      grJ->SetMarkerColor(mCol);
      grJ->SetLineColor(mCol);
      grJ->SetLineStyle(4);
      grJ->SetMarkerSize(1.5);
      
      can->cd(1);
      grN->Draw("plsame");
      can->cd(2);
      grA->Draw("plsame");
      can->cd(3);
      grRT->Draw("plsame");
      can->cd(4);
      grJ->Draw("plsame");

      leg1->AddEntry(grJ, legend, "p");
    }

  can->cd(6);
  leg1->Draw();
  leg9->Draw();
  leg12->Draw();

  /*
  //================================================Gain=======================================================

  //Wafer 11
  TString gainFilesW11[4]={"../data/slapp_1MIP_W_11_pad_B_type_1_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_11_pad_B_type_0_gain_0_VS_IR_a.tct",
    "../data/slapp_1MIP_W_11_pad_B_type_2_gain_1_VS_IR_a.tct", "../data/slapp_1MIP_W_11_pad_B_type_0_gain_0_VS_IR_a.tct"};
    
  for(Int_t k=0;k<4;k+=2)
    {
      info = gainFilesW11[k].Tokenize("_");
      wafer = ((TObjString*)(info->At(3)))->String();
      pad   = ((TObjString*)(info->At(5)))->String();
      type  = ((TObjString*)(info->At(7)))->String();
      sensor= ((TObjString*)(info->At(9)))->String();
        
      legend = pad+", Type-"+ type;
                
      cout<<"Processing File: "<<gainFilesW11[k]<<endl;
      //Read LGAD data
      AnalyzeTCTData lgad(gainFilesW11[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      cout<<"Processing File: "<<gainFilesW11[k+1]<<endl;
      //Read PIN data
      AnalyzeTCTData pin(gainFilesW11[k+1], 1.5);
      pin.CorrectBaseline();
      pin.CalcNoise();
      pin.CalculateWaveformProperties();
        
      Int_t nVL = lgad._nV1; //Number of Voltage points for LGAD
      Int_t nVP = pin._nV1; //Number of Voltage points for PIN
        
      vector<Float_t> chargeL(nVL);
      vector<Float_t> errL(nVL);
      vector<Float_t> chargeP(nVP);
      vector<Float_t> errP(nVP);

      Int_t arrS = (nVL-10)/5 + 1;
      //cout<<"==============================================>"<<arrS<<endl;
      Float_t voltage[arrS];
      Float_t gain[arrS];
      Float_t errG[arrS];
        
      for(Int_t i=0; i< nVP; ++i)
        {
	  chargeP[i] = pin._sigNormCharge[0][i];
	  errP[i] = pin._sigChargeError[0][i];
	}
        
      //Average Charge in Pin
      Float_t avgQ = TMath::Mean(chargeP.begin()+5,chargeP.end());
      //Error in Average Charge
      Float_t errAvg = TMath::StdDev(chargeP.begin()+5,chargeP.end())/TMath::Sqrt(nVP-5);
        
      Float_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
        
      for(Int_t i=10; i<nVL; i+=5)
        {
	  voltage[(i-10)/5] = i*step;
	  chargeL[i] = lgad._sigNormCharge[0][i];
	  errL[i] = lgad._sigChargeError[0][i];
	  
	  gain[(i-10)/5] = lgad._sigNormCharge[0][i]/avgQ;
	  errG[(i-10)/5] = gain[(i-10)/5]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
        }
        
      TGraphErrors *gr = new TGraphErrors(arrS,voltage,gain,0,errG);
        
      if(type == "0")
	gr->SetMarkerStyle(25);
      else if(type =="1")
	gr->SetMarkerStyle(21);
      else
	gr->SetMarkerStyle(49);

      gr->SetMarkerColor(kRed);
      gr->SetLineColor(kRed);
      gr->SetLineStyle(4);
      gr->SetMarkerSize(1.5);
        
      if(k==0)
        {
	  can->cd(5);
	  gr->Draw("apl");
	  gr->GetXaxis()->SetTitle("V_{bias} (V)");
	  gr->GetYaxis()->SetTitle("Gain");
        }
      else
	{
	  can->cd(5);
	  gr->Draw("plsame");
	}
      gr->GetYaxis()->SetRangeUser(0,100);
      legG11->AddEntry(gr, legend, "p");
    }
  
  //Wafer 12
  TString gainFilesW12[4]={"../w12Data/w_12_type_2_area_25_av_256", "../w12Data/w_12_type_2_area_25_NoGain_av_256",
    "../w12Data/w_12_type_3_area_25_av_256", "../w12Data/w_12_type_3_area_25_NoGain_av_256"};
      
  for(Int_t k=0;k<4;k+=2)
    {
      info = gainFilesW12[k].Tokenize("_");
      wafer = ((TObjString*)(info->At(1)))->String();
      type  = ((TObjString*)(info->At(3)))->String();
      pad   = ((TObjString*)(info->At(5)))->String();
	
      if(type == "1")
	type = "0";
      else if(type =="2")
	type = "1";
      else
	type = "2";
      
      if(pad == "6")
	pad = "A";
      else if(pad =="25")
	pad = "B";
      else
	pad = "C";

      legend = pad+", Type-"+ type;
        
      //Read data
        
      cout<<"Processing File: "<<gainFilesW12[k]<<endl;
      //Read LGAD data
      AnalyzeTCTData lgad(gainFilesW12[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      cout<<"Processing File: "<<gainFilesW12[k+1]<<endl;
      //Read PIN data
      AnalyzeTCTData pin(gainFilesW12[k+1], 1.5);
      pin.CorrectBaseline();
      pin.CalcNoise();
      pin.CalculateWaveformProperties();
        
      Int_t nVL = lgad._nV1; //Number of Voltage points for LGAD
      Int_t nVP = pin._nV1; //Number of Voltage points for PIN

      Int_t arrS = (nVL-10)/5 + 1;
      vector<Float_t> chargeL(nVL);
      vector<Float_t> errL(nVL);
      vector<Float_t> chargeP(nVP);
      vector<Float_t> errP(nVP);
        
      Float_t voltage[arrS];
      Float_t gain[arrS];
      Float_t errG[arrS];
        
      for(Int_t i=0; i< nVP; ++i)
        {
	  chargeP[i] = pin._sigNormCharge[0][i];
	  errP[i] = pin._sigChargeError[0][i];
        }
        
      //Average Charge in Pin
      Float_t avgQ = TMath::Mean(chargeP.begin()+5,chargeP.end());
      //Error in Average Charge
      Float_t errAvg = TMath::StdDev(chargeP.begin()+5,chargeP.end())/TMath::Sqrt(nVP-5);
        
      Float_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
        
      for(Int_t i=10; i<nVL; i+=5)
        {
	  voltage[(i-10)/5] = i*step;
	  chargeL[i] = lgad._sigNormCharge[0][i];
	  errL[i] = lgad._sigChargeError[0][i];
            
	  gain[(i-10)/5] = chargeL[i]/avgQ;
	  errG[(i-10)/5] = gain[(i-10)/5]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
        }
        
      TGraphErrors *gr = new TGraphErrors(arrS,voltage,gain,0,errG);

      if(type == "0")
	gr->SetMarkerStyle(25);
      else if(type =="1")
	gr->SetMarkerStyle(21);
      else
	gr->SetMarkerStyle(49);
	
      gr->SetMarkerColor(kGreen+2);
      gr->SetLineColor(kGreen+2);
      gr->SetLineStyle(4);
      gr->SetMarkerSize(1.5);
	
      can->cd(5);
      gr->Draw("plsame");
      legG12->AddEntry(gr, legend, "p");
    }

  //Wafer 9
  TString gainFilesW9[4]={"../w9Data/w_9_type_2_area_25_av_256", "../w9Data/w_9_type_2_area_25_NoGain_av_256",
    "../w9Data/w_9_type_3_area_25_av_256", "../w9Data/w_9_type_3_area_25_NoGain_av_256"};
      
  for(Int_t k=0;k<4;k+=2)
    {
      info = gainFilesW12[k].Tokenize("_");
      wafer = ((TObjString*)(info->At(1)))->String();
      type  = ((TObjString*)(info->At(3)))->String();
      pad   = ((TObjString*)(info->At(5)))->String();
	
      if(type == "1")
	type = "0";
      else if(type =="2")
	type = "1";
      else
	type = "2";
      
      if(pad == "6")
	pad = "A";
      else if(pad =="25")
	pad = "B";
      else
	pad = "C";

      legend = pad+", Type-"+ type;
        
      //Read data
        
      cout<<"Processing File: "<<gainFilesW9[k]<<endl;
      //Read LGAD data
      AnalyzeTCTData lgad(gainFilesW9[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      cout<<"Processing File: "<<gainFilesW9[k+1]<<endl;
      //Read PIN data
      AnalyzeTCTData pin(gainFilesW9[k+1], 1.5);
      pin.CorrectBaseline();
      pin.CalcNoise();
      pin.CalculateWaveformProperties();
        
      Int_t nVL = lgad._nV1; //Number of Voltage points for LGAD
      Int_t nVP = pin._nV1; //Number of Voltage points for PIN
        
      vector<Float_t> chargeL(nVL);
      vector<Float_t> errL(nVL);
      vector<Float_t> chargeP(nVP);
      vector<Float_t> errP(nVP);

      Int_t arrS = (nVL-10)/5 + 1;
      Float_t voltage[nVL];
      Float_t gain[nVL];
      Float_t errG[nVL];
        
      for(Int_t i=0; i< nVP; ++i)
        {
	  chargeP[i] = pin._sigNormCharge[0][i];
	  errP[i] = pin._sigChargeError[0][i];
        }
        
      //Average Charge in Pin
      Float_t avgQ = TMath::Mean(chargeP.begin()+5,chargeP.end());
      //Error in Average Charge
      Float_t errAvg = TMath::StdDev(chargeP.begin()+5,chargeP.end())/TMath::Sqrt(nVP-5);
        
      Float_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
        
      for(Int_t i=10; i<nVL; i+=5)
        {
	  voltage[(i-10)/5] = i*step;
	  chargeL[i] = lgad._sigNormCharge[0][i];
	  errL[i] = lgad._sigChargeError[0][i];
            
	  gain[(i-10)/5] = chargeL[i]/avgQ;
	  errG[(i-10)/5] = gain[(i-10)/5]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
        }
        
      TGraphErrors *gr = new TGraphErrors(arrS,voltage,gain,0,errG);

      if(type == "0")
	gr->SetMarkerStyle(25);
      else if(type =="1")
	gr->SetMarkerStyle(21);
      else
	gr->SetMarkerStyle(49);
	
      gr->SetMarkerColor(kMagenta);
      gr->SetLineColor(kMagenta);
      gr->SetLineStyle(4);
      gr->SetMarkerSize(1.5);
	
      can->cd(5);
      gr->Draw("plsame");
      legG9->AddEntry(gr, legend, "p");
    }
  */
  can->SaveAs("../figures/compareSlapWafersPadB.png","png");
  return 0;
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
