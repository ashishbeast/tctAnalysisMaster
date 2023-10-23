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
  TString files[8]={"../dataLeo/W_14_Type_3_Area_100_Hole_5_Av_256.tct", "../dataLeo/W_14_Type_3_Area_100_Hole_5.tct",
    "../dataLeo/W_14_Type_2_Area_100_Hole_5_Av_256.tct", "../dataLeo/W_14_Type_2_Area_100_Hole_5.tct",
    "../dataLeo/W_14_Type_1_Area_100_Av_256.tct", "../dataLeo/W_14_Type_1_Area_100.tct",
    "../dataLeo/W_14_Type_1_Area_25_Av_256.tct", "../dataLeo/W_14_Type_1_Area_25.tct"};
  
  TCanvas *can = new TCanvas("can","",1800,1200);
  can->Divide(2,2);
  TLegend *leg = new TLegend(0.6078464,0.625,0.1482471,0.914959,NULL,"brNDC");
  leg->SetHeader("Wafer 14","C");
  Int_t col[12] = {kRed,kBlue,kGreen,kBlack,13,6,7,8,9,28,34,49};
  Int_t ms[12] = {20,21,22,23,29,33,34,39,41,43,45,47};
    
  //gPad->SetLogy();
  TObjArray *info;
  TString wafer, pad, type, padType, legend, sensor;
    
  for(Int_t k=0;k<8;k+=2)
    {
      cout<<"Processing File: "<<files[k]<<endl;
      info = files[k].Tokenize("_");
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
        
      //slewRate File
      AnalyzeTCTData lgad(files[k]);
      lgad.CorrectBaseline();
      lgad.CalcNoise();
      lgad.CalculateWaveformProperties();
        
      //noise File
      AnalyzeTCTData noise(files[k+1]);
      noise.CorrectBaseline();
      noise.CalcNoise();
      noise.CalculateWaveformProperties();
        
      Int_t nV = lgad._nV1; //Number of Voltage points for LGAD
      Int_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
      Int_t arrS = nV;
        
      Float_t jitter[arrS];
      Float_t vBias[arrS];
    
      for(Int_t i=0; i<nV; i++)
        {
	  vBias[i] = i*step;
	  if(lgad._sigSlewRate[0][i] == 0)
	    jitter[i] = 0;
	  else
	    jitter[i] =  1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
        }
        
      TGraphErrors *grN = new TGraphErrors(arrS,vBias,lgad._sigNoise[0],0,0);
      TGraphErrors *grA = new TGraphErrors(arrS,vBias,lgad._sigAmplitude[0],0,0);
      TGraphErrors *grRT = new TGraphErrors(arrS,vBias,lgad._sigRiseTime[0],0,0);
      TGraphErrors *grJ = new TGraphErrors(arrS,vBias,jitter,0,0);

      Int_t mStyle, mCol;
      
      if(type == "0")
	mStyle = 25;
      else if(type =="1")
	mStyle = 21;
      else
	mStyle = 49;

      if(pad == "A")
	mCol = kRed;
      else if(pad =="B")
	mCol = kBlue;
      else
	mCol = kGreen+2;
            	  
      grN->SetMarkerStyle(mStyle);
      grN->SetMarkerColor(mCol);
      grN->SetLineColor(mCol);
      grN->SetMarkerSize(2);
      
      grA->SetMarkerStyle(mStyle);
      grA->SetMarkerColor(mCol);
      grA->SetLineColor(mCol);
      grA->SetMarkerSize(2);
      
      grRT->SetMarkerStyle(mStyle);
      grRT->SetMarkerColor(mCol);
      grRT->SetLineColor(mCol);
      grRT->SetMarkerSize(2);
      
      grJ->SetMarkerStyle(mStyle);
      grJ->SetMarkerColor(mCol);
      grJ->SetLineColor(mCol);
      grJ->SetMarkerSize(2);
      
      if(k==0)
        {
	  can->cd(1);
	  grN->Draw("ap");
	  grN->GetXaxis()->SetTitle("V_{bias} (V)");
	  grN->GetYaxis()->SetTitle("Noise (mV)");

	  can->cd(2);
	  grA->Draw("ap");
	  grA->GetXaxis()->SetTitle("V_{bias} (V)");
	  grA->GetYaxis()->SetTitle("Amplitude (mV)");

	  can->cd(3);
	  grRT->Draw("ap");
	  grRT->GetXaxis()->SetTitle("V_{bias} (V)");
	  grRT->GetYaxis()->SetTitle("Rise-Time (ns)");

	  can->cd(4);
	  grJ->Draw("ap");
	  grJ->GetXaxis()->SetTitle("V_{bias} (V)");
	  grJ->GetYaxis()->SetTitle("Jitter (ps)");
	  
        }
      else
	{
	  can->cd(1);
	  grN->Draw("psame");
	  can->cd(2);
	  grA->Draw("psame");
	  can->cd(3);
	  grRT->Draw("psame");
	  can->cd(4);
	  grJ->Draw("psame");
	  
	}

      grN->GetYaxis()->SetRangeUser(0,10);
      grN->GetXaxis()->SetRangeUser(300,610);

      grA->GetYaxis()->SetRangeUser(0,130);
      grA->GetXaxis()->SetRangeUser(300,610);

      grRT->GetYaxis()->SetRangeUser(1,4);
      grRT->GetXaxis()->SetRangeUser(300,610);
      
      grJ->GetYaxis()->SetRangeUser(0,300);
      grJ->GetXaxis()->SetRangeUser(300,610);

      leg->AddEntry(grJ, legend, "p");
    }
  for(Int_t i=1; i<=4;++i)
    {
      can->cd(i);  
      leg->Draw();
    }
  can->SaveAs("../figures/probeW14.png","png");
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
