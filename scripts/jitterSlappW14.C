//This macro calculates gain for multiples files
//Gain is defined as the ratio of charge collected by LGAD to charge collected by PIN sensor

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
  
  TCanvas *can = new TCanvas("can","",1200,1000);
  TLegend *leg = new TLegend(0.6078464,0.625,0.8482471,0.914959,NULL,"brNDC");
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

      legend = "Type-"+ type;
        
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
        
      TGraphErrors *gr = new TGraphErrors(arrS,vBias,jitter,0,0);
        
      if(type == "1")
	gr->SetMarkerStyle(25);
      else if(type =="2")
	gr->SetMarkerStyle(21);
      else
	gr->SetMarkerStyle(49);
        
      gr->SetMarkerColor(col[k/2]);
      gr->SetLineColor(col[k/2]);
      gr->SetMarkerSize(2);
        
      if(k==0)
        {
	  gr->Draw("ap");
	  gr->GetXaxis()->SetTitle("V_{bias} (V)");
	  gr->GetYaxis()->SetTitle("Jitter (ps)");
        }
      else
	gr->Draw("psame");
      gr->GetYaxis()->SetRangeUser(0,500);
      gr->GetXaxis()->SetRangeUser(300,500);
      leg->AddEntry(gr, legend, "p");
    }
    
  leg->Draw();
  can->SaveAs("../figures/jitterSLAPPW14.png","png");
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
