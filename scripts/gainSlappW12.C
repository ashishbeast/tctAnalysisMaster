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

    TString files[8]={"../w12Data/w_12_type_1_area_6_av_256", "../w12Data/w_12_type_1_area_6_NoGain_av_256",
      "../w12Data/w_12_type_1_area_100_av_256", "../w12Data/w_12_type_1_area_100_NoGain_av_256",
      "../w12Data/w_12_type_2_area_6_av_256", "../w12Data/w_12_type_2_area_6_NoGain_av_256",
      "../w12Data/w_12_type_3_area_6_av_256", "../w12Data/w_12_type_3_area_6_NoGain_av_256"};
      
    TCanvas *can = new TCanvas("can","",1200,1000);
    //gPad->SetLogy();
    TLegend *leg = new TLegend(0.6078464,0.625,0.2482471,0.914959,NULL,"brNDC");
    leg->SetHeader("Wafer 12","C");
    Int_t col[12] = {kRed,kBlue,kGreen,kBlack,13,6,7,8,9,28,34,49};
    Int_t ms[12] = {20,21,22,23,29,33,34,39,41,43,45,47};
    
    //gPad->SetLogy();
    TObjArray *info;
    TString wafer, pad, type, padType, legend, sensor;
    
    for(Int_t k=0;k<8;k+=2)
    {
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
        
        cout<<"Processing File: "<<files[k]<<endl;
        //Read LGAD data
        AnalyzeTCTData lgad(files[k]);
        lgad.CorrectBaseline();
        lgad.CalcNoise();
        lgad.CalculateWaveformProperties();
        
        cout<<"Processing File: "<<files[k+1]<<endl;
        //Read PIN data
        AnalyzeTCTData pin(files[k+1], 1.5);
        pin.CorrectBaseline();
        pin.CalcNoise();
        pin.CalculateWaveformProperties();
        
        Int_t nVL = lgad._nV1; //Number of Voltage points for LGAD
        Int_t nVP = pin._nV1; //Number of Voltage points for PIN
        
        vector<Float_t> chargeL(nVL);
        vector<Float_t> errL(nVL);
        vector<Float_t> chargeP(nVP);
        vector<Float_t> errP(nVP);
        
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
        
        for(Int_t i=0; i<nVL; ++i)
        {
            voltage[i] = i*step;
            //chargeL[i] = lgad._sigCharge[0][i];
            //errL[i] = lgad._sigChargeError[0][i];
            
            chargeL[i] = lgad._sigNormCharge[0][i];
            errL[i] = lgad._sigChargeError[0][i];
            
            gain[i] = chargeL[i]/avgQ;
            errG[i] = gain[i]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
        }
        
        TGraphErrors *gr = new TGraphErrors(nVL,voltage,gain,0,errG);

        if(type == "0")
	  gr->SetMarkerStyle(25);
	else if(type =="1")
	  gr->SetMarkerStyle(21);
	else
	  gr->SetMarkerStyle(49);

	if(pad == "A")
	  {
            gr->SetMarkerColor(kRed);
            gr->SetLineColor(kRed);
	  }
        else if(pad == "B")
	  {
            gr->SetMarkerColor(kBlue);
            gr->SetLineColor(kBlue);
	  }
        else
	  {
            gr->SetMarkerColor(kGreen+2);
            gr->SetLineColor(kGreen+2);
	  }

	gr->SetMarkerSize(2);
      
	if(k==0)
        {
            gr->Draw("apl");
            gr->GetXaxis()->SetTitle("V_{bias} (V)");
            gr->GetYaxis()->SetTitle("Gain");
        }
        else
            gr->Draw("plsame");
        gr->GetYaxis()->SetRangeUser(0,30);
	gr->GetXaxis()->SetRangeUser(0,400);
	
        leg->AddEntry(gr, legend, "p");
    }
    
    leg->Draw();
    can->SaveAs("../figures/gainSLAPPW12.png","png");
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
