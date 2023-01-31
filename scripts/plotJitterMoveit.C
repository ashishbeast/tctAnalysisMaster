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
    TString files[3]={"../data/MoVeIT_W13L_A_VS_IR_a.tct", "../data/MoVeIT_W13L_B_VS_IR_a.tct",
        "../data/MoVeIT_W13L_AreaC_VS_IR_d.tct"};
        
    TCanvas *canvas = new TCanvas("canvas","",1200,1000);
    TLegend *leg = new TLegend(0.4078464,0.625,0.4482471,0.914959,NULL,"brNDC");
    leg->SetHeader("Space LGAD","C");
    Int_t col[12] = {kRed,kBlue,kGreen,kBlack,13,6,7,8,9,28,34,49};
    Int_t ms[12] = {20,21,22,23,29,33,34,39,41,43,45,47};
    
    //gPad->SetLogy();
    TObjArray *info;
    TString wafer, pad, type, padType, legend, sensor;
    
    for(Int_t k=0;k<3;++k)
    {
        cout<<"Processing File: "<<files[k]<<endl;
        info = files[k].Tokenize("_");
        wafer = ((TObjString*)(info->At(1)))->String();
        pad   = ((TObjString*)(info->At(2)))->String();
        
        legend = "Wafer-"+wafer+", pad-"+pad;
        
        //Read data
        AnalyzeTCTData lgad(files[k]);
        lgad.CorrectBaseline();
        lgad.CalcNoise();
        lgad.CalculateSignalProperties();
        
        Int_t nV = lgad._nV1; //Number of Voltage points for LGAD
        Int_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
        //Int_t arrS = ((nV * step)/10) - 12;
        Int_t arrS = nV-20;
        Float_t jitter[arrS];
        Float_t err[arrS];
        Float_t voltage[arrS];
    
        for(Int_t i=20; i<nV; ++i)
        {
            voltage[i-20] = i*step;
            //jitter[i-20] = 1000*(lgad._noise * TMath::Sqrt(256) * lgad._sigRiseTime[0][i])/(0.8 * lgad._sigAmplitude[0][i]);
            //jitter[i-20] = 1000*(lgad._noise * TMath::Sqrt(256) / lgad._sigMaxSlope[0][i]);
            //jitter[i-20] = 1000*(lgad._noise * TMath::Sqrt(256) / lgad._sigSlope[0][i]);
        }
        
        TGraphErrors *gr = new TGraphErrors(arrS,voltage,jitter,0,0);
        gr->SetMarkerStyle(ms[k]);
        gr->SetMarkerColor(col[k]);
        gr->SetMarkerSize(1.5);
        gr->SetLineColor(col[k]);
        if(k==0)
        {
            gr->Draw("apl");
            gr->GetXaxis()->SetTitle("V_{bias} (V)");
            gr->GetYaxis()->SetTitle("Jitter (ps)");
        }
        else
            gr->Draw("plsame");
        gr->GetYaxis()->SetRangeUser(0,1000);
        gr->GetXaxis()->SetRangeUser(0,520);
        leg->AddEntry(gr, legend, "epl");
    }
    leg->Draw();
    canvas->SaveAs("../figures/jitterMoveit.png","png");
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
