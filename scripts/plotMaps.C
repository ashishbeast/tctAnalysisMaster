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
    SetStyle(true);
    TString files[3]={"../data/slapp_1MIP_W_8_pad_C_type_2_gain_1_XY_200_Hole_1_IR.tct",
        "../data/slapp_1MIP_W_8_pad_C_type_2_gain_1_XY_200_Hole_4_IR.tct",
        "../data/slapp_1MIP_W_8_pad_C_type_2_gain_1_XY_200_Hole_5_IR.tct"};
        
    TCanvas *canC = new TCanvas("","",2000,600);
    canC->Divide(3,1);
    TCanvas *canA = new TCanvas("","",2000,600);
    canA->Divide(3,1);
    
    TObjArray *info;
    TString wafer, pad, type, padType, legend, sensor, hole;
    for(Int_t k=0;k<3;++k)
    {
        cout<<"Processing File: "<<files[k]<<endl;
        info = files[k].Tokenize("_");
        wafer = ((TObjString*)(info->At(3)))->String();
        pad   = ((TObjString*)(info->At(5)))->String();
        hole  = ((TObjString*)(info->At(13)))->String();
  
        legend = "Hole-"+hole;
        //Read data
        AnalyzeTCTData lgad(files[k]);
        lgad.CorrectBaseline();
        lgad.CalcNoise();
        lgad.CalculateSignalProperties();
        lgad.PlotMaps(0);
        
        canC->cd(k+1);
        lgad._map[0]->SetTitle(legend);
        lgad._map[0]->Draw("COLZ");
        
        canA->cd(k+1);
        lgad._map[1]->SetTitle(legend);
        lgad._map[1]->Draw("COLZ");
    }
    canA->SaveAs("../figures/amplitudeUniformity.png","png");
    canC->SaveAs("../figures/chargeUniformity.png","png");
    return 0;
}

void SetStyle(Bool_t threeD)
{
    gErrorIgnoreLevel=kError; //Removes annoying Potential memory leak warnings
    gStyle->Reset("Plain");
    gStyle->SetOptTitle(1);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
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
