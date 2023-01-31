#include "AnalyzeTCTData.h"

#include <iostream>
#include <math.h>
#include <vector>
#include "TCanvas.h"

using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout<<"[WARNING] : \tUsage: ./bin/plotJitter data/signal.tct data/noise.tct \n";
        return -1;
    }
    
    //Slew Rate File
    AnalyzeTCTData lgad(argv[1]);
    lgad.CorrectBaseline();
    lgad.CalculateWaveformProperties();
    lgad.SetCanvasSettings(false);
    
    //Noise File
    AnalyzeTCTData noise(argv[2]);
    noise.CorrectBaseline();
    noise.CalculateWaveformProperties();
    
    
    Int_t nV = lgad._nV1;
    
    Int_t arrS = 9;
    
    Int_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
    
    Float_t jitter[arrS];
    Float_t vBias[arrS];
    
    for(Int_t i = 20; i<nV; i+=10)
    {
        vBias[(i-20)/10] = i*step;
        if(lgad._sigSlewRate[0][i] == 0)
            jitter[(i-20)/10] = 0;
        else
            jitter[(i-20)/10] = 1000 * noise._sigNoise[0][i] / lgad._sigSlewRate[0][i];
    }
    
    TCanvas *can = new TCanvas("can","",1200,1000);
    TGraphErrors *gr = new TGraphErrors(arrS,vBias,jitter,0,0);
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.5);
    gr->SetLineColor(kRed);
    gr->Draw("apl");
    gr->GetXaxis()->SetTitle("V_{bias} (V)");
    gr->GetYaxis()->SetTitle("Jitter (ps)");
    
    lgad.Save(can, "_Jitter", ".png");
    
    //can->SaveAs(lgad._outFile+"Jitter.png");
    return 0;
}
