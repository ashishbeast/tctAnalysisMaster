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
        cout<<"[WARNING] : \tUsage: ./bin/plotGain data/lgad.tct data/pin.tct \n";
        return -1;
    }
    
    //LGAD File
    AnalyzeTCTData lgad(argv[1]);
    lgad.CorrectBaseline();
    lgad.CalcNoise();
    lgad.CalculateWaveformProperties();
    lgad.SetCanvasSettings(false);
    
    //PIN File
    AnalyzeTCTData pin(argv[2]);
    pin.CorrectBaseline();
    pin.CalcNoise();
    pin.CalculateWaveformProperties();
    
    Int_t nVL = lgad._nV1; //Number of Voltage points for LGAD
    Int_t nVP = pin._nV1; //Number of Voltage points for PIN
    
    vector<Float_t> chargeL(nVL);
    vector<Float_t> errL(nVL);
    vector<Float_t> chargeP(nVP);
    vector<Float_t> errP(nVP);
    
    Float_t vBias[nVL];
    Float_t gain[nVL];
    Float_t errG[nVL];
    
    for(Int_t i=0; i< nVP; ++i)
    {
        chargeP[i] = pin._sigCharge[0][i];
        errP[i] = pin._sigChargeError[0][i];
    }
    
    //Average Charge in Pin
    Float_t avgQ = TMath::Mean(chargeP.begin()+5,chargeP.end());
    
    //Error in Average Charge
    Float_t errAvg = TMath::StdDev(chargeP.begin()+5,chargeP.end())/TMath::Sqrt(nVP-5);
    
    Float_t step = TMath::Abs(lgad._tct->V1[1]-lgad._tct->V1[0]);
    
    for(Int_t i=0; i<nVL; ++i)
    {
        vBias[i] = i*step;
        chargeL[i] = lgad._sigCharge[0][i];
        errL[i] = lgad._sigChargeError[0][i];
        gain[i] = chargeL[i]/avgQ;
        errG[i] = gain[i]*TMath::Sqrt(TMath::Power(errL[i]/chargeL[i],2)+TMath::Power(errAvg/avgQ,2));
    }
    
    TCanvas *can = new TCanvas("can","",1200,1000);
    TGraphErrors *gr = new TGraphErrors(nVL,vBias,gain,0,errG);
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.5);
    gr->SetLineColor(kRed);
    gr->Draw("apl");
    gr->GetXaxis()->SetTitle("V_{bias} (V)");
    gr->GetYaxis()->SetTitle("Gain");
    
    lgad.Save(can, "_Gain", ".png");
    return 0;
}
