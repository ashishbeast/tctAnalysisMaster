//This macro contains all the functions necessary to analyse TCT Data
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TClonesArray.h"
#include "TObjArray.h"

#define cI -1
#define tmin 0
#define tmax 30
#define wShift 40

Char_t TCTFile[100],TCTFile1[100],TCTFile2[100];
TObjArray *SensorInfo;
TString Sensor,inPath,outPath;

Bool_t isBM_ON;

Int_t col[12] = {kRed,kBlue,kBlack,13, kGreen,6,7,8,9,28,34,49};
Int_t ms[12] = {20,21,22,23,29,33,34,39,41,43,45,47};
Int_t nAC; //Number of Active Channels
Int_t *activeChannel;
Int_t nX, nY, nZ, nV, nV1, nV2;
Int_t stepX, stepY, stepZ, stepV, stepV1=0, stepV2=0;
TString titleCharge,titleAmplitude;

Float_t SignalShift(PSTCT readTCT)
{
    //Getting the active channels from Oscilloscope
    //Oscilloscope Stuff--------------------------------------------------------------
    Int_t scopeChannels[4] = {readTCT.WFOnOff[0],readTCT.WFOnOff[1],readTCT.WFOnOff[2],readTCT.WFOnOff[3]};
    nAC = 0;
    for(Int_t i = 0;i<4;++i)
        if(scopeChannels[i] == 1)
            nAC++;
    if(nAC == 0)
    {
        cout<<"ERROR: No Active Channels!!!"<<endl;
        exit(0);
    }
    else
        cout<<"Number of Active Channels: "<<nAC<<endl;
    activeChannel = new Int_t [nAC];
    Int_t countAC=0;
    for(Int_t i = 0;i<4;++i)
    {
        if(scopeChannels[i] != 1)
            continue;
        activeChannel[countAC] = i;
        countAC++;
    }
    
    //Get the number of Steps for each axis and the step
    nX = readTCT.Nx-1;
    nY = readTCT.Ny-1;
    nZ = readTCT.Nz-1;
    stepX = readTCT.dx;
    stepY = readTCT.dy;
    stepZ = readTCT.dz;
    
    nV1 = readTCT.NU1-1;
    if(nV1>1)
        stepV1 = readTCT.U1[1]-readTCT.U1[0];
    
    nV2 = readTCT.NU2-1;
    if(nV2>1)
        stepV2 = readTCT.U2[1]-readTCT.U2[0];
    
    if(stepV1 != 0)
    {
        nV = nV1;
        stepV = stepV1;
    }
    else
    {
        nV = nV2;
        stepV = stepV2;
    }
    
    //cout<<"(nX,nY,nZ) : "<<"("<<nX<<","<<nY<<","<<nZ<<")"<<endl;
    
    //Check for Beam Monitor
    Int_t bmIndex = readTCT.indx(rand()%(nX+1),rand()%(nY+1),rand()%(nZ+1),rand()%(nV+1),0);
    Float_t bmIntensity = readTCT.xyz[7][bmIndex];
    if(bmIntensity != 0.)
    {
        isBM_ON = kTRUE;
        cout<<"Beam Monitor Status: ON"<<endl;
    }
    else
    {
        isBM_ON = kFALSE;
        cout<<"Beam Monitor Status: OFF"<<endl;
    }
    
    //Check for Signal above the threshold
    TH1F *signal, *signalCopy;
    Float_t threshold = 5, amplitude;
    do
    {
        signal = readTCT.GetHA(activeChannel[0],rand()%(nX+1),rand()%(nY+1),rand()%(nZ+1),rand()%(nV+1),0);
        signalCopy = (TH1F*)signal->Clone();
        signalCopy->Scale(cI);
        amplitude = signalCopy->GetMaximum();
    }
    while(amplitude<threshold);
    
    Int_t nBins = signalCopy->GetNbinsX(); //Total bins
    Float_t time0 = readTCT.t0;
    Float_t timeBin = readTCT.dt;
    Float_t amp[nBins];
    Float_t time[nBins];
    
    Float_t ampMax=0;
    for(Int_t i=0;i<nBins;++i)
    {
        amp[i]= signalCopy->GetBinContent(i);
        time[i] = time0+i*timeBin;
        if(amp[i]>ampMax)
            ampMax = amp[i];
    }
    Float_t amp10 = ampMax*0.1; //10% of maximum signal
    Float_t time10;
    
    for(Int_t i=0;i<nBins-1;++i)
    {
        if(amp[i]>=amp10)
        {
            time10 = time[i];
            break;
        }
    }
    Float_t shift = (time10/1e-9)-1;//1 ns before the 10% of Max Signal
    cout<<"Shift in Waveform: "<<shift<<endl;
    return shift;
    //return (time10/1e-9)-1;
}

Float_t FindSignalShift(TString file)
{
    //------Setting the input and output Path
    SensorInfo = file.Tokenize("_");
    Sensor = ((TObjString*)(SensorInfo->At(0)))->String();
    inPath = "../data/";
    outPath = "../figures/";
    
    //Reading file data
    cout<<"=====Finding the shift in the Signal====="<<endl;
    strcpy(TCTFile,inPath+file);
    PSTCT readTCT((Char_t*)TCTFile,0,2);
    
    Double_t signalShift = SignalShift(readTCT);
    return signalShift;
}

//calculates the noise before the signal and sets it as the error in the histogram
TH1F* GetSetNoise(TH1F *his)
{
    Int_t left=1;                                   //First Bin
    Int_t right = his->GetXaxis()->FindBin(0.);     //The  last Bin before Trigger (signal)
    
    TH1D *hist = new TH1D("hist","",10000,-0.3,0.3);
    
    for(Int_t k = left; k<=right;++k)
        hist->Fill(his->GetBinContent(k));
    
    Double_t noise = hist->GetStdDev();
    
    for(Int_t k = left; k<=his->GetNbinsX();++k)
        his->SetBinError(k,noise);
    
    delete hist;
    return his;
}

Double_t GetNoise(TH1F *his)
{
    Int_t left=1;                                   //First Bin
    Int_t right = his->GetXaxis()->FindBin(0.);     //The  last Bin before Trigger (signal)
    
    TH1D *hist = new TH1D("hist","",10000,-0.3,0.3);
    
    for(Int_t k = left; k<=right;++k)
        hist->Fill(his->GetBinContent(k));
    
    Double_t noise = hist->GetStdDev();
    return noise;
}



//This function calculates the noise before signal for all the measurements
//The calculated noise is set as error in histogram bins
void CorrectNoise(PSTCT fRead)
{
    TH1D *nHist = new TH1D("nHist","", 100, -1,1);
    Int_t right, left = 1;
    
    Int_t nX = fRead.Nx-1;
    Int_t nY = fRead.Ny-1;
    Int_t nZ = fRead.Nz-1;
    Int_t nV = fRead.NU1-1;
    
    Double_t nWf = (nX+1)*(nY+1)*(nZ+1)*(nV+1); //Total Number of Acquisition
    
    TObjArray *his = new TObjArray(nWf);
    
    for(Int_t j = 0;j <= nV;++j)
        for(Int_t k = 0;k <= nZ;++k)
            for(Int_t l = 0;l <= nY;++l)
                for(Int_t m = 0;m <= nX;++m)
                {
                    TH1F *h1 = fRead.GetHA(2,m,l,k,j,0);
                    his->Add(h1);
                }
    
    for(Int_t i = 0; i < nWf; i++)
    {
        TH1F *hist = (TH1F*)his->At(i);
        right = hist->GetXaxis()->FindBin(0.);
        for(Int_t j = left;j<right;++j)
            nHist->Fill(hist->GetBinContent(j));
    }
    
    Double_t noise = nHist->GetStdDev();
    
    for(Int_t i = 0; i < nWf; i++)
    {
        TH1F *hist = (TH1F*)his->At(i);
        right = hist->GetNbinsX();
        for(Int_t j = left;j<right;++j)
            hist->SetBinError(j,noise);
    }
    // TCanvas *can = new TCanvas("can","",1200,1000);
    // nHist->Draw();
    // nHist->GetXaxis()->SetTitle("noise (mV)");
    // nHist->GetYaxis()->SetTitle("Count");
} 

//Return the normalised Amplitude
Double_t GetNormalisedAmplitude(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Double_t *error)
{
    Int_t Index = fRead.indx(x,y,z,v1,v2);
    Double_t beamMon = fRead.xyz[7][Index];  //BeamMonitor Intensity
    TH1F *h1 = fRead.GetHA(2,x,y,z,v1,v2);
    TH1F *hisA = (TH1F*)h1->Clone("hisA");
    hisA->Scale(cI);
    
    Double_t maxBin = hisA->GetMaximumBin();
    Double_t Amplitude = hisA->GetBinContent(maxBin);
    Double_t normAmp = Amplitude*(1000/beamMon);  //Dividing by Beam Monitor Intensity
    *error = hisA->GetBinError(maxBin);
    return normAmp;
}

//Return the normalised Charge
Double_t GetNormalisedCharge(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Double_t* error)
{
    Int_t Index = fRead.indx(x,y,z,v1,v2);
    Double_t beamMon = fRead.xyz[7][Index];  //BeamMonitor Intensity
    TH1F *h1 = fRead.GetHA(2,x,y,z,v1,v2);
    TH1F *hisC = (TH1F*)h1->Clone("hisC");
    hisC->Scale(cI);
    
    Double_t errC;
    Double_t ti = tmin;
    Double_t tf = tmax;
    
    Int_t mintime;
    Int_t maxtime;
    if(tf==-1111)
        maxtime=hisC->GetNbinsX()-1;
    else
        maxtime=hisC->GetXaxis()->FindBin(tf);
    if(ti==-1111)
        mintime=1;
    else
        mintime=hisC->GetXaxis()->FindBin(ti);
    
    Double_t charge = hisC->IntegralAndError(mintime,maxtime,errC,"width");
    Double_t normCharge = charge*(1000/beamMon);
    *error = errC;
    return normCharge;
}

Double_t GetAmplitude(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Double_t *error)
{
    TH1F *hA = fRead.GetHA(2,x,y,z,v1,v2);
    TH1F *hisA = (TH1F*)hA->Clone("hisA");
    hisA->Scale(cI);
    
    Double_t maxBin = hisA->GetMaximumBin();
    Double_t Amplitude = hisA->GetBinContent(maxBin);
    //Double_t normAmp = Amplitude*(1000/beamMon);  //Dividing by Beam Monitor Intensity
    *error = hisA->GetBinError(maxBin);
    return Amplitude;
}

Double_t GetCharge(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Double_t* error)
{
    TH1F *hC = fRead.GetHA(2,x,y,z,v1,v2);
    TH1F *hisC = (TH1F*)hC->Clone("hisC");
    hisC->Scale(cI);
    
    Double_t errC;
    Double_t ti = tmin;
    Double_t tf = tmax;
    
    Int_t mintime;
    Int_t maxtime;
    if(tf==-1111)
        maxtime=hisC->GetNbinsX()-1;
    else
        maxtime=hisC->GetXaxis()->FindBin(tf);
    if(ti==-1111)
        mintime=1;
    else
        mintime=hisC->GetXaxis()->FindBin(ti);
    
    Double_t charge = hisC->IntegralAndError(mintime,maxtime,errC,"width");
    *error = errC;
    return charge;
}

//Return the Ratio of Amplitude over Charge
Double_t GetRatio(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2)
{
    TH1F *hC = fRead.GetHA(2,x,y,z,v1,v2);
    hC->Scale(cI);
    
    Double_t ti = tmin;
    Double_t tf = tmax;
    
    Double_t Amp = hC->GetBinContent(hC->GetMaximumBin());
    Double_t Charge = hC->Integral(hC->GetXaxis()->FindBin(ti),hC->GetXaxis()->FindBin(tf));
    
    Double_t ratio = Amp/Charge;
    return ratio;
}


//Calculate Depletion voltage by fitting the parabola in amp and charge ratio
void GetDepletionVoltage(TGraphErrors* gr)
{
    Double_t yMin = TMath::MinElement(gr->GetN(),gr->GetY()); //First Change of Slope
    
    cout<<"Minimum Value of Ratio is: "<<yMin<<endl;
    Double_t xMin ; // Voltage Value corresponding to change of Slope
    Double_t dV = 0; //Depletion Voltage
    
    for(Int_t k = 0;k<gr->GetN();k++)
        if(gr->GetPointY(k) == yMin)
            xMin = gr->GetPointX(k);
    
    cout<<"V_{Bias} corresponding to minimum ratio is: "<<xMin<<endl;
    TCanvas *can = new TCanvas("can","",1200,800);
    
    TGraph *hr = (TGraph*)gr->Clone();
    
    Double_t xlow = 0;
    Double_t xhigh = -40;
    TF1 *fP = new TF1("fP","[2]*x*x+[1]*x+[0]",xMin-10,xMin+10);  //Pol 2 fitting
    fP->SetLineColor(kRed);
    hr->Fit(fP);
    hr->GetXaxis()->SetTitle("V_{bias} (V)");
    hr->GetYaxis()->SetTitle("Amplitude (mV)/Charge");
    hr->SetMarkerStyle(20);
    hr->SetMarkerColor(kBlue);
    hr->Draw("apl");
    hr->GetXaxis()->SetRangeUser(xMin-10,xMin+10);
    
    dV =  -1*fP->GetParameter(1)/(2*fP->GetParameter(2));    //x coordinate of the vertex
    cout<<"Depletion Voltage of Gain Layer is: "<<dV<<endl;
    TString legend = "V_{dep}^{GL}";
    legend += std::to_string(dV);
    TLegend *leg = new TLegend(0.3422371,0.6134021,0.5834725,0.9033505,NULL,"brNDC");
    leg->SetHeader(legend);
    leg->Draw();
}

//This function calculates the rise time (90/10) of the signal.
//GetRiseTime(Histogram, Threshold in Percentage)
//GetRiseTime(xAxis[],yAxis[], sizeArray, Threshold Value)
Double_t GetRiseTime(TH1F* his, Double_t low, Double_t high)
{
    Int_t n = his->GetNbinsX();
    Double_t sigMax = his->GetMaximum();
    
    Double_t sigLow  = sigMax * (low/100); //low limit of signal to find rise time
    Double_t sigHigh  = sigMax * (high/100); //high limit of signal
    Double_t tcrossLow = 0; //The time at which signal is low% of its max value
    Double_t tcrossHigh = 0;//The time at which signal is high% of its max value
    Bool_t gotMinimum = kFALSE;
    Bool_t gotMaximum = kFALSE;
    for(Int_t i=0;i<n-1;++i)
    {
        Double_t signal = his->GetBinContent(i);
        Double_t signalAbove = his->GetBinContent(i+1);
        Double_t time = his->GetXaxis()->GetBinCenter(i);
        Double_t timeAbove = his->GetXaxis()->GetBinCenter(i+1);
        
        if(gotMinimum==kFALSE)
        {
            if(sigLow-signal == 0)
            {
                tcrossLow = time;
                gotMinimum = kTRUE;
            }
            
            if(signal <sigLow && signalAbove >sigLow)
            {
                tcrossLow = time + (sigLow-signal) * (timeAbove-time)/(signalAbove-signal);
                gotMinimum = kTRUE;
            }
        }
        if(gotMaximum==kFALSE)
        {
            if(sigHigh-signal == 0)
            {
                tcrossHigh = time;
                gotMaximum = kTRUE;
                break;
            }
            
            if(signal <sigHigh && signalAbove >sigHigh)
            {
                tcrossHigh = time + (sigHigh-signal) * (timeAbove-time)/(signalAbove-signal);
                gotMaximum = kTRUE;
                break;
            }
        }
    }
    return tcrossHigh-tcrossLow;
}

//This is the function for Constant Fraction Discriminator
//CFD(Histogram, Threshold in Percentage)
//CFD(xAxis[],yAxis[], sizeArray, Threshold Value)
Double_t CFD(TH1F* his, Double_t threshold)
{
    Int_t n = his->GetNbinsX();
    Double_t sigMax = his->GetMaximum();
    if(threshold > 1)
        threshold = sigMax * (threshold/100);
    Double_t TCross;  //Crossing time for the threshold
    
    for(Int_t i=0;i<n-1;++i)
    {
        Double_t signal = his->GetBinContent(i);
        Double_t signalAbove = his->GetBinContent(i+1);
        Double_t time = his->GetXaxis()->GetBinCenter(i);
        Double_t timeAbove = his->GetXaxis()->GetBinCenter(i+1);
        if(threshold-signal == 0)
        {
            TCross = time;
            return TCross;
        }
        if(signal <= threshold && signalAbove > threshold)
        {
            TCross = time + (threshold-signal) * (timeAbove-time)/(signalAbove-signal);
            return TCross;
        }
    }
    return 0;
}



// Double_t GetToACFD(TH1F* his, Double_t threshold)
// {
//   Double_t sigMax = his->GetMaximum();
//   Bool_t gotMinimum = kFALSE;
//   Bool_t gotMaximum = kFALSE;
//   Int_t n = his->GetNbinsX();
//   Double_t xAxis[n];
//   Double_t yAxis[n];
//   Double_t errX[n];
//   Double_t errY[n];
//   for(int i = 0;i<n;++i)
//     {
//       xAxis[i]= his->GetXaxis()->GetBinCenter(i);
//       yAxis[i]= his->GetBinContent(i);
//       errX[i]=his->GetBinWidth(i);
//       errY[i]=his->GetBinError(i);
//     }
//   Double_t cfd = sigMax * (threshold/100);
//   Double_t ToAThres =0;

//   //The ToA is the time at threshold
//   //  TGraphErrors *hr = new TGraphErrors(n,xAxis,yAxis,errX,errY);
//   TGraphErrors *hr = new TGraphErrors(n,yAxis,xAxis,errY,errX);
//   for(Int_t i=0;i<n;++i)
//     {
//       if(gotMaximum==kFALSE)
// 	{
// 	  if(yAxis[i] >= cfd)
// 	    {
// 	      ToAThres =hr->Eval(yAxis[i],0,"");//The time at which signal crosses Threshold
// 	      gotMaximum = kTRUE;
// 	      break;
// 	    }
// 	  else
// 	    {
// 	      gotMaximum = kFALSE;
// 	    }
// 	}
//     }
//   return ToAThres;
// }

//Draws two axis left and right
void DrawPadWithTwoAxis(TGraphErrors *gr, TGraphErrors *hr,TString RightTitle, TLegend *leg)
{
    TCanvas *canvas = new TCanvas("canvas","",50,114,1132,987);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000); //will be transparent
    pad1->Draw();
    pad1->cd();
    gr->Draw("apl");
    pad1->Update();
    
    //compute the pad range with suitable margins
    Double_t ymin = TMath::MinElement(hr->GetN(),hr->GetY());    //Y axis min for 2nd graph
    Double_t ymax = TMath::MaxElement(hr->GetN(),hr->GetY())+5;  //Y axis max for 2nd graph
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    
    Double_t xmin = TMath::MinElement(gr->GetN(),gr->GetX());    //X axis min for 1st graph (becomes the same for 2nd graph)
    Double_t xmax = TMath::MaxElement(gr->GetN(),gr->GetX())+20; //X axis max for 2nd graph (becomes the maximum for 2nd graph as well)
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad2->Draw();
    pad2->cd();
    hr->Draw("][ plsame");
    pad2->Update();
    
    // draw axis on the right side of the pad
    TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
    axis->SetTitle(RightTitle);
    axis->SetTitleColor(kRed);
    axis->SetTitleFont(42);
    axis->SetLabelColor(kRed);
    axis->SetNdivisions(507);
    axis->Draw();
    
    leg->Draw();
    return canvas;
}

void RatioAmpCharge(PSTCT fRead)
{
    TString Title= "V_{bias} (V)";
    Int_t nWf=fRead.NU1-1;
    Double_t step= (fRead.U1[nWf]-fRead.U1[0])/nWf;
    Double_t *QS = new Double_t [nWf+1];
    Double_t *Err = new Double_t [nWf+1];
    Double_t errX[nWf+1];
    Double_t XAxis[nWf+1];
    
    for(Int_t i=0;i<=nWf;++i)
    {
        TH1F *his = fRead.GetHA(2,0,0,0,i,0);
        his->Scale(cI);
        TH1F *hA = (TH1F*)his->Clone("hA");
        TH1F *hC = (TH1F*)his->Clone("hC");
        errX[i] = 0;
        XAxis[i] = step*i;
        Double_t errC;
        Double_t errA;
        Double_t charge =  GetNormalisedCharge(fRead,0,0,0,i,0, &errC);
        Double_t amplitude = GetNormalisedAmplitude(fRead,0,0,0,i,0, &errA);
        QS[i] = amplitude/charge;
        Err[i] = QS[i]*TMath::Sqrt(TMath::Power(errC/charge,2)+TMath::Power(errA/amplitude,2));
    }
    TCanvas *canvas = new TCanvas("canvas","",1200,1000);
    TGraphErrors *gr = new TGraphErrors(nWf+1,XAxis,QS,errX,Err);
    gr->Draw("apl");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->SetMarkerSize(1);
    gr->SetLineColor(4);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle(Title);
    gr->GetYaxis()->SetTitle("Amplitude (mV)/Charge");
    //GetDepletionVoltage(gr);
    return canvas;
}
void FindFocus(PSTCT fRead)
{
    Int_t nO = fRead.Nz-1;
    Int_t stepO = fRead.dz;
    Double_t optDis[nO+1];
    Double_t fwhm[nO+1];
    Double_t errfwhm[nO+1];
    Int_t nS,stepS;
    Int_t proj;
    
    Int_t nX = fRead.Nx-1;
    Int_t nY = fRead.Ny-1;
    
    Double_t Xi;         //   low limit for fit  =0
    Double_t Xf;         //   high limit for fit =100
    Double_t FWHM = 10;  //   expected FWHM
    Double_t Level ;     //   Middle of S-curve
    
    //--------------------------------------------------------------------------------
    //Determine the scanning axis and set Parameters
    //________________________________________________________________________________
    if(nX>1 && nY>1)
    {
        cout<<"Error: Focus can not be calculated for XY Scan."<<endl;
        exit(0);
    }
    else if(nX>1)
    {
        cout<<"Scanning Axis is along X-Axis"<<endl;
        nS = nX;
        stepS = fRead.dx;
        if(stepS<0)
            stepS = -1 *stepS;
        Xi = 0;
        Xf = nS*stepS;
        Level = (Xf-Xi)/2;
        proj = 0;
    }
    else
    {
        cout<<"Scanning Axis is along Y-Axis"<<endl;
        nS = nY;
        stepS = fRead.dy;
        if(stepS<0)
            stepS = -1 *stepS;
        Xi = 0;
        Xf = nS*stepS;
        Level = (Xf-Xi)/2;
        proj = 1;
    }
    cout<<"Focus Point of the Z axis will be Calculated"<<endl;
    
    //-------------------------------------------------------------------------------
    //Fiting of the S-curve is done by error function which is defined as:
    // f(x) = erf((x-x0)/r)*M+y0
    // (x0,y0) are the coordinate of the middle point in the S-curve
    // y0 = Ymax/2
    // r is the radius
    // M = Ymax*,  is the maximum value of the S-curve (Saturation Value)
    // *if curve start from 0 and reaches a maximum value then M = Ymax/2
    //_______________________________________________________________________________
    
    TF1 *scurve=new TF1("scurve","TMath::Erf((x-[0])/[1])*[2]+[3]",Xi,Xf);
    scurve->SetParNames("L","SS","C1","C2");
    scurve->SetParameters(Level,FWHM,8,4);
    scurve->SetLineColor(kBlue);
    
    //scan over the measurements and select waveforms
    TCanvas *canvas = new TCanvas("canvas","",1200,1000);
    Double_t xMin = 0;
    for(Int_t i=0;i<=nO;++i)
    {
        Double_t charge[nS];
        Double_t errCh[nS];
        Double_t scanDis[nS];
        Double_t errX[nS];
        for(Int_t j=0;j<nS;++j)
        {
            scanDis[j] = j*stepS;
            Double_t errC;
            switch(proj)
            {
                case 0:
                    //charge[j] = GetNormalisedCharge(fRead,j,0,i,0,0,&errC);
                    charge[j] = GetCharge(fRead,j,0,i,0,0,&errC);
                    errCh[j] = errC;
                    break;
                case 1:
                    //charge[j] = GetNormalisedCharge(fRead,0,j,i,0,0,&errC);
                    charge[j] = GetCharge(fRead,0,j,i,0,0,&errC);
                    errCh[j] = errC;
                    break;
            }
            errX[j]=0;
        }
        TGraph *gr = new TGraph(nS,scanDis,charge);
        gr = new TGraph(nS,scanDis,charge);
        gr->SetLineColor((i+1)%8+1);
        
        gr->Fit(scurve,"R");
        if(i==0)
        {
            gr->Draw("apl");
            gr->GetXaxis()->SetTitle("Scanning distance [#mum]");
            gr->GetYaxis()->SetTitle("Charge (arb.)");
        }
        else
            gr->Draw("plsame");
        
        fwhm[i]= scurve->GetParameter(1)*2.35/TMath::Sqrt(2);
        errfwhm[i] = scurve->GetParError(1);
        optDis[i] = i*stepO;
        if(i==0) continue;
        if(fwhm[i]<fwhm[i-1])
            xMin = optDis[i];
    }
    
    //This part fits the FWHM graph and Calculates the focus spot
    TCanvas *canfocus = new TCanvas("canfocus","",1200,1000);
    
    //Function to fit FWHM
    //TF1 *fP = new TF1("fP","[2]*x*x+[1]*x+[0]",optDis[0],optDis[nO+1]);
    TF1 *fP = new TF1("fP","[2]*x*x+[1]*x+[0]",xMin-1000,xMin+1000);
    fP->SetLineColor(kBlue);
    
    //  TGraph *grFWHM=new TGraph(nO+1,optDis,fwhm);
    TGraphErrors *grFWHM=new TGraphErrors(nO+1,optDis,fwhm,0,errfwhm);
    grFWHM->SetMarkerStyle(20);
    grFWHM->SetMarkerSize(1);
    grFWHM->SetMarkerColor(kRed);
    grFWHM->SetLineColor(kRed);
    
    grFWHM->Fit("fP","RMS");
    grFWHM->GetHistogram()->GetXaxis()->SetTitle("Optical distance [#mum]");
    grFWHM->GetHistogram()->GetYaxis()->SetTitle("FWHM [#mum]");
    grFWHM->Draw("APL");
    
    
    Double_t foSpot = -1*fP->GetParameter(1)/(2*fP->GetParameter(2));
    Double_t SpotSize = fP->GetParameter(2)*TMath::Power(foSpot,2)+fP->GetParameter(1)*foSpot+fP->GetParameter(0);
    
    TString legend1 = "Z_{focus}: ";
    legend1 += std::to_string(foSpot);
    TString legend2 = "Spot Size = ";
    legend2 += std::to_string(SpotSize);
    
    //TLegend *leg = new TLegend(0.3422371,0.6134021,0.5834725,0.9033505,NULL,"brNDC");
    TPaveText *pt = new TPaveText(0.345576,0.7510684,0.7103506,0.8237179,"nbNDC");
    pt->SetTextColor(kBlue);
    pt->SetFillColor(0);
    //pt->SetBorderSize(0);
    pt->AddText(legend1);
    pt->AddText(legend2);
    pt->Draw("");
}

//This function draws a 2D profile of Laser Light
void DrawLaserProfile(PSTCT fRead)
{
    Int_t nZ = fRead.Nz-1;
    Int_t stepZ = fRead.dz;
    
    Int_t nX = fRead.Nx-1;
    Int_t nY = fRead.Ny-1;
    Int_t nS,stepS,proj;
    TString title;
    //--------------------------------------------------------------------------------
    //Determine the scanning axis
    //________________________________________________________________________________
    if(nX>1 && nY>1)
    {
        cout<<"Error: Laser Profile cannot be drawn for XY Scan."<<endl;
        exit(0);
    }
    else if(nX>1)
    {
        cout<<"Scanning Axis is along X-Axis"<<endl;
        nS = nX;
        stepS = fRead.dx;
        if(stepS<0)
            stepS = -1 *stepS;
        proj = 0;
        title = "X-Axis (#mum)";
    }
    else
    {
        cout<<"Scanning Axis is along Y-Axis"<<endl;
        nS = nY;
        stepS = fRead.dy;
        if(stepS<0)
            stepS = -1 *stepS;
        proj = 1;
        title = "Y-Axis (#mum)";
    }
    
    TH2D *laserProfile = new TH2D("laserProfile","",nS,0,nS*stepS,nZ,0,nZ*stepZ);
    
    //Charge Map
    TCanvas *canvasLaser = new TCanvas("canvasLaser","",1200,1000);
    for(Int_t i=0;i<=nZ;++i)
        for(Int_t j=0;j<=nS;++j)
        {
            Double_t errC;
            Double_t charge;
            if(proj == 0)
            {
                //charge = GetNormalisedCharge(fRead,j,0,i,0,0,&errC);
                charge = GetCharge(fRead,j,0,i,0,0,&errC);
            }
            else if(proj == 1)
            {
                //charge = GetNormalisedCharge(fRead,0,j,i,0,0,&errC);
                charge = GetCharge(fRead,0,j,i,0,0,&errC);
            }
            
            laserProfile->Fill(j*stepS,i*stepZ,charge);
        }
    laserProfile->Draw("COLZ");
    laserProfile->GetXaxis()->SetTitle(title);
    laserProfile->GetYaxis()->SetTitle("Z Axis [#mum]");
    laserProfile->GetZaxis()->SetTitle("Charge (arb.)");
    return canvasLaser;
}

void DrawChargeMap(PSTCT fRead)
{
    Int_t nX = fRead.Nx-1;
    Int_t nY = fRead.Ny-1;
    Int_t stepX = fRead.dx;
    Int_t stepY = fRead.dy;
    TH2D *xyC = new TH2D("xyC","",nX,0,nX*stepX,nY,0,nY*stepY);
    
    TCanvas *cMap = new TCanvas("cMap","",1200,1000);
    //Charge Map
    for(Int_t i=0;i<=nY;++i)
        for(Int_t j=0;j<=nX;++j)
        {
            Double_t errC;
            Double_t charge = GetNormalisedCharge(fRead,j,i,0,0,0,&errC);
            xyC->Fill(j*stepX,i*stepY,charge);
        }
    
    xyC->Draw("COLZ");
    xyC->GetXaxis()->SetTitle("X Axis [#mum]");
    xyC->GetYaxis()->SetTitle("Y Axis [#mum]");
    xyC->GetZaxis()->SetTitle("Norm. Charge (arb.)");
    return cMap;
}
void DrawAmplitudeMap(PSTCT fRead)
{
    Int_t nX = fRead.Nx-1;
    Int_t nY = fRead.Ny-1;
    Int_t stepX = fRead.dx;
    Int_t stepY = fRead.dy;
    TH2D *xyA = new TH2D("xyA","",nX,0,nX*stepX,nY,0,nY*stepY);
    
    TCanvas *aMap = new TCanvas("aMap","",1200,1000);
    //Amplitude Map
    for(Int_t i=0;i<=nY;++i)
        for(Int_t j=0;j<=nX;++j)
        {
            Double_t errA;
            Double_t amplitude = GetNormalisedAmplitude(fRead,j,i,0,0,0,&errA);
            xyA->Fill(j*stepX,i*stepY,amplitude);
        }
    xyA->Draw("COLZ");
    xyA->GetXaxis()->SetTitle("X Axis [#mum]");
    xyA->GetYaxis()->SetTitle("Y Axis [#mum]");
    xyA->GetZaxis()->SetTitle("Norm. Amplitude (mV)");
    return aMap;
}
//Set Style for canvas
void SetStyle(Bool_t threeDimension = kFALSE)
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
    //gStyle->SetPadRightMargin(0.03025943);
    gStyle->SetPadBottomMargin(0.1353861);//0.015
    gStyle->SetPadLeftMargin(0.1418293);//0.21
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kBlue);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.04,"xyz");
    gStyle->SetLabelSize(0.04,"y");
    gStyle->SetLabelOffset(0.005,"y");
    gStyle->SetLabelOffset(0.015,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.06,"xyz");
    //  gStyle->SetTitleOffset(1.14,"y");//1.81
    if(threeDimension)
    {
        gStyle->SetPadRightMargin(0.18);
        gStyle->SetTitleOffset(1.15,"y");//1.81
    }
    else
        gStyle->SetTitleOffset(0.81,"y");//1.81
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


//Return the normalised Amplitude
Double_t GetNormalizedAmplitude(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Double_t *error)
{
    Int_t Index = fRead.indx(x,y,z,v1,v2);
    Double_t beamMon = fRead.xyz[7][Index];  //BeamMonitor Intensity
    TH1F *h1 = fRead.GetHA(2,x,y,z,v1,v2);
    TH1F *hisA = (TH1F*)h1->Clone("hisA");
    hisA->Scale(cI);
    
    hisA = GetSetNoise(hisA);
    
    Double_t maxBin = hisA->GetMaximumBin();
    Double_t Amplitude = hisA->GetBinContent(maxBin);
    Double_t normAmp = Amplitude*(1000/beamMon);  //Dividing by Beam Monitor Intensity
    *error = hisA->GetBinError(maxBin);
    return normAmp;
}

//Return the normalised Charge
Double_t GetNormalizedCharge(PSTCT fRead, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Double_t* error)
{
    Int_t Index = fRead.indx(x,y,z,v1,v2);
    Double_t beamMon = fRead.xyz[7][Index];  //BeamMonitor Intensity
    TH1F *h1 = fRead.GetHA(2,x,y,z,v1,v2);
    TH1F *hisC = (TH1F*)h1->Clone("hisC");
    hisC->Scale(cI);
    
    hisC = GetSetNoise(hisC);
    
    Double_t errC;
    Double_t ti = tmin;
    Double_t tf = tmax;
    
    Int_t mintime;
    Int_t maxtime;
    if(tf==-1111)
        maxtime=hisC->GetNbinsX()-1;
    else
        maxtime=hisC->GetXaxis()->FindBin(tf);
    if(ti==-1111)
        mintime=1;
    else
        mintime=hisC->GetXaxis()->FindBin(ti);
    
    Double_t charge = hisC->IntegralAndError(mintime,maxtime,errC,"width");
    Double_t normCharge = charge*(1000/beamMon);
    *error = errC;
    return normCharge;
}



