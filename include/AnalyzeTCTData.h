#ifndef ANALYZETCTDATA_H
#define ANALYZETCTDATA_H


#include "ReadTCTFile.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TVectorF.h"

using namespace std;

class AnalyzeTCTData
{
public:
    AnalyzeTCTData(const char*  inFile); //Constructor
    ~AnalyzeTCTData(); //Destructor
    
    //Set Variables
    void SetPolarity(Int_t );
    void SetIntegralLimits(Float_t, Float_t);
    void SetTermination(Float_t);
    
    void GetActiveChannels();    //Find Active Channels
    
    Bool_t FindWFShift();        //Find Signal Starting Point
    
    void EstimateWFShift();        //Find Signal Starting Point
    
    void CorrectBaseline();      //Corrects the y shift in all Waveforms
        
    void CalculateWaveformProperties();
    
    void SetCanvasSettings(Bool_t threeD);     // style for plotting
    
    void AnalysisAction(); //Determines what to do with the file and perform Operations
    
    Int_t GetOperation(Int_t *binary); // Converts Binary to decimal
    
    void CalcNoise();  // calculate baseline and noise for all the events
    
    void CalcBaselineNoise();
    
    void SaveSignalShape();
    
    Float_t CalcSignalNoise(TH1F* his); // Calculate noise in a given Waveform
    
    Float_t CalcSlewRate(TH1F* his); //Calculate max slope in the leading edge of the signal
    
    Float_t CalcSlope(TH1F* his); //Calculate slope of the leading edge
    
    Float_t CalcSlope(TH1F* his, Float_t threshold); //Calculate slope of the leading edge for a given threshold
    
    Float_t CalcRiseTime(TH1F* his); // Calculate 20 to 80% rise time
    
    Float_t CalcCFD(TH1F *his, Float_t thr);  // Calculation of the CFD threhsold crossing
    
    Float_t CalcTOA(TH1F *his);  // Calculation of the Time of Arrival of signal
    
    Float_t CalcTOT(TH1F *his, Float_t thr);  // Calculation of the time over threshold
    
    Float_t CalcCharge(TH1F*, Float_t*);   //Calculation of the integral of the signal
    
    Float_t CalcCharge(Int_t channel, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Float_t *error);
    
    Float_t CalcAmplitude(TH1F*);   //Calculation of the integral of the signal
    
    Float_t CalcAmplitude(Int_t channel, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Float_t *error);
    
    Float_t CalcRatio(TH1F*, Float_t*); //Calculate Amplitude/Charge
    
    Float_t LinearInterpolation(vector<Float_t> x, vector<Float_t> y, Float_t thr);   //Calculate the interpolation between two points
    
    void PlotAxisSweep(Int_t operation);
    
    void PlotAxisSweep(Int_t channel, Int_t operation);
    
    //TH2F* GetMap(Int_t channel, Int_t mapType, Int_t mode);
    void PlotMaps(Int_t mapType);
    
    void PlotMaps(Int_t mapType, Int_t mode);
    
    void PlotMaps(Int_t channel, Int_t mapType, Int_t mode);
    
    void CalcFDV();   //Calculate the full depletion voltage of the sensor
    
    void GetFocus();  //Give back the focus of the laser
    
    void GetFocus(int channel);
    
    void Save(TCanvas* canvas, TString canvasType, TString ext);   //Save the canvas
    
public:
    ReadTCTFile *_tct;
    Float_t _tA;
    Int_t _polarity; //polarity of Signal
    Float_t _tmin;   //Signal Start
    Float_t _tmax;   //Signal Stop
    Float_t _noise;  //Baseline Noise (Average of all the measurements)
    Int_t _nAC = 0.; //Number of active channels
    Int_t* _aCH;     //Active Channels
    TString _fileName;
    TString _outFile;
    Int_t _indexMaxSignal;
    Int_t _chMaxSignal;
    Int_t _termination;
    TCanvas* _canvas;
    Int_t _events;
    
    Int_t _nX, _nY, _nZ, _nV1, _nV2;
    Float_t _threshold;
    Bool_t _isNoise;
    Bool_t _isAveraged;
    
    Bool_t _bmON = false;
    
    //Waveform Properties
    Float_t** _wavNoise;
    
    //Signal Properties
    Int_t _index;
    Float_t* _bmValue;
    TH1F*** _histo;
    Float_t** _sigNoise;
    Float_t** _sigCharge;
    Float_t** _sigChargeError;
    Float_t** _sigAmplitude;
    Float_t** _sigAmplitudeError;
    Float_t** _sigRatio;
    Float_t** _sigRatioError;
    Float_t** _sigRiseTime;
    Float_t** _sigSlope;
    Float_t** _sigSlewRate;
    Float_t** _sigJitter;
    Float_t** _sigTOA;
    Float_t** _sigTOT;
    Float_t*** _sigCFD;
    Float_t** _sigNormCharge;
    Float_t** _sigNormAmplitude;
    
    
    //2D Array
    TH2F** _map;
    TCanvas** _canvasMap;
    
    //3D Array
    Float_t*** _zValue;
};

#endif
