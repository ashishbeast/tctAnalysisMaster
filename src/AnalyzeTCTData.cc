#include "ReadTCTFile.h"
#include "AnalyzeTCTData.h"

#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TVectorF.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

using namespace std;

AnalyzeTCTData::AnalyzeTCTData(const char* inFile)
{
  printf("Reading File: %s \n", inFile);
  //Set Parameters
  _deltaT = 0.05; //OScilloscope time bin: 50 ps
  _nAvOsc = 256;
  _termination = 50;
  _polarity = -1;
  _tmin = 0.;
  _tmax = 20.;
  _noise = 0;
  _nAC = 0;
  //Noise level in TCT setup is about 1.8mV, a threshold of 4 Sigma should be used
  _threshold = 7.2; // The maximum signal should be over 7.2 mV for a given data to consider the data consists a signal from the sensor
    
  _outFile += inFile;
  _outFile = _outFile.ReplaceAll(".tct", "");
  _outFile = _outFile.ReplaceAll("data/", "figures/");
    
  SetCanvasSettings(false);
    
  //Read the TCT File without any shift
  _tct = new ReadTCTFile(inFile, 0.);
    
  _tct->Print();
  _events = _tct->_totEvents;
  if(_tct->xCords[7][0] != 0) //Check Beam Monitor Status
    _bmON = true;

  GetActiveChannels();   //Get Active Channels
    
  if(FindWFShift())   //Find the shift in the waveforms and Correct the Baseline using tA
    {
      delete _tct;
      _tct = new ReadTCTFile(inFile, _tA); //Read file setting 0 at _tA
    }
  else
    {
      cout<<"[  ERROR] Start of the signal can not be determined!! \n";
      cout<<"[WARNING] Start of the signal set to default Value (0.)!! \n";
      cout<<"[ STATUS] Noise will be estimated for the full waveform.\n";
    }
    
  //Set DAQ paremeters
  _nX = _tct->nX;
  _nY = _tct->nY;
  _nZ = _tct->nZ;
  _nV1 = _tct->nV1;
  _nV2 = _tct->nV2;
  
  if(_isNoise)
    {
      _histo         = new TH1F**[_nAC];
      _sigNoise   = new Float_t*[_nAC];
      for(Int_t iCh = 0; iCh < _nAC; ++iCh)
        {
	  _histo[iCh]         = new TH1F*[_events];
	  _sigNoise[iCh]   = new Float_t[_events];
        }
    }
  else
    {
      //Allocate Memory for variables
      _histo         = new TH1F**[_nAC];
      _sigNoise   = new Float_t*[_nAC];
      _sigCharge     = new Float_t*[_nAC];
      _sigChargeError= new Float_t*[_nAC];
      _sigAmplitude  = new Float_t*[_nAC];
      _SNR = new Float_t*[_nAC];
      _sigAmplitudeError = new Float_t*[_nAC];
      _sigRatio         = new Float_t*[_nAC];
      _sigRatioError    = new Float_t*[_nAC];
      _sigRiseTime   = new Float_t*[_nAC];
      _sigSlope   = new Float_t*[_nAC];
      _sigSlewRate   = new Float_t*[_nAC];
      _sigJitter   = new Float_t*[_nAC];
      _sigTOA        = new Float_t*[_nAC];
      _sigTOT        = new Float_t*[_nAC];
      _sigCFD        = new Float_t**[_nAC];
      if(_bmON)
        {
	  _sigNormCharge = new Float_t*[_nAC];
	  _sigNormAmplitude   = new Float_t*[_nAC];
        }
        
      for(Int_t iCh = 0; iCh < _nAC; ++iCh)
        {
	  _histo[iCh]         = new TH1F*[_events];
	  _sigNoise[iCh]   = new Float_t[_events];
	  _sigCharge[iCh]     = new Float_t[_events];
	  _sigChargeError[iCh]= new Float_t[_events];
	  _sigAmplitude[iCh]  = new Float_t[_events];
	  _SNR[iCh]  = new Float_t[_events];
	  _sigAmplitudeError[iCh] = new Float_t[_events];
	  _sigRatio[iCh]         = new Float_t[_events];
	  _sigRatioError[iCh]    = new Float_t[_events];
	  _sigRiseTime[iCh]   = new Float_t[_events];
	  _sigSlope[iCh]   = new Float_t[_events];
	  _sigSlewRate[iCh]   = new Float_t[_events];
	  _sigJitter[iCh]   = new Float_t[_events];
	  _sigTOA[iCh]        = new Float_t[_events];
	  _sigTOT[iCh]        = new Float_t[_events];
	  _sigCFD[iCh]        = new Float_t*[_events];
	  for(Int_t evt = 0; evt < _events; ++evt)
	    _sigCFD[iCh][evt] = new Float_t[9];     // 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%
	  if(_bmON)
            {
	      if(iCh == 0)
		_bmValue       = new Float_t[_events];
	      _sigNormCharge[iCh] = new Float_t[_events];
	      _sigNormAmplitude[iCh]   = new Float_t[_events];
            }
        }
    }
    
  //Allocate memory block for Map types (charge, amplitude, riseTime, CFD, toA, ratio(amplitude/charge))
  _zValue = new Float_t**[7];
  _map = new TH2F*[7];
  _canvasMap = new TCanvas*[7];
}

AnalyzeTCTData::AnalyzeTCTData(const char* inFile, Float_t threshold)
{
  //Set Parameters
  _termination = 50;
  _polarity = -1;
  _tmin = 0.;
  _tmax = 20.;
  _noise = 0;
  _nAC = 0;
  //Noise level in TCT setup is about 1.8mV, a threshold of 4 Sigma should be used
  _threshold = threshold; // The maximum signal should be over 7.2 mV for a given data to consider the data consists a signal from the sensor
    
  _outFile += inFile;
  _outFile = _outFile.ReplaceAll(".tct", "");
  _outFile = _outFile.ReplaceAll("data/", "figures/");
    
  SetCanvasSettings(false);
    
  //Read the TCT File without any shift
  _tct = new ReadTCTFile(inFile, 0.);
    
  _tct->Print();
  _events = _tct->_totEvents;
  if(_tct->xCords[7][0] != 0) //Check Beam Monitor Status
    _bmON = true;

  GetActiveChannels();   //Get Active Channels
    
  if(FindWFShift())   //Find the shift in the waveforms and Correct the Baseline using tA
    {
      delete _tct;
      _tct = new ReadTCTFile(inFile, _tA); //Read file setting 0 at _tA
    }
  else
    {
      cout<<"[  ERROR] Start of the signal can not be determined!! \n";
      cout<<"[WARNING] Start of the signal set to default Value (0.)!! \n";
      cout<<"[ STATUS] Noise will be estimated for the full waveform.\n";
    }
    
  //Set DAQ paremeters
  _nX = _tct->nX;
  _nY = _tct->nY;
  _nZ = _tct->nZ;
  _nV1 = _tct->nV1;
  _nV2 = _tct->nV2;
  
  if(_isNoise)
    {
      _histo         = new TH1F**[_nAC];
      _sigNoise   = new Float_t*[_nAC];
      for(Int_t iCh = 0; iCh < _nAC; ++iCh)
        {
	  _histo[iCh]         = new TH1F*[_events];
	  _sigNoise[iCh]   = new Float_t[_events];
        }
    }
  else
    {
      //Allocate Memory for variables
      _histo         = new TH1F**[_nAC];
      _sigNoise   = new Float_t*[_nAC];
      _sigCharge     = new Float_t*[_nAC];
      _sigChargeError= new Float_t*[_nAC];
      _sigAmplitude  = new Float_t*[_nAC];
      _SNR = new Float_t*[_nAC];
      _sigAmplitudeError = new Float_t*[_nAC];
      _sigRatio         = new Float_t*[_nAC];
      _sigRatioError    = new Float_t*[_nAC];
      _sigRiseTime   = new Float_t*[_nAC];
      _sigSlope   = new Float_t*[_nAC];
      _sigSlewRate   = new Float_t*[_nAC];
      _sigJitter   = new Float_t*[_nAC];
      _sigTOA        = new Float_t*[_nAC];
      _sigTOT        = new Float_t*[_nAC];
      _sigCFD        = new Float_t**[_nAC];
      if(_bmON)
        {
	  _sigNormCharge = new Float_t*[_nAC];
	  _sigNormAmplitude   = new Float_t*[_nAC];
        }
        
      for(Int_t iCh = 0; iCh < _nAC; ++iCh)
        {
	  _histo[iCh]         = new TH1F*[_events];
	  _sigNoise[iCh]   = new Float_t[_events];
	  _sigCharge[iCh]     = new Float_t[_events];
	  _sigChargeError[iCh]= new Float_t[_events];
	  _sigAmplitude[iCh]  = new Float_t[_events];
	  _SNR[iCh]  = new Float_t[_events];
	  _sigAmplitudeError[iCh] = new Float_t[_events];
	  _sigRatio[iCh]         = new Float_t[_events];
	  _sigRatioError[iCh]    = new Float_t[_events];
	  _sigRiseTime[iCh]   = new Float_t[_events];
	  _sigSlope[iCh]   = new Float_t[_events];
	  _sigSlewRate[iCh]   = new Float_t[_events];
	  _sigJitter[iCh]   = new Float_t[_events];
	  _sigTOA[iCh]        = new Float_t[_events];
	  _sigTOT[iCh]        = new Float_t[_events];
	  _sigCFD[iCh]        = new Float_t*[_events];
	  for(Int_t evt = 0; evt < _events; ++evt)
	    _sigCFD[iCh][evt] = new Float_t[9];     // 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%
	  if(_bmON)
            {
	      if(iCh == 0)
		_bmValue       = new Float_t[_events];
	      _sigNormCharge[iCh] = new Float_t[_events];
	      _sigNormAmplitude[iCh]   = new Float_t[_events];
            }
        }
    }
    
  //Allocate memory block for Map types (charge, amplitude, riseTime, CFD, toA, ratio(amplitude/charge))
  _zValue = new Float_t**[7];
  _map = new TH2F*[7];
  _canvasMap = new TCanvas*[7];
}

AnalyzeTCTData::~AnalyzeTCTData()
{
  delete _tct;
  delete _aCH;
  delete[] _map;
  delete[] _canvasMap;
  delete[] _zValue;
  delete[] _histo;
    
  if(_isNoise)
    delete[] _sigNoise;
  else
    {
      delete[] _sigNoise;
      delete[] _sigCharge;
      delete[] _sigChargeError;
      delete[] _sigAmplitude;
      delete[] _sigAmplitudeError;
      delete[] _sigRatio;
      delete[] _sigRatioError;
      delete[] _sigRiseTime;
        
      delete[] _sigJitter;
      delete[] _sigSlope;
      delete[] _sigSlewRate;
      delete[] _sigTOA;
      delete[] _sigTOT;
      delete[] _sigCFD;
      if(_bmON)
        {
	  delete[] _bmValue;
	  delete[] _sigNormCharge;
	  delete[] _sigNormAmplitude;
        }
    }
}

void AnalyzeTCTData::SetPolarity(Int_t pol)
{
  _polarity = pol;
  return;
}

void AnalyzeTCTData::SetAveragesInOscilloscope(Int_t nAv)
{
  _nAvOsc = nAv;
  return;
}

void AnalyzeTCTData::SetTermination(Float_t impedence)
{
  _termination = impedence;
  return;
}

void AnalyzeTCTData::SetIntegralLimits(Float_t ti, Float_t tf)
{
  _tmin = ti;
  _tmax = tf;
  return;
}

void AnalyzeTCTData::GetActiveChannels()
{
  for(Int_t iCh = 0; iCh < _tct->_nCH; ++iCh)
    if(_tct->_chStatus[iCh])
      _nAC++;
    
  if(_nAC == 0)
    {
      cout<<"[  ERROR] No Active Channels Found!! \n";
      cout<<"[ STATUS] Exiting Program... \n";
      exit(0);
    }
  else
    printf("[MESSAGE] %d Active Channel Found \n", _nAC);
    
  _aCH = new Int_t [_nAC];
    
  Int_t countAC = 0;
  for(Int_t iCh = 0; iCh < _tct->_nCH; ++iCh)
    if(_tct->_chStatus[iCh])
      {
	_aCH[countAC] = iCh;
	countAC++;
      }
    
  printf("[MESSAGE] Active Channels: ");
  for(Int_t iCh = 0; iCh < _nAC; ++iCh)
    printf("CH%d \n", _aCH[iCh]);
}

Bool_t AnalyzeTCTData::FindWFShift()
{
  Float_t maxAmp=0, amplitude, time;
  TH1F *signal, *maxSignal;
  TH1F **his = new TH1F*[_tct->_nCH];
    
  for(Int_t i=0;i<_events; ++i)
    {
      if(i==0)
	printf("[ STATUS] Estimating signal start Time (%d Waveforms).... \n",_events);
        
      for(Int_t j=0; j<_tct->_nCH; ++j)
	if(_tct->_chStatus[j])
	  his[j] = ((TH1F*) _tct->_ch[j]->At(i));
        
      for(Int_t j=0; j<_tct->_nCH; ++j)
	if(_tct->_chStatus[j])
	  {
	    signal = (TH1F*)his[j]->Clone();
	    signal->Scale(_polarity);
	    if(signal->GetMaximum() > maxAmp)
	      {
		maxAmp = signal->GetMaximum();
		maxSignal = (TH1F*)signal->Clone();
		_indexMaxSignal = i;
		_chMaxSignal = j;
	      }
	  }
    }
    
  delete[] his;
    
  if(maxAmp < _threshold)
    {
      printf("[   INFO] No Signal Found!!\n");
      _tA = _tct->t0 + (_tct->nPoints * _tct->dt);
      _isNoise = true;
      return false;
    }
    
  _isNoise = false;
    
  for(Int_t i=0; i< _tct->nPoints; ++i)
    {
      amplitude = maxSignal->GetBinContent(i+1);
      time = _tct->t0 + i*_tct->dt;
      if(amplitude >= maxAmp*0.5)
        {
	  _tA = (time/1e-9)-1; //2ns before the 50% of Max Signal
	  printf("[   INFO] Start of the Signal ==> %0.2f ns\n", _tA);
	  return true;
        }
    }
  _tA = 0.;
  return false;
}

void AnalyzeTCTData::CorrectBaseline()
{
  // Function corrects the baseline (DC offset) of all wafeforms
  // Float_t tA ; time denoting the start of the pulse
  // correction factor is calculated from all the bins before tA
    
  printf("[ STATUS] Initiating Baseline correction (%d Waveforms).... \n",_events);
    
  Int_t right, left=1;
  Float_t corr;
    
  for(Int_t i=0; i<_nAC; ++i)
    for(Int_t j=0;j<_events; ++j)
      {
	_histo[i][j] = ((TH1F*) _tct->_ch[_aCH[i]]->At(j));
	//right = _histo[i][j]->GetXaxis()->FindBin(_tA);
	right = _histo[i][j]->GetXaxis()->FindBin(0.);
	//_histo[i][j]->Integral(left,right);
	corr = _histo[i][j]->Integral(left,right)/(right-left);
	for(Int_t k=1; k < _histo[i][j]->GetNbinsX(); ++k)
	  _histo[i][j]->SetBinContent(k, _histo[i][j]->GetBinContent(k)-corr);
      }
  printf("[MESSAGE] Baseline Correction Finished!!! \n");
}

void AnalyzeTCTData::CalcNoise()
{
  // Function calculates the noise in the baseline of all wafeforms
  printf("[ STATUS] Calculating the noise in the Baseline.... \n");
    
  Int_t N = 0;
  Float_t mean = 0.;
  Float_t sum = 0;
  Float_t sumD2 = 0.;
    
  Int_t right, left=1;
    
  TH1F* signal;

  for(Int_t i=0; i<_nAC; ++i)
    for(Int_t j=0;j<_events; ++j)
      {
	signal = (TH1F*) _histo[i][j]->Clone();
	signal->Scale(_polarity);
	right = signal->GetXaxis()->FindBin(0.);
	for(Int_t k = left; k< right; ++k)
	  {
	    sum += signal->GetBinContent(k);
	    N++;
	  }
      }
    
  mean = sum / N ;
    
  for(Int_t i=0; i<_nAC; ++i)
    for(Int_t j=0;j<_events; ++j)
      {
	signal = (TH1F*) _histo[i][j]->Clone();
	signal->Scale(_polarity);
	right = signal->GetXaxis()->FindBin(0.);
            
	for(Int_t k = left; k< right; ++k)
	  sumD2 += pow(signal->GetBinContent(k) - mean, 2);
      }
  Float_t mu2 = sumD2 / (N - 1); // second central moment (variance)
  Float_t stdDev = sqrt(mu2); //std Dev
  _noise = stdDev;
    
  if(_noise < 1)
    _isAveraged = true;
  else
    _isAveraged = false;
    
  if(_isAveraged)
    {
      _noise *= sqrt(_nAvOsc);
      printf("[MESSAGE] Noise Estimation Finished ==> %f mV!!! \n", _noise);
      return;
    }
    
  printf("[MESSAGE] Noise Estimation Finished ==> %f mV!!! \n", _noise);
    
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  _canvas = new TCanvas("","Noise", 1200, 1000);
  TH1F *noiseHist = new TH1F("noiseHist","",200,-10,10); //0.1 mV binWidth
  for(Int_t i=0; i<_nAC; ++i)
    for(Int_t j=0;j<_events; ++j)
      {
	signal = (TH1F*) _histo[i][j]->Clone();
	signal->Scale(_polarity);
	right = signal->GetXaxis()->FindBin(0.);
	for(Int_t k = left; k< right; ++k)
	  noiseHist->Fill(signal->GetBinContent(k));
      }
    
  Float_t meanNoise = noiseHist->GetMean();
  Float_t height = noiseHist->GetBinContent(noiseHist->GetMaximumBin());
  Float_t stdD = noiseHist->GetStdDev();
    
  TF1 *Gaus = new TF1("Gaus", "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", -10, 10);
  Gaus->SetLineColor(kBlue);
  Gaus->SetLineWidth(3);
  Gaus->SetParameters(height, meanNoise, stdD);
    
  noiseHist->Draw();
  noiseHist->Fit("Gaus", "QRMS0");
  Gaus->Draw("lsame");
  noiseHist->GetXaxis()->SetTitle("Noise (mV)");
  noiseHist->GetYaxis()->SetTitle("Counts");
  noiseHist->SetLineColor(kRed);
  noiseHist->SetLineWidth(3);
    
  Save(_canvas, "_NoiseHistogram", ".png");
  _canvas = NULL;
  delete _canvas;
  return;
}

void AnalyzeTCTData::CalculateWaveformProperties()
{
  // Function calculates the signal properties of all events
  printf("[ STATUS] Calculating Waveform Properties...\n");
  TH1F *signal;
    
  if(_isNoise)
    {
      for(Int_t i=0;i<_nAC; ++i)
	for(Int_t j=0; j<_events; ++j)
	  {
	    signal = (TH1F*) _histo[i][j]->Clone();
	    signal->Scale(_polarity);
	    _sigNoise[i][j]    = CalcSignalNoise(signal);
	  }
      return;
    }
    
  for(Int_t i=0;i<_nAC; ++i)
    for(Int_t j=0; j<_events; ++j)
      {
	signal = (TH1F*) _histo[i][j]->Clone();
	signal->Scale(_polarity);
	_sigNoise[i][j]    = CalcSignalNoise(signal);
	_sigCharge[i][j]    = CalcCharge(signal, &_sigChargeError[i][j]);
	_sigAmplitude[i][j] = CalcAmplitude(signal);
	_SNR[i][j] = _sigAmplitude[i][j]/_sigNoise[i][j];
	_sigAmplitudeError[i][j] = _sigNoise[i][j];
	if(_sigCharge[i][j] == 0.)
	  {
	    _sigRatio[i][j] = 0.;
	    _sigRatioError[i][j] = 0.;
	  }
	else
	  {
	    _sigRatio[i][j] = _sigAmplitude[i][j]/_sigCharge[i][j];
	    _sigRatioError[i][j] = _sigRatio[i][j]*TMath::Sqrt(TMath::Power(_sigChargeError[i][j]/_sigCharge[i][j], 2)+ TMath::Power(_sigAmplitudeError[i][j]/_sigAmplitude[i][j], 2));
	  }
	_sigRiseTime[i][j] = CalcRiseTime(signal);
	_sigSlewRate[i][j] = CalcSlewRate(signal);
            
	if(_sigSlewRate[i][j] != 0)
	  _sigJitter[i][j] = (_sigNoise[i][j]/_sigSlewRate[i][j])*1000; //ps
	else
	  _sigJitter[i][j] = 0; //ps

	//_sigJitter[i][j] = (_sigNoise[i][j]/_sigSlope[i][j])*1000; //ps
	//_sigJitter[i][j] = 1000 * _sigNoise[i][j] * _sigRiseTime[i][j]/(0.8 * _sigAmplitude[i][j]);
            
	_sigTOA[i][j] = CalcCFD(signal, 30);
	_sigTOT[i][j] = CalcTOT(signal, 30);
	for(Int_t k =0; k<9; ++k)
	  _sigCFD[i][j][k] = CalcCFD(signal, 10*(k+1));
	if(_bmON)
	  {
	    _bmValue[j] = (_tct->xCords[7][j])/3000;
	    _sigNormCharge[i][j] = _sigCharge[i][j]/_bmValue[j];
	    _sigNormAmplitude[i][j] = _sigAmplitude[i][j]/_bmValue[j];
	  }
      }
  printf("[MESSAGE] Finished Calculating Signal Properties!!! \n");
}


Int_t AnalyzeTCTData::GetOperation(Int_t *binary)
{
  Int_t sum=0;
  for(Int_t i=0; i<5; ++i)
    sum += binary[i]*TMath::Power(2,i);
  return sum;
}

void AnalyzeTCTData::AnalysisAction()
{  
  Bool_t xON = false, yON = false, zON = false, v1ON = false, v2ON = false;
    
  if(_tct->nX > 1)     xON = true;
  if(_tct->nY > 1)     yON = true;
  if(_tct->nZ > 1)     zON = true;
  if(_tct->nV1 > 1)    v1ON = true;
  if(_tct->nV2 > 1)    v2ON = true;
    
  Int_t binary[5] = {xON, yON, zON, v1ON, v2ON};
    
  Int_t operation = GetOperation(binary);
    
  /*
    operation = 1 -> X Sweep
    operation = 2 -> Y Sweep
    operation = 3 -> XY Map
    operation = 4 -> Z Sweep
    operation = 5 -> (XZ) Focus
    operation = 6 -> (YZ) Focus
    operation = 8 -> V1 Sweep
    operation = 11 -> XY Map Evolution with Voltage
    operation = 16 -> V2 Sweep
  */
    
  switch(operation)
    {
    case 1:      PlotAxisSweep(operation);      break;
    case 2:      PlotAxisSweep(operation);      break;
    case 3:      PlotMaps(0);                   break;
    case 4:      PlotAxisSweep(operation);      break;
    case 5:      GetFocus();                    break;
    case 6:      GetFocus();                    break;
    case 8:      PlotAxisSweep(operation);      break;
    case 11:
      printf("[404 ERROR] Working on it!!\n");
      break;
    case 16:     PlotAxisSweep(operation);      break;
    default:
      printf("[  ERROR] Could not Resolve an Operation to perform...\n");
      printf("[WARNING] Check the File. Exiting Program!!\n");
      exit(0);
    }
  return;
}

void AnalyzeTCTData::SaveSignalShape()
{
  //TH1F *his =  _tct->GetHisto(_chMaxSignal, _indexMaxSignal);
    
  TH1F *maxSignal = (TH1F*) _histo[0][_indexMaxSignal]->Clone();
  maxSignal->Scale(_polarity);
    
  _canvas = new TCanvas("","Signal", 1200, 1000);
  TH1F *signal;
  for(Int_t j=0;j<_events; ++j)
    {
      //if(j!= _events-1)
      //    continue;
      signal = (TH1F*) _histo[0][j]->Clone();
      signal->Scale(_polarity);
      signal->SetLineColor(j%49+1);
      //signal->SetLineColor(kBlue);
      signal->SetLineWidth(2);
      //if(j==_events-1)
      if(j==0)
        {
	  signal->Draw("HIST L");
	  signal->GetXaxis()->SetTitle("Time (ns)");
	  signal->GetYaxis()->SetTitle("Voltage (mV)");
	  signal->GetYaxis()->SetRangeUser(maxSignal->GetMinimum()-5, maxSignal->GetMaximum()+5);
	  //signal->GetXaxis()->SetRangeUser(0, 10);
        }
      else
        {
	  signal->Draw("HIST LSAME");
	  signal->GetYaxis()->SetRangeUser(maxSignal->GetMinimum()-5, maxSignal->GetMaximum()+5);
	  //signal->GetXaxis()->SetRangeUser(0, 10);
        }
    }
  cout<<"[ STATUS] Saving Signal...\n";
  Save(_canvas, "_Signal", ".png");
    
  _canvas = NULL;
  delete _canvas;
  return;
}

Float_t AnalyzeTCTData::CalcSignalNoise(TH1F* his)
{
  Int_t N = 0;
  Float_t mean = 0.;
  Float_t sum = 0;
  Float_t sumD2 = 0.;
    
  Int_t right, left=1;
    
  right = his->GetXaxis()->FindBin(0.);
            
  for(Int_t k = left; k< right; ++k)
    {
      sum += his->GetBinContent(k);
      N++;
    }
    
  mean = sum / N ;
    
  for(Int_t k = left; k< right; ++k)
    sumD2 += pow(his->GetBinContent(k) - mean, 2);
    
  Float_t mu2 = sumD2 / (N - 1); // second central moment (variance)
  Float_t stdDev = sqrt(mu2); //std Dev

  if(_isAveraged)
    stdDev *= sqrt(_nAvOsc);

  return stdDev;
}

Float_t AnalyzeTCTData::CalcAmplitude(TH1F *his)
{
    
  return his->GetBinContent(his->GetMaximumBin());
}

Float_t AnalyzeTCTData::CalcCharge(TH1F *his, Float_t* error)
{
  Double_t err;
  Float_t charge;
    
  Int_t minTime = his->GetXaxis()->FindBin(_tmin);
  Int_t maxTime = his->GetXaxis()->FindBin(_tmax);
    
  charge = his->IntegralAndError(minTime, maxTime, err, "width");
  charge /= _termination;
  *error = (err/_termination);
  //if(charge < 0)
  //    return 0;
  return charge;
}

Float_t AnalyzeTCTData::CalcCharge(Int_t channel, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2, Float_t *error)
{
  TH1F *his, *signal;
  his = _tct->GetHisto(channel, x, y, z, v1, v2);
  signal = (TH1F*)his->Clone();
  signal->Scale(_polarity);
    
  Float_t err;
  Float_t charge = CalcCharge(signal, &err);
  *error = err;
  return charge;
}

Float_t AnalyzeTCTData::CalcRatio(TH1F* his, Float_t* error)
{
  Float_t errC, ratio;
  Float_t charge = CalcCharge(his, &errC);
  Float_t amplitude = CalcAmplitude(his);
  Float_t errA = CalcSignalNoise(his);
  if(charge==0)
    {
      *error = 0;
      return 0;
    }
  ratio = amplitude/charge;
  *error = ratio*TMath::Sqrt(TMath::Power(errC/charge, 2)+ TMath::Power(errA/amplitude, 2));
  return ratio;
}

Float_t AnalyzeTCTData::CalcRiseTime(TH1F* his)
{
  Float_t maxAmp = his->GetMaximum();
  Float_t t1, t2; //Used for RiseTime
    
  vector<Float_t> xData;
  vector<Float_t> yData;
    
  for(Int_t i=0; i<_tct->nPoints; ++i)
    {
      xData.push_back(his->GetXaxis()->GetBinCenter(i+1));
      yData.push_back(his->GetBinContent(i+1));
    }
    
  t1 = LinearInterpolation(xData, yData, 0.1*maxAmp);
  t2 = LinearInterpolation(xData, yData, 0.9*maxAmp);
    
  if((t2-t1) < 0)
    return 0;
    
  xData.clear();
  yData.clear();
    
  return t2-t1;
}

//Float_t AnalyzeTCTData::CalcSlope(TH1F* his, Float_t threshold) //Calculate slope of the leading edge
//{
//    Float_t sigMax = his->GetMaximum(); //Maximum Signal
//    Float_t thres1 = sigMax * (threshold-2)/100;
//    Float_t thres2 = sigMax * (threshold+2)/100;
//
//    Float_t slope, t1, t2;
//
//    vector<Float_t> xData;
//    vector<Float_t> yData;
//
//    for(Int_t i=0; i<_tct->nPoints; ++i)
//    {
//        xData.push_back(his->GetXaxis()->GetBinCenter(i+1));
//        yData.push_back(his->GetBinContent(i+1));
//    }
//    t1 = LinearInterpolation(xData, yData, thres1);
//    t2 = LinearInterpolation(xData, yData, thres2);
//
//    slope = (thres2-thres1)/(t2-t1);
//    return slope;
//}

//Float_t AnalyzeTCTData::CalcSlope(TH1F* his) //Calculate slope of the leading edge
//{
//    Float_t threshold = 40;
//    Float_t sigMax = his->GetMaximum(); //Maximum Signal
//    Float_t thres1 = sigMax * threshold/100;
//    Float_t thres2 = sigMax * (threshold-1)/100;
//
//    Float_t slope, t1, t2;
//
//    vector<Float_t> xData;
//    vector<Float_t> yData;
//
//    for(Int_t i=0; i<_tct->nPoints; ++i)
//    {
//        xData.push_back(his->GetXaxis()->GetBinCenter(i+1));
//        yData.push_back(his->GetBinContent(i+1));
//    }
//    t1 = LinearInterpolation(xData, yData, thres1);
//    t2 = LinearInterpolation(xData, yData, thres2);
//
//    slope = (thres2-thres1)/(t2-t1);
//    return slope;
//}

Float_t AnalyzeTCTData::CalcSlope(TH1F* his, Float_t threshold) //Calculate slope of the leading edge
{
  Float_t sigMax = his->GetMaximum(); //Maximum Signal
  Float_t thres = sigMax * threshold/100;
  Int_t cross = 0;
    
  for(Int_t i=1; i<_tct->nPoints; ++i)
    {
      if(his->GetBinContent(i) > thres)
        {
	  cross = i;
	  break;
        }
    }
    
  if(cross == 1)
    return 0;
    
  Float_t y1 = his->GetBinContent(cross);
  Float_t y2 = his->GetBinContent(cross-1);
  Float_t x1 = his->GetXaxis()->GetBinCenter(cross);
  Float_t x2 = his->GetXaxis()->GetBinCenter(cross-1);
    
  Float_t slope = (y2-y1)/(x2-x1);
  return slope;
}

Float_t AnalyzeTCTData::CalcSlope(TH1F* his) //Calculate slope of the leading edge
{
  Float_t threshold = 40;
  Float_t sigMax = his->GetMaximum(); //Maximum Signal
  Float_t thres = sigMax * threshold/100;
  Int_t cross = 0;
    
  for(Int_t i=1; i<_tct->nPoints; ++i)
    {
      if(his->GetBinContent(i) > thres)
        {
	  cross = i;
	  break;
        }
    }
    
  if(cross == 1)
    return 0;
    
  Float_t y1 = his->GetBinContent(cross);
  Float_t y2 = his->GetBinContent(cross-1);
  Float_t x1 = his->GetXaxis()->GetBinCenter(cross);
  Float_t x2 = his->GetXaxis()->GetBinCenter(cross-1);
    
  Float_t slope = (y2-y1)/(x2-x1);
  return slope;
}

Float_t AnalyzeTCTData::CalcSlewRate(TH1F* his) //Calculate the Max Slope in leading edge of signal
{
  Float_t slewRate = 0, slope = 0;
  Float_t thrLow = 20, thrHigh = 80;
    
  for(Int_t i=thrLow; i<=thrHigh; i+=10)
    {
      slope = CalcSlope(his, i);
      if(slope >= slewRate)
	slewRate = slope;
    }
  return slewRate;
}

Float_t AnalyzeTCTData::CalcCFD(TH1F *his, Float_t thr)
{
  //Check Threshold
  if(thr > 1)
    thr /= 100;
    
  Float_t maxAmp = his->GetMaximum();
  Float_t t;
    
  vector<Float_t> xData;
  vector<Float_t> yData;
    
  for(Int_t i=0; i<_tct->nPoints; ++i)
    {
      xData.push_back(his->GetXaxis()->GetBinCenter(i+1));
      yData.push_back(his->GetBinContent(i+1));
    }
    
  t = LinearInterpolation(xData, yData, thr*maxAmp);
    
  if(t<0 || t>5)
    return 0;
    
  xData.clear();
  yData.clear();
  return t;
}

Float_t AnalyzeTCTData::CalcTOA(TH1F *his)
{
  Float_t t = CalcCFD(his, 20);
    
  if(t<0 || t>2)
    return 0;
  return t;
}

Float_t AnalyzeTCTData::CalcTOT(TH1F *his, Float_t thr)
{
  //Check Threshold
  if(thr > 1)
    thr /= 100;
    
  Float_t t1 = CalcCFD(his, thr);
  Float_t t2;
    
  vector<Float_t> xData;
  vector<Float_t> yData;
    
  for(Int_t i=0; i<_tct->nPoints; ++i)
    {
      xData.push_back(his->GetXaxis()->GetBinCenter(i+1));
      yData.push_back(his->GetBinContent(i+1));
    }
    
  //search -> no point above the threshold found
  //check -> found point above threshold and checking for TOT
  enum state {search, check};
    
  state status = search;
    
  unsigned int abovePos = 0; //Index of the Point above Threshold
    
  for(unsigned int i = 0; i< yData.size(); ++i)
    {
      if(status == search)
        {
	  if(yData[i] > thr)
            {
	      abovePos = i;
	      status = check;
            }
        }
        
      if(status == check)
        {
	  if(yData[i] >  thr)
	    continue;
	  else
	    status = search;
            
	  if(yData[i] < thr)
	    break;
        }
    }
    
  if(abovePos == 0 || abovePos == yData.size()-1)
    return t1;
  Float_t y1 = yData[abovePos];
  Float_t y2 = yData[abovePos - 1];
  Float_t x1 = xData[abovePos];
  Float_t x2 = xData[abovePos - 1];
  xData.clear();
  yData.clear();
    
  Float_t a = (y2-y1)/(x2-x1);
  Float_t b = y1 - a*x1;
  t2 = (thr - b)/a;
    
  if((t2 - t1) < 0 || (t2 - t1) > 5)
    return 0;
  return t2-t1;
}

Float_t AnalyzeTCTData::LinearInterpolation(vector<Float_t> x, vector<Float_t> y, Float_t thr)
{
  if(y.size() != x.size() || y.size() == 0)
    return 100;
    
  UInt_t i = 0;
  for(; i < y.size(); ++i)
    if(y[i] > thr)
      break;
  if(i == 0)
    return x[0];
    
  //linear Interpolatin between points and above Threshold
  //y = a x + b
  Float_t y1 = y[i];
  Float_t y2 = y[i-1];
  Float_t x1 = x[i];
  Float_t x2 = x[i-1];
    
  Float_t a = (y2-y1)/(x2-x1);
  Float_t b = y1 - a*x1;
    
  return (thr - b)/a;
}

void AnalyzeTCTData::PlotAxisSweep(Int_t operation)
{
  cout<<"[ STATUS] Plotting Axis(X, Y, X, V1, V2) Sweep...\n";
    
  Int_t i, v=0, *s[5];
  Int_t pX, dpX;
  TString xTitle;
  TString yTitle[8] = { "Noise (mV)", "Charge (fC)", "Amplitude (mV)", "Rise Time (ns)", "CFD_{50%} (ns)", "Jitter (ps)", "Slew Rate(mV/ns)", "SNR"};
  //TString mode[8] = {"Noise", "Charge", "Amplitude", "RiseTime", "CFD", "Jitter", "SlewRate", "SNR"};
  const char* mode[8] = {"Noise", "Charge", "Amplitude", "RiseTime", "CFD", "Jitter", "SlewRate", "SNR"};
    
  switch(operation)
    {
    case 1://X
      pX = _tct->nX-1, dpX = TMath::Abs(_tct->stepX);
      s[0] = &i, s[1] = &v, s[2] = &v, s[3] = &v, s[4] = &v;
      xTitle = "X (#mum)";
      break;
    case 2://Y
      pX = _tct->nY-1, dpX = TMath::Abs(_tct->stepY);
      s[0] = &v, s[1] = &i, s[2] = &v, s[3] = &v, s[4] = &v;
      xTitle = "Y (#mum)";
      break;
    case 4://Z
      pX = _tct->nZ-1, dpX = TMath::Abs(_tct->stepZ);
      s[0] = &v, s[1] = &v, s[2] = &i, s[3] = &v, s[4] = &v;
      xTitle = "Z (#mum)";
      break;
    case 8://V1
      pX = _tct->nV1-1, dpX = TMath::Abs(_tct->V1[1]-_tct->V1[0]);
      s[0] = &v, s[1] = &v, s[2] = &v, s[3] = &i, s[4] = &v;
      xTitle = "V_{bias} (V)";
      break;
    case 16://V2
      pX = _tct->nV2-1, dpX = TMath::Abs(_tct->V2[1]-_tct->V2[0]);
      s[0] = &v, s[1] = &v, s[2] = &v, s[3] = &v, s[4] = &i;
      xTitle = "V_{bias} (V)";
      break;
    }
    
  Float_t **yAxis =  new Float_t*[8];
  Float_t **yAxisError =  new Float_t*[8];
    
  for(Int_t k =0;k < 8;++k)
    {
      yAxis[k] = new Float_t[pX];
      yAxisError[k] = new Float_t[pX];
    }
  TGraphErrors *gr[8];
  Float_t steps[pX];
    
  for(Int_t iCh = 0; iCh<_nAC; ++iCh)
    {
      for(i=0; i< pX; ++i)
        {
	  _index = _tct->GetIndex(*s[0], *s[1], *s[2], *s[3], *s[4]);
            
	  if(_isNoise)
            {
	      yAxis[0][i] = _sigNoise[iCh][_index];
	      yAxisError[0][i] = 0;
            }
	  else
            {
	      yAxis[0][i] = _sigNoise[iCh][_index];
	      yAxis[1][i] = _sigCharge[iCh][_index];
	      yAxis[2][i] = _sigAmplitude[iCh][_index];
	      yAxis[3][i] = _sigRiseTime[iCh][_index];
	      yAxis[4][i] = _sigCFD[iCh][_index][5];
	      yAxis[5][i] = _sigJitter[iCh][_index];
	      yAxis[6][i] = _sigSlewRate[iCh][_index];
	      yAxis[7][i] = _SNR[iCh][_index];
                
	      yAxisError[0][i] = 0;
	      yAxisError[1][i] = _sigChargeError[iCh][_index];
	      yAxisError[2][i] = _sigAmplitudeError[iCh][_index];
	      yAxisError[3][i] = 0;
	      yAxisError[4][i] = 0;
	      yAxisError[5][i] = 0;
	      yAxisError[6][i] = 0;
	      yAxisError[7][i] = 0;
                
            }
	  steps[i] = dpX*i;
        }
        
      // TGraphErrors *gr = new TGraphErrors(pX, steps, charge, 0, chargeError);
      // TGraphErrors *hr = new TGraphErrors(pX, steps, amplitude, 0, ampError);
        
      Int_t graphNumbers;
      if(_isNoise)
	graphNumbers = 1;
      else
	graphNumbers = 8;
        
      for(Int_t k =0; k< graphNumbers;++k)
        {
	  _canvas = new TCanvas("","Sweep Measurement", 1200, 1000);
	  gr[k] = new TGraphErrors(pX, steps, &yAxis[k][0], 0, &yAxisError[k][0]);
	  gr[k]->Draw("apl");
	  gr[k]->SetMarkerSize(1.5);
	  gr[k]->SetMarkerStyle(20);
	  gr[k]->SetMarkerColor(kBlack);
	  gr[k]->SetLineColor(kRed);
	  gr[k]->GetXaxis()->SetTitle(xTitle);
	  gr[k]->GetYaxis()->SetTitle(yTitle[k]);
	  if(k==3)
	    gr[k]->GetYaxis()->SetRangeUser(0,5);
	  //Save(_canvas, Form("_SweepMode%d_CH%d", k, _aCH[iCh]), ".png" );
	  Save(_canvas, Form("_CH%d_%s", _aCH[iCh], mode[k]), ".png" );
	  _canvas = NULL;
	  delete _canvas;
        }
    }
  cout<<"[MESSAGE] Finished Plotting Axis(X, Y, X, V1, V2) Sweep!!!\n";
  return;
}

void AnalyzeTCTData::PlotMaps(Int_t mapType)
{
  cout<<"[ STATUS] Plotting Position Sensitive Maps...\n";
  SetCanvasSettings(true);
  Int_t pX, pY, dpX, dpY, v=0;
  TString xTitle, yTitle;
  Int_t *s[5], i, j;
    
  TString zTitle[7] = {"Charge (fC)", "Amplitude (mV)", "Rise Time (ns)", "TOA_{50%} (ns)", "CFD_{30%} (ns)", "Ratio (Amplitude/Charge)", "TOT_{30%} (ns)"};
  const char* msg[7] = {"Charge", "Amplitude", "RiseTime", "TOA", "CFD", "Ratio", "TOT"};
    
  switch(mapType)
    {
    case 0://XY
      pX = _tct->nX, pY = _tct->nY, TMath::Abs(dpX = _tct->stepX), dpY = TMath::Abs(_tct->stepY);
      s[0] = &i, s[1] = &j, s[2] = &v, s[3] = &v, s[4] = &v;
      xTitle = "X (#mum)";
      yTitle = "Y (#mum)";
      break;
    case 1://XZ (LASER PROFILE)
      pX = _tct->nX, pY = _tct->nZ, TMath::Abs(dpX = _tct->stepX), dpY = TMath::Abs(_tct->stepZ);
      s[0] = &i, s[1] = &v, s[2] = &j, s[3] = &v, s[4] = &v;
      xTitle = "X (#mum)";
      yTitle = "Z (#mum)";
      break;
    case 2://XV1
      pX = _tct->nX, pY = _tct->nV1, dpX =TMath::Abs(_tct->stepX), dpY = TMath::Abs(_tct->V1[1]-_tct->V1[0]);
      s[0] = &i, s[1] = &v, s[2] = &v, s[3] = &j, s[4] = &v;
      xTitle = "X (#mum)";
      yTitle = "V_{bias} (V)";
      break;
    case 3://YZ (LASER PROFILE)
      pX = _tct->nY, pY = _tct->nZ, dpX = TMath::Abs(_tct->stepY), dpY = TMath::Abs(_tct->stepZ);
      s[0] = &v, s[1] = &i, s[2] = &j, s[3] = &v, s[4] = &v;
      xTitle = "Y (#mum)";
      yTitle = "Z (#mum)";
      break;
    case 4://YV1
      pX = _tct->nY, pY = _tct->nV1, dpX = TMath::Abs(_tct->stepY), dpY = TMath::Abs(_tct->V1[1]-_tct->V1[0]);
      s[0] = &v, s[1] = &i, s[2] = &v, s[3] = &j, s[4] = &v;
      xTitle = "Y (#mum)";
      yTitle = "V_{bias} (V)";
      break;
    case 5://ZV1
      pX = _tct->nZ,pY = _tct->nV1, dpX = TMath::Abs(_tct->stepZ), dpY = TMath::Abs(_tct->V1[1]-_tct->V1[0]);
      s[0] = &v, s[1] = &v, s[2] = &i, s[3] = &j, s[4] = &v;
      xTitle = "Z (#mum)";
      yTitle = "V_{bias} (V)";
      break;
    case 6://YX
      pX = _tct->nY, pY = _tct->nX, dpX = TMath::Abs(_tct->stepY), dpY = TMath::Abs(_tct->stepX);
      s[0] = &j, s[1] = &i, s[2] = &v, s[3] = &v, s[4] = &v;
      xTitle = "Y (#mum)";
      yTitle = "X (#mum)";
      break;
    }
    
  //Allocate memory
  for(Int_t mode=0; mode<7; ++mode)
    {
      _zValue[mode] = new Float_t*[pY];
        
      _map[mode] = new TH2F("","Position Map",pX,0,pX*dpX, pY,0,pY*dpY);
      _map[mode]->GetXaxis()->SetTitle(xTitle);
      _map[mode]->GetYaxis()->SetTitle(yTitle);
      _map[mode]->GetZaxis()->SetTitle(zTitle[mode]);
        
      _canvasMap[mode] = new TCanvas("", "Position Sensitive Map", 1200, 1000);
        
      for(j=0; j< pY; ++j)
	_zValue[mode][j] = new Float_t[pX];
    }
    
  for(Int_t iCh = 0; iCh<_nAC; ++iCh)
    {
      cout<<"[ STATUS] Processing CH "<<_aCH[iCh]<<std::flush;
      for(j = 0; j<pY; ++j)
	for(i = 0; i<pX; ++i)
	  {
	    _index = _tct->GetIndex(*s[0], *s[1], *s[2], *s[3], *s[4]);
	    _zValue[0][j][i] = _sigCharge[iCh][_index];
	    _zValue[1][j][i] = _sigAmplitude[iCh][_index];
	    _zValue[2][j][i] = _sigRiseTime[iCh][_index];
	    _zValue[3][j][i] = _sigCFD[iCh][_index][5];
	    _zValue[4][j][i] = _sigCFD[iCh][_index][3];
	    _zValue[5][j][i] = _sigRatio[iCh][_index];
	    _zValue[6][j][i] = _sigTOT[iCh][_index];
	  }
        
      for(Int_t mode=0; mode<2; ++mode)
        {
	  for(j = 0; j<pY; ++j)
	    for(i = 0; i<pX; ++i)
	      _map[mode]->Fill(i*dpX, j*dpY, _zValue[mode][j][i]);
            
	  _canvasMap[mode]->cd();
	  _map[mode]->Draw("COLZ");
	  Save(_canvasMap[mode], Form("_Map_CH%d_%s", _aCH[iCh], msg[mode]), ".png" );
	  cout<<"\r[ STATUS] Finished: "<<msg[mode]<<" Map!!   "<< std::flush;
        }
    }
  cout<<"\n[MESSAGE] Finished Plotting Position Sensitive Maps!!!\n";
  return;
}

void AnalyzeTCTData::GetFocus()
{
  cout<<"[ STATUS] Finding Focus of the LASER...\n";
  int pS, stepS, pO = _tct->nZ, stepO = TMath::Abs(_tct->stepZ);
  Float_t optDis[pO];
  Float_t fwhm[pO];
  Float_t fwhmErr[pO];
  TString xTitle;
  Int_t *s[5], i, j, v=0;
  //Int_t operation;
    
  if(_tct->nX > 1)
    {
      pS = _tct->nX, stepS = TMath::Abs(_tct->stepX);
      s[0] = &j, s[1] = &v, s[2] = &i, s[3] = &v, s[4] = &v;
      xTitle = "X (#mum)";
      //operation = 1;
    }
  else
    {
      pS = _tct->nY, stepS = TMath::Abs(_tct->stepY);
      s[0] = &v, s[1] = &j, s[2] = &i, s[3] = &v, s[4] = &v;
      xTitle = "Y (#mum)";
      //operation = 3;
    }
    
  Float_t charge[pO][pS];
  Float_t xData[pO][pS];
  Float_t chargeErr[pO][pS];
    
  //    Float_t maxCharge, minCharge = 100;
  //    for(i=0; i<pO; ++i)
  //        for(j=0; j<pS; ++j)
  //        {
  //            _index = _tct->GetIndex(*s[0], *s[1], *s[2], *s[3], *s[4]);
  //            if(_sigCharge[0][_index] > maxCharge)
  //                maxCharge = _sigCharge[0][_index];
  //            if(_sigCharge[0][_index] < minCharge)
  //                minCharge = _sigCharge[0][_index];
  //        }
  //
  Float_t xLow = 0;         // low limit for fit
  Float_t xHigh = pS*stepS; // high limit for fit
  Float_t FWHM = 10;        // expected FWHM
  Float_t mid = 0.5*(xHigh-xLow);
    
  /*-------------------------------------------------------------------------------
    Fiting of the S-curve is done by error function which is defined as:
    f(x) = erf((x-x0)/r)*M+y0
    (x0,y0) are the coordinate of the middle point in the S-curve
    y0 = Ymax/2
    r is the radius
    M = Ymax*,  is the maximum value of the S-curve (Saturation Value)
    *if curve start from 0 and reaches a maximum value then M = Ymax/2
    _______________________________________________________________________________*/
    
  TF1 *sCurve =  new TF1("sCurve", "TMath::Erf((x-[0])/[1])*[2]+[3]", xLow, xHigh);
  sCurve->SetParNames("L", "Spot Size", "const1", "const2");
  sCurve->SetParameters(mid, FWHM, 8, 4);
  sCurve->SetLineColor(kRed);
    
  TGraphErrors *gr[pO];
  _canvas = new TCanvas("", "Focus Measurement", 1800, 800);
  _canvas->Divide(2,1);
    
  _canvas->cd(1);
  Float_t maxCharge, minCharge;
  for(i=0; i<pO; ++i)
    {
      for(j=0; j<pS; ++j)
        {
	  xData[i][j] = j*stepS;
	  _index = _tct->GetIndex(*s[0], *s[1], *s[2], *s[3], *s[4]);
	  charge[i][j] = _sigCharge[0][_index];
	  chargeErr[i][j] = _sigChargeError[0][_index];
        }
      gr[i] = new TGraphErrors(pS, &xData[i][0], &charge[i][0], 0, &chargeErr[i][0]);
      gr[i]->SetMarkerStyle(21);
      gr[i]->SetMarkerSize(1);
      gr[i]->Fit(sCurve,"QMR");
        
      if(i==0)
        {
	  gr[i]->Draw("AP");
	  gr[i]->GetXaxis()->SetTitle("Scanning distance [#mum]");
	  gr[i]->GetYaxis()->SetTitle("Charge [fC]");
        }
      else
	gr[i]->Draw("PSAME");
        
      maxCharge = *max_element(charge[i], charge[i]+pS);
      minCharge = *min_element(charge[i], charge[i]+pS);
        
      gr[i]->GetYaxis()->SetRangeUser(minCharge-2, maxCharge+2);
        
      fwhm[i]= sCurve->GetParameter(1)*2.35/TMath::Sqrt(2);
      fwhmErr[i] = sCurve->GetParError(1)*2.35/TMath::Sqrt(2);
      optDis[i] = i*stepO;
    }
    
  //This part fits the FWHM graph and Calculates the focus spot
  _canvas->cd(2);
    
  //Function to fit FWHM
  Float_t min = *min_element(fwhm, fwhm+pO);
  Int_t minIndex = find(fwhm, fwhm+pO, min) - fwhm;
    
  TF1 *fP = new TF1("fP","[2]*x*x+[1]*x+[0]",optDis[minIndex]-1000,optDis[minIndex]+1000);
  fP->SetLineColor(kRed);
    
  TGraphErrors *grFWHM = new TGraphErrors(pO, optDis, fwhm, 0, fwhmErr);
  grFWHM->SetMarkerStyle(21);
  grFWHM->SetMarkerSize(1);
  grFWHM->SetMarkerColor(kBlue);
  grFWHM->SetLineColor(kBlue);
    
  grFWHM->GetHistogram()->GetXaxis()->SetTitle("Optical distance [#mum]");
  grFWHM->GetHistogram()->GetYaxis()->SetTitle("FWHM [#mum]");
  grFWHM->Draw("AP");
  grFWHM->Fit("fP","QRMS+");
    
  Float_t foSpot = -1 * fP->GetParameter(1)/(2*fP->GetParameter(2));
  Float_t SpotSize = fP->GetParameter(2)*TMath::Power(foSpot,2)+fP->GetParameter(1)*foSpot+fP->GetParameter(0);
    
  TString legend1 = Form("Z_{focus} = %0.0f #mum ", foSpot);
  TString legend2 = Form("Spot Size = %0.2f #mum ", SpotSize);
    
  TPaveText *pt = new TPaveText(0.345576,0.7510684,0.7103506,0.8237179,"nbNDC");
  pt->SetTextColor(kRed);
  pt->SetFillColor(0);
  pt->AddText(legend1);
  pt->AddText(legend2);
  pt->Draw("");
    
  Save(_canvas, Form("_Focus_CH%d", _aCH[0]), ".png" );
  _canvas = NULL;
  delete _canvas;
  //PlotMaps(operation);
  cout<<"[MESSAGE] Finished finding focus of the LASER!!!\n";
  return;
}

void AnalyzeTCTData::Save(TCanvas* canvas, TString canvasType, TString ext)
{
  TString fileName = _outFile+canvasType+ext;
  canvas->SaveAs(fileName);
}

void AnalyzeTCTData::SetCanvasSettings(Bool_t threeD = false)
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



/*
  Bool_t AnalyzeTCTData::FindWFShift()
  {
  printf("[ STATUS] Finding signal start Time (%d Waveforms).... \n",_tct->_totEvents);
 
  Float_t maxAmp=0, amplitude, time;
  TH1F *signal, *maxSignal;
 
  for(Int_t i=0;i<_tct->_totEvents; ++i)
  {
  for(Int_t j=0; j<_nAC; ++j)
  {
  _histo[j] = ((TH1F*) _tct->_ch[_aCH[j]]->At(i));
 
  signal = (TH1F*)_histo[j]->Clone();
  signal->Scale(_polarity);
 
  if(signal->GetMaximum() > maxAmp)
  {
  maxAmp = signal->GetMaximum();
  maxSignal = (TH1F*)signal->Clone();
  _indexMaxSignal = i;
  _chMaxSignal = j;
  }
  }
  }
 
  for(Int_t i=0; i< _tct->nPoints; ++i)
  {
  amplitude = maxSignal->GetBinContent(i+1);
  time = _tct->t0 + i*_tct->dt;
  if(amplitude > maxAmp*0.2)
  {
  _tA = (time/1e-9)-2; //2ns before the 20% of Max Signal
  printf("[   INFO] : Start of the Signal at %f ns\n", _tA);
  return true;
  }
  }
  _tA = 0.;
  return false;
  }
*/

/*
  void AnalyzeTCTData::CorrectBaseline()
  {
 
  // Function corrects the baseline (DC offset) of all wafeforms
  // Float_T_t tA ; time denoting the start of the pulse
  // correction factor is calculated from all the bins before tA
 
 
  Int_t *right = new Int_t [_tct->_nCH];
  Int_t *left = new Int_t [_tct->_nCH];
  Int_t *corr = new Int_t [_tct->_nCH];
  TH1F **his = new TH1F*[_tct->_nCH];
 
  for(Int_t i=0; i<_tct->_totEvents; ++i)
  {
  if(i==0)
  printf("[ STATUS] Initiating Baseline correction (%d Waveforms).... \n",_tct->_totEvents);
 
  for(Int_t j = 0; j<_tct->_nCH; ++j)
  if(_tct->_chStatus[j])
  his[j]=((TH1F *)_tct->_ch[j]->At(i));
 
  for(Int_t j=0;j<_tct->_nCH; ++j)
  {
  if(_tct->_chStatus[j])
  {
  right[j] = his[j]->GetXaxis()->FindBin(_tA);
  left[j] = 1;
  his[j]->Integral(left[j],right[j]);
  corr[j]=his[j]->Integral(left[j],right[j])/(right[j]-left[j]);
  }
  }
 
  // if(i%100==0)
  // 	printf("[ STATUS] Total Events Processed: %d\n", i);
 
  for(Int_t j=0; j < _tct->_nCH; ++j)
  if(_tct->_chStatus[j])
  for(Int_t k=1; k < his[j]->GetNbinsX(); ++k)
  his[j]->SetBinContent(k,his[j]->GetBinContent(k)-corr[j]);
  }
  delete[] right;
  delete[] left;
  delete[] corr;
  delete[] his;
  printf("[MESSAGE] Baseline Correction Finished!!! \n");
  }
 
  void AnalyzeTCTData::CalcNoise()
  {
  Int_t N = 0;
  Float_t mean = 0.;
  Float_t sum = 0;
  Float_t sumD2 = 0.;
  Float_t sumD4 = 0.;
 
  Int_t *right = new Int_t [_tct->_nCH];
  Int_t *left = new Int_t [_tct->_nCH];
  TH1F *signal;
  TH1F **his = new TH1F*[_tct->_nCH];
 
  for(Int_t i=0; i<_tct->_totEvents-1; ++i)
  {
  if(i==0)
  printf("[ STATUS] Calculating the noise in the Baseline.... \n");
 
  for(Int_t j = 0; j<_tct->_nCH; ++j)
  if(_tct->_chStatus[j])
  his[j]=((TH1F *)_tct->_ch[j]->At(i));
 
  for(Int_t j=0;j<_tct->_nCH; ++j)
  {
  if(_tct->_chStatus[j])
  {
  signal = (TH1F*)his[j]->Clone();
  signal->Scale(_polarity);
  right[j] = signal->GetXaxis()->FindBin(_tA);
  left[j] = 1;
  for(Int_t k = left[j]; k< right[j]; ++k)
  {
  sum += his[j]->GetBinContent(k);
  N++;
  }
  }
  }
  }
  mean = sum / N ;
  for(Int_t i=0; i<_tct->_totEvents; ++i)
  {
  for(Int_t j = 0; j<_tct->_nCH; ++j)
  if(_tct->_chStatus[j])
  his[j]=((TH1F *)_tct->_ch[j]->At(i));
 
  for(Int_t j=0;j<_tct->_nCH; ++j)
  {
  if(_tct->_chStatus[j])
  {
  right[j] = his[j]->GetXaxis()->FindBin(_tA);
  left[j] = 1;
  for(Int_t k = left[j]; k< right[j]; ++k)
  {
  sumD2 += pow(his[j]->GetBinContent(i) - mean, 2);
  sumD4 += pow(his[j]->GetBinContent(i) - mean, 4);
  }
  }
  }
  }
 
  Float_t mu2 = sumD2 / (N - 1); // second central moment (variance)
 
  Float_t stdDev = sqrt(mu2); //std Dev
 
  _noise = stdDev;
 
  printf("[MESSAGE] Noise Estimation Finished (%f mv)!!! \n", _noise);
 
  delete[] right;
  delete[] left;
  delete[] his;
  return;
  }
*/

/*
  void AnalyzeTCTData::EstimateWFShift()
  {
  Float_t maxAmp, minAmp;
  Int_t maxAmpIndex;
  Float_t *maxAmplitude = new Float_t[_events];
  Float_t *minAmplitude = new Float_t[_events];
 
  for(Int_t j=0;j<_events; ++j)
  {
  _histo[0][j] = ((TH1F*) _tct->_ch[_aCH[0]]->At(j));
  maxAmplitude[j] = _histo[0][j]->GetMaximum();
  minAmplitude[j] = _histo[0][j]->GetMinimum();
  }
 
  maxAmp = *max_element(maxAmplitude, maxAmplitude + _events);
  minAmp = *min_element(minAmplitude, minAmplitude + _events);
 
  if(TMath::Abs(minAmp) > TMath::Abs(maxAmp))
  {
  _polarity = -1;
  maxAmpIndex = find(minAmplitude, minAmplitude + _events, minAmp) - minAmplitude;
  }
  else
  {
  _polarity = 1;
  maxAmpIndex = find(maxAmplitude, maxAmplitude + _events, maxAmp) - maxAmplitude;
  }
 
  _histo[0][maxAmpIndex]->Scale(_polarity);
 
  _tA = _tct->t0 + _histo[0][maxAmpIndex]->GetMaximumBin() *_tct->dt;
 
  delete[] maxAmplitude;
  delete[] minAmplitude;
  return;
  }
*/
