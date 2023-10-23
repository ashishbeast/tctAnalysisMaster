#ifndef READTCTFILE_H
#define READTCTFILE_H

/*
TCT Data Format:
Line 0: Type of Measurement (11= waveform measurements, 12 = waveform maximum only, 13 = waveform maximum and integrals)
Line 1: Date
Line 2: Absolute time in seconds at start
Line 3: x0, dx, nX
Line 4: y0, dy, nY
Line 5: z0, dz, nZ
Line 6: channel status (ch1, ch2, ch3, ch4)
Line 7: nV1, V1[nV1]
Line 8: nV2, V2[nV2]
Line 9: t0, dt, nPoints
*/

#include <stdio.h>
#include <stdlib.h>

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TArray.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TText.h"

using namespace std;

class ReadTCTFile
{
private:
  Bool_t _swap = false;
  //Bool_t _swap = true;
  void byteSwap(Char_t *, Char_t*);
  void byteSwapIP(Float_t *, Int_t);
  
public:
  //Constructor
  //ReadTCTFile(const Char_t* inFileName);
  ReadTCTFile(const Char_t* inFileName, Float_t tA);
  //Destructor
  ~ReadTCTFile();

  //Functions
  void ReadWfs(Float_t tA);
  void ReadWfsBin(Float_t tA);
  TH1F* GetHisto(Int_t ch, Int_t index);
  TH1F* GetHisto(Int_t ch, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2);
  Int_t GetIndex(Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2);
  void Print();

  //Variables
  static const Int_t _nCH = 4;
  FILE *inFile;
  Int_t date[6];
  Int_t fileType, absTime;
  Int_t nX, nY, nZ, nV1, nV2;
  Float_t stepX, stepY, stepZ, stepV1, stepV2;
  Float_t x0, y0, z0;
  TArrayF V1, V2, I1, I2; //Array of Voltages
  Float_t t0, dt;
  Int_t nPoints;
  Int_t _nV;
  Int_t _nSteps;
  Int_t _chStatus[_nCH];
  Int_t _totEvents;
  Float_t _tA = 0.; //start of the signal time
  Int_t _polarity = -1;
  
  //Header 
  Float_t header[200];
  Float_t Temp;                    // temperature
  Float_t Source;               // type of e-h generation
  Char_t *User;                 // user taking the measurements
  Char_t *Sample;               // Sample name
  Char_t *Comment;              // Comment

  TClonesArray* _ch[_nCH];
  
  Float_t *xCords[22];
};

#endif
