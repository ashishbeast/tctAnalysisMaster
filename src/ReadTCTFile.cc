#include "ReadTCTFile.h"

#include <iostream>

ReadTCTFile::ReadTCTFile(const Char_t* inFileName, Float_t tA)
{
  Bool_t binFile = 1;
  Int_t userSize, sampleSize, commentSize, offset; 

  //Set All channels to OFF
  for(Int_t i =0; i < _nCH; ++i)
    _chStatus[i] = 0;

  //Initialize Measurement Preamble
  User = NULL;
  Comment = NULL;
  Sample = NULL;

  //Open File
  if(binFile)
    inFile = fopen(inFileName, "rb+");
  else
    inFile = fopen(inFileName, "r+");
      
  if(inFile == NULL)
    {
      printf("  [ERROR] : Unable to open File for Reading \n");
      printf("[Warning] : Exiting program!!! \n");
      exit(0);
    }
  
  if(!binFile) //Read ASCII File
    {
      fscanf(inFile, "%d", &fileType);
      if(!(fileType==11 || fileType==22 || fileType==33 || fileType==51 || fileType==81 || fileType==82))  // if it is something else exit
	{
	  printf("  [ERROR] : Can not read other formats than waveform: %d!\n", fileType);
	  printf("[Warning] : Exiting program!!! \n");
	  exit(0);
	}
      else
	{
	  fscanf(inFile,"%d %d %d %d %d %d\n", &date[0], &date[1], &date[2], &date[3], &date[4], &date[5]);
	  fscanf(inFile,"%d\n", &absTime);
	  fscanf(inFile,"%f %f %d\n", &x0, &stepX, &nX);
	  fscanf(inFile,"%f %f %d\n", &y0, &stepY, &nY);
	  fscanf(inFile,"%f %f %d\n", &z0, &stepZ, &nX);
	  _nSteps = nX*nY*nZ;
	  if(fileType==22)
	    fscanf(inFile,"%d %d %d\n", &_chStatus[0], &_chStatus[1], &_chStatus[2]);
	  if(fileType==33 || fileType==51 || fileType==81 || fileType==82)
	    fscanf(inFile,"%d %d %d %d\n", &_chStatus[0], &_chStatus[1], &_chStatus[2], &_chStatus[3]);
	  fscanf(inFile,"%d", &nV1);
	  V1=TArrayF(nV1);
	  for(Int_t i=0; i<nV1; ++i)
	    fscanf(inFile,"%f", &V1[i]);    
	  if(fileType!=11)
	    fscanf(inFile,"%d", &nV2);
	  else
	    nV2 = 1;
	  V2=TArrayF(nV2);
	  if(fileType!=11)
	    for(Int_t i=0; i<nV2; ++i)
	      fscanf(inFile,"%f", &V2[i]);

	  _nV = nV1+nV2;
	  _totEvents = _nSteps*nV1*nV2; //number of all waveforms

	  I1=TArrayF(nV2*nV1); 
	  I2=TArrayF(nV2*nV1); 

	  fscanf(inFile,"%f %f %d\n",&t0,&dt,&nPoints);

	  //Book Arrays for x, y, z, V1, V2, I1, I2, BM, time
	  for(Int_t i =0; i<9; ++i)
	    xCords[i] = new Float_t [_nSteps*_nV];

	  //initialize Histograms
	  for(Int_t i=0; i<_nCH; ++i)
	    {
	      if(_chStatus[i])
		{
		  _ch[i] = new TClonesArray("TH1F", _totEvents);
		  _ch[i]->BypassStreamer(kFALSE);
		}
	    }
	  ReadWfs(tA); //Read the eaveforms
	}
    }
  else //Read Binary File
    {
      Int_t read = fread((void *)header, sizeof(Float_t), 1, inFile);

      if(_swap)
	byteSwapIP(header, read);
      
      fileType = (Int_t) header[0];
      cout<<Form("Reading File type = %d \n",fileType);

      //Read in the Buffer
      rewind(inFile);
      read = fread((void *)header, sizeof(Float_t), 200, inFile);

      if(_swap)
	byteSwapIP(header, read);
      
      for(Int_t i=0;i<6; ++i)
	date[i]=(Int_t) header[i+1];

      absTime = (Int_t) header[7];

      //get Movement Matrix
      x0 = header[8];  stepX = header[9];  nX = (Int_t) header[10];
      y0 = header[11]; stepY = header[12]; nY = (Int_t) header[13];
      z0 = header[14]; stepZ = header[15]; nZ = (Int_t) header[16];

      _nSteps = nX*nY*nZ;
      
      //Adjust the reading of the header
      if(fileType==33 || fileType==51 || fileType==81 || fileType==82)
	offset = 1;
      else
	offset = 0;
      
      for(Int_t i =0; i < _nCH-1+offset; ++i)
	_chStatus[i] = (Int_t) header[17+i];

      //get Number of voltage steps for Power Supply 1
      nV1 = (Int_t) header[20+offset];
      V1 = TArrayF(nV1);
      for(Int_t i =0;i <nV1;++i)
	V1[i] = header[21+offset+i];

      //get Number of voltage steps for Power Supply 2
      nV2 = (Int_t) header[21+offset+nV1];
      V2 = TArrayF(nV2);
      for(Int_t i =0;i <nV2;++i)
	V2[i] = header[22+offset+nV1+i];

      _nV = nV1+nV2;
      _totEvents = _nSteps*nV1*nV2; //number of all waveforms
      
      //Initialize Current Arrays
      I1 = TArrayF(nV1*nV2);
      I2 = TArrayF(nV1*nV2);

      //get Time Scale
      t0 = header[22+offset+_nV];
      if(TMath::Abs(t0)>1e-3)
	t0 *= 1e-9;
      
      dt = header[23+offset+_nV];
      if(TMath::Abs(dt)>1e-3)
	dt *= 1e-9;

      nPoints = (Int_t) header[24+offset+_nV];

      //Header Information
      Temp = header[25+offset+_nV];
      Source = (Int_t) header[26+offset+_nV];

      //rewind to the apropriate position
      fseek(inFile, (28+_nV)*sizeof(Float_t), SEEK_SET);
      
      fread(&userSize, sizeof(Int_t), 1, inFile);
      if(_swap)
	byteSwapIP((Float_t *) &userSize, 1);
      User = new Char_t[userSize+1];
      fread(User, sizeof(Char_t), userSize, inFile);
      User[userSize] = '\0';
      
      fread(&sampleSize, sizeof(Int_t), 1, inFile);
      if(_swap)
	byteSwapIP((Float_t *) &sampleSize, 1);
      Sample = new Char_t[sampleSize+1];
      fread(Sample, sizeof(Char_t), sampleSize, inFile);
      Sample[sampleSize] = '\0';

      fread(&commentSize, sizeof(Int_t), 1, inFile);
      if(_swap)
	byteSwapIP((Float_t *) &commentSize, 1);
      Comment = new Char_t[commentSize+1];
      fread(Comment, sizeof(Char_t), commentSize, inFile);
      Comment[commentSize] = '\0';

      //Book Arrays for x, y, z, V1, V2, I1, I2, BM, time
      for(Int_t i =0; i<22; ++i)
	xCords[i] = new Float_t [_nSteps*_nV];

      //initialize Histograms
      for(Int_t i=0; i<_nCH; ++i)
	{
	  if(_chStatus[i])
	    {
	      _ch[i] = new TClonesArray("TH1F", _totEvents);
	      _ch[i]->BypassStreamer(kFALSE);
	    }
	}	
      ReadWfsBin(tA);
    }
}

//Destructor
ReadTCTFile::~ReadTCTFile()
{
}

//Read the waveform from Binary Format
void ReadTCTFile::ReadWfsBin(Float_t tA)
{
  Int_t index;
  Float_t tV1, tV2, tI1, tI2;

  //Char_t **hisname = new Char_t*[_nCH];
  //for(Int_t i=0; i< _nCH; ++i)
  //  histo[i] = *ch[i];

  Char_t hisname1[100];
  Char_t hisname2[100];
  Char_t hisname3[100];
  Char_t hisname4[100];

  TClonesArray &entryP0 = *_ch[0];
  TClonesArray &entryP1 = *_ch[1];
  TClonesArray &entryP2 = *_ch[2];
  TClonesArray &entryP3 = *_ch[3];
  
  Float_t buf[33000];
  
  for(Int_t i=0; i<nV1; ++i)
    {
      for(Int_t j=0; j<nV2; ++j)
	{
	  fread((void *)buf, sizeof(Float_t), 4, inFile);

	  if(_swap)
	    byteSwapIP(buf, 4);

	  tV1 = buf[0]; tV2 = buf[1];
	  tI1 = buf[2]; tI2 = buf[3];
	  
	  V1[i] = tV1;
	  I1[j+ i*nV2] = tI1;
	  V2[j] = tV2;
	  I2[j+ i*nV2] = tI2;

	  for(Int_t k=0; k<_nSteps; ++k)
	    {
	      index = k + _nSteps*j + nV2*_nSteps*i;

	      if(fileType==33 || fileType==22)
		fread((void *)buf, sizeof(Float_t), 4, inFile);

	      if(fileType==51)
		fread((void *)buf, sizeof(Float_t), 5, inFile);

	      if(fileType==81)
		fread((void *)buf, sizeof(Float_t), 8, inFile);

	      if(fileType==82)
		fread((void *)buf, sizeof(Float_t), 18, inFile);

	      if(_swap)
		byteSwapIP(buf, 18);

	      for(Int_t m=0; m<3; ++m)
		xCords[m][index] = buf[m];

	      xCords[7][index] = buf[3]; //Beam Monitor

	      if(fileType == 82)
		for(Int_t m=0; m<14; ++m)
		  xCords[m+8][index] = buf[m+4];

	      xCords[3][index] = tV1;
	      xCords[4][index] = tV2;
	      xCords[5][index] = tI1;
	      xCords[6][index] = tI2;
	      

	      // //Creating Histograms for each Index
	      // for(Int_t m=0; m<_nCH; ++m)
	      // 	{
	      // 	  if(_chStatus[m])
	      // 	    {
	      // 	      sprintf(hisname[m],"Channel %d: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", m, xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);

	      // 	      new(_ch[m][index]) TH1F((const Char_t *)(hisname[m]), (const Char_t *)(hisname[m]), nPoints, t0*1e-9, (t0+nPoints*dt)*1e-9);
	      // 	    }
	      // 	}

	      // //Filling Histograms
	      // for(Int_t m=0; m<_nCH; ++m)
	      // 	{
	      // 	  if(_chStatus[m])
	      // 	    {
	      // 	      fread(buf, sizeof(Float_t), nPoints, inFile);
	      // 	      for(Int_t p=0; p<nPoints;++p)
	      // 		((TH1F*) _ch[m][index])->SetBinContent(p+1, buf[p]);
	      // 	    }
	      // 	}

	       if(_chStatus[0])
		{
		  sprintf(hisname1,"Channel 1: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP0[index]) TH1F((const Char_t *)(hisname1), (const Char_t *)(hisname1), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      if(_chStatus[1])
		{
		  sprintf(hisname2,"Channel 2: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP1[index]) TH1F((const Char_t *) hisname2, (const Char_t *)(hisname2), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      if(_chStatus[2])
		{
		  sprintf(hisname3,"Channel 3: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP2[index]) TH1F((const Char_t *)(hisname3), (const Char_t *)(hisname3), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      if(_chStatus[3])
		{
		  sprintf(hisname4,"Channel 4: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP3[index]) TH1F((const Char_t *)(hisname4), (const Char_t *)(hisname4), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      
	      for(Int_t m=0; m<_nCH; ++m)
		{
		  if(_chStatus[m])
		    {
		      fread(buf, sizeof(Float_t), nPoints, inFile);
		      for(Int_t p=0; p<nPoints;++p)
			{
			  if(_swap)
			    byteSwapIP(buf, nPoints);

			  switch(m)
			    {
			    case 0:
			      ((TH1F*)entryP0[index])->SetBinContent(p+1, buf[p]);
			      break;
			    case 1:
			      ((TH1F*)entryP1[index])->SetBinContent(p+1, buf[p]);
			      break;
			    case 2:
			      ((TH1F*)entryP2[index])->SetBinContent(p+1, buf[p]);
			      break;
			    case 3:
			      ((TH1F*)entryP3[index])->SetBinContent(p+1, buf[p]);
			      break;
			    }
			}
		    }
		}
	    }
	}
    }
}

//Read the waveforms from ASCII Format
void ReadTCTFile::ReadWfs(Float_t tA)
{
  Int_t index;
  Float_t data, tV1, tV2, tI1, tI2;
  Char_t **hisname = new char*[_nCH];

  TClonesArray &entryP0 = *_ch[0];
  TClonesArray &entryP1 = *_ch[1];
  TClonesArray &entryP2 = *_ch[2];
  TClonesArray &entryP3 = *_ch[3];

  // TClonesArray histo[4];
  // for(Int_t i=0; i< _nCH; ++i)
  //   histo[i] = *ch[i];

  for(Int_t i=0; i<nV1; ++i)
    {
      for(Int_t j=0; j<nV2; ++j)
	{
	  fscanf(inFile, "%f %f %f %f", &tV1, &tV2, &tI1, &tI2);
	  
	  V1[i] = tV1;
	  I1[j+ i*nV2] = tI1;
	  V2[j] = tV2;
	  I2[j+ i*nV2] = tI2;

	  for(Int_t k=0; k<_nSteps; ++k)
	    {
	      index = k + _nSteps*j + nV2*_nSteps*i;

	      if(fileType == 22 || fileType == 33)
		{
		  for(Int_t m=0; m<4; ++m)
		    fscanf(inFile, "%f", &xCords[m][index]);

		  xCords[7][index] = xCords[3][index]; //Beam Monitor
		  xCords[3][index] = tV1;
		  xCords[4][index] = tV2;
		  xCords[5][index] = tI1;
		  xCords[6][index] = tI2;
		}
	      if(fileType == 11)
		for(Int_t m=0; m<5; ++m)
		  fscanf(inFile, "%f", &xCords[m][index]);

	      // //Creating Histograms for each Index
	      // for(Int_t m=0; m<_nCH; ++m)
	      // 	{
	      // 	  if(_chStatus[m])
	      // 	    {
	      // 	      sprintf(hisname[m],"Channel %d: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", m, xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
	      // 	      new(_ch[m][index]) TH1F((const Char_t *)(hisname[m]), (const Char_t *)(hisname[m]), nPoints, t0*1e-9-tA*1e-9, (t0+nPoints*dt)*1e-9-tA*1e-9);
	      // 	    }
	      // 	}

	      // //Filling Histograms
	      // for(Int_t m=0; m<_nCH; ++m)
	      // 	{
	      // 	  if(_chStatus[m])
	      // 	    {
	      // 	      for(Int_t p=0; p<nPoints;++p)
	      // 		{
	      // 		  fscanf(inFile, "%e", &data);
	      // 		  ((TH1F*) _ch[m][index])->SetBinContent(p+1, data);
	      // 		}
	      // 	    }
	      // 	}
	      
	      if(_chStatus[0])
		{
		  sprintf(hisname[0],"Channel 1: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP0[index]) TH1F((const Char_t *)(hisname[0]), (const Char_t *)(hisname[0]), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      if(_chStatus[1])
		{
		  sprintf(hisname[1],"Channel 2: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP1[index]) TH1F((const Char_t *)(hisname[1]), (const Char_t *)(hisname[1]), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      if(_chStatus[2])
		{
		  sprintf(hisname[2],"Channel 3: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP2[index]) TH1F((const Char_t *)(hisname[2]), (const Char_t *)(hisname[2]), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      if(_chStatus[3])
		{
		  sprintf(hisname[3],"Channel 4: x=%5.3e,y=%5.3e,z=%5.3e,V1=%4.2f, V2=%4.2f ", xCords[0][index], xCords[1][index], xCords[2][index], xCords[3][index], xCords[4][index]);
		  new(entryP3[index]) TH1F((const Char_t *)(hisname[3]), (const Char_t *)(hisname[3]), nPoints, t0*1e9-tA, (t0+nPoints*dt)*1e9-tA);
		}
	      
	      for(Int_t m=0; m<_nCH; ++m)
		{
		  if(_chStatus[m])
		    {
		      for(Int_t p=0; p<nPoints;++p)
			{
			  fscanf(inFile, "%e", &data);
			  switch(m)
			    {
			    case 0:
			      ((TH1F*)entryP0[index])->SetBinContent(p+1, data);
			      break;
			    case 1:
			      ((TH1F*)entryP1[index])->SetBinContent(p+1, data);
			      break;
			    case 2:
			      ((TH1F*)entryP2[index])->SetBinContent(p+1, data);
			      break;
			    case 3:
			      ((TH1F*)entryP3[index])->SetBinContent(p+1, data);
			      break;
			    }
			}
		    }
		}
	    }
	}
    }
}


TH1F* ReadTCTFile::GetHisto(Int_t ch, Int_t index)
{
  TH1F *his;
  his = (TH1F*)_ch[ch]->At(index);
  return his;
}

TH1F* ReadTCTFile::GetHisto(Int_t ch, Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2)
{
  return (GetHisto(ch, GetIndex(x, y, z, v1, v2)));
}

// Function returns the index of the position of the waveform corresponding to x,y,z,v1,v2 in the linear array of waveforms.
Int_t ReadTCTFile::GetIndex(Int_t x, Int_t y, Int_t z, Int_t v1, Int_t v2)
{
  if(x > nX-1 || x<0)
    {
      printf("[  ERROR] Index x out of Range!! \n");
      return 0;
    }
  if(y > nY-1 || y<0)
    {
      printf("[  ERROR] Index y out of Range!! \n");
      return 0;
    }
  if(z > nZ-1 || z<0)
    {
      printf("[  ERROR] Index z out of Range!! \n");
      return 0;
    }
  if(v1 > nV1-1 || v1<0)
    {
      printf("[  ERROR] Index v1 out of Range!! \n");
      return 0;
    }
  if(v2 > nV2-1 || v2<0)
    {
      printf("[  ERROR] Index v2 out of Range!! \n");
      return 0;
    }
  return ((x+nX*y+(nX*nY)*z) + _nSteps*v2 +(nV2*_nSteps)*v1);
}


void ReadTCTFile::Print()
{
  //Prints the information about the class and its members
  printf("*************************************\n");
  printf("[            Date]: %d.%d.%d \n", date[0], date[1], date[2]);
  printf("[            Time]: %d:%d:%d \n", date[3], date[4], date[5]);
  if(User!=NULL)
    printf("[            User]: %s \n",User);
  if(Sample!=NULL)
    printf("[          Sample]: %s \n",Sample);
  printf("[ Active Channels]: Ch1=%d Ch2=%d Ch3=%d Ch4=%d \n",_chStatus[0], _chStatus[1], _chStatus[2], _chStatus[3]);
  printf("[ Position Points]: %d (X=%d, Y=%d, Z=%d) \n", _nSteps, nX, nY, nZ);
  printf("[       Positions]: r0=(%f,%f,%f) dr=(%f,%f,%f) \n", x0, y0, z0, stepX, stepY, stepZ);
  printf("[      Time scale]: points=%d, t0=%e, dt=%e \n", nPoints, t0, dt);
  printf("[     Temperature]: %f \n", Temp);
  printf("[ Generation Type]: %4.0f \n", Source);
  printf("[ Voltage Sources]: NU1=%d , NU2=%d \n", nV1, nV2);
  
  if(Comment!=NULL)
    printf("[         Comment]: %s \n", Comment);
  printf("*************************************\n");
  if(nV1>1 and nV2>1)
    {
      for(Int_t i=0; i<nV1; ++i) 
	for(Int_t j=0; j<nV2; ++j) 
	  printf("V1, V2(%f,%f)::I1,I2(%e,%e)\n",V1[i],V2[j],I1[j+nV2*i],I2[j+nV2*i]);
    }
  printf("\n");
}

// byte swaping (LABVIEW,HPUX g++)<->(LINUX g++) 
void ReadTCTFile::byteSwap(Char_t *a, Char_t *b)
{
  Char_t c = *a;
  *a = *b;
  *b = c;
}

void ReadTCTFile::byteSwapIP(Float_t *in, Int_t size)
{
  Char_t *sr;
  while(size--)
    {
      sr = (Char_t *)in;
      byteSwap(&sr[0], &sr[3]);
      byteSwap(&sr[1], &sr[2]);
      in++;
    }
}

