#include "AnalyzeTCTData.h"

#include <iostream>
#include <math.h>
#include <vector>
#include "TCanvas.h"
#include "TString.h"
#include "TLine.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc < 2)
    {
      cout<<"[WARNING] : \tUsage: ./analyseTCT fileName.tct \n";
      return -1;
    }
    
  //Read and Analyze the file
  for(Int_t i = 1; i < argc; ++i)
    {
      AnalyzeTCTData *pin = new AnalyzeTCTData(argv[i], 0.9);
      pin->CorrectBaseline();
      
      Int_t nV = pin->_nV1;
      Int_t step = TMath::Abs(pin->_tct->V1[1]-pin->_tct->V1[0]);

      TH1* tempHist = (TH1*) pin->_histo[0][pin->_indexMaxSignal]->Clone();
      tempHist->Scale(-1);
      Float_t lowLim = tempHist->GetMinimum()-5;
      Float_t highLim =  tempHist->GetMaximum()+5;
      Float_t thres[2] = {0., 0.};

      TCanvas *can = new TCanvas("can", pin->_fileName, 1200, 1000);

      TPaveText *bias = new TPaveText(0.528798,0.397541,0.8384808,0.4579918,"nbNDC");
      bias->SetFillColor(kWhite);
      for(Int_t j=0;j<nV; ++j)
	{
	
	  TH1* tempSig = (TH1*) pin->_histo[0][j]->Clone();
	  tempSig->Scale(-1);
	  thres[0] = 0.2 * tempSig->GetMaximum();
	  thres[1] = 0.8 * tempSig->GetMaximum();
	  tempSig->Draw("HIST L");
	  tempSig->SetLineColor(kRed);
	  tempSig->GetXaxis()->SetTitle("Time (ns)");
	  tempSig->GetYaxis()->SetTitle("Voltage (mV)");
	  tempSig->GetYaxis()->SetRangeUser(lowLim, highLim);
	  tempSig->GetXaxis()->SetRangeUser(0, 10);
	  bias->Clear();
	  bias->AddText(Form("V_{bias} = %d V", j*step))->SetTextColor(kBlue);
	  bias->Draw();
	  TLine *lineThres[2];
	  for(Int_t i=0; i<2;++i)
	    {
	      lineThres[i] = new TLine(0, thres[i], 10, thres[i]);
	      lineThres[i]->SetLineColor(kBlack);
	      lineThres[i]->SetLineWidth(4);
	      lineThres[i]->SetLineStyle(9);
	      lineThres[i]->Draw("lsame");
	    }
	  TLine *lineCross[2];

	  vector<Float_t> xData;
	  vector<Float_t> yData;
	  
	  for(Int_t i=0; i<pin->_tct->nPoints; ++i)
	    {
	      // Push the signal part only into the vector
	      // 30ns of the signal will be used to find the riseTime
	      if(tempSig->GetXaxis()->GetBinCenter(i+1) < 0. || tempSig->GetXaxis()->GetBinCenter(i+1) > 30.)
		continue;
	      xData.push_back(tempSig->GetXaxis()->GetBinCenter(i+1));
	      yData.push_back(tempSig->GetBinContent(i+1));
	    }
	  
	  Float_t timeCross[2];
	  timeCross[0] = pin->LinearInterpolation(xData, yData, thres[0]);
	  timeCross[1] = pin->LinearInterpolation(xData, yData, thres[1]);
	  
	  for(Int_t i=0; i<2;++i)
	    {
	      lineCross[i] = new TLine(timeCross[i], lowLim, timeCross[i], highLim);
	      lineCross[i]->SetLineColor(kBlue);
	      lineCross[i]->SetLineWidth(4);
	      lineCross[i]->SetLineStyle(9);
	      lineCross[i]->Draw("lsame");
	    }

	  usleep(5000);
	  gPad->Modified();
	  gPad->Update();
	  can->Print(TString::Format("%s.gif+", pin->_fileName.Data()), ".gif");
	}
    }
  return 0;
}
