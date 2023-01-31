#include "AnalyzeTCTData.h"

#include <iostream>
#include <math.h>
#include <vector>

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
        AnalyzeTCTData tct(argv[i]);
        tct.CorrectBaseline();
        tct.CalcNoise();
        tct.CalculateWaveformProperties();
        tct.SaveSignalShape();
        tct.AnalysisAction();
    }
    return 0;
}
