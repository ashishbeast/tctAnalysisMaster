#include "ReadTCTFile.h"

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout<<"[WARNING] : \tUsage reaTCTFile fileName.tct\n";
        return -1;
    }
    
    for(Int_t i = 1; i < argc; ++i)
    {
        ReadTCTFile *tct = new ReadTCTFile(argv[1], 0.);
        //Print Headers
        tct->Print();
    }
    return 0;
}
