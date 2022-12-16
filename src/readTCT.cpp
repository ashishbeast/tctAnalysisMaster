#include "ReadTCTFile.h"

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      cout<<"[WARNING] : \tUsage reaTCTFile fileName.tct\n";
      return -1;
    }

  ReadTCTFile tct(argv[1], 0.);

  //Print Headers
  tct.Print();
  
  return 0;
}
