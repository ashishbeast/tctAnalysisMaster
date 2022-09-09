#include "AnalyzeTCTData.h"

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      cout<<"[WARNING] : \tUsage: ./analyseTCT fileName.tct \n";
      return -1;
    }

  //Read the file
  AnalyzeTCTData tct(argv[1]);
  return 0;
}
