g++ chargeErrorLgads.C AnalyzeTCTData.cc ReadTCTFile.cc -o chargeErrors `root-config --libs --cflags`
./chargeErrors
rm chargeErrors


