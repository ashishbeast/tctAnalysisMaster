g++ plotCharge.C AnalyzeTCTData.cc ReadTCTFile.cc -o chargeSlapp `root-config --libs --cflags --ldflags`
./chargeSlapp
rm chargeSlapp

