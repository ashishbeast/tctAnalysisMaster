g++ plotJitter.C AnalyzeTCTData.cc ReadTCTFile.cc -o jitterSlapp `root-config --libs --cflags`
./jitterSlapp
rm jitterSlapp

