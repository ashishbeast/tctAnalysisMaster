g++ jitterMoveit.C AnalyzeTCTData.cc ReadTCTFile.cc -o jitterMoveit `root-config --libs --cflags`
./jitterMoveit
rm jitterMoveit

