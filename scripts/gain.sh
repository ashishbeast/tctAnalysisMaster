g++ plotGain.C AnalyzeTCTData.cc ReadTCTFile.cc -o gainSlapp `root-config --libs --cflags --ldflags`
./gainSlapp
rm gainSlapp
