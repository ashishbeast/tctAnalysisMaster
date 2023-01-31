g++ plotMaps.C AnalyzeTCTData.cc ReadTCTFile.cc -o mapSlapp `root-config --libs --cflags --ldflags`
./mapSlapp
rm mapSlapp

