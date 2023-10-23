g++ gainThickness.C AnalyzeTCTData.cc ReadTCTFile.cc -o gainThick `root-config --libs --cflags`
./gainThick
rm gainThick



