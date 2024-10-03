g++ plotAmplitudeVariation.C AnalyzeTCTData.cc ReadTCTFile.cc -o ampVar `root-config --libs --cflags --ldflags`
./ampVar
rm ampVar

