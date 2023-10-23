g++ plotAllWaveforms.C AnalyzeTCTData.cc ReadTCTFile.cc -o plotWaveforms `root-config --libs --cflags`
./plotWaveforms $1
rm plotWaveforms
