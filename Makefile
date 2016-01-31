all: TreeReader.run

TreeReader.run: TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o TopTools/ttbar_Categorization/ttbar_category.o TopTools/SF_btag/BTagCalibrationStandalone.o TopTools/SF_btag/BTagSFUtil.o TopTools/SF_Lumi/SFLumi.o TreeReader.o 
	g++ -o TreeReader.run TreeReader.o TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o TopTools/ttbar_Categorization/ttbar_category.o TopTools/SF_btag/BTagSFUtil.o TopTools/SF_btag/BTagCalibrationStandalone.o TopTools/SF_Lumi/SFLumi.o `root-config --libs`

TreeReader.o: TreeReader.C
	g++ -std=c++11 -static -I`root-config --incdir` -ITopTools/SF_ID-ISO-Trigger -ITopTools/ttbar_Categorization -ITopTools/SF_btag -ITopTools/SF_Lumi -c -g TreeReader.C

TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o: TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.C
	g++ -std=c++11 -static -I `root-config --incdir` -c TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.C -o TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o

TopTools/ttbar_Categorization/ttbar_category.o: TopTools/ttbar_Categorization/ttbar_category.C
	g++ -std=c++11 -static -I `root-config --incdir` -c TopTools/ttbar_Categorization/ttbar_category.C -o TopTools/ttbar_Categorization/ttbar_category.o

TopTools/SF_btag/BTagSFUtil.o: TopTools/SF_btag/BTagSFUtil.C
	g++ -std=c++11 -static -I `root-config --incdir` -c TopTools/SF_btag/BTagSFUtil.C -o TopTools/SF_btag/BTagSFUtil.o

TopTools/SF_btag/BTagCalibrationStandalone.o: TopTools/SF_btag/BTagCalibrationStandalone.C
	g++ -std=c++11 -static -I `root-config --incdir` -c TopTools/SF_btag/BTagCalibrationStandalone.C -o TopTools/SF_btag/BTagCalibrationStandalone.o

TopTools/SF_Lumi/SFLumi.o: TopTools/SF_Lumi/SFLumi.C
	g++ -std=c++11 -static -I `root-config --incdir` -c TopTools/SF_Lumi/SFLumi.C -o TopTools/SF_Lumi/SFLumi.o

clean:
	rm TreeReader.run TreeReader.o TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o TopTools/ttbar_Categorization/ttbar_category.o TopTools/SF_btag/BTagCalibrationStandalone.o TopTools/SF_btag/BTagSFUtil.o TopTools/SF_Lumi/SFLumi.o
