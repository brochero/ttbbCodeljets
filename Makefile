all: TreeReader.run

TreeReader.run: TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o TreeReader.o 
	g++ -o TreeReader.run TreeReader.o TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o `root-config --libs`

TreeReader.o: TreeReader.C
	g++ -I`root-config --incdir` -ITopTools/SF_ID-ISO-Trigger -c -g TreeReader.C

TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o: TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.C
	g++ -I `root-config --incdir` -c TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.C -o TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o

clean:
	rm TreeReader.run TreeReader.o TopTools/SF_ID-ISO-Trigger/SFIDISOTrigger.o
