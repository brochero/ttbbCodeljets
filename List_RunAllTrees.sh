#!/bin/sh
file="Tree_LepJets_v7-6-1_NewSF_btag_Spring15-bx25_2170pb-1"

#######################################
# $2 should be '-o OutputName... etc' #
#######################################

./TreeReader.run -i ${file}_DataSingleEG $1
./TreeReader.run -i ${file}_DataSingleMu $1

./TreeReader.run -i ${file}_ttbar_PowhegPythia   $1
./TreeReader.run -i ${file}_ttbar_MadgraphPythia $1
./TreeReader.run -i ${file}_ttbar_MCatNLOPythia  $1

./TreeReader.run -i ${file}_ttbar_PowhegPythia $1 -cat tt
./TreeReader.run -i ${file}_ttbar_PowhegPythia $1 -cat ttbb
./TreeReader.run -i ${file}_ttbar_PowhegPythia $1 -cat ttb
./TreeReader.run -i ${file}_ttbar_PowhegPythia $1 -cat ttcc
./TreeReader.run -i ${file}_ttbar_PowhegPythia $1 -cat ttLF
./TreeReader.run -i ${file}_ttbar_PowhegPythia $1 -cat ttjj

./TreeReader.run -i ${file}_ttbar_PowhegPythiaBkg $1

./TreeReader.run -i ${file}_tW            $1
./TreeReader.run -i ${file}_tbarW         $1
./TreeReader.run -i ${file}_tbar_tchannel $1
./TreeReader.run -i ${file}_t_tchannel    $1

./TreeReader.run -i ${file}_WW $1
./TreeReader.run -i ${file}_WZ $1
./TreeReader.run -i ${file}_ZZ $1 

./TreeReader.run -i ${file}_WJets_MCatNLO         $1 
./TreeReader.run -i ${file}_ZJets_M50_MCatNLO     $1
./TreeReader.run -i ${file}_ZJets_M10to50_MCatNLO $1

./TreeReader.run -i ${file}_QCD_MuEnr_20to30    $1
./TreeReader.run -i ${file}_QCD_MuEnr_30to50    $1
./TreeReader.run -i ${file}_QCD_MuEnr_50to80    $1
./TreeReader.run -i ${file}_QCD_MuEnr_80to120   $1
./TreeReader.run -i ${file}_QCD_MuEnr_120to170  $1
./TreeReader.run -i ${file}_QCD_MuEnr_170to300  $1
./TreeReader.run -i ${file}_QCD_MuEnr_300to470  $1
./TreeReader.run -i ${file}_QCD_MuEnr_470to600  $1
./TreeReader.run -i ${file}_QCD_MuEnr_800to1000 $1
#./TreeReader.run -i ${file}_QCD_MuEnr_1000toInf $1

./TreeReader.run -i ${file}_QCD_EGEnr_15to20      $1
./TreeReader.run -i ${file}_QCD_EGEnr_20to30      $1
./TreeReader.run -i ${file}_QCD_EGEnr_30to50      $1
./TreeReader.run -i ${file}_QCD_EGEnr_50to80      $1
./TreeReader.run -i ${file}_QCD_EGEnr_80to120     $1
#./TreeReader.run -i ${file}_QCD_EGEnr_120to170    $1
./TreeReader.run -i ${file}_QCD_EGEnr_170to300    $1
./TreeReader.run -i ${file}_QCD_EGEnr_300toInf    $1




