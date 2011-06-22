#!/bin/bash -f
# usage createGoodList.sh CERT.json

certjson=$1

# do the V1 OR V2 for all the tasks
compareJSON.py --or DoubleMu_v1/res/lumiSummary.json DoubleMu_v2/res/lumiSummary.json DoubleMu.json
compareJSON.py --or DoubleElectron_v1/res/lumiSummary.json DoubleElectron_v2/res/lumiSummary.json DoubleElectron.json
compareJSON.py --or MuEG_v1/res/lumiSummary.json MuEG_v2/res/lumiSummary.json MuEG.json
compareJSON.py --or SingleMu_v1/res/lumiSummary.json SingleMu_v2/res/lumiSummary.json SingleMu.json
compareJSON.py --or SingleElectron_v1/res/lumiSummary.json SingleElectron_v2/res/lumiSummary.json SingleElectron.json

# and now the AND of the important double lepton PDs
compareJSON.py --and DoubleMu.json DoubleElectron.json output1.json
compareJSON.py --and MuEG.json output1.json DoubleLeptonPDs.json

# and now do the OR with the certification json and calc the lumi
lumiCalc.py -c frontier://LumiCalc/CMS_LUMI_PROD -i DoubleLeptonPDs.json --nowarning overview >& DoubleLeptonPD.lumi &
compareJSON.py --and $certjson DoubleLeptonPDs.json DoubleLeptonPDsGood.json
lumiCalc.py -c frontier://LumiCalc/CMS_LUMI_PROD -i DoubleLeptonPDsGood.json --nowarning overview >& DoubleLeptonPDGood.lumi &

echo "END. Final jsons are:"
echo "---> from crab reportDoubleLeptonPDs.json (lumi = DoubleLeptonPD.lumi)"
echo "---> from crab report in AND with DQM json DoubleLeptonPDsGood.json (lumi = DoubleLeptonPDGood.lumi)"
