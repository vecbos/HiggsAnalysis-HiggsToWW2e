#include "HiggsAnalysis/HiggsToWW2e/plugins/LumiAnalyzer.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/Algorithms.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

LumiAnalyzer::LumiAnalyzer(const edm::ParameterSet&){
  edm::LogWarning("")<<"LumiAnalyzer::LumiAnalyzer()"<<std::endl;
  nEvent=0;
  nLumiSec=0;
  nRun=0;
  totalLumi=0;

  //lumiOutput = new TProofOutputFile("lumiOutput.root","RECREATE"); //should be the same output file of the tree... how?
  //store_ = new fwlite::TFileService((filePath_+"/PFAnalysis_"+sampleName_+".root").c_str());
  edm::Service<TFileService> store_;
  //edm::Service<TFileService> store_;

  // tree = new TTree("lumiTree","tree with lumi info");
  tree = store_->mkdir("lumiInfo").make<TTree>("lumiTree","tree with lumi info");
  tree->Branch("run",&run,"run/I");
  tree->Branch("lumiSection",&lumiSection,"lumiSection/I");

  tree->Branch("totalLumiByLS",&totalLumiByLS,"totalLumiByLS/D");
  tree->Branch("nEventByLS",&nEventByLS,"nEventByLS/I");

  tree->Branch("totalLumiByRun",&totalLumiByRun,"totalLumiByRun/D");
  tree->Branch("nLumiSecByRun",&nLumiSecByRun,"nLumiSecByRun/I");
  tree->Branch("nEventByRun",&nEventByRun,"nEventByRun/I");

  tree->Branch("totalLumi",&totalLumi,"totalLumi/D");
  tree->Branch("nLumiSec",&nLumiSec,"nLumiSec/I");
  tree->Branch("nEvent",&nEvent,"nEvent/I");
  tree->Branch("nRun",&nRun,"nRun/I");
}


LumiAnalyzer::~LumiAnalyzer(){
  edm::LogWarning("") <<"total number of events "<<nEvent<<std::endl;
  edm::LogWarning("") <<"total number of lumiSec "<<nLumiSec<<std::endl;
  edm::LogWarning("") <<"total number of runs "<<nRun<<std::endl;
  edm::LogWarning("") <<"total lumi "<<totalLumi<<std::endl;

  for(unsigned int i=0; i<runList.size(); i++){
    edm::LogWarning("") <<"run "<<runList[i]<<" with "<<lumiSecByRun[i]<<" lumisec;  "<<eventByRun[i]<<" events;  "<<lumiByRun[i]<<" luminosity"<<std::endl;
  }
  //lumiOutput->cd();
  // tree->Write();
  edm::LogWarning("") <<"LumiAnalyzer::~LumiAnalyzer()"<<std::endl;

}


void LumiAnalyzer::beginRun(edm::Run const& iRun , edm::EventSetup const& eventSetup){
  totalLumiByRun=0;
  nLumiSecByRun=0;
  nEventByRun=0;

  run = iRun.run();
  nRun++;
}

void LumiAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& eventSetup) {

  edm::LogWarning("")<<"LumiAnalyzer::beginLuminosityBlock"<<std::endl;
  edm::Handle<LumiSummary> lumiSummary;
  lumiBlock.getByLabel("lumiProducer", lumiSummary);
  if(lumiSummary->isValid()){
    edm::LogWarning("")<<"LumiAnalyzer::lumiSummary is valid"<<std::endl;

    lumiSection = lumiSummary->lsNumber();
    nLumiSec++;
    nLumiSecByRun++;
  }  

  nEventByLS=0;
  totalLumiByLS=0;
}

void LumiAnalyzer::analyze(edm::Event const&, edm::EventSetup const&){

  nEventByLS++;
  nEventByRun++;
  nEvent++;

  //edm::LogWarning("") <<"nEventByLS "<<nEventByLS<<std::endl;
}

void LumiAnalyzer::endLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& eventSetup) {

  edm::LogWarning("")<<"LumiAnalyzer::endLuminosityBlock"<<std::endl;
  edm::Handle<LumiSummary> lumiSummary;
  lumiBlock.getByLabel("lumiProducer", lumiSummary);
  //This is not working properly in 384, it will work in 386
  /*if(lumiSummary->isValid()){
    float liveFraction = (*lumiSummary).liveFrac();
    float deliveredLumi = (*lumiSummary).avgInsDelLumiErr();
    totalLumi = totalLumi+ liveFraction*deliveredLumi;
    totalLumiByRun = totalLumiByRun + liveFraction*deliveredLumi;
    }*/
  //This is the recipe from zhen for 384
  if(lumiSummary->isValid()){
    edm::LogWarning("")<<"LumiAnalyzer::lumiSummary is valid"<<std::endl;
    float integratedRecorded;
    if(lumiSummary->nTriggerLine()) {
      float integratedDelivered=lumiSummary->avgInsDelLumi()*lumiSummary->lumiSectionLength();
      float deadFrac=float (lumiSummary->deadcount())/float(lumiSummary->l1info(0).ratecount*lumiSummary->l1info(0).prescale);
      integratedRecorded=integratedDelivered*(1.0-deadFrac); //this is in unit of ub-1
    }
    else
      integratedRecorded=0;

    totalLumi = totalLumi+ integratedRecorded;
    totalLumiByRun = totalLumiByRun + integratedRecorded;
    totalLumiByLS = integratedRecorded;
  }  
  /*edm::LogWarning("") <<"totalLumi "<<totalLumi<<std::endl;
   edm::LogWarning("") <<"totalLumiByRun "<<totalLumiByRun<<std::endl;
   edm::LogWarning("") <<"totalLumiByLS "<<totalLumiByLS<<std::endl;

   edm::LogWarning("") <<"nEventByLS "<<nEventByLS<<std::endl;*/

  tree->Fill();
}

void LumiAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const&){
  runList.push_back(run);

  lumiByRun.push_back(totalLumiByRun);
  lumiSecByRun.push_back(nLumiSecByRun);
  eventByRun.push_back(nEventByRun);
}



DEFINE_FWK_MODULE(LumiAnalyzer);

