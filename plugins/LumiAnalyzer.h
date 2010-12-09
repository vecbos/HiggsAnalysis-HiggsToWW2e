#ifndef __LumiAnalyzer__
#define __LumiAnalyzer__

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
//#include "TProofOutputFile.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

class LumiAnalyzer : public edm::EDAnalyzer {

  public:


    ///Default contructor
  LumiAnalyzer(const edm::ParameterSet&);
  LumiAnalyzer(fwlite::TFileService* store_=0);

    ///Default destructor
    virtual ~LumiAnalyzer();

    ///Method where the analysis is done.
  virtual void analyze(edm::Event const&, edm::EventSetup const&);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& eventSetup);
  virtual void endLuminosityBlock(const edm::LuminosityBlock& , const edm::EventSetup&);
      
  private:

  int run;
  int lumiSection;

  int nEvent;
  int nLumiSec;
  double totalLumi;
  int nRun;

  int nLumiSecByRun;
  int nEventByRun;
  double totalLumiByRun;

  int nEventByLS;
  double totalLumiByLS;

  std::vector<int> runList;
  std::vector<int> lumiSecByRun;
  std::vector<int> eventByRun;
  std::vector<double> lumiByRun;

  TTree* tree;
  //TFile*  lumiOutput;
  //fwlite::TFileService *store_;
  };


#endif
