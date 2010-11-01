#ifndef CmsHLTObjectFiller_h
#define CmsHLTObjectFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TTree.h>
#include <string>

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class CmsHLTObjectFiller {
      
 public:

  /// Constructors
  CmsHLTObjectFiller(CmsTree * tree ,edm::ParameterSet& iConfig);

  /// Destructor
  virtual ~CmsHLTObjectFiller();

  /// Write the trigger bits to the tree
  void writeHLTObjectToTree (  const edm::Event & iEvent) ;

  /// Initialize config service
  void beginRun( const edm::Run & iRun, const edm::EventSetup & iSetup ) ;

 protected:
  HLTConfigProvider         hltConfig_;
  std::string               nameProcess_; 
  edm::InputTag             tagTriggerResults_;   
  edm::InputTag             tagTriggerEvent_;       

  CmsTree* tree_;

};

#endif // CmsHLTObjectFiller_h
