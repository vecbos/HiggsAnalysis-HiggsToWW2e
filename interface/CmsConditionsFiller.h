#ifndef CmsConditionsFiller_h
#define CmsConditionsFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include <TTree.h>

#include <string>
#include <vector>

struct CmsConditionsFillerData {
  int *nHLT;
  std::vector<int> *trgIndices;
  std::vector<std::string> *trgNames;
  
  public:
  void initialise(std::vector<std::string>* trgnames);
  void clearVectors();
};

class CmsConditionsFiller {

 public:

  /// Constructor
  CmsConditionsFiller(CmsTree *, edm::ParameterSet&, std::vector<std::string>* trgNames );

  /// Destructor
  virtual ~CmsConditionsFiller();

  /// dump the "conditions", i.e. trigger mask right now.
  virtual void writeConditionsToTree(const edm::Event&, bool firstEvent);

  /// create the branches at the first event
  void setHLTResults();

 protected:
  
  /// the data container of this filler
  CmsConditionsFillerData *privateData_;
  
  /// the physical tree
  CmsTree *cmstree;

  std::string               nameProcess_; 
  edm::InputTag             tagTriggerResults_;   
  edm::InputTag             tagTriggerEvent_;

};

#endif /// CmsConditionsFiller_h
