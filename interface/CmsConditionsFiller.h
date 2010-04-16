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
};

class CmsConditionsFiller {

 public:

  /// Constructor
  CmsConditionsFiller(TTree * );

  /// Destructor
  virtual ~CmsConditionsFiller();

  /// dump the "conditions", i.e. trigger mask right now.
  virtual void writeConditionsToTree(const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);

  /// set the HLT results handle
  virtual void setHLTResults(edm::InputTag triggerTag, const edm::Event& iEvent);

 protected:
  
  void writeConditionsInfo(const edm::EventSetup&);
  void treeConditionsInfo(const std::string colPrefix, const std::string colSuffix);

  /// the data container of this filler
  CmsConditionsFillerData *privateData_;
  
  /// the physical tree
  TTree *tree;

  /// the trigger names container
  edm::TriggerNames hltNames_;
  /// to be stored in the ROOT tree
  int nHLT_;
  std::vector< std::string > trgNames_;
  std::vector<unsigned int> trgIndices_;

  /// used to initialize HLT only once if needed
  bool HLTinitialised_;

};

#endif /// CmsConditionsFiller_h
