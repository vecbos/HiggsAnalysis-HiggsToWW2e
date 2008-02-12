#ifndef CmsConditionsFiller_h
#define CmsConditionsFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include "FWCore/Framework/interface/TriggerNames.h"

#include <string>
#include <vector>

struct CmsConditionsFillerData {
};

class CmsConditionsFiller {

 public:

  /// Constructor
  CmsConditionsFiller(CmsTree *, int maxLenght=500, bool noOutputIfLimitsReached=false );

  /// Destructor
  virtual ~CmsConditionsFiller();

  /// dump the "conditions", i.e. trigger mask right now.
  virtual void writeConditionsToTree(const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);

  /// set the HLT results handle
  virtual void setHLTResults(edm::Handle<edm::TriggerResults> & trh);

 protected:
  
  void writeConditionsInfo(const edm::EventSetup&);
  void treeConditionsInfo(const std::string colPrefix, const std::string colSuffix);

  /// the data container of this filler
  CmsConditionsFillerData *privateData_;
  
  /// the physical tree
  CmsTree *cmstree;

  /// do not write array elements above this threshold
  bool hitLimitsMeansNoOutput_;
  int maxLenght_;

  /// the trigger names container
  edm::TriggerNames hltNames_;
  
  /// used to initialize HLT only once if needed
  bool HLTinitialised_;

};

#endif /// CmsConditionsFiller_h
