#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"


/// Constructor
CmsConditionsFiller::CmsConditionsFiller(CmsTree *cmsTree, int maxLenght, bool noOutputIfLimitsReached) :
  privateData_(new CmsConditionsFillerData)
{

  cmstree = cmsTree;
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  HLTinitialised_ = false;

}

/// Destructor
CmsConditionsFiller::~CmsConditionsFiller() {
}



void CmsConditionsFiller::writeConditionsToTree(const edm::EventSetup& iSetup,
						const std::string &columnPrefix, const std::string &columnSuffix,
						bool dumpData) {

  writeConditionsInfo(iSetup);
  treeConditionsInfo(columnPrefix, columnPrefix);

  if ( dumpData ) cmstree->dumpData();

}



void CmsConditionsFiller::setHLTResults(edm::Handle<edm::TriggerResults> & trh) {
  
  hltNames_.init(*trh);
  HLTinitialised_ = true;

}



void CmsConditionsFiller::writeConditionsInfo(const edm::EventSetup& iSetup) {
}



void CmsConditionsFiller::treeConditionsInfo(const std::string colPrefix, const std::string colSuffix) {
  
  
  if ( HLTinitialised_ ) {

    std::vector< std::string > hltStringNames; 
    hltStringNames = hltNames_.triggerNames();

    /// put 1 column for each path with the name and the corresponding index
    std::vector< std::string >::const_iterator hltPath;
    for( hltPath = hltStringNames.begin(); hltPath != hltStringNames.end(); ++hltPath) {

      unsigned int trgIndex = hltNames_.triggerIndex( *hltPath );
      cmstree->column((colPrefix+(*hltPath)+colSuffix).c_str(), trgIndex, 0, "Reco");

    }

  }

  else {

    cmstree->column((colPrefix+"UNKNOWN"+colSuffix).c_str(), -1, 0, "Reco");
    
  }
  
}
