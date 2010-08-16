#include "DataFormats/Common/interface/TriggerResults.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"

/// Constructor
CmsConditionsFiller::CmsConditionsFiller(CmsTree *cmsTree, std::vector<std::string>* trgNames) :
  privateData_(new CmsConditionsFillerData)
{
  cmstree = cmsTree;
  privateData_->initialise(trgNames);
}

/// Destructor
CmsConditionsFiller::~CmsConditionsFiller() {
  delete privateData_->nHLT;
  delete privateData_->trgIndices;
}

void CmsConditionsFiller::writeConditionsToTree(edm::InputTag triggerTag, const edm::Event& iEvent, bool firstEvent) {

  privateData_->clearVectors();

  edm::Handle<edm::TriggerResults> trh;
  try {iEvent.getByLabel(triggerTag,trh);} 
  catch( cms::Exception& ex ) { std::cout << "Trigger results: " << triggerTag << " not found"; }

  /// the trigger names container
  edm::TriggerNames hltNames_ = iEvent.triggerNames(*trh);

  std::vector< std::string > hltStringNames; 
  hltStringNames = hltNames_.triggerNames();
  int TrSize = hltStringNames.size();
  
  *(privateData_->nHLT) = TrSize;
  
  /// put 1 column for each path with the name and the corresponding index
  std::vector< std::string >::const_iterator hltPath;
  for( hltPath = hltStringNames.begin(); hltPath != hltStringNames.end(); ++hltPath) {
    
    unsigned int trgIndex = hltNames_.triggerIndex( *hltPath );
    
    privateData_->trgNames->push_back( *hltPath );
    privateData_->trgIndices->push_back( trgIndex );
    
  }

  cmstree->column("nHLT", TrSize, 0, "HLT");
  cmstree->column("indexHLT", *privateData_->trgIndices, "nHLT", 0, "HLT");
  if(firstEvent) cmstree->getTree()->Branch("nameHLT", &(*privateData_->trgNames));

}



void CmsConditionsFillerData::initialise(std::vector<std::string>* trgnames) {
  nHLT = new int;
  trgNames = trgnames;
  trgIndices = new vector<int>;
}

void CmsConditionsFillerData::clearVectors() {
  trgNames->clear();
  trgIndices->clear();
}
