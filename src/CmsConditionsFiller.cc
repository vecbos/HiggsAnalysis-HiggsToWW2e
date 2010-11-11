#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"

/// Constructor
CmsConditionsFiller::CmsConditionsFiller(CmsTree *cmsTree, edm::ParameterSet& iConfig, std::vector<std::string>* trgNames) :
  privateData_(new CmsConditionsFillerData)
{
  cmstree = cmsTree;
  nameProcess_        = iConfig.getParameter<std::string>("processName");
  tagTriggerResults_  = iConfig.getParameter<edm::InputTag>("triggerResults");
  tagTriggerEvent_    = iConfig.getParameter<edm::InputTag>("triggerSummaryAOD");
  privateData_->initialise(trgNames);
}

/// Destructor
CmsConditionsFiller::~CmsConditionsFiller() {
  delete privateData_->nHLT;
  delete privateData_->trgIndices;
}

void CmsConditionsFiller::writeConditionsToTree(const edm::Event& iEvent, bool firstEvent) {

  privateData_->clearVectors();

  edm::Handle< trigger::TriggerEvent > handleTriggerEvent;
  if((tagTriggerResults_.process()).compare("AUTO") == 0) {
    iEvent.getByLabel( edm::InputTag(tagTriggerEvent_.label(), tagTriggerEvent_.instance()), handleTriggerEvent );
    if (!handleTriggerEvent.failedToGet()) {
      const edm::Provenance *meta = handleTriggerEvent.provenance();
      if (meta->processName() != nameProcess_) {
        nameProcess_ = meta->processName();
        if(firstEvent) std::cout << "CmsConditionsFiller: using Automatic process name for HLT results. Discovered the following process: " << nameProcess_ << std::endl;
        tagTriggerResults_ = edm::InputTag( tagTriggerResults_.label(), tagTriggerResults_.instance(), nameProcess_ );
        tagTriggerEvent_   = edm::InputTag( tagTriggerEvent_.label(),   tagTriggerEvent_.instance(),   nameProcess_ );
      }
    }
  } else {
    iEvent.getByLabel( tagTriggerEvent_, handleTriggerEvent );
  }

  edm::Handle<edm::TriggerResults> trh;
  try {iEvent.getByLabel(tagTriggerResults_,trh);} 
  catch( cms::Exception& ex ) { std::cout << "Trigger results: " << tagTriggerResults_ << " not found"; }

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
