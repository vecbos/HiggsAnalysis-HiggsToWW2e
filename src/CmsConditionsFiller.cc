#include "DataFormats/Common/interface/TriggerResults.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"

/// Constructor
CmsConditionsFiller::CmsConditionsFiller(TTree *Tree) :
  privateData_(new CmsConditionsFillerData)
{

  tree = Tree;
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

  if ( dumpData ) tree->Fill();

}



void CmsConditionsFiller::setHLTResults(edm::InputTag triggerTag,
					const edm::Event& iEvent) {
  
  edm::Handle<edm::TriggerResults> trh;
  try {iEvent.getByLabel(triggerTag,trh);} 
  catch( cms::Exception& ex ) { std::cout << "Trigger results: " << triggerTag << " not found"; }

  hltNames_.init(*trh);
  HLTinitialised_ = true;

  tree->Branch("nHLT", &nHLT_, "nHLT/I");
  tree->Branch("nameHLT", &trgNames_);
  tree->Branch("indexHLT", &trgIndices_);

}



void CmsConditionsFiller::writeConditionsInfo(const edm::EventSetup& iSetup) {
}



void CmsConditionsFiller::treeConditionsInfo(const std::string colPrefix, const std::string colSuffix) {
  
  
  if ( HLTinitialised_ ) {

    std::vector< std::string > hltStringNames; 
    hltStringNames = hltNames_.triggerNames();
    int TrSize = hltStringNames.size();

    nHLT_ = TrSize;

    /// put 1 column for each path with the name and the corresponding index
    std::vector< std::string >::const_iterator hltPath;
    for( hltPath = hltStringNames.begin(); hltPath != hltStringNames.end(); ++hltPath) {

      unsigned int trgIndex = hltNames_.triggerIndex( *hltPath );
      
      trgNames_.push_back( *hltPath );
      trgIndices_.push_back( trgIndex );

    }

  }

}
