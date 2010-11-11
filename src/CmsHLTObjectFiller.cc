#include "HiggsAnalysis/HiggsToWW2e/interface/CmsHLTObjectFiller.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"


CmsHLTObjectFiller::CmsHLTObjectFiller(CmsTree * tree,edm::ParameterSet& iConfig){
  
  tree_ = tree;
  nameProcess_        = iConfig.getParameter<std::string>("processName"); 
  tagTriggerResults_  = iConfig.getParameter<edm::InputTag>("triggerResults");   
  tagTriggerEvent_    = iConfig.getParameter<edm::InputTag>("triggerSummaryAOD") ;       
  
}

CmsHLTObjectFiller::~CmsHLTObjectFiller(){}


void CmsHLTObjectFiller::beginRun(const edm::Event & iEvent,  const edm::Run & iRun, const edm::EventSetup & iSetup ){
  // Initialize HLTConfigProvider
  bool changed( true );

  edm::Handle< trigger::TriggerEvent > handleTriggerEvent;
  if((tagTriggerResults_.process()).compare("AUTO") == 0) {
    iEvent.getByLabel( edm::InputTag(tagTriggerEvent_.label(), tagTriggerEvent_.instance()), handleTriggerEvent );
    if (!handleTriggerEvent.failedToGet()) {
      const edm::Provenance *meta = handleTriggerEvent.provenance();
      if (meta->processName() != nameProcess_) {
        nameProcess_ = meta->processName();
        tagTriggerResults_ = edm::InputTag( tagTriggerResults_.label(), tagTriggerResults_.instance(), nameProcess_ );
        tagTriggerEvent_   = edm::InputTag( tagTriggerEvent_.label(),   tagTriggerEvent_.instance(),   nameProcess_ );
      }
    }
  } else {
    iEvent.getByLabel( tagTriggerEvent_, handleTriggerEvent );
  }

  if ( ! hltConfig_.init( iRun, iSetup, nameProcess_, changed ) ) {
    edm::LogError( "errorHltConfigExtraction" ) << "HLT config extraction error with process name " << nameProcess_;
  } else if ( hltConfig_.size() <= 0 ) {
    edm::LogError( "hltConfigSize" ) << "HLT config size error";
  }
  
}



void CmsHLTObjectFiller::writeHLTObjectToTree(const edm::Event & iEvent){
  //4-momenta of trigger objects
  std::vector<float> hltPt,hltEta,hltPhi,hltMass;
  //indices of passing objects
  std::vector<int> passingIndex;
  //association of paths to passing objects
  std::vector<int> passingIndexPerPath;
  std::vector<int> passingIndexRanges;

  //read trigger results
  edm::Handle< edm::TriggerResults > handleTriggerResults;
  iEvent.getByLabel( tagTriggerResults_, handleTriggerResults );
  edm::Handle< trigger::TriggerEvent > handleTriggerEvent;
  iEvent.getByLabel( tagTriggerEvent_, handleTriggerEvent );
 
  //Load info used in other trigger dump modules to check consitency
  edm::TriggerNames hltNam;
  hltNam = iEvent.triggerNames(*handleTriggerResults);
  std::vector<std::string> hltNames;
  hltNames = hltNam.triggerNames();
  unsigned TrSize = hltNames.size();
  
  //prepare tables
  const unsigned sizePaths( hltConfig_.size() );
  const unsigned sizeFilters( handleTriggerEvent->sizeFilters() );
  const unsigned sizeObjects( handleTriggerEvent->sizeObjects() );

  if(sizePaths != TrSize) 
    edm::LogError( "HLTObjectFiller" ) << "HLT config size mismatch";

  hltPt.resize(sizeObjects);
  hltEta.resize(sizeObjects);
  hltPhi.resize(sizeObjects);
  hltMass.resize(sizeObjects);
  passingIndexPerPath.resize(sizePaths);
  passingIndexRanges.resize(sizePaths);

  //loop and copy all trigger objects
  for ( size_t iO = 0; iO < sizeObjects ; ++iO ) {
    trigger::TriggerObject triggerObject( handleTriggerEvent->getObjects().at( iO ) );
    hltPt[iO]=triggerObject.pt();
    hltEta[iO]=triggerObject.eta();
    hltPhi[iO]=triggerObject.phi();
    hltMass[iO]=triggerObject.mass();
  }


  int totalpassed=0;
  //std::cout << "X" <<std::endl;

  //loop over paths (iP = trigger index)
  for ( size_t iP = 0; iP < sizePaths; ++iP ) {

    if(iP != hltNam.triggerIndex(hltConfig_.triggerName( iP ) ))
       edm::LogError( "HLTObjectFiller" ) << "Trigger index";


    //only save info for successful paths
    if(!handleTriggerResults->wasrun(iP) || !handleTriggerResults->accept(iP)){
      passingIndexPerPath[iP]=-1;
      passingIndexRanges[iP]=-1;
      //std::cout << "Path: " << hltConfig_.triggerName( iP ) << " failed or was not run" <<std::endl;
      continue;
    }

    //find the last filter in the path
    const std::string namePath( hltConfig_.triggerName( iP ) );
    const unsigned indexLastFilterPath( handleTriggerResults->index( iP ) );
    unsigned indexLastFilterPathModules( indexLastFilterPath + 1 );
    unsigned indexLastFilterFilters( sizeFilters );

    //std::cout << "finding last filter for: " << namePath << " " << iP <<" starting at module " << indexLastFilterPathModules << std::endl;

    while ( indexLastFilterPathModules > 0 ) {
      --indexLastFilterPathModules;
      const std::string labelLastFilterModules( hltConfig_.moduleLabel( iP, indexLastFilterPathModules ) );      
      //std::cout << "looking for module: " << labelLastFilterModules << " (" << indexLastFilterPathModules <<")" <<std::endl;
      indexLastFilterFilters = handleTriggerEvent->filterIndex( edm::InputTag( labelLastFilterModules, "", nameProcess_ ) );
      //std::cout <<"number " << indexLastFilterFilters <<" of" << sizeFilters <<std::endl;
      if ( indexLastFilterFilters < sizeFilters ) break;
    }
    if(indexLastFilterFilters == sizeFilters){// no last filter? L1 Passthrough?
      // no last filter => no passing objects
      passingIndexPerPath[iP]=-1;
      passingIndexRanges[iP]=0;
      //std::cout << "Path: " << hltConfig_.triggerName( iP ) << "has no last filter" <<std::endl;
      continue;
    }

    //retrieve indeces of passing objects
    const trigger::Keys& KEYS(handleTriggerEvent->filterKeys(indexLastFilterFilters));
    passingIndexRanges[iP]=KEYS.size();
    if(KEYS.size()==0){// no passing objects => nothing to store; 
      passingIndexPerPath[iP]=-1;
      //std::cout << "Path: " << hltConfig_.triggerName( iP ) << " has no passing objects" <<std::endl;
      continue;
    }
    passingIndexPerPath[iP]=totalpassed;
    totalpassed+=KEYS.size();
    //std::cout << "Path: " << hltConfig_.triggerName( iP ) << "has passing objects" <<std::endl;
    for(trigger::Keys::const_iterator key=KEYS.begin();key!=KEYS.end();key++){
      passingIndex.push_back(*key);
      //std::cout << "  eta:" << hltEta[*key] <<"  phi: " << hltPhi[*key] << "  pt: " << hltPt[*key] <<" index:" << *key <<  std::endl;
    }
  }

  //save index collections
  //std::cout <<"nTriggerPaths" << std::endl;
  tree_->column("nTriggerPaths", sizePaths, 0, "Reco" );
  //std::cout << "nTriggerObsPassing"<< std::endl;
  tree_->column("nTriggerObsPassing",totalpassed, 0, "Reco" );
  //std::cout << "sizePassing"<< std::endl;
  tree_->column("sizePassing",passingIndexRanges,"nTriggerPaths" , 0, "Reco");
  //std::cout <<"indexPassing" << std::endl;
  tree_->column("indexPassing",passingIndex,"nTriggerObsPassing" , 0, "Reco");
  //std::cout << "indexPassingPerPath"<< std::endl;
  tree_->column("indexPassingPerPath",passingIndexPerPath,"nTriggerPaths" , 0, "Reco");

  //save trigger objects
  //std::cout <<"nTriggerObs" << std::endl;
  tree_->column("nTriggerObs", sizeObjects, 0, "Reco" );
  //std::cout <<"triggerObsPt" << std::endl;
  tree_->column("triggerObsPt",hltPt,"nTriggerObs" , 0, "Reco");
  //std::cout <<"triggerObsPhi" << std::endl;
  tree_->column("triggerObsPhi",hltPhi,"nTriggerObs" , 0, "Reco");
  //std::cout <<"triggerObsEta" << std::endl;
  tree_->column("triggerObsEta",hltEta,"nTriggerObs" , 0, "Reco");
  //std::cout <<"triggerObsMass" << std::endl;
  tree_->column("triggerObsMass",hltMass,"nTriggerObs" , 0, "Reco");

}
