//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsTriggerTreeFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco;

struct CmsTriggerTreeFillerData {
  CmsTree *cmstree;
};


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsTriggerTreeFiller::CmsTriggerTreeFiller(CmsTree *cmsTree,edm::ParameterSet& iConfig)
{
  cmstree=cmsTree;
  nameProcess_        = iConfig.getParameter<std::string>("processName");
  tagTriggerResults_  = iConfig.getParameter<edm::InputTag>("triggerResults");
  tagTriggerEvent_    = iConfig.getParameter<edm::InputTag>("triggerSummaryAOD");
}



CmsTriggerTreeFiller::~CmsTriggerTreeFiller() {
}



void
CmsTriggerTreeFiller::writeTriggerToTree (const edm::Event& iEvent,
					  const std::string & columnPrefix, const std::string & columnSuffix) 
{

  Handle< trigger::TriggerEvent > handleTriggerEvent;
  if((tagTriggerResults_.process()).compare("AUTO") == 0) {
    iEvent.getByLabel( edm::InputTag(tagTriggerEvent_.label(), tagTriggerEvent_.instance()), handleTriggerEvent );
    if (!handleTriggerEvent.failedToGet()) {
      const edm::Provenance *meta = handleTriggerEvent.provenance();
      if (meta->processName() != nameProcess_) {
        nameProcess_ = meta->processName();
        tagTriggerResults_ = InputTag( tagTriggerResults_.label(), tagTriggerResults_.instance(), nameProcess_ );
        tagTriggerEvent_   = InputTag( tagTriggerEvent_.label(),   tagTriggerEvent_.instance(),   nameProcess_ );
      }
    }
  } else {
    iEvent.getByLabel( tagTriggerEvent_, handleTriggerEvent );
  }

  edm::Handle<edm::TriggerResults> trh;
  try {iEvent.getByLabel(tagTriggerResults_, trh);} 
  catch( cms::Exception& ex ) { edm::LogWarning("CmsTriggerTreeFiller") << "Trigger results: " << tagTriggerResults_ << " not found"; }
  if (!trh.isValid())
      throw cms::Exception("ProductNotValid") << "TriggerResults product not valid";


  vector<int> Trfired;
  edm::TriggerNames hltNam;
  hltNam = iEvent.triggerNames(*trh);
  std::vector<std::string> hltNames;
  hltNames = hltNam.triggerNames();
  int TrSize = hltNames.size();

  // use words in 32 bits: use 30/32 bits per word
  int nWords = (TrSize-1)/30+1;
  Trfired.resize(nWords);
  for(int block=0; block<nWords; block++) Trfired[block]=0;
  
  for ( int tr=0; tr<TrSize; ++tr) {
    int block = tr/30;
    int pos = tr%30;
    int passed = ( trh->accept(tr) ) ? 1 : 0;
    Trfired[block] = Trfired[block] | (passed << pos);
  }

  std::string nTrgString = columnPrefix+"n"+columnSuffix;
  cmstree->column(nTrgString.c_str(), nWords, 0, "Reco" );
  cmstree->column((columnPrefix+"fired"+columnSuffix).c_str(), Trfired, nTrgString.c_str(), 0, "Reco");

  return ;
}



