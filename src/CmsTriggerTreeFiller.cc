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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"

#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"
#include "CLHEP/HepMC/GenParticle.h"

#include "FWCore/Framework/interface/TriggerNames.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco
;
struct CmsTriggerTreeFillerData {
  CmsTree *cmstree;
};


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsTriggerTreeFiller::CmsTriggerTreeFiller(CmsTree *cmsTree)
{
  cmstree=cmsTree;
}



CmsTriggerTreeFiller::~CmsTriggerTreeFiller() {
}



void
CmsTriggerTreeFiller::writeTriggerToTree (edm::Handle<edm::TriggerResults> & trh,
					  const std::string & columnPrefix, const std::string & columnSuffix) 
{

  vector<bool> Trfired;
  edm::TriggerNames hltNam;
  hltNam.init(*trh);
  std::vector<std::string> hltNames;
  hltNames = hltNam.triggerNames();
  int TrSize = hltNames.size();
  Trfired.resize(TrSize);
  
  for ( int tr=0; tr<TrSize; ++tr) {
    Trfired[tr] = false;
    if ( trh->accept(tr) ) Trfired[tr] = true ;
  }
  
  std::string nTrgString = columnPrefix+"n"+columnSuffix;
  cmstree->column(nTrgString.c_str(), TrSize, 0, "Reco" );
  cmstree->column((columnPrefix+"fired"+columnSuffix).c_str(), Trfired, nTrgString.c_str(), 0, "Reco");

  return ;
}



