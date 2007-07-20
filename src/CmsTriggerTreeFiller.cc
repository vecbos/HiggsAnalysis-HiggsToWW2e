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


CmsTriggerTreeFiller::CmsTriggerTreeFiller(CmsTree *cmstree):
  privateData_(new CmsTriggerTreeFillerData)
{
  privateData_->cmstree=cmstree;

  m_TrigNames.push_back ("SingleElectron"       ) ;
  m_TrigNames.push_back ("RelaxedSingleElectron") ;
  m_TrigNames.push_back ("DoubleElectron"       ) ;
  m_TrigNames.push_back ("RelaxedDoubleElectron") ;
  m_TrigNames.push_back ("SinglePhoton"         ) ;
  m_TrigNames.push_back ("RelaxedSinglePhoton"  ) ;
  m_TrigNames.push_back ("DoublePhoton"         ) ;
  m_TrigNames.push_back ("RelaxedDoublePhoton"  ) ;
  m_TrigNames.push_back ("SingleEMHighEt"       ) ;
  m_TrigNames.push_back ("SingleEMVeryHighEt"   ) ;

}


// ---------------------------------------------------------------


CmsTriggerTreeFiller::~CmsTriggerTreeFiller() {
}


// ---------------------------------------------------------------


void
CmsTriggerTreeFiller::writeTriggerToTree (edm::Handle<edm::TriggerResults> & trh,
					  const std::string & columnPrefix, const std::string & columnSuffix) 
{
  int TrSize = trh->size();
  bool Trfired [10] ;
  // loop over trigger paths
  for (int tr=0 ; tr < 10 ; ++tr) {
    int ind = trh->find (m_TrigNames[tr]) ; // sistema di default per dire che nn ha trovato nulla:
    Trfired[tr] = false;                  // ritorna la size
    if ((ind < TrSize) && (trh->accept (ind)))
      Trfired[tr] = true ;
  } // loop over triggers
 
  privateData_->cmstree->column ((columnPrefix+"singleElePassed"+columnSuffix).c_str(),      Trfired[0], "Reco");
  privateData_->cmstree->column ((columnPrefix+"singleEleRelaxPassed"+columnSuffix).c_str(), Trfired[1], "Reco");
  privateData_->cmstree->column ((columnPrefix+"doubleElePassed"+columnSuffix).c_str(),      Trfired[2], "Reco");
  privateData_->cmstree->column ((columnPrefix+"doubleEleRelaxPassed"+columnSuffix).c_str(), Trfired[3], "Reco");
  privateData_->cmstree->column ((columnPrefix+"singlePhoPassed"+columnSuffix).c_str(),      Trfired[4], "Reco");
  privateData_->cmstree->column ((columnPrefix+"singlePhoRelaxPassed"+columnSuffix).c_str(), Trfired[5], "Reco");
  privateData_->cmstree->column ((columnPrefix+"doublePhoPassed"+columnSuffix).c_str(),      Trfired[6], "Reco");
  privateData_->cmstree->column ((columnPrefix+"doublePhoRelaxPassed"+columnSuffix).c_str(), Trfired[7], "Reco");
  privateData_->cmstree->column ((columnPrefix+"highEMPassed"+columnSuffix).c_str(),         Trfired[8], "Reco");
  privateData_->cmstree->column ((columnPrefix+"veryHighEMPassed"+columnSuffix).c_str(),     Trfired[9], "Reco");
  
  return ;
}
