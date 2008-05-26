//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsGenInfoFiller
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
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGenInfoFiller.h"

#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco
;
struct CmsGenInfoFillerData {
  CmsTree *cmstree;
  double processID;
  double ptHat;
  double genFilterEff;
  double genCrossXsec;
  double weight;
};


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsGenInfoFiller::CmsGenInfoFiller(CmsTree *cmstree):
  privateData_(new CmsGenInfoFillerData)
{
  privateData_->cmstree=cmstree;

}


// ---------------------------------------------------------------


CmsGenInfoFiller::~CmsGenInfoFiller() {
}


// ---------------------------------------------------------------
void
CmsGenInfoFiller::writeGenInfoToTree (double processID, double ptHat, double genFilterEff,  double genXsec,  double weight, double AlpgenID)
{
  
  std::cout << processID << " " << ptHat << " " << genFilterEff << " " << genXsec << " " << weight << std::endl;
  privateData_->cmstree->column ("genProcessId", processID, 0., "Gen");
  privateData_->cmstree->column ("genPtHat", ptHat, 0., "Gen");
  privateData_->cmstree->column ("genFilterEff", genFilterEff, 0., "Gen");
  privateData_->cmstree->column ("genXsec", genXsec, 0., "Gen");
  privateData_->cmstree->column ("genWeight", weight, 0., "Gen");
  privateData_->cmstree->column ("genAlpgenID", AlpgenID, 0., "Gen");
  return;
}

// ---------------------------------------------------------------
void 
CmsGenInfoFiller::writeGenInfoToTree (edm::Handle<edm::GenInfoProduct> & gi, edm::Handle<edm::HepMCProduct>& mc, double weight )
{
  const HepMC::GenEvent *genEvt = mc->GetEvent();
    
  double processID = genEvt->signal_process_id();
  double pthat = genEvt->event_scale(); 

  double external_cross_section = gi->external_cross_section(); // is the precalculated one written in the cfg file -- units is pb
  double filter_eff = gi->filter_efficiency();

  writeGenInfoToTree(processID,pthat,filter_eff, external_cross_section, weight, 0);
}
