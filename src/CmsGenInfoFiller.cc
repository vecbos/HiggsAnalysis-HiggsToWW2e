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
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGenInfoFiller.h"

#include "FWCore/Framework/interface/TriggerNames.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco;



//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsGenInfoFiller::CmsGenInfoFiller(CmsTree *cmsTree):
  privateData_(new CmsGenInfoFillerData)
{
  cmstree=cmsTree;
  privateData_->initialise();
}


// ---------------------------------------------------------------


CmsGenInfoFiller::~CmsGenInfoFiller() {
  // delete here the privateData_ members
  delete privateData_->ptHat;
  delete privateData_->processID;
  delete privateData_->weight;
  delete privateData_->alphaQCD;
  delete privateData_->alphaQED;
}


void CmsGenInfoFiller::writeGenInfoToTree ( edm::Handle<GenEventInfoProduct> & gei )
{
  *(privateData_->ptHat)     = gei->qScale();
  *(privateData_->processID) = gei->signalProcessID();
  *(privateData_->weight)    = gei->weight();
  *(privateData_->alphaQCD)  = gei->alphaQCD();
  *(privateData_->alphaQED)  = gei->alphaQED();

  double thePtHat     = gei->qScale();
  double theProcessID = gei->signalProcessID();
  double theWeight    = gei->weight();
  double theAlphaQCD  = gei->alphaQCD();
  double theAlphaQED  = gei->alphaQED();

  cmstree->column("genPtHat",     thePtHat,     0., "Gen");
  cmstree->column("genProcessId", theProcessID, 0., "Gen");
  cmstree->column("genWeight",    theWeight,    0., "Gen");
  cmstree->column("genAlphaQCD",  theAlphaQCD,  0., "Gen");
  cmstree->column("genAlphaQED",  theAlphaQED,  0., "Gen");
}

void CmsGenInfoFillerData::initialise() {

  ptHat = new float;
  processID = new float;
  weight = new float;
  alphaQCD = new float;
  alphaQED = new float;

}

