// -*- C++ -*-
//
// Package:    HtoWWElectrons
// Class:      HtoWWEleCloner
// 

// user include files
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h" 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/HtoWWEleCloner.h"

#include <iostream>
#include <TMath.h>


HtoWWEleCloner::HtoWWEleCloner (const edm::ParameterSet& conf) :
  electronProducer_       (conf.getParameter<edm::InputTag>("src")) ,
  hwwEleCandidCollection_ (conf.getParameter<std::string> ("HWWgsfEleCandCollection"))
{
  // to be put in the event
  produces<reco::CandidateCollection> (hwwEleCandidCollection_) ;
}  


// --------------------------------------------------------
  

HtoWWEleCloner::~HtoWWEleCloner(){}


// --------------------------------------------------------
  

void
HtoWWEleCloner::produce (edm::Event& e, const edm::EventSetup& iSetup)
{
  // get reconstructed electrons
  edm::Handle<reco::GsfElectronCollection> electrons;
  try { e.getByLabel (electronProducer_,electrons) ; }
  catch ( cms::Exception& ex ) 
    { 
      edm::LogWarning ("catching") << "Can't get collection " ;
      return ;
    }

  std::auto_ptr<reco::CandidateCollection> recoEleCandColl (
    new reco::CandidateCollection );

  typedef reco::GsfElectronCollection::const_iterator citer ;
  for (citer MyEle = electrons->begin () ; 
       MyEle != electrons->end () ; 
       ++MyEle) 
    {
      // filling the CandidateCollection
      reco::Candidate* MyEleCand = (*MyEle).clone () ;
      recoEleCandColl->push_back (MyEleCand) ;
    }
    
  // put selected information in the event
  e.put ( recoEleCandColl, hwwEleCandidCollection_) ;
}


// --------------------------------------------------------
  

void 
HtoWWEleCloner::beginJob (const edm::EventSetup&)
{}


// --------------------------------------------------------
  

void 
HtoWWEleCloner::endJob() 
{}


