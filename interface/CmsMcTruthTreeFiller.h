// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HtoWWTreeDumper
// Description:
//      Class CmsMcTruthTreeFiller
//      Simple class for dumping RECO (or AOD) contents to an ntuple
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  21 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsMcTruthTreeFiller_h
#define CmsMcTruthTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include "CLHEP/HepMC/GenEvent.h"
#include <string>
#include <TTree.h>


class CmsMcTruthTreeFillerData;

class CmsMcTruthTreeFiller {

 public:

  // Constructor
  CmsMcTruthTreeFiller(CmsTree *);

  // Destructor
  virtual ~CmsMcTruthTreeFiller();

  // Write the content of the collection
  void writeCollectionToTree( edm::InputTag mcTruthCollection, const edm::Event& iEvent, int range=100 );

  // Modifiers

 private:

  CmsMcTruthTreeFillerData *privateData_;

};

#endif // CmsMcTruthTreeFiller_h
