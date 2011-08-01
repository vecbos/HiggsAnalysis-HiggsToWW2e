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
  void writeCollectionToTree( edm::InputTag mcTruthCollection, std::vector<std::string>* lheComments, const edm::Event& iEvent, int range=100, bool firstEvent=false );

  // Modifiers
  void saveLHEComments(bool what) {saveLHE_ = what; }

 private:

  CmsMcTruthTreeFillerData *privateData_;
  bool saveLHE_;

};

#endif // CmsMcTruthTreeFiller_h
