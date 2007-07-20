// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HtoWWTreeDumper
// Description:
//      Class CmsTriggerTreeFiller
//      Simple class for dumping RECO (or AOD) contents to an ntuple
//      
// Original Author:  Alessio Ghezzi, Pietro Govoni
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsTriggerTreeFiller_h
#define CmsTriggerTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

class CmsTriggerTreeFillerData;

class CmsTriggerTreeFiller {

 public:

  //! Constructors
  CmsTriggerTreeFiller(CmsTree *);

  //! Destructor
  virtual ~CmsTriggerTreeFiller();

  void writeTriggerToTree (edm::Handle<edm::TriggerResults> & trh,
                           const std::string & columnPrefix, const std::string & columnSuffix) ;

 private:
  
  CmsTriggerTreeFillerData *privateData_;
  std::vector<std::string> m_TrigNames ;

};

#endif // CmsTriggerTreeFiller_h
