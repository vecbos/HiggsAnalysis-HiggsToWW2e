// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsPdfWeightFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Maurizio Pierini
//         Created:  Mon Nov  14 11:15:34 CEST 2008
//
//-----------------------------------------------------------------------

#ifndef CmsPdfWeightFiller_h
#define CmsPdfWeightFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

class CmsPdfWeightFiller {

public:
  
  CmsPdfWeightFiller(CmsTree *tree);
  virtual ~CmsPdfWeightFiller();

  void writePdfWeightToTree(edm::InputTag pdfSet, const edm::Event&, const edm::EventSetup&, 
			    const std::string &columnPrefix, const std::string &columnSuffix,
			    bool dumpData=false);

protected:

  CmsTree *cmstree;
  
  std::string *trkIndexName_;

};

#endif // CmsPdfWeightFiller_h
