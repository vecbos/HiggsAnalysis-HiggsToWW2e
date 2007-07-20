//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsMcTruthTreeFiller
//      Simple class for dumping MC truth info to an ntuple. 
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  21 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "CLHEP/HepMC/GenParticle.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"


#include <string>

struct CmsMcTruthTreeFillerData {
  CmsTree *cmstree;
  
};

//----------------
// Constructors --
//----------------
CmsMcTruthTreeFiller::CmsMcTruthTreeFiller(CmsTree *cmstree):
  privateData_(new CmsMcTruthTreeFillerData)
{
  privateData_->cmstree=cmstree; 
}

//--------------
// Destructor --
//--------------
CmsMcTruthTreeFiller::~CmsMcTruthTreeFiller() {
  delete privateData_;
}

//-------------
// Methods   --
//-------------
void CmsMcTruthTreeFiller::writeCollectionToTree(const reco::CandidateCollection *collection, int range) {

  vector<float> pMC,massMC,thetaMC,etaMC,phiMC,energyMC;
  vector<float> xMC,yMC,zMC;
  vector<int> idMC,mothMC,nDauMC;

  reco::CandidateCollection::const_iterator cand;
  for(cand=collection->begin(); cand!=collection->end(); cand++) {
    if((int)pMC.size()>range) break;
    pMC.push_back(cand->p());
    massMC.push_back(cand->mass());
    thetaMC.push_back(cand->theta());
    etaMC.push_back(cand->eta());
    phiMC.push_back(cand->phi());
    energyMC.push_back(cand->energy());
    idMC.push_back(cand->pdgId());
    nDauMC.push_back(cand->numberOfDaughters());
    
    int indMom=-1;
    int idx=0;
    reco::CandidateCollection::const_iterator candIter;
    for(candIter=collection->begin(); candIter!=collection->end(); candIter++) {
      if(&(*candIter)==cand->mother()) {
	indMom=idx;
	break;
      }
      idx++;
    }
    mothMC.push_back(indMom);
    
    // decay vertex
    xMC.push_back(cand->vx());
    yMC.push_back(cand->vy());
    zMC.push_back(cand->vz());
      
  }
  
  std::string theName = "Mc";
  std::string indName = "n"+theName;
  
  privateData_->cmstree->column(indName.c_str(),pMC.size(),0,theName.c_str());
  privateData_->cmstree->column(("p"+theName).c_str(), pMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("mass"+theName).c_str(), massMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("theta"+theName).c_str(), thetaMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("eta"+theName).c_str(), etaMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("phi"+theName).c_str(), phiMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("id"+theName).c_str(), idMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("moth"+theName).c_str(), mothMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("nDau"+theName).c_str(), nDauMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("x"+theName).c_str(), xMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("y"+theName).c_str(), yMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("z"+theName).c_str(), zMC, indName.c_str(), 0, theName.c_str());
    
}
  
