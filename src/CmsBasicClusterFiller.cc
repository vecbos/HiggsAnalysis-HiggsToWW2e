//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsBasicClusterFiller
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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
// FC: added SC 
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

//#include "DataFormats/EcalRecHit/interface/EcalSeverityLevel.h"
//#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
//#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsBasicClusterFiller.h"

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


CmsBasicClusterFiller::CmsBasicClusterFiller(CmsTree *cmsTree, int maxBC):  privateData_(new CmsBasicClusterFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");
  maxBC_=maxBC;
  privateData_->initialiseCandidate();
  //closestProb_ = DetId(0);
  //severityClosestProb_ = -1;

}

//--------------
// Destructor --
//--------------

CmsBasicClusterFiller::~CmsBasicClusterFiller() 
{
  // delete here the vector ptr's  delete privateData_->nBC;
  delete privateData_->eta;
  delete privateData_->theta;
  delete privateData_->phi;
  delete privateData_->etaCrystal;
  delete privateData_->phiCrystal;
  delete privateData_->iEta;
  delete privateData_->iPhi;
  delete privateData_->thetaTilt;
  delete privateData_->phiTilt;
  delete privateData_->indexSC;
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out


void CmsBasicClusterFiller::writeCollectionToTree(edm::InputTag collectionTag,
						  const edm::Event& iEvent, const edm::EventSetup& iSetup,
						  const std::string &columnPrefix, const std::string &columnSuffix,
						  bool dumpData) 
{
  
  Handle<BasicClusterCollection> collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsBasicClusterFiller") << "Can't get BC Collection: " << collectionTag; }
  const BasicClusterCollection *collection = collectionHandle.product();

  privateData_->clear();
  
  if(collection) 
    {
      if((int)collection->size() > maxBC_)
	{
	  edm::LogError("CmsBasicClusterFiller") << "Track length " << collection->size() 
						 << " is too long for declared max length for tree "
						 << maxBC_ 
						 << ". Collection will be truncated ";
	}
      
      *(privateData_->nBC) = collection->size();

      
      // for super cluster link 
      Handle<SuperClusterCollection> EcalSC;
      try { iEvent.getByLabel(EcalSC_, EcalSC); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsBasicClusterFiller") << "Can't get SC Collection: " << EcalSC_; }
      const SuperClusterCollection *ESCCollection = EcalSC.product();
      
      
      // for cluster shape variables
      Handle< EcalRecHitCollection > EcalBarrelRecHits;
      try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsBasicClusterFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
      const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();
      
      Handle< EcalRecHitCollection > EcalEndcapRecHits;
      try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsBasicClusterFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
      const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();

      try { iEvent.getByLabel(Calotowers_, calotowers_); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get primary calotowers collection" << Calotowers_; }
      
      nBCCandOfSC = 0;
      BasicClusterCollection::const_iterator cand;
      for(cand=collection->begin(); cand!=collection->end(); cand++) 
	{
	  // fill basic kinematics
	  writeBCInfo(&(*cand),iEvent,iSetup,EBRecHits,EERecHits,ESCCollection);
	}
    }
  else 
    {
      *(privateData_->nBC) = 0;
    }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  
  //  int blockSize = (collection) ? collection->size() : 0;
  
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),nBCCandOfSC,0,"Reco");
  
//   if(collection) 
//     {
      treeBCInfo(columnPrefix,columnSuffix);
//     }

  if(dumpData) cmstree->dumpData();

}






void CmsBasicClusterFiller::writeBCInfo(const BasicCluster *cand, 
                                        const edm::Event& iEvent, 
                                        const edm::EventSetup& iSetup,
                                        const EcalRecHitCollection *EBRecHits,
                                        const EcalRecHitCollection *EERecHits,
					const SuperClusterCollection *ESCCollection )

{
  
  // int limit_to_print=1;
  
  int index_SC=-1; 
  int index_SC_ref=-1; 
  
  SuperClusterCollection::const_iterator iSC;
  for(iSC=ESCCollection->begin(); iSC!=ESCCollection->end(); iSC++) {
    index_SC++;
    // find SC constituents 
    for(CaloClusterPtrVector::const_iterator bcItr = iSC->clustersBegin(); bcItr != iSC->clustersEnd(); bcItr++) {
      if( cand->energy()== (*bcItr)->energy() 
          && cand->eta()== (*bcItr)->eta() 
          && cand->phi()== (*bcItr)->phi()  ) {
        index_SC_ref=index_SC;
      }
    }
  }
  
  nBCCandOfSC++;
  privateData_->indexSC->push_back(index_SC_ref);
  
  privateData_->eta->push_back((float)cand->position().eta());
  privateData_->theta->push_back((float)cand->position().theta());
  privateData_->phi->push_back((float)cand->position().phi());
  
  EcalClusterLocal _ecalLocal;
  float etacry,phicry,thetatilt,phitilt;
  int ieta,iphi;
  if(cand->hitsAndFractions().at(0).first.subdetId() == EcalBarrel){
    _ecalLocal.localCoordsEB(*cand,iSetup,etacry,phicry,ieta,iphi,thetatilt,phitilt);
  }else{
    _ecalLocal.localCoordsEE(*cand,iSetup,etacry,phicry,ieta,iphi,thetatilt,phitilt);
  }
  privateData_->etaCrystal->push_back(etacry);
  privateData_->phiCrystal->push_back(phicry);
  privateData_->iEta->push_back(ieta);
  privateData_->iPhi->push_back(iphi);
  privateData_->thetaTilt->push_back(thetatilt);
  privateData_->phiTilt->push_back(phitilt);

}


void CmsBasicClusterFiller::treeBCInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"etaCrystal"+colSuffix).c_str(), *privateData_->etaCrystal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiCrystal"+colSuffix).c_str(), *privateData_->phiCrystal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"iEta"+colSuffix).c_str(), *privateData_->iEta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"iPhi"+colSuffix).c_str(), *privateData_->iPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thetaTilt"+colSuffix).c_str(), *privateData_->thetaTilt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiTilt"+colSuffix).c_str(), *privateData_->phiTilt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"indexSC"+colSuffix).c_str(), *privateData_->indexSC, nCandString.c_str(), 0, "Reco");
}

void CmsBasicClusterFillerData::initialiseCandidate() 
{
  eta = new vector<float>; 
  theta = new vector<float>; 
  phi = new vector<float>;
  etaCrystal = new vector<float>;
  phiCrystal = new vector<float>;
  iEta = new vector<int>;
  iPhi = new vector<int>;
  thetaTilt = new vector<float>;
  phiTilt= new vector<float>;
  indexSC = new vector<int>;
  nBC =  new int;
}

void CmsBasicClusterFillerData::clear() 
{
  eta->clear(); 
  theta->clear();
  phi->clear();
  etaCrystal->clear();
  phiCrystal->clear();
  iEta->clear();
  iPhi->clear();
  thetaTilt->clear();
  phiTilt->clear();
  indexSC->clear();
}
