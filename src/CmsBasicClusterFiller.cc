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
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

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
}

//--------------
// Destructor --
//--------------

CmsBasicClusterFiller::~CmsBasicClusterFiller() 
{
  // delete here the vector ptr's
  delete privateData_->nCrystals;
  delete privateData_->nOverlap3x3;
  delete privateData_->energy;
  delete privateData_->eta;
  delete privateData_->theta;
  delete privateData_->phi;
  delete privateData_->seedEnergy;
  delete privateData_->eMax;
  delete privateData_->e3x3;
  delete privateData_->e5x5;
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
      
      // for cluster shape variables
      Handle< EcalRecHitCollection > EcalBarrelRecHits;
      try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsBasicClusterFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
      const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();
      
      Handle< EcalRecHitCollection > EcalEndcapRecHits;
      try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsBasicClusterFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
      const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();
      
      BasicClusterCollection::const_iterator cand;
      for(cand=collection->begin(); cand!=collection->end(); cand++) 
	{
	  // fill basic kinematics
	  writeBCInfo(&(*cand),iEvent,iSetup,EBRecHits,EERecHits);
	}
    }
  else 
    {
      *(privateData_->nBC) = 0;
    }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  
  int blockSize = (collection) ? collection->size() : 0;
  
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
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
                                        const EcalRecHitCollection *EERecHits) 
{

  std::vector< std::pair<DetId, float> > ids = cand->hitsAndFractions();
  
  privateData_->nCrystals->push_back((int)ids.size());
  privateData_->eta->push_back((float)cand->position().eta());
  privateData_->theta->push_back((float)cand->position().theta());
  privateData_->phi->push_back((float)cand->position().phi());
  
  const EcalRecHitCollection *rechits = 0;
  float seedEta = cand->position().eta();
  if( fabs(seedEta) < 1.479 ) rechits = EBRecHits;
  else rechits = EERecHits; 

  // find the seed energy
  // find the eventual noisy channels inside the cluster and record their energy
  double noisyChanEnergy = 0.0;
  DetId seedId;
  EcalRecHitCollection::const_iterator seedItr = rechits->begin();

  for(std::vector< std::pair<DetId,float> >::const_iterator idItr = ids.begin(); idItr != ids.end(); ++idItr) {
    DetId id = idItr->first;
    if(id.det() != DetId::Ecal) { continue; }
    EcalRecHitCollection::const_iterator hitItr = rechits->find(id);
    if(hitItr == rechits->end()) { continue; }
    if(hitItr->energy() > seedItr->energy()) {
      seedItr = hitItr;
      seedId = id;
    }
    
    if( hitItr->recoFlag() == EcalRecHit::kDead || EcalRecHit::kFaultyHardware ) noisyChanEnergy += hitItr->energy();

  }

  double energyToBeRemoved = (removeBadChannels_) ? noisyChanEnergy : 0.0;

  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);

  bool validTopologyAndGeometry = false;
  
  if ( pTopology.isValid() && pGeometry.isValid() ) {
    
    validTopologyAndGeometry = true;
    
    const CaloTopology *topology = pTopology.product();
    
    float eMax = EcalClusterTools::eMax( *cand, &(*rechits) ) - energyToBeRemoved;    
    float e3x3 = EcalClusterTools::e3x3( *cand, &(*rechits), topology ) - energyToBeRemoved;
    float e5x5 = EcalClusterTools::e5x5( *cand, &(*rechits), topology ) - energyToBeRemoved;

    privateData_->e3x3->push_back(e3x3);
    privateData_->e5x5->push_back(e5x5);
    privateData_->eMax->push_back(eMax);
    
  } else {
    privateData_->e3x3->push_back(-1.0);
    privateData_->e5x5->push_back(-1.0);
    privateData_->eMax->push_back(-1.0);
  }

  privateData_->energy->push_back((float)cand->energy() - energyToBeRemoved);
  privateData_->seedEnergy->push_back(seedItr->energy() - energyToBeRemoved);

  //  DetId seedId = cand->seed(); // it seems it doesn't work
  
  int nCry3x3=0;

  std::vector< std::pair<DetId, float> >::const_iterator clusItr = ids.begin();
  
  if(seedId.subdetId() == EcalBarrel) {
    EBDetId ebId = (EBDetId)seedId;
    for(int icry=0; icry<9; ++icry) {
      unsigned int row    = icry/3;
      unsigned int column = icry%3;
      int icryEta = ebId.ieta()+column-1;
      int icryPhi = ebId.iphi()+row-1;
      if ( EBDetId::validDetId(icryEta, icryPhi) ) {
        EBDetId id3x3 = EBDetId(icryEta, icryPhi, EBDetId::ETAPHIMODE);
        for( clusItr=ids.begin(); clusItr!=ids.end(); ++clusItr ) {
          EBDetId idInCluster = (EBDetId)clusItr->first;
          if ( idInCluster == id3x3 ) {
            nCry3x3++;
            break;
          }
        }
      }
    }
  } else if(seedId.subdetId() == EcalEndcap) {
    EEDetId eeId = (EEDetId)seedId;
    for(int icry=0; icry<9; ++icry) {
      unsigned int row    = icry/3;
      unsigned int column = icry%3;
      int icryX = eeId.ix()+column-1;
      int icryY = eeId.iy()+row-1;
      int iz = (eeId.positiveZ()) ? 1 : -1;
      if ( EEDetId::validDetId(icryX, icryY, iz) ) {
        EEDetId id3x3 = EEDetId(icryX, icryY, iz, EEDetId::XYMODE);
        for( clusItr=ids.begin(); clusItr!=ids.end(); ++clusItr ) {
          if ( clusItr->first == id3x3 ) {
            nCry3x3++;
            break;
          }
        }
      }
    }
  }

  privateData_->nOverlap3x3->push_back(nCry3x3);

}







void CmsBasicClusterFiller::treeBCInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"nCrystals"+colSuffix).c_str(), *privateData_->nCrystals, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nOverlap3x3"+colSuffix).c_str(), *privateData_->nOverlap3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedEnergy"+colSuffix).c_str(), *privateData_->seedEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
}



void CmsBasicClusterFillerData::initialiseCandidate() 
{
  nCrystals = new vector<int>; 
  nOverlap3x3 = new vector<int>;
  energy = new vector<float>; 
  eta = new vector<float>; 
  theta = new vector<float>; 
  phi = new vector<float>;;
  seedEnergy = new vector<float>;
  eMax = new vector<float>;
  e3x3 = new vector<float>;
  e5x5 = new vector<float>;
  nBC =  new int;
}

void CmsBasicClusterFillerData::clear() 
{
  nCrystals->clear();
  nOverlap3x3->clear();
  energy->clear();
  eta->clear();
  theta->clear();
  phi->clear();
  seedEnergy->clear();
  eMax->clear();
  e3x3->clear();
  e5x5->clear();
}
