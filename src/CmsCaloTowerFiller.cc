//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsCaloTowerFiller
//
// Original Author:  Christopher Rogan
//         Created:  Thur Feb  14 11:01:00 CEST 2008
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
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCaloTowerFiller.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoLocalCalo/CaloTowersCreator/interface/HcalMaterials.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoLocalCalo/CaloTowersCreator/interface/CaloTowersCreationAlgo.h"



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


CmsCaloTowerFiller::CmsCaloTowerFiller(CmsTree *cmsTree, 
				       edm::InputTag hbheLabel,
				       edm::InputTag hoLabel,
				       edm::InputTag hfLabel,
				       std::vector<edm::InputTag> ecalLabels,
				       int maxTracks, int maxMCTracks,
				       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsCaloTowerFillerData)
{
  cmstree=cmsTree;

  saveCaloTowerExtras_=true;

  trkIndexName_ = new std::string("n");

  hbheLabel_ = hbheLabel;
  hoLabel_ = hoLabel;
  hfLabel_ = hfLabel;
  ecalLabels_ = ecalLabels;

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsCaloTowerFiller::CmsCaloTowerFiller(CmsTree *cmsTree, 
				       edm::InputTag hbheLabel,
				       edm::InputTag hoLabel,
				       edm::InputTag hfLabel,
				       std::vector<edm::InputTag> ecalLabels,
				       bool fatTree, 
				       int maxTracks, int maxMCTracks,
				       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,fatTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsCaloTowerFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");

  hbheLabel_ = hbheLabel;
  hoLabel_ = hoLabel;
  hfLabel_ = hfLabel;
  ecalLabels_ = ecalLabels;
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsCaloTowerFiller::~CmsCaloTowerFiller() {

  // delete here the vector ptr's
  delete privateData_->energy;
  delete privateData_->x;
  delete privateData_->y;
  delete privateData_->z;
  delete privateData_->CALO;
  delete privateData_->CaloIndex;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsCaloTowerFiller::saveCaloTowerExtras(bool what) { saveCaloTowerExtras_=what; }

void CmsCaloTowerFiller::writeCollectionToTree(edm::InputTag caloTowersLabel,
					 const edm::Event& iEvent, const edm::EventSetup& iSetup,
					 const std::string &columnPrefix, const std::string &columnSuffix,
					 bool dumpData) {

  privateData_->clearTrkVectors();
  int icell = 0;

  Handle<CaloTowerCollection> collection;
  try { iEvent.getByLabel(caloTowersLabel, collection); }
  catch ( cms::Exception& ex ) {
    std::cout << "Did not found " << caloTowersLabel << std::endl;
  }

  if(collection->size()) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      LogInfo("CmsCaloTowerFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      LogInfo("CmsCaloTowerFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ 
			       << ". Collection will be truncated ";
    }
  
    //*(privateData_->ncand) = collection->size();

    edm::Handle<HBHERecHitCollection> hbhe;
    edm::Handle<HORecHitCollection> ho;
    edm::Handle<HFRecHitCollection> hf;
    edm::Handle<EcalRecHitCollection> eb;
    edm::Handle<EcalRecHitCollection> ee;
    edm::ESHandle<HcalTopology> htopo;
    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<IdealGeometryRecord>().get(htopo);
    iSetup.get<CaloGeometryRecord>().get(pG);
    const HcalTopology *theHcalTopology = htopo.product();
    const CaloGeometry *geo = pG.product();
    const CaloSubdetectorGeometry* theTowerGeometry = geo->getSubdetectorGeometry(DetId::Calo,CaloTowerDetId::SubdetId);
    
    if(saveCaloTowerExtras_) {
      try { iEvent.getByLabel(hbheLabel_,hbhe); }
      catch ( cms::Exception& ex ) { 
// 	LogWarning("CmsCaloTowerFiller") << "Can't get RecHit collection " 
// 				   << hbheLabel_;   
 	std::cout << "Can't get RecHit collection " 
		  << hbheLabel_ << std::endl;   

      }     
      try { iEvent.getByLabel(hoLabel_,ho); }
      catch ( cms::Exception& ex ) { 
	std::cout << "Can't get RecHit collection " 
		  << hoLabel_ << std::endl;   
      }
      try { iEvent.getByLabel(hfLabel_,hf); }
      catch ( cms::Exception& ex ) { 
	std::cout << "Can't get RecHit collection " 
		  << hfLabel_ << std::endl;   
      }
      std::vector<edm::InputTag>::const_iterator i;
      i=ecalLabels_.begin();
      try { iEvent.getByLabel(*i,eb); }
      catch ( cms::Exception& ex ) { 
	std::cout << "Can't get RecHit collection " 
		  << (*i) << std::endl;   
      }
      i++;
      try { iEvent.getByLabel(*i,ee); }
      catch ( cms::Exception& ex ) { 
	std::cout << "Can't get RecHit collection " 
		  << (*i) << std::endl;   
      }
    }


    int itower = -1;
    
    edm::View<reco::Candidate>::const_iterator cand;
    for(unsigned idx = 0; idx < collection->size(); idx++) {
      itower++;
      // fill basic kinematics
      //if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      // fill CaloTower extra informations
      if(saveCand_ && saveCaloTowerExtras_) { 
	//const CaloTower *thisCaloTower = dynamic_cast< const CaloTower * > ( &(*cand) );
	const CaloTower *thisCaloTower = &((*collection)[idx]);
	size_t nConst = thisCaloTower->constituentsSize();
	
	//privateData_->alpha->push_back(*jetVtxAlphaItr);
	GlobalPoint p = theTowerGeometry->getGeometry(thisCaloTower->id())->getPosition();
	privateData_->energy->push_back(thisCaloTower->energy());
	privateData_->x->push_back(p.x());
	privateData_->y->push_back(p.y());
	privateData_->z->push_back(p.z());
	privateData_->CALO->push_back(0);
	privateData_->CaloIndex->push_back(itower);
	icell++;
	for(unsigned int i = 0; i < nConst; i++){
	  DetId id = thisCaloTower->constituent((size_t)i);
	  DetId::Detector det = id.det();
	  
	  // ECAL CELL 
	  if(det == DetId::Ecal){
	    EcalSubdetector subdet = (EcalSubdetector)(id.subdetId());
	    if(subdet == EcalBarrel){
	      EcalRecHitCollection::const_iterator shit = eb->find(id);
	      GlobalPoint q = geo->getPosition(id);
	      privateData_->energy->push_back(shit->energy());
	      privateData_->x->push_back(q.x());
	      privateData_->y->push_back(q.y());
	      privateData_->z->push_back(q.z());
	      privateData_->CALO->push_back(100);
	      privateData_->CaloIndex->push_back(itower);
	      icell++;
	    }
	    if(subdet == EcalEndcap){
	      EcalRecHitCollection::const_iterator shit = ee->find(id);
	      GlobalPoint q = geo->getPosition(id);
	      privateData_->energy->push_back(shit->energy());
	      privateData_->x->push_back(q.x());
	      privateData_->y->push_back(q.y());
	      privateData_->z->push_back(q.z());
	      privateData_->CALO->push_back(200);
	      privateData_->CaloIndex->push_back(itower);
	      icell++;
	    }
	    continue;
	  }
	  // HCAL CELL
	  if(det == DetId::Hcal){
	    HcalDetId hcalDetId(id);
	    HcalSubdetector subdet = hcalDetId.subdet();

	    if(subdet == HcalBarrel){
	      HBHERecHitCollection::const_iterator shit = hbhe->find(id);
	      GlobalPoint q = geo->getPosition(id);
	      privateData_->energy->push_back(shit->energy());
	      privateData_->x->push_back(q.x());
	      privateData_->y->push_back(q.y());
	      privateData_->z->push_back(q.z());
	      privateData_->CALO->push_back(10+hcalDetId.depth());
	      privateData_->CaloIndex->push_back(itower);
	      icell++;
	    }
	    if(subdet == HcalEndcap){
	      //single or double tower?
	      int d = hcalDetId.depth();
	      int ieta = hcalDetId.ietaAbs();
	      if(hcalDetId.ietaAbs() < theHcalTopology->firstHEDoublePhiRing()){
		HBHERecHitCollection::const_iterator shit = hbhe->find(id);
		GlobalPoint q = geo->getPosition(id);
		if(d == 3 && ieta == 28){
		  privateData_->energy->push_back(shit->energy()/2.0);
		} else {
		  privateData_->energy->push_back(shit->energy());
		}
		privateData_->x->push_back(q.x());
		privateData_->y->push_back(q.y());
		privateData_->z->push_back(q.z());
		privateData_->CALO->push_back(20+hcalDetId.depth());
		privateData_->CaloIndex->push_back(itower);
		icell++;
	      } else {
		HBHERecHitCollection::const_iterator shit = hbhe->find(id);
		GlobalPoint q = geo->getPosition(id);
		if(d == 3 && ieta == 28){
		  privateData_->energy->push_back(shit->energy()/2.0);
		} else {
		  privateData_->energy->push_back(shit->energy());
		}
		privateData_->x->push_back(q.x());
		privateData_->y->push_back(q.y());
		privateData_->z->push_back(q.z());
		privateData_->CALO->push_back(30+hcalDetId.depth());
		privateData_->CaloIndex->push_back(itower);
		icell++;
	      }
	      
	    }
	    if(subdet == HcalOuter){
	      HORecHitCollection::const_iterator shit = ho->find(id);
	      GlobalPoint q = geo->getPosition(id);
	      privateData_->energy->push_back(shit->energy());
	      privateData_->x->push_back(q.x());
	      privateData_->y->push_back(q.y());
	      privateData_->z->push_back(q.z());
	      privateData_->CALO->push_back(40+hcalDetId.depth());
	      privateData_->CaloIndex->push_back(itower);
	      icell++;
	    }
	    if(subdet == HcalForward){
	      HFRecHitCollection::const_iterator shit = hf->find(id);
	      GlobalPoint q = geo->getPosition(id);
	      privateData_->energy->push_back(shit->energy());
	      privateData_->x->push_back(q.x());
	      privateData_->y->push_back(q.y());
	      privateData_->z->push_back(q.z());
	      privateData_->CALO->push_back(50+hcalDetId.depth());
	      privateData_->CaloIndex->push_back(itower);
	      icell++;
	      
	    }
	    
	    continue;
	  }
	  std::cout << "BAD CELL det " << det << endl;
	}

	// privateData_->alpha->push_back(*jetVtxAlphaItr);
// 	jetVtxAlphaItr++;
	
// 	// em, had fractions
// 	const CaloJet *thisRecoJet = dynamic_cast< const CaloJet * > ( &(*cand) );
// 	if( thisRecoJet != 0 ) { 
// 	  privateData_->emFrac->push_back( thisRecoJet->emEnergyFraction() );
// 	  privateData_->hadFrac->push_back( thisRecoJet->energyFractionHadronic() );
// 	}
// 	else {
// 	  privateData_->emFrac->push_back( -1. );
// 	  privateData_->hadFrac->push_back( -1. );
// 	}
      }
      else {
	// privateData_->alpha->push_back( -1. );
// 	privateData_->emFrac->push_back( -1. );
// 	privateData_->hadFrac->push_back( -1. );
      }
    }
  }
  else {
  *(privateData_->ncand) = 0;
  }
  *(privateData_->ncand) = icell;
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  //int blockSize = (collection) ? collection->size() : 0;
  int blockSize = (collection->size()) ? icell : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  //  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  if(saveCaloTowerExtras_) treeCaloTowerInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}









void CmsCaloTowerFiller::treeCaloTowerInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"x"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"y"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"z"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"CALO"+colSuffix).c_str(), *privateData_->CALO, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"CaloIndex"+colSuffix).c_str(), *privateData_->CaloIndex, nCandString.c_str(), 0, "Reco");

}







void CmsCaloTowerFillerData::initialise() {

  initialiseCandidate();
  energy = new vector<float>;
  x = new vector<float>;
  y = new vector<float>;
  z = new vector<float>;
  CALO = new vector<int>;
  CaloIndex = new vector<int>;

}

void CmsCaloTowerFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  energy->clear();
  x->clear();
  y->clear();
  z->clear();
  CALO->clear();
  CaloIndex->clear();

}
