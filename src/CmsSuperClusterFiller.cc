//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsSuperClusterFiller
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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsSuperClusterFiller.h"

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


CmsSuperClusterFiller::CmsSuperClusterFiller(CmsTree *cmsTree, int maxSC):  privateData_(new CmsSuperClusterFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");
  maxSC_=maxSC;
  privateData_->initialiseCandidate();

}

//--------------
// Destructor --
//--------------

CmsSuperClusterFiller::~CmsSuperClusterFiller() 
{
  // delete here the vector ptr's
  delete privateData_->nBC;
  delete privateData_->nCrystals;
  delete privateData_->rawEnergy;
  delete privateData_->energy;
  delete privateData_->esEnergy;
  delete privateData_->seedEnergy;
  delete privateData_->seedX;
  delete privateData_->seedY;
  delete privateData_->eta;
  delete privateData_->theta;
  delete privateData_->phi;
  delete privateData_->phiWidth;
  delete privateData_->etaWidth;
  delete privateData_->e3x3;
  delete privateData_->e5x5;
  delete privateData_->eMax;
  delete privateData_->e2x2;
  delete privateData_->e2nd;
  delete privateData_->hOverE;
  delete privateData_->covIEtaIEta;
  delete privateData_->covIEtaIPhi;
  delete privateData_->covIPhiIPhi;
  delete privateData_->sMaj;
  delete privateData_->sMin;
  delete privateData_->alpha;
  delete privateData_->e1x5;
  delete privateData_->e2x5Max;
  delete privateData_->e4SwissCross;
  delete privateData_->time;
  delete privateData_->chi2;
  delete privateData_->recoFlag;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out


void CmsSuperClusterFiller::writeCollectionToTree(edm::InputTag collectionTag,
						  const edm::Event& iEvent, const edm::EventSetup& iSetup,
						  const std::string &columnPrefix, const std::string &columnSuffix,
						  bool dumpData)
{
  
  Handle<SuperClusterCollection> collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get SC Collection: " << collectionTag; }
  const SuperClusterCollection *collection = collectionHandle.product();

  privateData_->clear();
  
  if(collection) 
    {
      if((int)collection->size() > maxSC_)
	{
	  edm::LogError("CmsSuperClusterFiller") << "Track length " << collection->size() 
						 << " is too long for declared max length for tree "
						 << maxSC_ 
						 << ". Collection will be truncated ";
	}
      
      *(privateData_->nSC) = collection->size();
  
      try { iEvent.getByLabel(Calotowers_, calotowers_); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get primary calotowers collection" << Calotowers_; }

      // for cluster shape variables
      Handle< EcalRecHitCollection > EcalBarrelRecHits;
      try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
      const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();
      
      Handle< EcalRecHitCollection > EcalEndcapRecHits;
      try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
      const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();
      
      SuperClusterCollection::const_iterator cand;
      for(cand=collection->begin(); cand!=collection->end(); cand++) {
        // fill basic kinematics
        writeSCInfo(&(*cand),iEvent,iSetup,EBRecHits,EERecHits);
      }
    }
  else {
    *(privateData_->nSC) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  
  int blockSize = (collection) ? collection->size() : 0;
    
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  treeSCInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

}






void CmsSuperClusterFiller::writeSCInfo(const SuperCluster *cand, 
                                        const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                        const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits) {

  // fill the SC infos
  privateData_->nBC->push_back((int)cand->clustersSize());
  privateData_->nCrystals->push_back((int)cand->hitsAndFractions().size());
  privateData_->rawEnergy->push_back((float)cand->rawEnergy());
  privateData_->energy->push_back((float)cand->energy());
  privateData_->esEnergy->push_back((float)cand->preshowerEnergy());
  privateData_->phiWidth->push_back((float)cand->phiWidth());
  privateData_->etaWidth->push_back((float)cand->etaWidth());
  privateData_->eta->push_back((float)cand->position().eta());
  privateData_->theta->push_back((float)cand->position().theta());
  privateData_->phi->push_back((float)cand->position().phi());
  
  // fill the seed basic cluster shapes
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);

  if ( pTopology.isValid() && pGeometry.isValid() ) {
    
    const CaloTopology *topology = pTopology.product();

    const EcalRecHitCollection *rechits = 0;

    // seed crystal properties
    const Ptr<CaloCluster> theSeed = cand->seed();

    float seedEta = theSeed->position().eta();

    if( fabs(seedEta) < 1.479 ) rechits = EBRecHits;
    else rechits = EERecHits; 

      float eMax = EcalClusterTools::eMax( *theSeed, &(*rechits) );
      float e3x3 = EcalClusterTools::e3x3( *theSeed, &(*rechits), topology );
      float e5x5 = EcalClusterTools::e5x5( *theSeed, &(*rechits), topology );
      float e2x2 = EcalClusterTools::e2x2( *theSeed, &(*rechits), topology );
      float e2nd = EcalClusterTools::e2nd( *theSeed, &(*rechits) );
      float e1x5 = EcalClusterTools::e1x5( *theSeed, &(*rechits), topology );
      float e2x5Max = EcalClusterTools::e2x5Max( *theSeed, &(*rechits), topology );
      float e4SwissCross = ( EcalClusterTools::eLeft( *theSeed, &(*rechits), topology ) +
                             EcalClusterTools::eRight( *theSeed, &(*rechits), topology ) +
                             EcalClusterTools::eTop( *theSeed, &(*rechits), topology ) +
                             EcalClusterTools::eBottom( *theSeed, &(*rechits), topology ) );

      privateData_->e3x3->push_back(e3x3);
      privateData_->e5x5->push_back(e5x5);
      privateData_->eMax->push_back(eMax);
      privateData_->e2x2->push_back(e2x2);
      privateData_->e2nd->push_back(e2nd);
      privateData_->e1x5->push_back(e1x5);
      privateData_->e2x5Max->push_back(e2x5Max);
      privateData_->e4SwissCross->push_back(e4SwissCross);

      // local covariances: instead of using absolute eta/phi it counts crystals normalised
      std::vector<float> vLocCov = EcalClusterTools::localCovariances( *theSeed, &(*rechits), topology );
      
      float covIEtaIEta = vLocCov[0];
      float covIEtaIPhi = vLocCov[1];
      float covIPhiIPhi = vLocCov[2];
      
      privateData_->covIEtaIEta->push_back(covIEtaIEta);
      privateData_->covIEtaIPhi->push_back(covIEtaIPhi);
      privateData_->covIPhiIPhi->push_back(covIPhiIPhi);

      // seed second moments wrt principal axes:
      Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*theSeed, *rechits );
      privateData_->sMaj->push_back(moments.sMaj);
      privateData_->sMin->push_back(moments.sMin);
      // angle between sMaj and phi direction:
      privateData_->alpha->push_back(moments.alpha);



      std::pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *theSeed, &(*rechits) );
      DetId seedCrystalId = maxRH.first;
      EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
      
      privateData_->time->push_back((float)seedRH->time());
      privateData_->chi2->push_back((float)seedRH->chi2());
      privateData_->recoFlag->push_back((int)seedRH->recoFlag());
      privateData_->seedEnergy->push_back((float)maxRH.second);

      if(EcalSubdetector(seedCrystalId.subdetId()) == EcalBarrel) {
        EBDetId id(seedCrystalId);
        privateData_->seedX->push_back(id.ieta());
        privateData_->seedY->push_back(id.iphi());
      } else {
        EEDetId id(seedCrystalId);
        privateData_->seedX->push_back(id.ix());
        privateData_->seedY->push_back(id.iy());        
      }

  } else {
    privateData_->e3x3->push_back(-1.);
    privateData_->e5x5->push_back(-1.);
    privateData_->eMax->push_back(-1.);
    privateData_->e2x2->push_back(-1.);
    privateData_->e2nd->push_back(-1.);
    privateData_->covIEtaIEta->push_back(-1.);
    privateData_->covIEtaIPhi->push_back(-1.);
    privateData_->covIPhiIPhi->push_back(-1.);
    privateData_->sMaj->push_back(-1.);
    privateData_->sMin->push_back(-1.);
    privateData_->alpha->push_back(-1.);
    privateData_->e1x5->push_back(-1);
    privateData_->e2x5Max->push_back(-1);
    privateData_->e4SwissCross->push_back(-1);
    privateData_->time->push_back(-999.);
    privateData_->chi2->push_back(-999.);
    privateData_->recoFlag->push_back(-1);
    privateData_->seedEnergy->push_back(-1.);
  }

  // calculate H/E
  float hOverEConeSize = 0.15;
  float hOverEPtMin = 0.;
  EgammaTowerIsolation *towerIso1 = new EgammaTowerIsolation(hOverEConeSize,0.,hOverEPtMin,1,calotowers_.product()) ;
  EgammaTowerIsolation *towerIso2 = new EgammaTowerIsolation(hOverEConeSize,0.,hOverEPtMin,2,calotowers_.product()) ;
  
  float TowerHcalESum1 = towerIso1->getTowerESum(cand);
  float TowerHcalESum2 = towerIso2->getTowerESum(cand);
  float hcalESum = TowerHcalESum1 + TowerHcalESum2;
  
  privateData_->hOverE->push_back(hcalESum/cand->energy());

  delete towerIso1;
  delete towerIso2;

}

void CmsSuperClusterFiller::treeSCInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"nBC"+colSuffix).c_str(), *privateData_->nBC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nCrystals"+colSuffix).c_str(), *privateData_->nCrystals, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"rawEnergy"+colSuffix).c_str(), *privateData_->rawEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esEnergy"+colSuffix).c_str(), *privateData_->esEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiWidth"+colSuffix).c_str(), *privateData_->phiWidth, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"etaWidth"+colSuffix).c_str(), *privateData_->etaWidth, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x2"+colSuffix).c_str(), *privateData_->e2x2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2nd"+colSuffix).c_str(), *privateData_->e2nd, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e1x5"+colSuffix).c_str(), *privateData_->e1x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Max"+colSuffix).c_str(), *privateData_->e2x5Max, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e4SwissCross"+colSuffix).c_str(), *privateData_->e4SwissCross, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIEtaIEta"+colSuffix).c_str(), *privateData_->covIEtaIEta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIEtaIPhi"+colSuffix).c_str(), *privateData_->covIEtaIPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIPhiIPhi"+colSuffix).c_str(), *privateData_->covIPhiIPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"sMaj"+colSuffix).c_str(), *privateData_->sMaj, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"sMin"+colSuffix).c_str(), *privateData_->sMin, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"alpha"+colSuffix).c_str(), *privateData_->alpha, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverE"+colSuffix).c_str(), *privateData_->hOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlag"+colSuffix).c_str(), *privateData_->recoFlag, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"time"+colSuffix).c_str(), *privateData_->time, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2"+colSuffix).c_str(), *privateData_->chi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedEnergy"+colSuffix).c_str(), *privateData_->seedEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedX"+colSuffix).c_str(), *privateData_->seedX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedY"+colSuffix).c_str(), *privateData_->seedY, nCandString.c_str(), 0, "Reco");
}


void CmsSuperClusterFillerData::initialiseCandidate() 
{
  nBC = new vector<int>;
  nCrystals = new vector<int>; 
  rawEnergy = new vector<float>; 
  energy = new vector<float>; 
  esEnergy = new vector<float>; 
  eta = new vector<float>; 
  theta = new vector<float>; 
  phi = new vector<float>;
  phiWidth = new vector<float>;
  etaWidth = new vector<float>;
  e3x3 = new vector<float>;
  e5x5 = new vector<float>;
  eMax = new vector<float>;
  e2x2 = new vector<float>;
  e2nd = new vector<float>;
  e1x5 = new vector<float>;
  e2x5Max = new vector<float>;
  e4SwissCross = new vector<float>;
  hOverE = new vector<float>;
  covIEtaIEta = new vector<float>;
  covIEtaIPhi = new vector<float>;
  covIPhiIPhi = new vector<float>;
  sMaj = new vector<float>;
  sMin = new vector<float>;
  alpha = new vector<float>;
  recoFlag = new vector<int>;
  time = new vector<float>;
  chi2 = new vector<float>;
  seedEnergy = new vector<float>;
  seedX = new vector<float>;
  seedY = new vector<float>;
  nSC =  new int;
}

void CmsSuperClusterFillerData::clear() 
{
  nBC->clear();
  nCrystals->clear();
  rawEnergy->clear();
  energy->clear();
  esEnergy->clear();
  eta->clear(); 
  theta->clear();
  phi->clear();
  phiWidth->clear();
  etaWidth->clear();
  e3x3->clear();
  e5x5->clear();
  eMax->clear();
  e2x2->clear();
  e2nd->clear();
  e1x5->clear();
  e2x5Max->clear();
  e4SwissCross->clear();
  hOverE->clear();
  covIEtaIEta->clear();
  covIEtaIPhi->clear();
  covIPhiIPhi->clear();
  sMaj->clear();
  sMin->clear();
  alpha->clear();
  recoFlag->clear();
  time->clear();
  chi2->clear();
  seedEnergy->clear();
  seedX->clear();
  seedY->clear();
}
