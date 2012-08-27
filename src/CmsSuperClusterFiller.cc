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
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsSuperClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/ESHelper.h"

#include "HiggsAnalysis/HiggsToGammaGamma/interface/GBRForest.h"

#include <TTree.h>

#include <string>
#include <iostream>

using namespace edm;
using namespace reco;
using namespace std;



//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsSuperClusterFiller::CmsSuperClusterFiller(CmsTree *cmsTree, int maxSC):  privateData_(new CmsSuperClusterFillerData), eCorrector_(0), pCorrector_(0), photonFixE_(0), photonFixP_(0)
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
  delete privateData_->seedClusterEnergy;
  delete privateData_->seedEnergy;
  delete privateData_->seedX;
  delete privateData_->seedY;
  delete privateData_->xPos;
  delete privateData_->yPos;
  delete privateData_->zPos;
  delete privateData_->eta;
  delete privateData_->theta;
  delete privateData_->phi;
  delete privateData_->phiWidth;
  delete privateData_->etaWidth;
  delete privateData_->e3x3;
  delete privateData_->e3x1;
  delete privateData_->e1x3;
  delete privateData_->e4x4;
  delete privateData_->e5x5;
  delete privateData_->eMax;
  delete privateData_->e2x2;
  delete privateData_->e2nd;
  delete privateData_->e2x5Left;
  delete privateData_->e2x5Right;
  delete privateData_->e2x5Top;
  delete privateData_->e2x5Bottom;
  delete privateData_->eLeft;
  delete privateData_->eRight;
  delete privateData_->eTop;
  delete privateData_->eBottom;
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
  delete privateData_->esEffsIxIx;
  delete privateData_->esEffsIyIy;
  delete privateData_->esL1Energy;
  delete privateData_->esL2Energy;
  delete privateData_->esL1Strips;
  delete privateData_->esL2Strips;
//   delete privateData_->etaC;
//   delete privateData_->etaS;
//   delete privateData_->etaM;
//   delete privateData_->phiC;
//   delete privateData_->phiS;
//   delete privateData_->phiM;
//   delete privateData_->xC;
//   delete privateData_->xS;
//   delete privateData_->xM;
//   delete privateData_->xZ;
//   delete privateData_->yC;
//   delete privateData_->yS;
//   delete privateData_->yM;
//   delete privateData_->yZ;
  delete privateData_->photonFix_phoE;
  delete privateData_->photonFix_phoSigma;
  delete privateData_->photonFix_eleE;
  delete privateData_->photonFix_eleSigma;
  delete privateData_->regrCorr_phoE;
  delete privateData_->regrCorr_phoSigma;
  delete privateData_->regrCorr_eleE;
  delete privateData_->regrCorr_eleSigma;
  delete privateData_;
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
  
      if(Calotowers_.label().size()!=0) {
        try { iEvent.getByLabel(Calotowers_, calotowers_); }
        catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get primary calotowers collection" << Calotowers_; }
      }
      // for cluster shape variables
      Handle< EcalRecHitCollection > EcalBarrelRecHits;
      try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
      const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();
      
      Handle< EcalRecHitCollection > EcalEndcapRecHits;
      try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
      const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();

      Handle< EcalRecHitCollection > EcalPreshowerRecHits;
      try { iEvent.getByLabel(ESRecHits_, EcalPreshowerRecHits); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get ECAL preshower rec hits Collection" << ESRecHits_; }
      const EcalRecHitCollection *ESRecHits = EcalPreshowerRecHits.product();
      
      SuperClusterCollection::const_iterator cand;
      for(cand=collection->begin(); cand!=collection->end(); cand++) {
        // fill basic kinematics
        writeSCInfo(&(*cand),iEvent,iSetup,EBRecHits,EERecHits,ESRecHits);
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

  delete trkIndexName_;

}






void CmsSuperClusterFiller::writeSCInfo(const SuperCluster *cand, 
                                        const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                        const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits,
                                        const EcalRecHitCollection *ESRecHits) {

  // fill the SC infos
  privateData_->nBC->push_back((int)cand->clustersSize());
  privateData_->nCrystals->push_back((int)cand->hitsAndFractions().size());
  privateData_->rawEnergy->push_back((float)cand->rawEnergy());
  privateData_->seedClusterEnergy->push_back((float)cand->seed()->energy());
  privateData_->energy->push_back((float)cand->energy());
  privateData_->esEnergy->push_back((float)cand->preshowerEnergy());
  privateData_->phiWidth->push_back((float)cand->phiWidth());
  privateData_->etaWidth->push_back((float)cand->etaWidth());
  privateData_->eta->push_back((float)cand->position().eta());
  privateData_->theta->push_back((float)cand->position().theta());
  privateData_->phi->push_back((float)cand->position().phi());
  privateData_->xPos->push_back((float)cand->position().x());
  privateData_->yPos->push_back((float)cand->position().y());
  privateData_->zPos->push_back((float)cand->position().z());
  
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
      float e3x1 = EcalClusterTools::e3x1( *theSeed, &(*rechits), topology );
      float e1x3 = EcalClusterTools::e1x3( *theSeed, &(*rechits), topology );
      float e4x4 = EcalClusterTools::e4x4( *theSeed, &(*rechits), topology );
      float e5x5 = EcalClusterTools::e5x5( *theSeed, &(*rechits), topology );
      float e2x2 = EcalClusterTools::e2x2( *theSeed, &(*rechits), topology );
      float e2nd = EcalClusterTools::e2nd( *theSeed, &(*rechits) );
      float e1x5 = EcalClusterTools::e1x5( *theSeed, &(*rechits), topology );
      float e2x5Max = EcalClusterTools::e2x5Max( *theSeed, &(*rechits), topology );
      float e2x5Left = EcalClusterTools::e2x5Left( *theSeed, &(*rechits), topology );
      float e2x5Right = EcalClusterTools::e2x5Right( *theSeed, &(*rechits), topology );
      float e2x5Top = EcalClusterTools::e2x5Top( *theSeed, &(*rechits), topology );
      float e2x5Bottom = EcalClusterTools::e2x5Bottom( *theSeed, &(*rechits), topology );
      float eLeft = EcalClusterTools::eLeft( *theSeed, &(*rechits), topology );
      float eRight = EcalClusterTools::eRight( *theSeed, &(*rechits), topology );
      float eTop = EcalClusterTools::eTop( *theSeed, &(*rechits), topology );
      float eBottom = EcalClusterTools::eBottom( *theSeed, &(*rechits), topology );
      float e4SwissCross = ( EcalClusterTools::eLeft( *theSeed, &(*rechits), topology ) +
                             EcalClusterTools::eRight( *theSeed, &(*rechits), topology ) +
                             EcalClusterTools::eTop( *theSeed, &(*rechits), topology ) +
                             EcalClusterTools::eBottom( *theSeed, &(*rechits), topology ) );

      privateData_->e3x3->push_back(e3x3);
      privateData_->e3x1->push_back(e3x1);
      privateData_->e1x3->push_back(e1x3);
      privateData_->e4x4->push_back(e4x4);
      privateData_->e5x5->push_back(e5x5);
      privateData_->eMax->push_back(eMax);
      privateData_->e2x2->push_back(e2x2);
      privateData_->e2nd->push_back(e2nd);
      privateData_->e1x5->push_back(e1x5);
      privateData_->e2x5Max->push_back(e2x5Max);
      privateData_->e2x5Left->push_back(e2x5Left);
      privateData_->e2x5Right->push_back(e2x5Right);
      privateData_->e2x5Top->push_back(e2x5Top);
      privateData_->e2x5Bottom->push_back(e2x5Bottom);
      privateData_->eLeft->push_back(eLeft);
      privateData_->eRight->push_back(eRight);
      privateData_->eTop->push_back(eTop);
      privateData_->eBottom->push_back(eBottom);
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

      // calculate H/E
      float hcalESum = 0.;
      if(Calotowers_.label().size()!=0) { 
        float hOverEConeSize = 0.15;
        float hOverEPtMin = 0.;
        EgammaTowerIsolation *towerIso1 = new EgammaTowerIsolation(hOverEConeSize,0.,hOverEPtMin,1,calotowers_.product()) ;
        EgammaTowerIsolation *towerIso2 = new EgammaTowerIsolation(hOverEConeSize,0.,hOverEPtMin,2,calotowers_.product()) ;
        
        float TowerHcalESum1 = towerIso1->getTowerESum(cand);
        float TowerHcalESum2 = towerIso2->getTowerESum(cand);
        hcalESum = TowerHcalESum1 + TowerHcalESum2;
        
        privateData_->hOverE->push_back(hcalESum/cand->energy());
        delete towerIso1;
        delete towerIso2;
      } else privateData_->hOverE->push_back(-999.);

      // === Preshower cluster shapes === //
      // ES geometry
      const CaloSubdetectorGeometry *es_geometry = pGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
      const CaloSubdetectorGeometry *& es_geometry_p = es_geometry;

      CaloSubdetectorTopology *es_topology_p = 0;
      if (es_geometry) es_topology_p = new EcalPreshowerTopology(pGeometry);

      // make the map of rechits
      es_rechits_map_.clear();
      if (ESRecHits) {
        for (EcalRecHitCollection::const_iterator it = ESRecHits->begin(); it != ESRecHits->end(); ++it) {
          // remove bad ES rechits
          if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
          //Make the map of DetID, EcalRecHit pairs
          es_rechits_map_.insert(std::make_pair(it->id(), *it));   
        }
      }

      float ixix, iyiy;
      ixix = iyiy = 0;
      // calculate preshower cluster shape
      ESHelper *gES = new ESHelper(es_rechits_map_, es_geometry_p, es_topology_p);
      if (fabs(cand->eta()) > 1.6 && fabs(cand->eta()) < 3) {
        std::vector<float> elESHits0 = gES->getESHits(cand->x(), cand->y(), cand->z(), 0);
        std::vector<float> elESShape = gES->getESShape(elESHits0);
        ixix = elESShape[0];
        iyiy = elESShape[1];
      }

      delete gES;

      privateData_->esEffsIxIx->push_back(ixix);
      privateData_->esEffsIyIy->push_back(iyiy);

      // the preshowerClustersBegin() iteration only works for PF superclusters...
      float esL1Energy, esL2Energy;
      esL1Energy = esL2Energy = 0.;
      int esL1Strips, esL2Strips;
      esL1Strips = esL2Strips = 0;
      for (reco::CaloCluster_iterator it=cand->preshowerClustersBegin(); it != cand->preshowerClustersEnd(); it++) {
        if (ESDetId(((*it)->hitsAndFractions())[0].first).plane() == 1) {
          esL1Energy += (float) (*it)->energy();
          esL1Strips += (*it)->hitsAndFractions().size();
        } else {
          esL2Energy += (float) (*it)->energy();
          esL2Strips += (*it)->hitsAndFractions().size();
        }
      }

      privateData_->esL1Energy->push_back(esL1Energy);
      privateData_->esL2Energy->push_back(esL2Energy);
      privateData_->esL1Strips->push_back(esL1Strips);
      privateData_->esL2Strips->push_back(esL2Strips);

      if (photonFixE_ || photonFixP_ || eCorrector_ || pCorrector_ )
	{
      const CaloSubdetectorGeometry* subDetGeometry =0 ;
      if ( fabs(seedEta) < 1.479 )
	subDetGeometry =  pGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
      else
	subDetGeometry =  pGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
      math::XYZPoint unconvPos = posCalculator_->Calculate_Location(cand->seed()->hitsAndFractions(), &(*rechits) ,subDetGeometry,pGeometry->getSubdetectorGeometry(DetId::Ecal,EcalPreshower));
      float r9 =e3x3/(cand->rawEnergy());
      
      float phot_energy;
      math::XYZPoint phot_pos;
      if ( fabs(seedEta) < 1.479 && r9>0.94 )
	{
	  const reco::SuperCluster* sc=&(*cand);
	  double deltaE = energyCorrectionF->getValue(*sc, 1);
	  //	  double deltaE=0;
	  phot_energy=e5x5*(1.0 +  deltaE/cand->rawEnergy() );
	  phot_pos=unconvPos;
	}
      else if ( fabs(seedEta) > 1.479 && r9>0.95 )
	{
	  phot_energy=e5x5 + cand->preshowerEnergy();
	  phot_pos=unconvPos;
	}
      else
	{
	  phot_pos = cand->position();
	  phot_energy=cand->energy();
	}
      
      if (photonFixE_)
	{
	  photonFixE_->setClusterParameters(phot_energy,(float)phot_pos.eta(),(float)phot_pos.eta(),e3x3/(float)cand->rawEnergy());
 	  privateData_->photonFix_eleE->push_back(photonFixE_->fixedEnergy());
 	  privateData_->photonFix_eleSigma->push_back(photonFixE_->sigmaEnergy());
 	}
      else
	{
 	  privateData_->photonFix_eleE->push_back(-1);
 	  privateData_->photonFix_eleSigma->push_back(-1);
	}

      if (photonFixP_)
	{

	  photonFixP_->setClusterParameters(phot_energy,(float)phot_pos.eta(),(float)phot_pos.eta(),e3x3/(float)cand->rawEnergy());
 	  privateData_->photonFix_phoE->push_back(photonFixP_->fixedEnergy());
 	  privateData_->photonFix_phoSigma->push_back(photonFixP_->sigmaEnergy());
	}
      else
	{
 	  privateData_->photonFix_phoE->push_back(-1);
 	  privateData_->photonFix_phoSigma->push_back(-1);
	}


      if (eCorrector_)
	{
	  photonFixE_->setClusterParameters(phot_energy,cand->position().eta(),cand->position().phi(),r9); 
	    
	  Bool_t isbarrel =  (fabs(cand->position().eta()) < 1.48);
	  
	  if (isbarrel) {
	    eCorrector_->fVals[0]  = cand->rawEnergy();
	    eCorrector_->fVals[1]  = e3x3/cand->rawEnergy(); //r9
	    eCorrector_->fVals[2]  = cand->position().eta();
	    eCorrector_->fVals[3]  = cand->position().phi();
	    eCorrector_->fVals[4]  = e5x5/cand->rawEnergy();
	    eCorrector_->fVals[5]  = photonFixE_->etaC();
	    eCorrector_->fVals[6]  = photonFixE_->etaS();
	    eCorrector_->fVals[7]  = photonFixE_->etaM();
	    eCorrector_->fVals[8]  = photonFixE_->phiC();
	    eCorrector_->fVals[9]  = photonFixE_->phiS();
	    eCorrector_->fVals[10] = photonFixE_->phiM();    
	    eCorrector_->fVals[11] = hcalESum/cand->energy();
	    eCorrector_->fVals[12] = cand->etaWidth();
	    eCorrector_->fVals[13] = cand->phiWidth();
	    eCorrector_->fVals[14] = sqrt(covIEtaIEta);
	  }
	  else {
	    eCorrector_->fVals[0]  = cand->rawEnergy();
	    eCorrector_->fVals[1]  = e3x3/cand->rawEnergy(); //r9
	    eCorrector_->fVals[2]  = cand->position().eta();
	    eCorrector_->fVals[3]  = cand->position().phi();
	    eCorrector_->fVals[4]  = e5x5/cand->rawEnergy();
	    eCorrector_->fVals[5]  = cand->preshowerEnergy()/cand->rawEnergy();
	    eCorrector_->fVals[6]  = photonFixE_->xZ();
	    eCorrector_->fVals[7]  = photonFixE_->xC();
	    eCorrector_->fVals[8]  = photonFixE_->xS();
	    eCorrector_->fVals[9]  = photonFixE_->xM();
	    eCorrector_->fVals[10] = photonFixE_->yZ();
	    eCorrector_->fVals[11] = photonFixE_->yC();
	    eCorrector_->fVals[12] = photonFixE_->yS();
	    eCorrector_->fVals[13] = photonFixE_->yM();
	    eCorrector_->fVals[14] = hcalESum/cand->energy();
	    eCorrector_->fVals[15] = cand->etaWidth();
	    eCorrector_->fVals[16] = cand->phiWidth();
	    eCorrector_->fVals[17] = sqrt(covIEtaIEta);
	  }
	  
	  const Double_t varscale = 1.253;
	  Double_t den;
	  const GBRForest *reader;
	  const GBRForest *readervar;
	  if (isbarrel) {
	    den = cand->rawEnergy();
	    reader = eCorrector_->fReadereb;
	    readervar = eCorrector_->fReaderebvariance;
	  }
	  else {
	    den = cand->rawEnergy() + cand->preshowerEnergy();
	    reader = eCorrector_->fReaderee;
	    readervar = eCorrector_->fReadereevariance;
	  }
	  
	  privateData_->regrCorr_eleE->push_back(reader->GetResponse(eCorrector_->fVals)*den);
	  privateData_->regrCorr_eleSigma->push_back(readervar->GetResponse(eCorrector_->fVals)*den*varscale);
	}
      else
	{
	  privateData_->regrCorr_eleE->push_back(-1);
	  privateData_->regrCorr_eleSigma->push_back(-1);
	}

      if (pCorrector_)
	{
	  photonFixP_->setClusterParameters(phot_energy,phot_pos.eta(),phot_pos.phi(),r9); 
	    
	  Bool_t isbarrel = (fabs(phot_pos.eta()) < 1.48); //to comply with PhotonFix asserts
	  
	  if (isbarrel) {
	    pCorrector_->fVals[0]  = cand->rawEnergy();
	    pCorrector_->fVals[1]  = e3x3/cand->rawEnergy(); //r9
	    pCorrector_->fVals[2]  = cand->position().eta();
	    pCorrector_->fVals[3]  = cand->position().phi();
	    pCorrector_->fVals[4]  = e5x5/cand->rawEnergy();
	    pCorrector_->fVals[5]  = photonFixP_->etaC();
	    pCorrector_->fVals[6]  = photonFixP_->etaS();
	    pCorrector_->fVals[7]  = photonFixP_->etaM();
	    pCorrector_->fVals[8]  = photonFixP_->phiC();
	    pCorrector_->fVals[9]  = photonFixP_->phiS();
	    pCorrector_->fVals[10] = photonFixP_->phiM();    
	    pCorrector_->fVals[11] = hcalESum/cand->energy();
	    pCorrector_->fVals[12] = cand->etaWidth();
	    pCorrector_->fVals[13] = cand->phiWidth();
	    pCorrector_->fVals[14] = sqrt(covIEtaIEta);
	  }
	  else {
	    pCorrector_->fVals[0]  = cand->rawEnergy();
	    pCorrector_->fVals[1]  = e3x3/cand->rawEnergy(); //r9
	    pCorrector_->fVals[2]  = cand->position().eta();
	    pCorrector_->fVals[3]  = cand->position().phi();
	    pCorrector_->fVals[4]  = e5x5/cand->rawEnergy();
	    pCorrector_->fVals[5]  = cand->preshowerEnergy()/cand->rawEnergy();
	    pCorrector_->fVals[6]  = photonFixP_->xZ();
	    pCorrector_->fVals[7]  = photonFixP_->xC();
	    pCorrector_->fVals[8]  = photonFixP_->xS();
	    pCorrector_->fVals[9]  = photonFixP_->xM();
	    pCorrector_->fVals[10] = photonFixP_->yZ();
	    pCorrector_->fVals[11] = photonFixP_->yC();
	    pCorrector_->fVals[12] = photonFixP_->yS();
	    pCorrector_->fVals[13] = photonFixP_->yM();
	    pCorrector_->fVals[14] = hcalESum/cand->energy();
	    pCorrector_->fVals[15] = cand->etaWidth();
	    pCorrector_->fVals[16] = cand->phiWidth();
	    pCorrector_->fVals[17] = sqrt(covIEtaIEta);
	  }
	  
	  const Double_t varscale = 1.253;
	  Double_t den;
	  const GBRForest *reader;
	  const GBRForest *readervar;
	  if (isbarrel) {
	    den = cand->rawEnergy();
	    reader = pCorrector_->fReadereb;
	    readervar = pCorrector_->fReaderebvariance;
	  }
	  else {
	    den = cand->rawEnergy() + cand->preshowerEnergy();
	    reader = pCorrector_->fReaderee;
	    readervar = pCorrector_->fReadereevariance;
	  }
	  
	  privateData_->regrCorr_phoE->push_back(reader->GetResponse(pCorrector_->fVals)*den);
	  privateData_->regrCorr_phoSigma->push_back(readervar->GetResponse(pCorrector_->fVals)*den*varscale);
	}
      else
	{
	  privateData_->regrCorr_phoE->push_back(-1);
	  privateData_->regrCorr_phoSigma->push_back(-1);
	}

        }

  } else {
    privateData_->e3x3->push_back(-1.);
    privateData_->e3x1->push_back(-1.);
    privateData_->e1x3->push_back(-1.);
    privateData_->e4x4->push_back(-1.);
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
    privateData_->photonFix_phoE->push_back(-1);
    privateData_->photonFix_phoSigma->push_back(-1);

    privateData_->photonFix_eleE->push_back(-1);
    privateData_->photonFix_eleSigma->push_back(-1);

    privateData_->regrCorr_eleE->push_back(-1);
    privateData_->regrCorr_eleSigma->push_back(-1);
    privateData_->regrCorr_phoE->push_back(-1);
    privateData_->regrCorr_phoSigma->push_back(-1);
  }

}

void CmsSuperClusterFiller::treeSCInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"nBC"+colSuffix).c_str(), *privateData_->nBC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nCrystals"+colSuffix).c_str(), *privateData_->nCrystals, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"rawEnergy"+colSuffix).c_str(), *privateData_->rawEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedClusterEnergy"+colSuffix).c_str(), *privateData_->seedClusterEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esEnergy"+colSuffix).c_str(), *privateData_->esEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"xPos"+colSuffix).c_str(), *privateData_->xPos, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"yPos"+colSuffix).c_str(), *privateData_->yPos, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"zPos"+colSuffix).c_str(), *privateData_->zPos, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiWidth"+colSuffix).c_str(), *privateData_->phiWidth, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"etaWidth"+colSuffix).c_str(), *privateData_->etaWidth, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x1"+colSuffix).c_str(), *privateData_->e3x1, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e1x3"+colSuffix).c_str(), *privateData_->e1x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e4x4"+colSuffix).c_str(), *privateData_->e4x4, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x2"+colSuffix).c_str(), *privateData_->e2x2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2nd"+colSuffix).c_str(), *privateData_->e2nd, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e1x5"+colSuffix).c_str(), *privateData_->e1x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Max"+colSuffix).c_str(), *privateData_->e2x5Max, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Left"+colSuffix).c_str(), *privateData_->e2x5Left, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Right"+colSuffix).c_str(), *privateData_->e2x5Right, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Top"+colSuffix).c_str(), *privateData_->e2x5Top, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Bottom"+colSuffix).c_str(), *privateData_->e2x5Bottom, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eLeft"+colSuffix).c_str(), *privateData_->eLeft, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eRight"+colSuffix).c_str(), *privateData_->eRight, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eTop"+colSuffix).c_str(), *privateData_->eTop, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eBottom"+colSuffix).c_str(), *privateData_->eBottom, nCandString.c_str(), 0, "Reco");
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
  cmstree->column((colPrefix+"esEffsIxIx"+colSuffix).c_str(), *privateData_->esEffsIxIx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esEffsIyIy"+colSuffix).c_str(), *privateData_->esEffsIyIy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esL1Energy"+colSuffix).c_str(), *privateData_->esL1Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esL2Energy"+colSuffix).c_str(), *privateData_->esL2Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esL1Strips"+colSuffix).c_str(), *privateData_->esL1Strips, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esL2Strips"+colSuffix).c_str(), *privateData_->esL2Strips, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"etaC"+colSuffix).c_str(), *privateData_->etaC, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"etaS"+colSuffix).c_str(), *privateData_->etaS, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"etaM"+colSuffix).c_str(), *privateData_->etaM, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"phiC"+colSuffix).c_str(), *privateData_->phiC, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"phiS"+colSuffix).c_str(), *privateData_->phiS, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"phiM"+colSuffix).c_str(), *privateData_->phiM, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"xC"+colSuffix).c_str(), *privateData_->xC, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"xS"+colSuffix).c_str(), *privateData_->xS, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"xM"+colSuffix).c_str(), *privateData_->xM, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"xZ"+colSuffix).c_str(), *privateData_->xZ, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"yC"+colSuffix).c_str(), *privateData_->yC, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"yS"+colSuffix).c_str(), *privateData_->yS, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"yM"+colSuffix).c_str(), *privateData_->yM, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"yZ"+colSuffix).c_str(), *privateData_->yZ, nCandString.c_str(), 0, "Reco");
      if (photonFixE_ || photonFixP_ || eCorrector_ || pCorrector_ )
	{
	  cmstree->column((colPrefix+"photonFix_phoE"+colSuffix).c_str(), *privateData_->photonFix_phoE, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"photonFix_phoSigma"+colSuffix).c_str(), *privateData_->photonFix_phoSigma, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"photonFix_eleE"+colSuffix).c_str(), *privateData_->photonFix_eleE, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"photonFix_eleSigma"+colSuffix).c_str(), *privateData_->photonFix_eleSigma, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"regrCorr_phoE"+colSuffix).c_str(), *privateData_->regrCorr_phoE, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"regrCorr_phoSigma"+colSuffix).c_str(), *privateData_->regrCorr_phoSigma, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"regrCorr_eleE"+colSuffix).c_str(), *privateData_->regrCorr_eleE, nCandString.c_str(), 0, "Reco");
	  cmstree->column((colPrefix+"regrCorr_eleSigma"+colSuffix).c_str(), *privateData_->regrCorr_eleSigma, nCandString.c_str(), 0, "Reco");
	}
}


void CmsSuperClusterFillerData::initialiseCandidate() 
{
  nBC = new vector<int>;
  nCrystals = new vector<int>; 
  rawEnergy = new vector<float>; 
  energy = new vector<float>; 
  seedClusterEnergy = new vector<float>;
  esEnergy = new vector<float>; 
  esEffsIxIx = new vector<float>;
  esEffsIyIy = new vector<float>;
  esL1Energy = new vector<float>;
  esL2Energy = new vector<float>;
  esL1Strips = new vector<int>;
  esL2Strips = new vector<int>;
  eta = new vector<float>; 
  theta = new vector<float>; 
  phi = new vector<float>;
  xPos = new vector<float>;
  yPos = new vector<float>;
  zPos = new vector<float>;
  phiWidth = new vector<float>;
  etaWidth = new vector<float>;
  e3x3 = new vector<float>;
  e3x1 = new vector<float>;
  e1x3 = new vector<float>;
  e4x4 = new vector<float>;
  e5x5 = new vector<float>;
  eMax = new vector<float>;
  e2x2 = new vector<float>;
  e2nd = new vector<float>;
  e1x5 = new vector<float>;
  e2x5Max = new vector<float>;
  e2x5Left = new vector<float>;
  e2x5Right = new vector<float>;
  e2x5Top = new vector<float>;
  e2x5Bottom = new vector<float>;
  eLeft = new vector<float>;
  eRight = new vector<float>;
  eTop = new vector<float>;
  eBottom = new vector<float>;
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
  photonFix_phoE= new vector<float>;
  photonFix_phoSigma= new vector<float>;
  photonFix_eleE= new vector<float>;
  photonFix_eleSigma= new vector<float>;
  regrCorr_phoE= new vector<float>;
  regrCorr_phoSigma= new vector<float>;
  regrCorr_eleE= new vector<float>;
  regrCorr_eleSigma= new vector<float>;  
  nSC =  new int;
//   etaC= new vector<float>;
//   etaS= new vector<float>;
//   etaM= new vector<float>;
//   phiC= new vector<float>;
//   phiS= new vector<float>;
//   phiM= new vector<float>;
//   xC= new vector<float>;
//   xS= new vector<float>;
//   xM= new vector<float>;
//   xZ= new vector<float>;
//   yC= new vector<float>;
//   yS= new vector<float>;
//   yM= new vector<float>;
//   yZ= new vector<float>;

}

void CmsSuperClusterFillerData::clear() 
{
  nBC->clear();
  nCrystals->clear();
  rawEnergy->clear();
  energy->clear();
  seedClusterEnergy->clear();
  esEnergy->clear();
  esEffsIxIx->clear();
  esEffsIyIy->clear();
  esL1Energy->clear();
  esL2Energy->clear();
  esL1Strips->clear();
  esL2Strips->clear();
  eta->clear(); 
  theta->clear();
  phi->clear();
  xPos->clear();
  yPos->clear();
  zPos->clear();
  phiWidth->clear();
  etaWidth->clear();
  e3x3->clear();
  e3x1->clear();
  e1x3->clear();
  e4x4->clear();
  e5x5->clear();
  eMax->clear();
  e2x2->clear();
  e2nd->clear();
  e1x5->clear();
  e2x5Max->clear();
  e2x5Left->clear();
  e2x5Right->clear();
  e2x5Top->clear();
  e2x5Bottom->clear();
  eLeft->clear();
  eRight->clear();
  eTop->clear();
  eBottom->clear();
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
//   etaC->clear();
//   etaS->clear();
//   etaM->clear();
//   phiC->clear();
//   phiS->clear();
//   phiM->clear();
//   xC->clear();
//   xS->clear();
//   xM->clear();
//   xZ->clear();
//   yC->clear();
//   yS->clear();
//   yM->clear();
//   yZ->clear();
  photonFix_phoE->clear();
   photonFix_phoSigma->clear();
   photonFix_eleE->clear();
   photonFix_eleSigma->clear();
   regrCorr_phoE->clear();
   regrCorr_phoSigma->clear();
   regrCorr_eleE->clear();
   regrCorr_eleSigma->clear();  
}
