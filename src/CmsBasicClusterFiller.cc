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
  delete privateData_->nCrystals;
  delete privateData_->energy;
  delete privateData_->seedEnergy;
  delete privateData_->eta;
  delete privateData_->theta;
  delete privateData_->phi;
  delete privateData_->e3x3;
  delete privateData_->e5x5;
  delete privateData_->eMax;
  delete privateData_->e2x2;
  delete privateData_->e2nd;
  delete privateData_->covIEtaIEta;
  delete privateData_->covIEtaIPhi;
  delete privateData_->covIPhiIPhi;
  delete privateData_->eTop;
  delete privateData_->eBottom;
  delete privateData_->eLeft;
  delete privateData_->eRight;
  delete privateData_->e2x5Max;
  delete privateData_->e2x5Top;
  delete privateData_->e2x5Bottom;
  delete privateData_->e2x5Left;
  delete privateData_->e2x5Right;
  delete privateData_->etaCrystal;
  delete privateData_->phiCrystal;
  delete privateData_->iEta;
  delete privateData_->iPhi;
  delete privateData_->thetaTilt;
  delete privateData_->phiTilt;
  delete privateData_->time;
  delete privateData_->chi2;
  delete privateData_->recoFlag;
  //delete privateData_->sevClosProbl;
  //delete privateData_->idClosProbl;
  //delete privateData_->fracClosProbl;
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
  for(iSC=ESCCollection->begin(); iSC!=ESCCollection->end(); iSC++)
    {
      index_SC++;
      // find SC constituents 
      for(CaloClusterPtrVector::const_iterator bcItr = iSC->clustersBegin(); bcItr != iSC->clustersEnd(); bcItr++)
        {
          if( cand->energy()== (*bcItr)->energy() 
              && cand->eta()== (*bcItr)->eta() 
              && cand->phi()== (*bcItr)->phi()  ) {
            index_SC_ref=index_SC;
          }
        }
    }

  if(  index_SC > -1 ) {

    nBCCandOfSC++;
    privateData_->indexSC->push_back(index_SC_ref);

      

    std::vector< std::pair<DetId, float> > ids = cand->hitsAndFractions();
  
    privateData_->nCrystals->push_back((int)cand->hitsAndFractions().size());
    privateData_->energy->push_back((float)cand->energy());
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

      float seedEta = cand->position().eta();

      if( fabs(seedEta) < 1.479 ) rechits = EBRecHits;
      else rechits = EERecHits; 

      float eMax = EcalClusterTools::eMax( *cand, &(*rechits) );
      float e3x3 = EcalClusterTools::e3x3( *cand, &(*rechits), topology );
      float e5x5 = EcalClusterTools::e5x5( *cand, &(*rechits), topology );
      float e2x2 = EcalClusterTools::e2x2( *cand, &(*rechits), topology );
      float e2nd = EcalClusterTools::e2nd( *cand, &(*rechits) );

      float eTop    = EcalClusterTools::eTop    (*cand, &(*rechits), topology );
      float eBottom = EcalClusterTools::eBottom (*cand, &(*rechits), topology );
      float eLeft   = EcalClusterTools::eLeft   (*cand, &(*rechits), topology );
      float eRight  = EcalClusterTools::eRight  (*cand, &(*rechits), topology );

      float e2x5Max    = EcalClusterTools::e2x5Max    (*cand, &(*rechits), topology );
      float e2x5Top    = EcalClusterTools::e2x5Top    (*cand, &(*rechits), topology );
      float e2x5Bottom = EcalClusterTools::e2x5Bottom (*cand, &(*rechits), topology );
      float e2x5Left   = EcalClusterTools::e2x5Left   (*cand, &(*rechits), topology );
      float e2x5Right  = EcalClusterTools::e2x5Right  (*cand, &(*rechits), topology );

      

      privateData_->e3x3->push_back(e3x3);
      privateData_->e5x5->push_back(e5x5);
      privateData_->eMax->push_back(eMax);
      privateData_->e2x2->push_back(e2x2);
      privateData_->e2nd->push_back(e2nd);

      privateData_->eTop->    push_back(eTop);
      privateData_->eBottom-> push_back(eBottom);
      privateData_->eLeft->   push_back(eLeft);
      privateData_->eRight->  push_back(eRight);

      privateData_->e2x5Max->    push_back(e2x5Max);
      privateData_->e2x5Top->    push_back(e2x5Top);
      privateData_->e2x5Bottom-> push_back(e2x5Bottom);
      privateData_->e2x5Left->   push_back(e2x5Left);
      privateData_->e2x5Right->  push_back(e2x5Right);

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

      // local covariances: instead of using absolute eta/phi it counts crystals normalised
      std::vector<float> vLocCov = EcalClusterTools::localCovariances( *cand, &(*rechits), topology );
      
      float covIEtaIEta = vLocCov[0];
      float covIEtaIPhi = vLocCov[1];
      float covIPhiIPhi = vLocCov[2];
      
      privateData_->covIEtaIEta->push_back(covIEtaIEta);
      privateData_->covIEtaIPhi->push_back(covIEtaIPhi);
      privateData_->covIPhiIPhi->push_back(covIPhiIPhi);

      std::pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *cand, &(*rechits) );
      DetId seedCrystalId = maxRH.first;
      EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
      
      privateData_->time->push_back((float)seedRH->time());
      privateData_->chi2->push_back((float)seedRH->chi2());
      privateData_->recoFlag->push_back((int)seedRH->recoFlag());
      privateData_->seedEnergy->push_back((float)maxRH.second);

      // channel status
      edm::ESHandle<EcalChannelStatus> pChannelStatus;
      iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
      const EcalChannelStatus *ch_status = pChannelStatus.product();

    } else {
      privateData_->e3x3->push_back(-1.);
      privateData_->e5x5->push_back(-1.);
      privateData_->eMax->push_back(-1.);
      privateData_->e2x2->push_back(-1.);
      privateData_->e2nd->push_back(-1.);
      privateData_->eTop->    push_back(-1.);
      privateData_->eBottom-> push_back(-1.);
      privateData_->eLeft->   push_back(-1.);
      privateData_->eRight->  push_back(-1.);
      privateData_->e2x5Max->    push_back(-1.);
      privateData_->e2x5Top->    push_back(-1.);
      privateData_->e2x5Bottom-> push_back(-1.);
      privateData_->e2x5Left->   push_back(-1.);
      privateData_->e2x5Right->  push_back(-1.);
      privateData_->covIEtaIEta->push_back(-1.);
      privateData_->covIEtaIPhi->push_back(-1.);
      privateData_->covIPhiIPhi->push_back(-1.);
      privateData_->time->push_back(-999.);
      privateData_->chi2->push_back(-999.);
      privateData_->recoFlag->push_back(-1);
      privateData_->seedEnergy->push_back(-1.);
      privateData_->indexSC->push_back(-1);
    }
  }

}


void CmsBasicClusterFiller::treeBCInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"nCrystals"+colSuffix).c_str(), *privateData_->nCrystals, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x2"+colSuffix).c_str(), *privateData_->e2x2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2nd"+colSuffix).c_str(), *privateData_->e2nd, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eTop"+colSuffix).c_str(), *privateData_->eTop, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eBottom"+colSuffix).c_str(), *privateData_->eBottom, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eLeft"+colSuffix).c_str(), *privateData_->eLeft, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eRight"+colSuffix).c_str(), *privateData_->eRight, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Max"+colSuffix).c_str(), *privateData_->e2x5Max, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Top"+colSuffix).c_str(), *privateData_->e2x5Top, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Bottom"+colSuffix).c_str(), *privateData_->e2x5Bottom, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Left"+colSuffix).c_str(), *privateData_->e2x5Left, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x5Right"+colSuffix).c_str(), *privateData_->e2x5Right, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"etaCrystal"+colSuffix).c_str(), *privateData_->etaCrystal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiCrystal"+colSuffix).c_str(), *privateData_->phiCrystal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"iEta"+colSuffix).c_str(), *privateData_->iEta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"iPhi"+colSuffix).c_str(), *privateData_->iPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thetaTilt"+colSuffix).c_str(), *privateData_->thetaTilt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiTilt"+colSuffix).c_str(), *privateData_->phiTilt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIEtaIEta"+colSuffix).c_str(), *privateData_->covIEtaIEta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIEtaIPhi"+colSuffix).c_str(), *privateData_->covIEtaIPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIPhiIPhi"+colSuffix).c_str(), *privateData_->covIPhiIPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlag"+colSuffix).c_str(), *privateData_->recoFlag, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"time"+colSuffix).c_str(), *privateData_->time, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2"+colSuffix).c_str(), *privateData_->chi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedEnergy"+colSuffix).c_str(), *privateData_->seedEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"indexSC"+colSuffix).c_str(), *privateData_->indexSC, nCandString.c_str(), 0, "Reco");
}
/*
// copied from: RecoEcal/EgammaCoreTools/src/EcalClusterSeverityLevelAlgo.cc
float CmsBasicClusterFiller::fractionAroundClosestProblematic( const edm::EventSetup& iSetup, const reco::CaloCluster & cluster,
                                                               const EcalRecHitCollection & recHits, const EcalChannelStatus & chStatus, const CaloTopology* topology )
{ 
  DetId closestProb = closestProblematic(iSetup, cluster , recHits, chStatus, topology).first;
  if (closestProb.null())
    return 0.;

  std::vector<DetId> neighbours = topology->getWindow(closestProb,3,3);
  std::vector<DetId>::const_iterator itn;

  std::vector< std::pair<DetId, float> > hitsAndFracs = cluster.hitsAndFractions();
  std::vector< std::pair<DetId, float> >::const_iterator it;

  float fraction = 0.;

  for ( itn = neighbours.begin(); itn != neighbours.end(); ++itn )
    { 
      for ( it = hitsAndFracs.begin(); it != hitsAndFracs.end(); ++it )
        { 
          DetId id = (*it).first;
          if ( id != (*itn) )
            continue;
          EcalRecHitCollection::const_iterator jrh = recHits.find( id );
          if ( jrh == recHits.end() )
            { 
              edm::LogError("EcalClusterSeverityLevelAlgo") << "The cluster DetId " << id.rawId() << " is not in the recHit collection!!";
              return -1;
            }

          fraction += (*jrh).energy() * (*it).second  / cluster.energy();
        }
    }
  return fraction;
}

std::pair <DetId,int> CmsBasicClusterFiller::closestProblematic(const edm::EventSetup& iSetup, const reco::CaloCluster & cluster, 
                                                                const EcalRecHitCollection & recHits, const EcalChannelStatus & chStatus , 
                                                                const CaloTopology* topology )
{
  DetId seed=EcalClusterTools::getMaximum(cluster,&recHits).first;
  if ( (seed.det() != DetId::Ecal) || 
       (EcalSubdetector(seed.subdetId()) != EcalBarrel) )
    {
      //method not supported if not in Barrel
      edm::LogError("EcalClusterSeverityLevelAlgo") << "The cluster seed is not in the BARREL";
      return std::make_pair<DetId,int>(DetId(0),-1);
    }

  int minDist=9999; DetId closestProb(0);   
  int severityClosestProb=-1;
  //Get a window of DetId around the seed crystal
  std::vector<DetId> neighbours = topology->getWindow(seed,51,11);

  for ( std::vector<DetId>::const_iterator it = neighbours.begin(); it != neighbours.end(); ++it ) 
    {
      EcalRecHitCollection::const_iterator jrh = recHits.find(*it);
      if ( jrh == recHits.end() ) 
        continue;
      //Now checking rh flag
      edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
      iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);

      uint32_t sev = sevlv->severityLevel( jrh->id(), recHits );
      if (sev == EcalSeverityLevel::kGood)
        continue;
      //Find the closest DetId in eta,phi space (distance defined by deta^2 + dphi^2)
      int deta=distanceEta(EBDetId(seed),EBDetId(*it));
      int dphi=distancePhi(EBDetId(seed),EBDetId(*it));
      if (sqrt(deta*deta + dphi*dphi) < minDist) {
        closestProb = *it;
        severityClosestProb = sev;
      }
    }
  
  closestProb_ = closestProb;
  severityClosestProb_ = severityClosestProb;
  
  return std::make_pair(closestProb,severityClosestProb);
}
*/
int CmsBasicClusterFiller::distanceEta(const EBDetId& a,const EBDetId& b)
{
  if (a.ieta() * b.ieta() > 0)
    return abs(a.ieta()-b.ieta());
  else
    return abs(a.ieta()-b.ieta())-1;
}

int CmsBasicClusterFiller::distancePhi(const EBDetId& a,const EBDetId& b)
{
  if (abs(a.iphi() -b.iphi()) > 180)
    return abs(a.iphi() - b.iphi()) - 180;
  else
    return abs(a.iphi()-b.iphi());
}


void CmsBasicClusterFillerData::initialiseCandidate() 
{
  nCrystals = new vector<int>; 
  energy = new vector<float>; 
  eta = new vector<float>; 
  theta = new vector<float>; 
  phi = new vector<float>;
  e3x3 = new vector<float>;
  e5x5 = new vector<float>;
  eMax = new vector<float>;
  e2x2 = new vector<float>;
  e2nd = new vector<float>;
  eTop    = new vector<float>;
  eBottom = new vector<float>;
  eLeft   = new vector<float>;
  eRight  = new vector<float>;
  e2x5Max    = new vector<float>;
  e2x5Top    = new vector<float>;
  e2x5Bottom = new vector<float>;
  e2x5Left   = new vector<float>;
  e2x5Right  = new vector<float>;
  etaCrystal = new vector<float>;
  phiCrystal = new vector<float>;
  iEta = new vector<int>;
  iPhi = new vector<int>;
  thetaTilt = new vector<float>;
  phiTilt= new vector<float>;
  covIEtaIEta = new vector<float>;
  covIEtaIPhi = new vector<float>;
  covIPhiIPhi = new vector<float>;
  recoFlag = new vector<int>;
  time = new vector<float>;
  chi2 = new vector<float>;
  seedEnergy = new vector<float>;
  //idClosProbl = new vector<int>;
  //sevClosProbl = new vector<int>;
  //fracClosProbl = new vector<float>;
  indexSC = new vector<int>;
  nBC =  new int;
}

void CmsBasicClusterFillerData::clear() 
{
  nCrystals->clear();
  energy->clear();
  eta->clear(); 
  theta->clear();
  phi->clear();
  e3x3->clear();
  e5x5->clear();
  eMax->clear();
  e2x2->clear();
  e2nd->clear();
  eTop->clear();
  eBottom->clear();
  eLeft->clear();
  eRight->clear();
  e2x5Max->clear();
  e2x5Top->clear();
  e2x5Bottom->clear();
  e2x5Left->clear();
  e2x5Right->clear();
  etaCrystal->clear();
  phiCrystal->clear();
  iEta->clear();
  iPhi->clear();
  thetaTilt->clear();
  phiTilt->clear();
  covIEtaIEta->clear();
  covIEtaIPhi->clear();
  covIPhiIPhi->clear();
  recoFlag->clear();
  time->clear();
  chi2->clear();
  seedEnergy->clear();
  //idClosProbl->clear();
  //sevClosProbl->clear();
  //fracClosProbl->clear();
  indexSC->clear();
}
