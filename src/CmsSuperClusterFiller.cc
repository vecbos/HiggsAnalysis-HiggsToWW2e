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

#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

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
  closestProb_ = DetId(0);
  severityClosestProb_ = -1;
  doTrackProp_ = false;

}

//--------------
// Destructor --
//--------------

CmsSuperClusterFiller::~CmsSuperClusterFiller() 
{
  // delete here the vector ptr's
  delete privateData_->nBC;
  delete privateData_->nCrystals;
  delete privateData_->iAlgo;
  delete privateData_->rawEnergy;
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
  delete privateData_->hOverE;
  delete privateData_->covIEtaIEta;
  delete privateData_->covIEtaIPhi;
  delete privateData_->covIPhiIPhi;
  delete privateData_->trackIndex;
  delete privateData_->trackDeltaR;
  delete privateData_->trackDeltaPhi;
  delete privateData_->trackDeltaEta;
  delete privateData_->gsfTrackIndex;
  delete privateData_->gsfTrackDeltaR;
  delete privateData_->gsfTrackDeltaPhi;
  delete privateData_->gsfTrackDeltaEta;
  delete privateData_->pxVtxPropagatedNegCharge;
  delete privateData_->pyVtxPropagatedNegCharge;
  delete privateData_->pzVtxPropagatedNegCharge;
  delete privateData_->pxVtxPropagatedPosCharge;
  delete privateData_->pyVtxPropagatedPosCharge;
  delete privateData_->pzVtxPropagatedPosCharge;
  delete privateData_->time;
  delete privateData_->chi2;
  delete privateData_->recoFlag;
  delete privateData_->channelStatus;
  delete privateData_->sevClosProbl;
  delete privateData_->idClosProbl;
  delete privateData_->fracClosProbl;
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
  
      // to match track-SC
      try { iEvent.getByLabel(Tracks_, tracks_); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get track collection" << Tracks_; }

      try { iEvent.getByLabel(GsfTracks_, gsfTracks_); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get GSF track collection" << GsfTracks_; }

      try { iEvent.getByType(theBeamSpot_); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get beam spot "; }

      try { iEvent.getByLabel("offlinePrimaryVertices", hVtx_); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsSuperClusterFiller") << "Can't get primary vertex collection: offlinePrimaryVertices"; }

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
        if ( doTrackProp_ ) {
          // fill CTF track - SC match
          writeTrackInfo(&(*cand),iEvent,iSetup,&(*tracks_),track);
          // fill GSF track - SC match
          writeTrackInfo(&(*cand),iEvent,iSetup,&(*gsfTracks_),gsftrack);
        }
        // fill trajectory propagation at vertex
        writeSCVtxPropagationInfo(&(*cand),iEvent,iSetup);
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
  if ( doTrackProp_ ) treeTrackInfo(columnPrefix,columnSuffix);
  treeSCVtxPropagationInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

}






void CmsSuperClusterFiller::writeSCInfo(const SuperCluster *cand, 
                                        const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                        const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits) {

  // fill the SC infos
  privateData_->nBC->push_back((int)cand->clustersSize());
  privateData_->nCrystals->push_back((int)cand->hitsAndFractions().size());
  privateData_->iAlgo->push_back((int)cand->seed()->algo());
  privateData_->rawEnergy->push_back((float)cand->rawEnergy());
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

      privateData_->e3x3->push_back(e3x3);
      privateData_->e5x5->push_back(e5x5);
      privateData_->eMax->push_back(eMax);
      privateData_->e2x2->push_back(e2x2);
      privateData_->e2nd->push_back(e2nd);

      // local covariances: instead of using absolute eta/phi it counts crystals normalised
      std::vector<float> vLocCov = EcalClusterTools::localCovariances( *theSeed, &(*rechits), topology );
      
      float covIEtaIEta = vLocCov[0];
      float covIEtaIPhi = vLocCov[1];
      float covIPhiIPhi = vLocCov[2];
      
      privateData_->covIEtaIEta->push_back(covIEtaIEta);
      privateData_->covIEtaIPhi->push_back(covIEtaIPhi);
      privateData_->covIPhiIPhi->push_back(covIPhiIPhi);

      std::pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *theSeed, &(*rechits) );
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

      EcalChannelStatusMap::const_iterator chit = pChannelStatus->find( seedCrystalId );
      EcalChannelStatusCode chStatusCode = 1;
      if ( chit != pChannelStatus->end() ) {
	chStatusCode = *chit;
      } else {
	edm::LogError("EcalRecHitProducerError") << "No channel status found for xtal "
						 << "! something wrong with EcalChannelStatus in your DB? ";
      }
      int cStatusFlag = (int)(chStatusCode.getStatusCode() & 0x001F);
      privateData_->channelStatus->push_back(cStatusFlag);

      if( fabs(seedEta) < 1.479 ) {
        float frac = fractionAroundClosestProblematic( *cand, *rechits, *ch_status, topology);
        privateData_->fracClosProbl->push_back(frac);
        if( closestProb_.null() ) {
          privateData_->idClosProbl->push_back(-1);
          privateData_->sevClosProbl->push_back(-1);
        } else {
          privateData_->idClosProbl->push_back(closestProb_.rawId());
          privateData_->sevClosProbl->push_back(severityClosestProb_);
        }
      } else {
        privateData_->fracClosProbl->push_back(-1);
        privateData_->idClosProbl->push_back(-1);
        privateData_->idClosProbl->push_back(-1);
        privateData_->sevClosProbl->push_back(-1);
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
    privateData_->time->push_back(-999.);
    privateData_->chi2->push_back(-999.);
    privateData_->recoFlag->push_back(-1);
    privateData_->channelStatus->push_back(-1);
    privateData_->seedEnergy->push_back(-1.);
    privateData_->fracClosProbl->push_back(-1);
    privateData_->idClosProbl->push_back(-1);
    privateData_->idClosProbl->push_back(-1);
    privateData_->sevClosProbl->push_back(-1);
  }

  // calculate H/E
  float hOverEConeSize = 0.15;
  float hOverEPtMin = 0.;
  EgammaTowerIsolation *towerIso1 = new EgammaTowerIsolation(hOverEConeSize,0.,hOverEPtMin,1,calotowers_.product()) ;
  EgammaTowerIsolation *towerIso2 = new EgammaTowerIsolation(hOverEConeSize,0.,hOverEPtMin,2,calotowers_.product()) ;
  
  float TowerHcalESum1 = towerIso1->getTowerESum(cand);
  float TowerHcalESum2 = towerIso2->getTowerESum(cand);
  float hcalESum = TowerHcalESum1 + TowerHcalESum2;
  
  privateData_->hOverE->push_back(hcalESum);

  delete towerIso1;
  delete towerIso2;

}

void CmsSuperClusterFiller::writeTrackInfo(const reco::SuperCluster *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                           const reco::TrackCollection *theTracks, int trackType) {

  math::XYZPoint xyzVertexPos;
  GlobalPoint gpVertexPos;
  GlobalPoint origin;

  if ( hVtx_->size()>0 ){ 
    float theMaxPt  = -999.;
    VertexCollection::const_iterator thisVertex;
    for(thisVertex = hVtx_->begin(); thisVertex != hVtx_->end(); ++thisVertex){      
      float SumPt = 0.0;
      if((*thisVertex).tracksSize() > 0){
        std::vector<TrackBaseRef >::const_iterator thisTrack;
        for( thisTrack=(*thisVertex).tracks_begin(); thisTrack!=(*thisVertex).tracks_end(); thisTrack++){
          // if((**thisTrack).charge()==-1 || (**thisTrack).charge()==1) SumPt += (**thisTrack).pt();
          SumPt += (**thisTrack).pt();
        }}
      if (SumPt>theMaxPt){ 
        theMaxPt = SumPt; 
        gpVertexPos  = GlobalPoint((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
        xyzVertexPos = math::XYZVector((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
      }}
  }
  else{
    gpVertexPos  = GlobalPoint(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());
    xyzVertexPos = math::XYZVector(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());
  }  
  origin = GlobalPoint(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());

  // magnetic field
  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  
  float bestDeltaR = 999.;
  float bestDeltaPhi = 999.;
  float bestDeltaEta = 999.;
  int bestTrack = -1;

  int trackIndex=0;
  TrackCollection::const_iterator trIter;      
  for (trIter=theTracks->begin(); trIter!=theTracks->end(); trIter++) {
    
    int trackQ          = trIter->charge();
    float trackPt       = trIter->p()*sin(trIter->theta());
    
    if (trackPt > 1) {

      // preso da HLTrigger/Egamma/src/HLTElectronDetaDphiFilter.cc
      const math::XYZVector trackMom = trIter->momentum();
      math::XYZPoint SCcorrPosition(cand->x()-gpVertexPos.x(), cand->y()-gpVertexPos.y(), cand->z()-gpVertexPos.z());
      float etaScCorr = SCcorrPosition.eta();                                                            // eta sc va corretto per il beam spot / vertex
      float deltaEta  = fabs(etaScCorr - trIter->eta());                                           // eta traccia al vtx non va corretto x beam spot (gia' incluso nel fit)
      // eta traccia non va propagato al calorimetro tanto non curva 
      ECALPositionCalculator posCalc;
      float phiTrCorr = posCalc.ecalPhi(&(*theMagField), trackMom, xyzVertexPos, trackQ);          // phi traccia al vtx non va corretto x beam spot (gia' incluso nel fit)
      // ma phi traccia va propagato al calo e qui serve il constraint del vtx 
      float deltaPhi  = fabs(cand->phi() - phiTrCorr);
      if(deltaPhi>6.283185308) deltaPhi -= 6.283185308;
      if(deltaPhi>3.141592654) deltaPhi = 6.283185308-deltaPhi;
      float deltaR = sqrt (deltaEta*deltaEta + deltaPhi*deltaPhi);

      if (deltaR < bestDeltaR){
        bestDeltaPhi    = deltaPhi;
        bestDeltaEta    = deltaEta;
        bestDeltaR      = deltaR;
        bestTrack = trackIndex;
      }
    }
    trackIndex++;
  }

  if ( trackType == track ) {
    privateData_->trackIndex->push_back(bestTrack);
    privateData_->trackDeltaR->push_back(bestDeltaR);
    privateData_->trackDeltaPhi->push_back(bestDeltaPhi);
    privateData_->trackDeltaEta->push_back(bestDeltaEta);
  } else if ( trackType == gsftrack ) {
    privateData_->gsfTrackIndex->push_back(bestTrack);
    privateData_->gsfTrackDeltaR->push_back(bestDeltaR);
    privateData_->gsfTrackDeltaPhi->push_back(bestDeltaPhi);
    privateData_->gsfTrackDeltaEta->push_back(bestDeltaEta);
  }

}

void CmsSuperClusterFiller::writeTrackInfo(const reco::SuperCluster *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                           const reco::GsfTrackCollection *theTracks, int trackType) {

  math::XYZPoint xyzVertexPos;
  GlobalPoint gpVertexPos;
  GlobalPoint origin;

  if ( hVtx_->size()>0 ){ 
    float theMaxPt  = -999.;
    VertexCollection::const_iterator thisVertex;
    for(thisVertex = hVtx_->begin(); thisVertex != hVtx_->end(); ++thisVertex){      
      float SumPt = 0.0;
      if((*thisVertex).tracksSize() > 0){
        std::vector<TrackBaseRef >::const_iterator thisTrack;
        for( thisTrack=(*thisVertex).tracks_begin(); thisTrack!=(*thisVertex).tracks_end(); thisTrack++){
          // if((**thisTrack).charge()==-1 || (**thisTrack).charge()==1) SumPt += (**thisTrack).pt();
          SumPt += (**thisTrack).pt();
        }}
      if (SumPt>theMaxPt){ 
        theMaxPt = SumPt; 
        gpVertexPos  = GlobalPoint((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
        xyzVertexPos = math::XYZVector((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
      }}
  }
  else{
    gpVertexPos  = GlobalPoint(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());
    xyzVertexPos = math::XYZVector(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());
  }  
  origin = GlobalPoint(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());

  // magnetic field
  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  
  float bestDeltaR = 999.;
  float bestDeltaPhi = 999.;
  float bestDeltaEta = 999.;
  int bestTrack = -1;

  int trackIndex=0;
  GsfTrackCollection::const_iterator trIter;      
  for (trIter=theTracks->begin(); trIter!=theTracks->end(); trIter++) {
    
    int trackQ          = trIter->charge();
    float trackPt       = trIter->p()*sin(trIter->theta());
    
    if (trackPt > 1) {

      // preso da HLTrigger/Egamma/src/HLTElectronDetaDphiFilter.cc
      const math::XYZVector trackMom = trIter->momentum();
      math::XYZPoint SCcorrPosition(cand->x()-gpVertexPos.x(), cand->y()-gpVertexPos.y(), cand->z()-gpVertexPos.z());
      float etaScCorr = SCcorrPosition.eta();                                                            // eta sc va corretto per il beam spot / vertex
      float deltaEta  = fabs(etaScCorr - trIter->eta());                                           // eta traccia al vtx non va corretto x beam spot (gia' incluso nel fit)
      // eta traccia non va propagato al calorimetro tanto non curva 
      ECALPositionCalculator posCalc;
      float phiTrCorr = posCalc.ecalPhi(&(*theMagField), trackMom, xyzVertexPos, trackQ);          // phi traccia al vtx non va corretto x beam spot (gia' incluso nel fit)
      // ma phi traccia va propagato al calo e qui serve il constraint del vtx 
      float deltaPhi  = fabs(cand->phi() - phiTrCorr);
      if(deltaPhi>6.283185308) deltaPhi -= 6.283185308;
      if(deltaPhi>3.141592654) deltaPhi = 6.283185308-deltaPhi;
      float deltaR = sqrt (deltaEta*deltaEta + deltaPhi*deltaPhi);

      if (deltaR < bestDeltaR){
        bestDeltaPhi    = deltaPhi;
        bestDeltaEta    = deltaEta;
        bestDeltaR      = deltaR;
        bestTrack = trackIndex;
      }
    }
    trackIndex++;
  }

  if ( trackType == track ) {
    privateData_->trackIndex->push_back(bestTrack);
    privateData_->trackDeltaR->push_back(bestDeltaR);
    privateData_->trackDeltaPhi->push_back(bestDeltaPhi);
    privateData_->trackDeltaEta->push_back(bestDeltaEta);
  } else if ( trackType == gsftrack ) {
    privateData_->gsfTrackIndex->push_back(bestTrack);
    privateData_->gsfTrackDeltaR->push_back(bestDeltaR);
    privateData_->gsfTrackDeltaPhi->push_back(bestDeltaPhi);
    privateData_->gsfTrackDeltaEta->push_back(bestDeltaEta);
  }

}

void CmsSuperClusterFiller::writeSCVtxPropagationInfo(const reco::SuperCluster *cand,
                                                      const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  GlobalPoint vertexPos;

  // choose the Highest pT primary vertex
  if ( hVtx_->size()>0 ){ 
    float theMaxPt  = -999.;
    VertexCollection::const_iterator thisVertex;
    for(thisVertex = hVtx_->begin(); thisVertex != hVtx_->end(); ++thisVertex){      
      float SumPt = 0.0;
      if((*thisVertex).tracksSize() > 0){
        std::vector<TrackBaseRef >::const_iterator thisTrack;
        for( thisTrack=(*thisVertex).tracks_begin(); thisTrack!=(*thisVertex).tracks_end(); thisTrack++){
          SumPt += (**thisTrack).pt();
        }}
      if (SumPt>theMaxPt){ 
        theMaxPt = SumPt; 
        vertexPos  = GlobalPoint((*thisVertex).x(), (*thisVertex).y(), (*thisVertex).z()); 
      }
    }
  } else {
    vertexPos  = GlobalPoint(theBeamSpot_->position().x(),theBeamSpot_->position().y(),theBeamSpot_->position().z());
  }

  // magnetic field
  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);

  // propagate the SC to the primary vertex
  const GlobalPoint clusterPos(cand->position().x(), cand->position().y(), cand->position().z());    

  float clusterEnergy = cand->energy();
  // electron hypothesis
  FreeTrajectoryState ftsElectron = myFTS(&(*theMagField), clusterPos, vertexPos, clusterEnergy, 1);    
  // positron hypothesis
  FreeTrajectoryState ftsPositron = myFTS(&(*theMagField), clusterPos, vertexPos, clusterEnergy, -1);    
  
  privateData_->pxVtxPropagatedNegCharge->push_back(ftsElectron.momentum().x());
  privateData_->pyVtxPropagatedNegCharge->push_back(ftsElectron.momentum().y());
  privateData_->pzVtxPropagatedNegCharge->push_back(ftsElectron.momentum().z());
  privateData_->pxVtxPropagatedPosCharge->push_back(ftsPositron.momentum().x());
  privateData_->pyVtxPropagatedPosCharge->push_back(ftsPositron.momentum().y());
  privateData_->pzVtxPropagatedPosCharge->push_back(ftsPositron.momentum().z());

}

void CmsSuperClusterFiller::treeSCInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"nBC"+colSuffix).c_str(), *privateData_->nBC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nCrystals"+colSuffix).c_str(), *privateData_->nCrystals, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"iAlgo"+colSuffix).c_str(), *privateData_->iAlgo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"rawEnergy"+colSuffix).c_str(), *privateData_->rawEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2x2"+colSuffix).c_str(), *privateData_->e2x2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e2nd"+colSuffix).c_str(), *privateData_->e2nd, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIEtaIEta"+colSuffix).c_str(), *privateData_->covIEtaIEta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIEtaIPhi"+colSuffix).c_str(), *privateData_->covIEtaIPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"covIPhiIPhi"+colSuffix).c_str(), *privateData_->covIPhiIPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverE"+colSuffix).c_str(), *privateData_->hOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlag"+colSuffix).c_str(), *privateData_->recoFlag, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"channelStatus"+colSuffix).c_str(), *privateData_->channelStatus, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"time"+colSuffix).c_str(), *privateData_->time, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2"+colSuffix).c_str(), *privateData_->chi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"seedEnergy"+colSuffix).c_str(), *privateData_->seedEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"idClosProbl"+colSuffix).c_str(), *privateData_->idClosProbl, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"sevClosProbl"+colSuffix).c_str(), *privateData_->sevClosProbl, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"fracClosProbl"+colSuffix).c_str(), *privateData_->fracClosProbl, nCandString.c_str(), 0, "Reco");
}


void CmsSuperClusterFiller::treeTrackInfo(const std::string colPrefix, const std::string colSuffix)
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(), *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDeltaR"+colSuffix).c_str(), *privateData_->trackDeltaR, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDeltaPhi"+colSuffix).c_str(), *privateData_->trackDeltaPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDeltaEta"+colSuffix).c_str(), *privateData_->trackDeltaEta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(), *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackDeltaR"+colSuffix).c_str(), *privateData_->gsfTrackDeltaR, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackDeltaPhi"+colSuffix).c_str(), *privateData_->gsfTrackDeltaPhi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackDeltaEta"+colSuffix).c_str(), *privateData_->gsfTrackDeltaEta, nCandString.c_str(), 0, "Reco");
}

void CmsSuperClusterFiller::treeSCVtxPropagationInfo(const std::string colPrefix, const std::string colSuffix) 
{
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"pxVtxPropagatedNegCharge"+colSuffix).c_str(), *privateData_->pxVtxPropagatedNegCharge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pyVtxPropagatedNegCharge"+colSuffix).c_str(), *privateData_->pyVtxPropagatedNegCharge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pzVtxPropagatedNegCharge"+colSuffix).c_str(), *privateData_->pzVtxPropagatedNegCharge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pxVtxPropagatedPosCharge"+colSuffix).c_str(), *privateData_->pxVtxPropagatedPosCharge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pyVtxPropagatedPosCharge"+colSuffix).c_str(), *privateData_->pyVtxPropagatedPosCharge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pzVtxPropagatedPosCharge"+colSuffix).c_str(), *privateData_->pzVtxPropagatedPosCharge, nCandString.c_str(), 0, "Reco");
}

// copied from: RecoEcal/EgammaCoreTools/src/EcalClusterSeverityLevelAlgo.cc
float CmsSuperClusterFiller::fractionAroundClosestProblematic( const reco::CaloCluster & cluster,
                                                               const EcalRecHitCollection & recHits, const EcalChannelStatus & chStatus, const CaloTopology* topology )
{ 
  DetId closestProb = closestProblematic(cluster , recHits, chStatus, topology).first;
  //  std::cout << "%%%%%%%%%%% Closest prob is " << EBDetId(closestProb) << std::endl;
  if (closestProb.null())
    return 0.;

  std::vector<DetId> neighbours = topology->getWindow(closestProb,3,3);
  std::vector<DetId>::const_iterator itn;

  std::vector< std::pair<DetId, float> > hitsAndFracs = cluster.hitsAndFractions();
  std::vector< std::pair<DetId, float> >::const_iterator it;

  float fraction = 0.;

  for ( itn = neighbours.begin(); itn != neighbours.end(); ++itn )
    { 
      //      std::cout << "Checking detId " << EBDetId((*itn)) << std::endl;
      for ( it = hitsAndFracs.begin(); it != hitsAndFracs.end(); ++it )
        { 
          DetId id = (*it).first;
          if ( id != (*itn) )
            continue;
          //      std::cout << "Is in cluster detId " << EBDetId(id) << std::endl;
          EcalRecHitCollection::const_iterator jrh = recHits.find( id );
          if ( jrh == recHits.end() )
            { 
              edm::LogError("EcalClusterSeverityLevelAlgo") << "The cluster DetId " << id.rawId() << " is not in the recHit collection!!";
              return -1;
            }

          fraction += (*jrh).energy() * (*it).second  / cluster.energy();
        }
    }
  //  std::cout << "%%%%%%%%%%% Fraction is " << fraction << std::endl;
  return fraction;
}

std::pair <DetId,int> CmsSuperClusterFiller::closestProblematic(const reco::CaloCluster & cluster, 
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
      uint32_t sev = EcalSeverityLevelAlgo::severityLevel( jrh->id(), recHits, chStatus );
      if (sev == EcalSeverityLevelAlgo::kGood)
        continue;
      //      std::cout << "[closestProblematic] Found a problematic channel " << EBDetId(*it) << " " << flag << std::endl;
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
  
  return std::make_pair<DetId,int>(closestProb,severityClosestProb);
}

int CmsSuperClusterFiller::distanceEta(const EBDetId& a,const EBDetId& b)
{
  if (a.ieta() * b.ieta() > 0)
    return abs(a.ieta()-b.ieta());
  else
    return abs(a.ieta()-b.ieta())-1;
}

int CmsSuperClusterFiller::distancePhi(const EBDetId& a,const EBDetId& b)
{
  if (abs(a.iphi() -b.iphi()) > 180)
    return abs(a.iphi() - b.iphi()) - 180;
  else
    return abs(a.iphi()-b.iphi());
}



void CmsSuperClusterFillerData::initialiseCandidate() 
{
  nBC = new vector<int>;
  nCrystals = new vector<int>; 
  iAlgo = new vector<int>;
  rawEnergy = new vector<float>; 
  energy = new vector<float>; 
  eta = new vector<float>; 
  theta = new vector<float>; 
  phi = new vector<float>;
  e3x3 = new vector<float>;
  e5x5 = new vector<float>;
  eMax = new vector<float>;
  e2x2 = new vector<float>;
  e2nd = new vector<float>;
  hOverE = new vector<float>;
  covIEtaIEta = new vector<float>;
  covIEtaIPhi = new vector<float>;
  covIPhiIPhi = new vector<float>;
  trackIndex = new vector<int>;
  trackDeltaR = new vector<float>;
  trackDeltaPhi = new vector<float>;
  trackDeltaEta = new vector<float>;
  gsfTrackIndex = new vector<int>;
  gsfTrackDeltaR = new vector<float>;
  gsfTrackDeltaPhi = new vector<float>;
  gsfTrackDeltaEta = new vector<float>;
  pxVtxPropagatedNegCharge = new vector<float>;
  pyVtxPropagatedNegCharge = new vector<float>;
  pzVtxPropagatedNegCharge = new vector<float>;
  pxVtxPropagatedPosCharge = new vector<float>;
  pyVtxPropagatedPosCharge = new vector<float>;
  pzVtxPropagatedPosCharge = new vector<float>;
  recoFlag = new vector<int>;
  channelStatus = new vector<int>;
  time = new vector<float>;
  chi2 = new vector<float>;
  seedEnergy = new vector<float>;
  idClosProbl = new vector<int>;
  sevClosProbl = new vector<int>;
  fracClosProbl = new vector<float>;
  nSC =  new int;
}

void CmsSuperClusterFillerData::clear() 
{
  nBC->clear();
  nCrystals->clear();
  iAlgo->clear();
  rawEnergy->clear();
  energy->clear();
  eta->clear(); 
  theta->clear();
  phi->clear();
  e3x3->clear();
  e5x5->clear();
  eMax->clear();
  e2x2->clear();
  e2nd->clear();
  hOverE->clear();
  covIEtaIEta->clear();
  covIEtaIPhi->clear();
  covIPhiIPhi->clear();
  trackIndex->clear();
  trackDeltaR->clear();
  trackDeltaPhi->clear();
  trackDeltaEta->clear();
  gsfTrackIndex->clear();
  gsfTrackDeltaR->clear();
  gsfTrackDeltaPhi->clear();
  gsfTrackDeltaEta->clear();
  pxVtxPropagatedNegCharge->clear();
  pyVtxPropagatedNegCharge->clear();
  pzVtxPropagatedNegCharge->clear();
  pxVtxPropagatedPosCharge->clear();
  pyVtxPropagatedPosCharge->clear();
  pzVtxPropagatedPosCharge->clear();
  recoFlag->clear();
  channelStatus->clear();
  time->clear();
  chi2->clear();
  seedEnergy->clear();
  idClosProbl->clear();
  sevClosProbl->clear();
  fracClosProbl->clear();
}
