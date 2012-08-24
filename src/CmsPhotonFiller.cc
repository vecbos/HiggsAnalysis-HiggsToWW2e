//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsTreeFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPhotonFiller.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"                                                                                      
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"                                                                             
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"              


#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <string>

using namespace edm;
using namespace reco;
using namespace reco::isodeposit;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPhotonFiller::CmsPhotonFiller(CmsTree *cmsTree, 
                                 int maxTracks, int maxMCTracks,
                                 bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPhotonFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

//--------------
// Destructor --
//--------------

CmsPhotonFiller::~CmsPhotonFiller() {

  // delete here the vector ptr's
  delete privateData_->fiducialFlags;
  delete privateData_->recoFlags;
  delete privateData_->superClusterIndex;
  delete privateData_->PFsuperClusterIndex;
  
  delete privateData_->hOverE;
  delete privateData_->hTowOverE;
  delete privateData_->dr03HollowTkSumPt;
  delete privateData_->dr03TkSumPt;
  delete privateData_->dr03EcalRecHitSumEt;
  delete privateData_->dr03HcalTowerSumEt;
  delete privateData_->dr04HollowTkSumPt;
  delete privateData_->dr04TkSumPt;
  delete privateData_->dr04EcalRecHitSumEt;
  delete privateData_->dr04HcalTowerSumEt;
  delete privateData_->chargedHadronIso;
  delete privateData_->neutralHadronIso;
  delete privateData_->photonIso;
  delete privateData_->hasPixelSeed;
  delete privateData_->hasMatchedConversion;

  delete privateData_->dr01chPFIso;
  delete privateData_->dr01nhPFIso;
  delete privateData_->dr01phPFIso;  

  delete privateData_->dr02chPFIso;
  delete privateData_->dr02nhPFIso;
  delete privateData_->dr02phPFIso;  

  delete privateData_->dr03chPFIso;
  delete privateData_->dr03nhPFIso;
  delete privateData_->dr03phPFIso;  

  delete privateData_->dr04chPFIso;
  delete privateData_->dr04nhPFIso;
  delete privateData_->dr04phPFIso;  

  delete privateData_->dr05chPFIso;
  delete privateData_->dr05nhPFIso;
  delete privateData_->dr05phPFIso;  

  delete privateData_->dr06chPFIso;
  delete privateData_->dr06nhPFIso;
  delete privateData_->dr06phPFIso;  

  delete privateData_->ncand;
  delete privateData_;
}


//-------------
// Methods   --
//-------------

void CmsPhotonFiller::writeCollectionToTree(edm::InputTag collectionTag,
					      const edm::Event& iEvent, const edm::EventSetup& iSetup,
					      const std::string &columnPrefix, const std::string &columnSuffix,
					      bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPhotonFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ << " and no output flag is set."
				   << " No tracks written to tuple for this event ";
      return;
    }
  
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPhotonFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ 
				   << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    // superclusters
    try { iEvent.getByLabel(EcalSuperClusters_, hSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get merged ECAL supercluster collection: " << EcalSuperClusters_; }

    // for conversions with full vertex fit
    iEvent.getByLabel("offlineBeamSpot", bsHandle);

    try { iEvent.getByLabel(conversionsProducer_, hConversions); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Can't get conversions collection " << conversionsProducer_; }

    try{
      iEvent.getByLabel(pfCandidates_,pfCands);
      iEvent.getByLabel(primaryVertices_,PVs);
    }catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Couldn't get PFCandidates/PVs:  " << pfCandidates_ << "/" << primaryVertices_;}

    for(int index = 0; index < (int)collection->size(); index++) {

      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);

      const PhotonRef photonRef = collection->refAt(index).castTo<PhotonRef>();

      if ( !(photonRef.isNull()) ) {

        // fill Cluster Adapter
        SuperClusterRef sclusRef = photonRef->superCluster();
        //        SuperClusterRef pfclusRef = photonRef->pflowSuperCluster();
        SuperClusterRef pfclusRef = photonRef->superCluster(); // placeholder
        writeEcalInfo(photonRef,iEvent,iSetup,sclusRef,pfclusRef);

      } else {
        edm::LogWarning("CmsPhotonFiller") << "Warning! The collection seems to be not made by "
                                             << "electrons, electron-specific infos will be set to default.";
      }

    }

  }
  else {
    *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  int blockSizePerPV = (PVs.isValid()) ? blockSize*PVs->size() : 0;
  std::string nCandPVString = columnPrefix+(*trkIndexName_)+"PerPV"+columnSuffix; 
  cmstree->column(nCandPVString.c_str(),blockSizePerPV,0,"Reco");
  
  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  treeEcalInfo(columnPrefix,columnSuffix);
  
  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}




void CmsPhotonFiller::writeEcalInfo(const PhotonRef photonRef, 
                                    const edm::Event& iEvent, const edm::EventSetup& iSetup, 
                                    SuperClusterRef sclusRef, SuperClusterRef pfclusRef){
				

  if(photonRef.isNonnull()) {

    // fiducial flags in ECAL
    int packed_sel = -1;
    int isEB = ( photonRef->isEB() ) ? 1 : 0;
    int isEE = ( photonRef->isEE() ) ? 1 : 0;
    int isGap = 0; // photon has not this method, but mantain the synch with ele
    int isEBEEGap = ( photonRef->isEBEEGap() ) ? 1 : 0;
    int isEBGap = ( photonRef->isEBGap() ) ? 1 : 0;
    int isEBEtaGap = ( photonRef->isEBEtaGap() ) ? 1 : 0;
    int isEBPhiGap = ( photonRef->isEBPhiGap() ) ? 1 : 0;
    int isEEGap = ( photonRef->isEEGap() ) ? 1 : 0;
    int isEEDeeGap = ( photonRef->isEEDeeGap() ) ? 1 : 0;
    int isEERingGap = ( photonRef->isEERingGap() ) ? 1 : 0;
    
    packed_sel = ( isEB << 9 ) | ( isEE << 8 ) | ( isGap << 7 ) |
      ( isEBEEGap << 6 ) | ( isEBGap << 5 ) | ( isEBEtaGap << 4 ) | ( isEBPhiGap << 3 ) |
      ( isEEGap << 2 ) | ( isEEDeeGap << 1 ) | isEERingGap;

    privateData_->fiducialFlags->push_back(packed_sel);

    //    int packed_reco; // the following available only in > 420
    //    int isStdPho = ( photonRef->isStandardPhoton() ) ? 1 : 0;
    //    int isPFPho = ( photonRef->isPFlowPhoton() ) ? 1 : 0;
    //    packed_reco = ( isStdPho << 1 ) | isPFPho;
    //    privateData_->recoFlags->push_back( packed_reco );
    privateData_->recoFlags->push_back( -1 );

    // link to the supercluster. We have to match manually because in the merged collection
    // the reference is broken.
    if ( sclusRef.isNonnull() ) {
      int scInd=-1;
      for(unsigned int isclu=0; isclu!=hSuperClusters->size(); ++isclu) {
        const reco::SuperClusterRef thisScRef(hSuperClusters,isclu);
        if(sclusRef->energy()==thisScRef->energy() && sclusRef->position()==thisScRef->position()) {
          scInd = thisScRef.key();
          break;
        }
      }
      privateData_->superClusterIndex->push_back( scInd );
    } else {
      privateData_->superClusterIndex->push_back( -1 );
    }

    //    if ( isPFPho && pfclusRef.isNonnull() ) {
//     if ( pfclusRef.isNonnull() ) { 
//       privateData_->PFsuperClusterIndex->push_back( pfclusRef.key() );
//     } else {
//       privateData_->PFsuperClusterIndex->push_back( -1 );
//     }
    privateData_->PFsuperClusterIndex->push_back( -1 );
    
    // isolations
    privateData_->hOverE->push_back(photonRef->hadronicOverEm());
    // privateData_->hTowOverE->push_back(photonRef->hadTowOverEm()); // available in CMSSW > 5.2.X
    privateData_->hTowOverE->push_back(-999.);
    privateData_->dr03TkSumPt->push_back(photonRef->trkSumPtSolidConeDR03());
    privateData_->dr03HollowTkSumPt->push_back(photonRef->trkSumPtHollowConeDR03());
    privateData_->dr03EcalRecHitSumEt->push_back(photonRef->ecalRecHitSumEtConeDR03());
    privateData_->dr03HcalTowerSumEt->push_back(photonRef->hcalTowerSumEtConeDR03());
    privateData_->dr04TkSumPt->push_back(photonRef->trkSumPtSolidConeDR04());
    privateData_->dr04HollowTkSumPt->push_back(photonRef->trkSumPtHollowConeDR04());
    privateData_->dr04EcalRecHitSumEt->push_back(photonRef->ecalRecHitSumEtConeDR04());
    privateData_->dr04HcalTowerSumEt->push_back(photonRef->hcalTowerSumEtConeDR04());
    privateData_->chargedHadronIso->push_back(photonRef->chargedHadronIso());
    privateData_->neutralHadronIso->push_back(photonRef->neutralHadronIso());
    privateData_->photonIso->push_back(photonRef->photonIso());
    privateData_->hasPixelSeed->push_back(int(photonRef->hasPixelSeed()));

    // conversion rejection
    bool matchesConv = ConversionTools::hasMatchedConversion(*sclusRef,hConversions,bsHandle->position());
    privateData_->hasMatchedConversion->push_back(matchesConv);

    //do the PF Isolation w.r.t. EVERY vertex

    PFCandidateCollection::const_iterator cand;
    VertexCollection::const_iterator vtx;
    AbsVetos emptyVetos;
    
    // neutral hadron and photon pf isolation are always computed with respect to 0,0,0
    // charged hadron isolation is computed w.r.t. each vertex
      
    TVector3 phoPos(photonRef->superCluster()->x(),photonRef->superCluster()->y(),photonRef->superCluster()->z());
    //vetoes
    //photon
    const float dRVetoBarrel_p = 0.;
    const float dRVetoEndcap_p = 0.07;
    const float etaStripBarrel_p = 0.015;
    const float etaStripEndcap_p = 0.;
    //charged
    const float dRVetoBarrel_c = 0.02;
    const float dRVetoEndcap_c = 0.02;
    const float dzMax_c = 0.2;
    const float dxyMax_c  = 0.1;

    float dRVeto_p = dRVetoEndcap_p;
    float etaStrip_p = etaStripEndcap_p;
    float dRVeto_c = dRVetoEndcap_c;
    if(fabs(phoPos.Eta()) < 1.48){ //barrel
      dRVeto_p = dRVetoBarrel_p;
      etaStrip_p = etaStripBarrel_p;
      dRVeto_c = dRVetoBarrel_c;
    }
    const int nCones=6;
    float neutralIsos[nCones], photonIsos[nCones]; 
    for(int i=0;i<nCones;i++){
      neutralIsos[i]=0; photonIsos[i]=0;
    }

    for(vtx = PVs->begin(); vtx != PVs->end(); vtx++){ 
      //calculate the vector w.r.t. this vtx
      TVector3 vtxPos(vtx->x(),vtx->y(),vtx->z());
      TVector3 phoDirFromVtx = (phoPos-vtxPos).Unit();

      float chargedIsos[nCones];
      for(int i=0;i<nCones;i++) chargedIsos[i] = 0;

      for(cand = pfCands->begin(); cand != pfCands->end(); cand++){
	if(cand->pt() < 0.0001) continue;
	if( fabs(cand->energy() - photonRef->energy()) < 1e-6 
	    && fabs(cand->eta() - photonRef->eta())    < 1e-6
	    && fabs(cand->phi() - photonRef->phi())    < 1e-6) continue; // this candidate is the photon!
	TVector3 candPos;
	candPos.SetPtEtaPhi(cand->pt(), cand->eta(), cand->phi());

	if( cand->particleId() == reco::PFCandidate::h) {   // charged
	  float dz = fabs(cand->trackRef()->dz(vtx->position()));
	  if (dz > dzMax_c) continue;  //vertex of this cand not compatible with this reco vtx
	  float dxy = fabs(cand->trackRef()->dxy(vtx->position()));
	  if(fabs(dxy) > dxyMax_c) continue;

	  float dR = phoDirFromVtx.DeltaR(candPos);
	  if(dR < dRVeto_c) continue;
	  for(int i=0;i<nCones;i++){
	    if(dR > 0.1*(i+1)) continue; // fill cones in increasing 0.1 radius
	    chargedIsos[i] += cand->pt();
	  }
	} // end charged 
      
	if(vtx!=PVs->begin()) continue; // only fill these for the first PV

	TVector3 candVtx(cand->vertex().x(),cand->vertex().y(),cand->vertex().z());
	TVector3 phoDirFromCandVtx = phoPos-candVtx;

	//neutral VETO
	if( cand->particleId() == reco::PFCandidate::h0){
	  float dR = phoDirFromCandVtx.DeltaR(candPos);
	  for(int i=0;i<nCones;i++){
	    if(dR > 0.1*(i+1))continue;
	    neutralIsos[i]+=cand->pt();
	  }	  
	}

	//PHOTON ISO
	if( cand->particleId() == reco::PFCandidate::gamma){
	  if(cand->superClusterRef().isNonnull() && photonRef->superCluster().isNonnull()){
	    if(cand->superClusterRef()== photonRef->superCluster()) continue; // veto the candidate matched if its a pfPhoton matched to the reco::Photon
	  }

	  float dEta = fabs(phoDirFromCandVtx.Eta() - cand->momentum().Eta());
	  if (dEta < etaStrip_p) continue;
	  float dR = phoDirFromCandVtx.DeltaR(candPos);
	  if(dR < dRVeto_p) continue;
	  for(int i=0;i<nCones;i++){
	    if(dR > 0.1*(i+1))continue;
	    photonIsos[i]+=cand->pt();
	  }
	}
      } // end pfCand loop
      privateData_->dr01chPFIso->push_back(chargedIsos[0]);
      privateData_->dr02chPFIso->push_back(chargedIsos[1]);
      privateData_->dr03chPFIso->push_back(chargedIsos[2]);
      privateData_->dr04chPFIso->push_back(chargedIsos[3]);
      privateData_->dr05chPFIso->push_back(chargedIsos[4]);
      privateData_->dr06chPFIso->push_back(chargedIsos[5]);

    } // end vertex loop
    privateData_->dr01nhPFIso->push_back(neutralIsos[0]);
    privateData_->dr01phPFIso->push_back(photonIsos[0]);
    privateData_->dr02nhPFIso->push_back(neutralIsos[1]);
    privateData_->dr02phPFIso->push_back(photonIsos[1]);
    privateData_->dr03nhPFIso->push_back(neutralIsos[2]);
    privateData_->dr03phPFIso->push_back(photonIsos[2]);
    privateData_->dr04nhPFIso->push_back(neutralIsos[3]);
    privateData_->dr04phPFIso->push_back(photonIsos[3]);
    privateData_->dr05nhPFIso->push_back(neutralIsos[4]);
    privateData_->dr05phPFIso->push_back(photonIsos[4]);
    privateData_->dr06nhPFIso->push_back(neutralIsos[5]);
    privateData_->dr06phPFIso->push_back(photonIsos[5]);

  } else {
    privateData_->fiducialFlags->push_back(-1);
    privateData_->recoFlags->push_back(-1);
    privateData_->superClusterIndex->push_back( -1 );
    privateData_->PFsuperClusterIndex->push_back( -1 );

    privateData_->hOverE->push_back(-999.);
    privateData_->hTowOverE->push_back(-999.);
    privateData_->dr03TkSumPt->push_back(-999.);
    privateData_->dr03HollowTkSumPt->push_back(-999.);
    privateData_->dr03EcalRecHitSumEt->push_back(-999.);
    privateData_->dr03HcalTowerSumEt->push_back(-999.);
    privateData_->dr04TkSumPt->push_back(-999.);
    privateData_->dr04HollowTkSumPt->push_back(-999.);
    privateData_->dr04EcalRecHitSumEt->push_back(-999.);
    privateData_->dr04HcalTowerSumEt->push_back(-999.);
    privateData_->chargedHadronIso->push_back(-999.);
    privateData_->neutralHadronIso->push_back(-999.);
    privateData_->photonIso->push_back(-999.);
    privateData_->hasPixelSeed->push_back(0);
    privateData_->hasMatchedConversion->push_back(false);

    privateData_->dr01nhPFIso->push_back(-999.);
    privateData_->dr01phPFIso->push_back(-999.);

    privateData_->dr02nhPFIso->push_back(-999.);
    privateData_->dr02phPFIso->push_back(-999.);

    privateData_->dr03nhPFIso->push_back(-999.);
    privateData_->dr03phPFIso->push_back(-999.);

    privateData_->dr04nhPFIso->push_back(-999.);
    privateData_->dr04phPFIso->push_back(-999.);

    privateData_->dr05nhPFIso->push_back(-999.);
    privateData_->dr05phPFIso->push_back(-999.);

    privateData_->dr06nhPFIso->push_back(-999.);
    privateData_->dr06phPFIso->push_back(-999.);

    for(unsigned int i=0;i<PVs->size();i++){
      privateData_->dr01chPFIso->push_back(-999.);
      privateData_->dr02chPFIso->push_back(-999.);
      privateData_->dr03chPFIso->push_back(-999.);
      privateData_->dr04chPFIso->push_back(-999.);
      privateData_->dr05chPFIso->push_back(-999.);
      privateData_->dr06chPFIso->push_back(-999.);
    }
  }

}

void CmsPhotonFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"fiducialFlags"+colSuffix).c_str(), *privateData_->fiducialFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlags"+colSuffix).c_str(), *privateData_->recoFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"superClusterIndex"+colSuffix).c_str(), *privateData_->superClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFsuperClusterIndex"+colSuffix).c_str(), *privateData_->PFsuperClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverE"+colSuffix).c_str(), *privateData_->hOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hTowOverE"+colSuffix).c_str(), *privateData_->hTowOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03TkSumPt"+colSuffix).c_str(), *privateData_->dr03TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HollowTkSumPt"+colSuffix).c_str(), *privateData_->dr03HollowTkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr03EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr03HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04TkSumPt"+colSuffix).c_str(), *privateData_->dr04TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HollowTkSumPt"+colSuffix).c_str(), *privateData_->dr04HollowTkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr04EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr04HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedHadronIso"+colSuffix).c_str(), *privateData_->chargedHadronIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronIso"+colSuffix).c_str(), *privateData_->neutralHadronIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"photonIso"+colSuffix).c_str(), *privateData_->photonIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hasPixelSeed"+colSuffix).c_str(), *privateData_->hasPixelSeed, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hasMatchedConversion"+colSuffix).c_str(),  *privateData_->hasMatchedConversion, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr01NeutralHadronPFIso"+colSuffix).c_str(),  *privateData_->dr01nhPFIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr01PhotonPFIso"+colSuffix).c_str(),  *privateData_->dr01phPFIso, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"dr02NeutralHadronPFIso"+colSuffix).c_str(),  *privateData_->dr02nhPFIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr02PhotonPFIso"+colSuffix).c_str(),  *privateData_->dr02phPFIso, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"dr03NeutralHadronPFIso"+colSuffix).c_str(),  *privateData_->dr03nhPFIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03PhotonPFIso"+colSuffix).c_str(),  *privateData_->dr03phPFIso, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"dr04NeutralHadronPFIso"+colSuffix).c_str(),  *privateData_->dr04nhPFIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04PhotonPFIso"+colSuffix).c_str(),  *privateData_->dr04phPFIso, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"dr05NeutralHadronPFIso"+colSuffix).c_str(),  *privateData_->dr05nhPFIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr05PhotonPFIso"+colSuffix).c_str(),  *privateData_->dr05phPFIso, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"dr06NeutralHadronPFIso"+colSuffix).c_str(),  *privateData_->dr06nhPFIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr06PhotonPFIso"+colSuffix).c_str(),  *privateData_->dr06phPFIso, nCandString.c_str(), 0, "Reco");

  std::string nCandPVString = colPrefix+(*trkIndexName_)+"PerPV"+colSuffix; 
  cmstree->column((colPrefix+"dr01ChargedHadronPFIso"+colSuffix).c_str(),  *privateData_->dr01chPFIso, nCandPVString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr02ChargedHadronPFIso"+colSuffix).c_str(),  *privateData_->dr02chPFIso, nCandPVString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03ChargedHadronPFIso"+colSuffix).c_str(),  *privateData_->dr03chPFIso, nCandPVString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04ChargedHadronPFIso"+colSuffix).c_str(),  *privateData_->dr04chPFIso, nCandPVString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr05ChargedHadronPFIso"+colSuffix).c_str(),  *privateData_->dr05chPFIso, nCandPVString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr06ChargedHadronPFIso"+colSuffix).c_str(),  *privateData_->dr06chPFIso, nCandPVString.c_str(), 0, "Reco");
  
}




void CmsPhotonFillerData::initialise() {
  
  initialiseCandidate();

  fiducialFlags = new vector<int>;
  recoFlags = new vector<int>;

  superClusterIndex = new vector<int>;
  PFsuperClusterIndex = new vector<int>;

  hOverE                   = new vector<float>;
  hTowOverE                = new vector<float>;
  dr03TkSumPt              = new vector<float>;
  dr03HollowTkSumPt        = new vector<float>;
  dr03EcalRecHitSumEt      = new vector<float>;
  dr03HcalTowerSumEt       = new vector<float>;
  dr04TkSumPt              = new vector<float>;
  dr04HollowTkSumPt        = new vector<float>;
  dr04EcalRecHitSumEt      = new vector<float>;
  dr04HcalTowerSumEt       = new vector<float>;
  chargedHadronIso         = new vector<float>;
  neutralHadronIso         = new vector<float>;
  photonIso                = new vector<float>;
  hasPixelSeed             = new vector<int>;
  hasMatchedConversion     = new vector<bool>;
  
  dr01chPFIso              = new vector<float>;
  dr01nhPFIso              = new vector<float>;
  dr01phPFIso              = new vector<float>;

  dr02chPFIso              = new vector<float>;
  dr02nhPFIso              = new vector<float>;
  dr02phPFIso              = new vector<float>;

  dr03chPFIso              = new vector<float>;
  dr03nhPFIso              = new vector<float>;
  dr03phPFIso              = new vector<float>;

  dr04chPFIso              = new vector<float>;
  dr04nhPFIso              = new vector<float>;
  dr04phPFIso              = new vector<float>;

  dr05chPFIso              = new vector<float>;
  dr05nhPFIso              = new vector<float>;
  dr05phPFIso              = new vector<float>;

  dr06chPFIso              = new vector<float>;
  dr06nhPFIso              = new vector<float>;
  dr06phPFIso              = new vector<float>;

}

void CmsPhotonFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  fiducialFlags->clear();
  recoFlags->clear();

  superClusterIndex->clear();
  PFsuperClusterIndex->clear();

  dr03TkSumPt              ->clear();
  dr03HollowTkSumPt        ->clear();
  dr03EcalRecHitSumEt      ->clear();
  dr03HcalTowerSumEt       ->clear();
  dr04TkSumPt              ->clear();
  dr04HollowTkSumPt        ->clear();
  dr04EcalRecHitSumEt      ->clear();
  dr04HcalTowerSumEt       ->clear();
  
  chargedHadronIso         ->clear();
  neutralHadronIso         ->clear();
  photonIso                ->clear();
  hasPixelSeed             ->clear();
  hasMatchedConversion     ->clear();

  dr01chPFIso              ->clear();
  dr01nhPFIso              ->clear();
  dr01phPFIso              ->clear();

  dr02chPFIso              ->clear();
  dr02nhPFIso              ->clear();
  dr02phPFIso              ->clear();

  dr03chPFIso              ->clear();
  dr03nhPFIso              ->clear();
  dr03phPFIso              ->clear();

  dr04chPFIso              ->clear();
  dr04nhPFIso              ->clear();
  dr04phPFIso              ->clear();

  dr05chPFIso              ->clear();
  dr05nhPFIso              ->clear();
  dr05phPFIso              ->clear();

  dr06chPFIso              ->clear();
  dr06nhPFIso              ->clear();
  dr06phPFIso              ->clear();

}
