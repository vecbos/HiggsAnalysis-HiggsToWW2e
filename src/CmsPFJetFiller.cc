//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsPFJetFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Sep  29 11:01:00 CEST 2008
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFJetFiller.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <TVector3.h>
#include <TLorentzVector.h>

#include <string>

using namespace edm;
using namespace reco;



QGLikelihoodVars computeQGLikelihoodVars( const PFJet* pfjet, float R=0., float ptratio=0.);


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPFJetFiller::CmsPFJetFiller(CmsTree *cmsTree, 
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFJetFillerData)
{
  cmstree=cmsTree;

  saveJetBTag_ = false;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsPFJetFiller::~CmsPFJetFiller() {
  
  // delete here the vector ptr's
  delete privateData_->charge;
  delete privateData_->energy;
  delete privateData_->et;
  delete privateData_->momentum;
  delete privateData_->theta;
  delete privateData_->eta;
  delete privateData_->phi;
  delete privateData_->x;
  delete privateData_->y;
  delete privateData_->z;
  delete privateData_->uncorrx;
  delete privateData_->uncorry;
  delete privateData_->uncorrz;
  delete privateData_->vertexX;
  delete privateData_->vertexY;
  delete privateData_->vertexZ;
  delete privateData_->ncand;

  delete privateData_->chargedHadronEnergy;
  delete privateData_->neutralHadronEnergy;
  delete privateData_->photonEnergy;
  delete privateData_->electronEnergy;
  delete privateData_->muonEnergy;
  delete privateData_->HFHadronEnergy;
  delete privateData_->HFEMEnergy;
  delete privateData_->chargedHadronMultiplicity;
  delete privateData_->neutralHadronMultiplicity;
  delete privateData_->photonMultiplicity;
  delete privateData_->electronMultiplicity;
  delete privateData_->muonMultiplicity;
  delete privateData_->HFHadronMultiplicity;
  delete privateData_->HFEMMultiplicity;
  delete privateData_->combinedSecondaryVertexBJetTags;
  delete privateData_->combinedSecondaryVertexMVABJetTags; 
  delete privateData_->jetBProbabilityBJetTags; 
  delete privateData_->jetProbabilityBJetTags;
  delete privateData_->simpleSecondaryVertexHighEffBJetTags;
  delete privateData_->simpleSecondaryVertexHighPurBJetTags;
  delete privateData_->trackCountingHighPurBJetTags;
  delete privateData_->trackCountingHighEffBJetTags;
  delete privateData_->trackCountingVeryHighEffBJetTags;
  delete privateData_->uncorrEnergy;
  delete privateData_->area;
  delete privateData_->weightedDz1;
  delete privateData_->weightedDz2;

  // for backward compatibility with existing trees
  delete privateData_->chargedEmEnergy;
  delete privateData_->neutralEmEnergy;

  // additional variables for Marini's likelihood:
  delete privateData_->ptD;
  delete privateData_->rmsCand;
  delete privateData_->axis1;
  delete privateData_->axis2;
  delete privateData_->pull;
  delete privateData_->r_ch;
  delete privateData_->tana;

  // additional variables for PU studies
  delete privateData_->jetIdMvaSimple;
  delete privateData_->jetIdMvaFull;
  delete privateData_->jetIdMvaPhilV1;
  delete privateData_->nChargedIdMva;
  delete privateData_->nNeutralsIdMva;
  delete privateData_->dZIdMva;
  delete privateData_->dR2MeanIdMva;
  delete privateData_->dRMeanIdMva;
  delete privateData_->frac01IdMva;
  delete privateData_->frac02IdMva;
  delete privateData_->frac03IdMva;
  delete privateData_->frac04IdMva;
  delete privateData_->frac05IdMva;
  delete privateData_->betaIdMva;
  delete privateData_->betastarIdMva;
  delete privateData_->betastarclassicIdMva;
  delete privateData_->betastar;
  delete privateData_->rmsCandsHand;
  delete privateData_;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsPFJetFiller::saveJetBTag(bool what) { saveJetBTag_=what; }

void CmsPFJetFiller::setBTags(edm::ParameterSet btagcollections) { BTagCollections_ = btagcollections; }

// Set boolean control options for quantities that are written out

void CmsPFJetFiller::writeCollectionToTree(edm::InputTag collectionTag,
                                           const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                           const std::string &columnPrefix, const std::string &columnSuffix,
                                           bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  // the same, but in PFJetCollection format
  edm::Handle<reco::PFJetCollection> jets;      
  try { iEvent.getByLabel(collectionTag, jets); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << collectionTag; }

  edm::Handle< reco::VertexCollection>  primaryVertexColl  ;
  try { iEvent.getByLabel(_verticesTag, primaryVertexColl); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get primary vertex collection: offlinePrimaryVertices"; }
  VertexCollection::const_iterator vMax = primaryVertexColl->begin();
  bestPrimaryVertex_ = *vMax;

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFJetFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFJetFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    edm::Handle<reco::JetTagCollection> combinedSecondaryVertexBJetTags,
      combinedSecondaryVertexMVABJetTags, 
      jetBProbabilityBJetTags, 
      jetProbabilityBJetTags,
      simpleSecondaryVertexHighEffBJetTags,
      simpleSecondaryVertexHighPurBJetTags,
      trackCountingHighPurBJetTags,
      trackCountingHighEffBJetTags,
      trackCountingVeryHighEffBJetTags;

    if(saveJetBTag_) {
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("combinedSecondaryVertexBJetTags"), combinedSecondaryVertexBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("combinedSecondaryVertexMVABJetTags"), combinedSecondaryVertexMVABJetTags); 
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("jetBProbabilityBJetTags"), jetBProbabilityBJetTags); 
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("jetProbabilityBJetTags"), jetProbabilityBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("simpleSecondaryVertexHighEffBJetTags"), simpleSecondaryVertexHighEffBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("simpleSecondaryVertexHighPurBJetTags"), simpleSecondaryVertexHighPurBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("trackCountingHighPurBJetTags"), trackCountingHighPurBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("trackCountingHighEffBJetTags"), trackCountingHighEffBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("trackCountingVeryHighEffBJetTags"), trackCountingVeryHighEffBJetTags);
    }

    const JetCorrector* corrector = JetCorrector::getJetCorrector (m_jcs, iSetup);

    int index = 0;
    edm::View<reco::Candidate>::const_iterator cand;
    reco::PFJetCollection::const_iterator jetit=jets->begin();
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill jet extra informations

      // em, had fractions
      const PFJet *thisPFJet = dynamic_cast< const PFJet * > ( &(*cand) );
      const PFJetRef thisPFJetRef = collection->refAt(index).castTo<PFJetRef>();
      if( thisPFJet != 0 ) { 
        privateData_->chargedHadronEnergy->push_back( thisPFJet->chargedHadronEnergy() );
        privateData_->neutralHadronEnergy->push_back( thisPFJet->neutralHadronEnergy() );
        privateData_->photonEnergy->push_back( thisPFJet->photonEnergy() );
        privateData_->electronEnergy->push_back( thisPFJet->electronEnergy() );
        privateData_->muonEnergy->push_back( thisPFJet->muonEnergy() );
        privateData_->HFHadronEnergy->push_back( thisPFJet->HFHadronEnergy() );
        privateData_->HFEMEnergy->push_back( thisPFJet->HFEMEnergy() );
        privateData_->chargedHadronMultiplicity->push_back( thisPFJet->chargedHadronMultiplicity() );
        privateData_->neutralHadronMultiplicity->push_back( thisPFJet->neutralHadronMultiplicity() );
        privateData_->photonMultiplicity->push_back( thisPFJet->photonMultiplicity() );
        privateData_->electronMultiplicity->push_back( thisPFJet->electronMultiplicity() );
        privateData_->muonMultiplicity->push_back( thisPFJet->muonMultiplicity() );
        privateData_->HFHadronMultiplicity->push_back( thisPFJet->HFHadronMultiplicity() );
        privateData_->HFEMMultiplicity->push_back( thisPFJet->HFEMMultiplicity() );
        
        // for backward compatibility with existing trees
        privateData_->chargedEmEnergy->push_back( thisPFJet->chargedEmEnergy() );
        privateData_->neutralEmEnergy->push_back( thisPFJet->neutralEmEnergy() );
        privateData_->area->push_back( thisPFJet->jetArea() );

        // compute marini's variables
        QGLikelihoodVars qgvars = computeQGLikelihoodVars(thisPFJet);
        privateData_->ptD->push_back( qgvars.ptD );
        privateData_->rmsCand->push_back( qgvars.rmsCand );
        privateData_->axis1->push_back( qgvars.axis1 );
        privateData_->axis2->push_back( qgvars.axis2 );
        privateData_->pull->push_back( qgvars.pull );
        privateData_->r_ch->push_back( qgvars.r_ch );
        privateData_->tana->push_back( qgvars.tana );

	float sumPt_cands,sumPt2_cands,rms_cands,sumTrkPtBetaStar,sumTrkPt,betastar_;
        sumPt_cands=sumPt2_cands=rms_cands=sumTrkPtBetaStar=sumTrkPt=betastar_=0.0;
        float dznum1, dznum2, dzdenom;
        dznum1 = dznum2 = dzdenom = 0.0;
        float dzjet1_, dzjet2_;
        dzjet1_ = dzjet2_ = -1.0;
        
	vector<reco::PFCandidatePtr> constituents = thisPFJet->getPFConstituents();
        for ( std::vector<reco::PFCandidatePtr>::const_iterator pfcand=constituents.begin(); pfcand != constituents.end(); pfcand++) {          
          // compute RMS
	  math::XYZTLorentzVectorD const& p4t = (*pfcand)->p4();
	  TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());   // of the candidate
	  TLorentzVector jetp4;
	  jetp4.SetPtEtaPhiE(thisPFJet->pt(), thisPFJet->eta(), thisPFJet->phi(), thisPFJet->energy());  // of the jet
	  if(p4.Pt()!=0){
	    sumPt_cands += p4.Pt();
	    sumPt2_cands += (p4.Pt()*p4.Pt());
	    float deltaR = jetp4.DeltaR(p4);
	    rms_cands += (p4.Pt()*p4.Pt()*deltaR*deltaR);
	  }
          
          // and dzjet, betastar
          reco::TrackRef i_trk = (*pfcand)->trackRef();
          if(i_trk.isNonnull()) {

            // calculating average dzjet. The two output differ only in the calculation of corrected deltaZ of the track wrt the PV
            // first is Guillelmo's deltaZ method (below). second is Track dsz method.
            float dzCorr1 = DzCorrected(i_trk,bestPrimaryVertex_);
            float dzCorr2 = i_trk->dsz(bestPrimaryVertex_.position());
            float pt = (*pfcand)->pt();
            dznum1 += pt*pt * dzCorr1;
            dznum2 += pt*pt * dzCorr2;
            dzdenom += pt*pt;

            sumTrkPt += pt;

            bool isFirstVtx = false;
            for(reco::Vertex::trackRef_iterator i_vtxTrk = bestPrimaryVertex_.tracks_begin(); i_vtxTrk != bestPrimaryVertex_.tracks_end(); ++i_vtxTrk) {
              reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>()); // match the jet track to the track from the vertex
              if (trkRef == i_trk) {
                isFirstVtx=true; 
                break;
              }
            }

            bool isOtherVtx = false;
            if (!isFirstVtx) {   // if not associated to PV check others... 
              for(unsigned iotherVtx=1; iotherVtx<primaryVertexColl->size();iotherVtx++) {
                if (!((*primaryVertexColl)[iotherVtx].isFake()) && 
                    (*primaryVertexColl)[iotherVtx].ndof() >= 4 && 
                    fabs((*primaryVertexColl)[iotherVtx].z()) <= 24) {
                  for(reco::Vertex::trackRef_iterator i_vtxTrk = (*primaryVertexColl)[iotherVtx].tracks_begin(); 
                      i_vtxTrk != (*primaryVertexColl)[iotherVtx].tracks_end(); ++i_vtxTrk) {
                    reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
                    if (trkRef == i_trk) {
                      isOtherVtx=true;
                      break;
                    }
                  }
                }}}
            if(!isFirstVtx && isOtherVtx) { sumTrkPtBetaStar += i_trk->pt(); } 
          }
        }  // loop overt tracks
        
        if (sumTrkPt > 0) {
          using namespace std;
          betastar_ = sumTrkPtBetaStar/sumTrkPt;   
          dzjet1_ = dznum1/dzdenom;
          dzjet2_ = dznum2/dzdenom;
        } else {
          betastar_ = -999;
          dzjet1_ = -999;
          dzjet2_ = -999.;
        }

        privateData_->weightedDz1->push_back( dzjet1_ );
        privateData_->weightedDz2->push_back( dzjet2_ );
        
        float forRms = rms_cands/sumPt2_cands;
	privateData_->betastar->push_back(betastar_);
	privateData_->rmsCandsHand->push_back(forRms);

      } else {
        privateData_->chargedHadronEnergy->push_back( -1. );
        privateData_->neutralHadronEnergy->push_back( -1. );
        privateData_->photonEnergy->push_back( -1. );
        privateData_->electronEnergy->push_back( -1. );
        privateData_->muonEnergy->push_back( -1. );
        privateData_->chargedHadronMultiplicity->push_back( -1. );
        privateData_->neutralHadronMultiplicity->push_back( -1. );
        privateData_->photonMultiplicity->push_back( -1. );
        privateData_->electronMultiplicity->push_back( -1. );
        privateData_->muonMultiplicity->push_back( -1. );
        privateData_->weightedDz1->push_back( -1. );
        privateData_->weightedDz2->push_back( -1. );

        // for backward compatibility with existing trees
        privateData_->chargedEmEnergy->push_back( -1. );
        privateData_->neutralEmEnergy->push_back( -1. );

	// for PU studies
	privateData_->betastar->push_back(-999.);            // chiara
	privateData_->rmsCandsHand->push_back(-999.);        // chiara
 	privateData_->jetIdMvaSimple->push_back(-1.);
 	privateData_->jetIdMvaFull->push_back(-1.);
 	privateData_->jetIdMvaPhilV1->push_back(-1.);
 	privateData_->nChargedIdMva->push_back(-1.);
 	privateData_->nNeutralsIdMva->push_back(-1.);
 	privateData_->dZIdMva->push_back(-1.);
 	privateData_->dR2MeanIdMva->push_back(-1.);
 	privateData_->dRMeanIdMva->push_back(-1.);
 	privateData_->frac01IdMva->push_back(-1.);
 	privateData_->frac02IdMva->push_back(-1.);
 	privateData_->frac03IdMva->push_back(-1.);
 	privateData_->frac04IdMva->push_back(-1.);
 	privateData_->frac05IdMva->push_back(-1.);
 	privateData_->betaIdMva->push_back(-1.);
 	privateData_->betastarIdMva->push_back(-1.);
 	privateData_->betastarclassicIdMva->push_back(-1.);
      }

      // fill the btag algorithms output
      if(saveJetBTag_) {
        privateData_->combinedSecondaryVertexBJetTags->push_back( (*combinedSecondaryVertexBJetTags)[index].second );
        privateData_->combinedSecondaryVertexMVABJetTags->push_back( (*combinedSecondaryVertexMVABJetTags)[index].second ); 
        privateData_->jetBProbabilityBJetTags->push_back( (*jetBProbabilityBJetTags)[index].second ); 
        privateData_->jetProbabilityBJetTags->push_back( (*jetProbabilityBJetTags)[index].second );
        privateData_->simpleSecondaryVertexHighEffBJetTags->push_back( (*simpleSecondaryVertexHighEffBJetTags)[index].second );
        privateData_->simpleSecondaryVertexHighPurBJetTags->push_back( (*simpleSecondaryVertexHighPurBJetTags)[index].second );
        privateData_->trackCountingHighPurBJetTags->push_back( (*trackCountingHighPurBJetTags)[index].second );
        privateData_->trackCountingHighEffBJetTags->push_back( (*trackCountingHighEffBJetTags)[index].second );
        privateData_->trackCountingVeryHighEffBJetTags->push_back( (*trackCountingVeryHighEffBJetTags)[index].second );
      } else {
        privateData_->combinedSecondaryVertexBJetTags->push_back( -1. );
        privateData_->combinedSecondaryVertexMVABJetTags->push_back( -1. ); 
        privateData_->jetBProbabilityBJetTags->push_back( -1. ); 
        privateData_->jetProbabilityBJetTags->push_back( -1. );
        privateData_->simpleSecondaryVertexHighEffBJetTags->push_back( -1. );
        privateData_->simpleSecondaryVertexHighPurBJetTags->push_back( -1. );
        privateData_->trackCountingHighPurBJetTags->push_back( -1. );
        privateData_->trackCountingHighEffBJetTags->push_back( -1. );
        privateData_->trackCountingVeryHighEffBJetTags->push_back( -1. );
      }

      // run the correction on the fly.
      PFJet  correctedJet = *jetit;
      float scale = corrector->correction(correctedJet,iEvent,iSetup);
      correctedJet.scaleEnergy(scale);

      privateData_->charge->push_back((int)correctedJet.charge());
      privateData_->energy->push_back(correctedJet.energy());
      privateData_->et->push_back(correctedJet.et());
      privateData_->momentum->push_back(correctedJet.p());
      privateData_->theta->push_back(correctedJet.theta());
      privateData_->eta->push_back(correctedJet.eta());
      privateData_->phi->push_back(correctedJet.phi());
      privateData_->x->push_back(correctedJet.momentum().x());
      privateData_->y->push_back(correctedJet.momentum().y());
      privateData_->z->push_back(correctedJet.momentum().z());
      privateData_->vertexX->push_back(correctedJet.vx());
      privateData_->vertexY->push_back(correctedJet.vy());
      privateData_->vertexZ->push_back(correctedJet.vz());
      privateData_->uncorrx->push_back(thisPFJet->momentum().x());
      privateData_->uncorry->push_back(thisPFJet->momentum().y());
      privateData_->uncorrz->push_back(thisPFJet->momentum().z());
      privateData_->uncorrEnergy->push_back(thisPFJet->energy());

      PileupJetIdentifier jetIdentifer_vars = _jetId_algos[0]->computeIdVariables(&correctedJet, scale, &bestPrimaryVertex_, *(primaryVertexColl.product()));
      privateData_->dRMeanIdMva->push_back( jetIdentifer_vars.dRMean() );
      privateData_->frac01IdMva->push_back( jetIdentifer_vars.frac01() );
      privateData_->frac02IdMva->push_back( jetIdentifer_vars.frac02() );
      privateData_->frac03IdMva->push_back( jetIdentifer_vars.frac03() );
      privateData_->frac04IdMva->push_back( jetIdentifer_vars.frac04() );
      privateData_->frac05IdMva->push_back( jetIdentifer_vars.frac05() );
      privateData_->nNeutralsIdMva->push_back( jetIdentifer_vars.nNeutrals() );
      privateData_->betaIdMva->push_back( jetIdentifer_vars.beta() );
      privateData_->betastarIdMva->push_back( jetIdentifer_vars.betaStar() );
      privateData_->dZIdMva->push_back( jetIdentifer_vars.dZ() );
      privateData_->nChargedIdMva->push_back( jetIdentifer_vars.nCharged() );
      privateData_->dR2MeanIdMva->push_back( jetIdentifer_vars.dR2Mean() );
      privateData_->betastarclassicIdMva->push_back( jetIdentifer_vars.betaStarClassic() );

      for(unsigned int imva=0; imva<_jetId_algos.size(); imva++){
        PileupJetIdAlgo* ialgo = (_jetId_algos[imva]);
        ialgo->set(jetIdentifer_vars);
        PileupJetIdentifier id = ialgo->computeMva();
        // mva values
        if (imva==0) privateData_->jetIdMvaSimple->push_back( id.mva() );
        if (imva==1) privateData_->jetIdMvaFull->push_back( id.mva() );
        if (imva==2) privateData_->jetIdMvaPhilV1->push_back( id.mva() );
      }

      index++;
      jetit++;
    }
  } else {
    *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  treeJetInfo(columnPrefix,columnSuffix);
  
  if(dumpData) cmstree->dumpData();
  
}

void CmsPFJetFiller::treeCandInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"uncorrenergy"+colSuffix).c_str(), *privateData_->uncorrEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"uncorrpx"+colSuffix).c_str(), *privateData_->uncorrx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"uncorrpy"+colSuffix).c_str(), *privateData_->uncorry, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"uncorrpz"+colSuffix).c_str(), *privateData_->uncorrz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexX"+colSuffix).c_str(), *privateData_->vertexX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexY"+colSuffix).c_str(), *privateData_->vertexY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexZ"+colSuffix).c_str(), *privateData_->vertexZ, nCandString.c_str(), 0, "Reco");
}

void CmsPFJetFiller::treeJetInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"chargedHadronEnergy"+colSuffix).c_str(), *privateData_->chargedHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronEnergy"+colSuffix).c_str(), *privateData_->neutralHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"photonEnergy"+colSuffix).c_str(), *privateData_->photonEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"electronEnergy"+colSuffix).c_str(), *privateData_->electronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muonEnergy"+colSuffix).c_str(), *privateData_->muonEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFHadronEnergy"+colSuffix).c_str(), *privateData_->HFHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFEMEnergy"+colSuffix).c_str(), *privateData_->HFEMEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedHadronMultiplicity"+colSuffix).c_str(), *privateData_->chargedHadronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronMultiplicity"+colSuffix).c_str(), *privateData_->neutralHadronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"photonMultiplicity"+colSuffix).c_str(), *privateData_->photonMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"electronMultiplicity"+colSuffix).c_str(), *privateData_->electronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muonMultiplicity"+colSuffix).c_str(), *privateData_->muonMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFHadronMultiplicity"+colSuffix).c_str(), *privateData_->HFHadronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFEMMultiplicity"+colSuffix).c_str(), *privateData_->HFEMMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"area"+colSuffix).c_str(), *privateData_->area, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"weightedDz1"+colSuffix).c_str(), *privateData_->weightedDz1, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"weightedDz2"+colSuffix).c_str(), *privateData_->weightedDz2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"betastar"+colSuffix).c_str(), *privateData_->betastar, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"rmsCandsHand"+colSuffix).c_str(), *privateData_->rmsCandsHand, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"jetIdMvaSimple"+colSuffix).c_str(), *privateData_->jetIdMvaSimple, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"jetIdMvaFull"+colSuffix).c_str(), *privateData_->jetIdMvaFull, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"jetIdMvaPhilV1"+colSuffix).c_str(), *privateData_->jetIdMvaPhilV1, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nChargedIdMva"+colSuffix).c_str(), *privateData_->nChargedIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nNeutralsIdMva"+colSuffix).c_str(), *privateData_->nNeutralsIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dZIdMva"+colSuffix).c_str(), *privateData_->dZIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dR2MeanIdMva"+colSuffix).c_str(), *privateData_->dR2MeanIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dRMeanIdMva"+colSuffix).c_str(), *privateData_->dRMeanIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"frac01IdMva"+colSuffix).c_str(), *privateData_->frac01IdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"frac02IdMva"+colSuffix).c_str(), *privateData_->frac02IdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"frac03IdMva"+colSuffix).c_str(), *privateData_->frac03IdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"frac04IdMva"+colSuffix).c_str(), *privateData_->frac04IdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"frac05IdMva"+colSuffix).c_str(), *privateData_->frac05IdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"betaIdMva"+colSuffix).c_str(), *privateData_->betaIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"betastarIdMva"+colSuffix).c_str(), *privateData_->betastarIdMva, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"betastarclassicIdMva"+colSuffix).c_str(), *privateData_->betastarclassicIdMva, nCandString.c_str(), 0, "Reco");

  // for backward compatibility with existing trees 
  cmstree->column((colPrefix+"chargedEmEnergy"+colSuffix).c_str(), *privateData_->chargedEmEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralEmEnergy"+colSuffix).c_str(), *privateData_->neutralEmEnergy, nCandString.c_str(), 0, "Reco");
  if(saveJetBTag_) {
    cmstree->column((colPrefix+"combinedSecondaryVertexBJetTags"+colSuffix).c_str(), *privateData_->combinedSecondaryVertexBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"combinedSecondaryVertexMVABJetTags"+colSuffix).c_str(), *privateData_->combinedSecondaryVertexMVABJetTags, nCandString.c_str(), 0, "Reco"); 
    cmstree->column((colPrefix+"jetBProbabilityBJetTags"+colSuffix).c_str(), *privateData_->jetBProbabilityBJetTags, nCandString.c_str(), 0, "Reco"); 
    cmstree->column((colPrefix+"jetProbabilityBJetTags"+colSuffix).c_str(), *privateData_->jetProbabilityBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"simpleSecondaryVertexHighEffBJetTags"+colSuffix).c_str(), *privateData_->simpleSecondaryVertexHighEffBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"simpleSecondaryVertexHighPurBJetTags"+colSuffix).c_str(), *privateData_->simpleSecondaryVertexHighPurBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighPurBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighPurBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighEffBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighEffBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingVeryHighEffBJetTags"+colSuffix).c_str(), *privateData_->trackCountingVeryHighEffBJetTags, nCandString.c_str(), 0, "Reco");
  }
  cmstree->column((colPrefix+"rmsCand"+colSuffix).c_str(), *privateData_->rmsCand, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ptD"+colSuffix).c_str(), *privateData_->ptD, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"axis1"+colSuffix).c_str(), *privateData_->axis1, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"axis2"+colSuffix).c_str(), *privateData_->axis2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pull"+colSuffix).c_str(), *privateData_->pull, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"r_ch"+colSuffix).c_str(), *privateData_->r_ch, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"tana"+colSuffix).c_str(), *privateData_->tana, nCandString.c_str(), 0, "Reco");
}







void CmsPFJetFillerData::initialise() {

  charge = new vector<int>;
  energy = new vector<float>;
  et = new vector<float>;
  momentum = new vector<float>;
  theta = new vector<float>;
  eta = new vector<float>;
  phi = new vector<float>;
  x = new vector<float>;
  y = new vector<float>;
  z = new vector<float>;
  uncorrx = new vector<float>;
  uncorry = new vector<float>;
  uncorrz = new vector<float>;
  vertexX = new vector<float>;
  vertexY = new vector<float>;
  vertexZ = new vector<float>;
  ncand = new int;

  chargedHadronEnergy = new vector<float>;
  neutralHadronEnergy = new vector<float>;
  chargedEmEnergy = new vector<float>;
  neutralEmEnergy = new vector<float>;
  photonEnergy = new vector<float>;
  electronEnergy = new vector<float>;
  muonEnergy = new vector<float>;
  HFHadronEnergy = new vector<float>;
  HFEMEnergy = new vector<float>;
  chargedHadronMultiplicity = new vector<int>;
  neutralHadronMultiplicity = new vector<int>;
  photonMultiplicity = new vector<int>;
  electronMultiplicity = new vector<int>;
  muonMultiplicity = new vector<int>;
  HFHadronMultiplicity = new vector<int>;
  HFEMMultiplicity = new vector<int>;
  combinedSecondaryVertexBJetTags = new vector<float>;
  combinedSecondaryVertexMVABJetTags = new vector<float>; 
  jetBProbabilityBJetTags = new vector<float>; 
  jetProbabilityBJetTags = new vector<float>;
  simpleSecondaryVertexHighEffBJetTags = new vector<float>;
  simpleSecondaryVertexHighPurBJetTags = new vector<float>;
  trackCountingHighPurBJetTags = new vector<float>;
  trackCountingHighEffBJetTags = new vector<float>;
  trackCountingVeryHighEffBJetTags = new vector<float>;
  uncorrEnergy = new vector<float>;
  area = new vector<float>;
  weightedDz1 = new vector<float>;
  weightedDz2 = new vector<float>;
  ptD = new vector<float>;
  rmsCand = new vector<float>;
  axis1 = new vector<float>;
  axis2 = new vector<float>;
  pull = new vector<float>;
  r_ch = new vector<float>;
  tana = new vector<float>;
  betastar = new vector<float>;
  jetIdMvaSimple = new vector<float>;
  jetIdMvaFull = new vector<float>;
  jetIdMvaPhilV1 = new vector<float>;
  nChargedIdMva = new vector<float>;
  nNeutralsIdMva = new vector<float>;
  dZIdMva = new vector<float>;
  dR2MeanIdMva = new vector<float>;
  dRMeanIdMva = new vector<float>;
  frac01IdMva = new vector<float>;
  frac02IdMva = new vector<float>;
  frac03IdMva = new vector<float>;
  frac04IdMva = new vector<float>;
  frac05IdMva = new vector<float>;
  betaIdMva = new vector<float>;
  betastarIdMva = new vector<float>;
  betastarclassicIdMva = new vector<float>;
  rmsCandsHand = new vector<float>;
}

void CmsPFJetFillerData::clearTrkVectors() {

  charge->clear();
  energy->clear();
  et->clear();
  momentum->clear();
  theta->clear();
  eta->clear();
  phi->clear();
  x->clear();
  y->clear();
  z->clear();
  uncorrx->clear();
  uncorry->clear();
  uncorrz->clear();
  vertexX->clear();
  vertexY->clear();
  vertexZ->clear();

  chargedHadronEnergy->clear();
  neutralHadronEnergy->clear();
  chargedEmEnergy->clear();
  neutralEmEnergy->clear();
  photonEnergy->clear();
  electronEnergy->clear();
  muonEnergy->clear();
  HFHadronEnergy->clear();
  HFEMEnergy->clear();
  chargedHadronMultiplicity->clear();
  neutralHadronMultiplicity->clear();
  photonMultiplicity->clear();
  electronMultiplicity->clear();
  muonMultiplicity->clear();
  HFHadronMultiplicity->clear();
  HFEMMultiplicity->clear();
  combinedSecondaryVertexBJetTags->clear();
  combinedSecondaryVertexMVABJetTags->clear(); 
  jetBProbabilityBJetTags->clear(); 
  jetProbabilityBJetTags->clear();
  simpleSecondaryVertexHighEffBJetTags->clear();
  simpleSecondaryVertexHighPurBJetTags->clear();
  trackCountingHighPurBJetTags->clear();
  trackCountingHighEffBJetTags->clear();
  trackCountingVeryHighEffBJetTags->clear();
  uncorrEnergy->clear();
  area->clear();
  weightedDz1->clear();
  weightedDz2->clear();
  ptD->clear();
  rmsCand->clear();
  axis1->clear();
  axis2->clear();
  pull->clear();
  r_ch->clear();
  tana->clear();
  betastar -> clear();
  jetIdMvaSimple -> clear();
  jetIdMvaFull -> clear();
  jetIdMvaPhilV1 -> clear();
  nChargedIdMva -> clear();
  nNeutralsIdMva -> clear();
  dZIdMva -> clear();
  dR2MeanIdMva -> clear();
  dRMeanIdMva -> clear();
  frac01IdMva -> clear();
  frac02IdMva -> clear();
  frac03IdMva -> clear();
  frac04IdMva -> clear();
  frac05IdMva -> clear();
  betaIdMva -> clear(); 
  betastarIdMva -> clear();
  betastarclassicIdMva -> clear();
  rmsCandsHand -> clear();
}



QGLikelihoodVars computeQGLikelihoodVars( const PFJet* pfjet, float R, float ptratio ) {

  std::vector<fastjet::PseudoJet> *input_particles = new std::vector<fastjet::PseudoJet>;

  for(long int i=0;i<pfjet->nConstituents();++i)
     {
     input_particles->push_back(fastjet::PseudoJet( pfjet->getJetConstituentsQuick()[i]->px(),
                          pfjet->getJetConstituentsQuick()[i]->py(),
                          pfjet->getJetConstituentsQuick()[i]->pz(),
                          pfjet->getJetConstituentsQuick()[i]->energy()
                      ) );
     }

  
  //Compute rmsCand
  float rmsCand=  -999.;
  float ptD=      -999.;
  float axis1=    -999.;
  float axis2=    -999.;
  float pull=     -999.;
  float r_ch=     -999.;
  float tana =    -999.;
	//compute the variables
  {
	float SumW=0;
	float SumW2=0;
	float SumDeta=0;
	float SumDeta2=0;
	float SumDphi=0;
	float SumDphi2=0;
	float SumDetaDphi=0;
	float SumDR=0;
	
	float Eta0=pfjet->eta();
	float Phi0=pfjet->phi();
	
	for(int i=0;i<pfjet->nConstituents();++i)
		{
		double pt=pfjet->getJetConstituentsQuick()[i]->pt();
		double eta=pfjet->getJetConstituentsQuick()[i]->eta();
		double phi=pfjet->getJetConstituentsQuick()[i]->phi();
		TLorentzVector cand;
		cand.SetPtEtaPhiE(pt, eta, phi, pt); //only interested in direction
		TLorentzVector jet;
		jet.SetPtEtaPhiE(pfjet->pt(), pfjet->eta(), pfjet->phi(), pfjet->pt()); //only interested in direction
		double dphi = 2*atan(tan((phi-Phi0)/2));      
		double deta = eta-Eta0;
		double dR = jet.DeltaR(cand);
		SumW+=pt;
		SumW2+=pt*pt;
		SumDeta+=pt*pt*deta;
		SumDeta2+=pt*pt*deta*deta;
		SumDphi+=pt*pt*dphi;
		SumDphi2+=pt*pt*dphi*dphi;
		SumDetaDphi+=pt*pt*deta*dphi;
		SumDR+=pt*pt*dR*dR;
		}
	float ave_deta = SumDeta/SumW2;
	float ave_dphi = SumDphi/SumW2;
	float  ave_deta2 = SumDeta2/SumW2;
	float  ave_dphi2 = SumDphi2/SumW2;
      float a = ave_deta2-ave_deta*ave_deta;
      float b = ave_dphi2-ave_dphi*ave_dphi;
      float c = -(SumDetaDphi/SumW2-ave_deta*ave_dphi);
      float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
      if (a+b+delta > 0) {
        axis1 = sqrt(0.5*(a+b+delta));
      }
      if (a+b-delta > 0) {  
        axis2 = sqrt(0.5*(a+b-delta));
      }
      if (c != 0) {
        tana = 0.5*(b-a+delta)/c;
      }	

      ptD =sqrt( SumW2/ (SumW*SumW));
      rmsCand = SumDR/SumW2;

      if(pfjet->nConstituents()>0)r_ch=pfjet->getJetConstituentsQuick()[0]->pt()/SumW;
	//-------calculate pull------
    	float ddetaR_sum(0.0), ddphiR_sum(0.0),ddetaR_sum_QC(0.0), ddphiR_sum_QC(0.0);
    	for(int i=0; i<pfjet->nConstituents(); ++i) {
			double pt=pfjet->getJetConstituentsQuick()[i]->pt();
			double eta=pfjet->getJetConstituentsQuick()[i]->eta();
			double phi=pfjet->getJetConstituentsQuick()[i]->phi();
			double dphi = 2*atan(tan((phi-Phi0)/2));      
			double deta = eta-Eta0;
  		    float weight = pt*pt;
  		    float ddeta, ddphi,ddR;
  		    ddeta = deta - ave_deta ;//jetPart_deta[i] - ave_deta ; 
  		    ddphi = 2*atan(tan(( dphi - ave_dphi)/2.)) ;
  		    ddR = sqrt(ddeta*ddeta + ddphi*ddphi);
  		    ddetaR_sum += ddR*ddeta*weight;
  		    ddphiR_sum += ddR*ddphi*weight;
  		  }//second loop over constituents  
  		  if (SumW2 > 0) {
  		    float ddetaR_ave = ddetaR_sum/SumW2;
  		    float ddphiR_ave = ddphiR_sum/SumW2;
  		    pull = sqrt(ddetaR_ave*ddetaR_ave+ddphiR_ave*ddphiR_ave);
  		  }
  } //close compute brackets
  
  //export variables
  QGLikelihoodVars vars;
  vars.ptD = ptD;
  vars.rmsCand = rmsCand;
  vars.axis1= axis1;
  vars.axis2= axis2;
  vars.pull=  pull;
  vars.r_ch=  r_ch;
  vars.tana = tana;

  return vars;

}

float CmsPFJetFiller::DzCorrected(reco::TrackRef trk, reco::Vertex vtx)
{
  // Compute Dxy with respect to a given position
  TVector3 momPerp(trk->px(),trk->py(),0);
  TVector3 posPerp(trk->vx()-vtx.position().x(),trk->vy()-vtx.position().y(),0);
  return trk->vz()-vtx.position().z() - posPerp.Dot(momPerp)/trk->pt() * (trk->pz()/trk->pt());
}
