//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsHcalNoiseFiller
//
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCaloFlagLabels.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsHcalNoiseFiller.h"

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


CmsHcalNoiseFiller::CmsHcalNoiseFiller(CmsTree *cmsTree)
{
   cmstree=cmsTree;
}

CmsHcalNoiseFiller::CmsHcalNoiseFiller(CmsTree *cmsTree, bool fatTree)
{
   cmstree=cmsTree;
}


//--------------
// Destructor --
//--------------

CmsHcalNoiseFiller::~CmsHcalNoiseFiller()
{
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsHcalNoiseFiller::writeHcalNoiseSummaryToTree(edm::InputTag rechitHBHETag,
      edm::InputTag rechitHFTag,
      edm::InputTag noiseSummaryTag,
      const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   using namespace edm;

   Handle<HBHERecHitCollection> hHBHERecHits;
   iEvent.getByLabel(rechitHBHETag, hHBHERecHits);
   
   Handle<HFRecHitCollection> hHFRecHits;
   iEvent.getByLabel(rechitHFTag, hHFRecHits);

   Handle<HcalNoiseSummary> hSummary;
   iEvent.getByLabel(noiseSummaryTag, hSummary);

   int NumberOfIsolatedHBHENoise = 0;
   double SumEOfIsolatedHBHENoise = 0;
   int NumberOfFlatHBHENoise = 0;
   double SumEOfFlatHBHENoise = 0;
   int NumberOfSpikeHBHENoise = 0;
   double SumEOfSpikeHBHENoise = 0;
   int NumberOfTriangleHBHENoise = 0;
   double SumEOfTriangleHBHENoise = 0;
   int NumberOfTimingHBHENoise = 0;
   double SumEOfTimingHBHENoise = 0;
   int NumberOfHBHENoise = 0;
   double SumEOfHBHENoise = 0;

   for(int i = 0; i < (int)hHBHERecHits->size(); i++)
   {
      bool IsNoise = false;
      double Energy = (*hHBHERecHits)[i].energy();

      if((*hHBHERecHits)[i].flagField(8) > 0)   // Out-of-time noise
      {
         IsNoise = true;
         NumberOfTimingHBHENoise = NumberOfTimingHBHENoise + 1;
         SumEOfTimingHBHENoise = SumEOfTimingHBHENoise + Energy;
      }
      if((*hHBHERecHits)[i].flagField(11) > 0)   // isolated noise
      {
         IsNoise = true;
         NumberOfIsolatedHBHENoise = NumberOfIsolatedHBHENoise + 1;
         SumEOfIsolatedHBHENoise = SumEOfIsolatedHBHENoise + Energy;
      }
      if((*hHBHERecHits)[i].flagField(12) > 0)   // flat noise
      {
         IsNoise = true;
         NumberOfFlatHBHENoise = NumberOfFlatHBHENoise + 1;
         SumEOfFlatHBHENoise = SumEOfFlatHBHENoise + Energy;
      }
      if((*hHBHERecHits)[i].flagField(13) > 0)   // spike-like noise
      {
         IsNoise = true;
         NumberOfSpikeHBHENoise = NumberOfSpikeHBHENoise + 1;
         SumEOfSpikeHBHENoise = SumEOfSpikeHBHENoise + Energy;
      }
      if((*hHBHERecHits)[i].flagField(14) > 0)   // triangle noise
      {
         IsNoise = true;
         NumberOfTriangleHBHENoise = NumberOfTriangleHBHENoise + 1;
         SumEOfTriangleHBHENoise = SumEOfTriangleHBHENoise + Energy;
      }
      if(IsNoise == true)
      {
         NumberOfHBHENoise = NumberOfHBHENoise + 1;
         SumEOfHBHENoise = SumEOfHBHENoise + Energy;
      }
   }

   int NumberOfHFPMTHits = 0;
   double SumEOfHFPMTHits = 0;

   for(int i = 0; i < (int)hHFRecHits->size(); i++)
   {
      double Energy = (*hHFRecHits)[i].energy();

      if((*hHFRecHits)[i].flagField(HcalCaloFlagLabels::HFLongShort) > 0)
      {
         NumberOfHFPMTHits = NumberOfHFPMTHits + 1;
         SumEOfHFPMTHits = SumEOfHFPMTHits + Energy;
      }
   }

   bool FailHPDHits = false;
   bool FailHPDNoOtherHits = false;
   bool FailMaxZeros = false;
   bool FailE2E10 = false;
   bool FailICHEPFilter = false;

   if(hSummary->maxHPDHits() >= 17)
      FailHPDHits = true;
   if(hSummary->maxHPDNoOtherHits() >= 10)
      FailHPDNoOtherHits = true;
   if(hSummary->maxZeros() >= 10)
      FailMaxZeros = true;
   if(hSummary->minE2Over10TS() < 0.70 || hSummary->maxE2Over10TS() > 0.96)
      FailE2E10 = true;
   FailICHEPFilter = FailHPDHits || FailHPDNoOtherHits || FailMaxZeros || FailE2E10;

   cmstree->column("nHBHENoise", NumberOfHBHENoise, int(0), "HcalNoise");
   cmstree->column("sumEHBHENoise",SumEOfHBHENoise, double(0), "HcalNoise");
   cmstree->column("nTimingHBHENoise", NumberOfTimingHBHENoise, int(0), "HcalNoise");
   cmstree->column("sumETimingHBHENoise", SumEOfTimingHBHENoise, double(0), "HcalNoise");
   cmstree->column("nIsolatedHBHENoise", NumberOfIsolatedHBHENoise, int(0), "HcalNoise");
   cmstree->column("sumEIsolatedHBHENoise", SumEOfIsolatedHBHENoise, double(0), "HcalNoise");
   cmstree->column("nFlatHBHENoise", NumberOfFlatHBHENoise, int(0), "HcalNoise");
   cmstree->column("sumEFlatHBHENoise", SumEOfFlatHBHENoise, double(0), "HcalNoise");
   cmstree->column("nSpikeHBHENoise", NumberOfSpikeHBHENoise, int(0), "HcalNoise");
   cmstree->column("sumESpikeHBHENoise", SumEOfSpikeHBHENoise, double(0), "HcalNoise");
   cmstree->column("nTriangleHBHENoise", NumberOfTriangleHBHENoise, int(0), "HcalNoise");
   cmstree->column("sumETriangleHBHENoise", SumEOfTriangleHBHENoise, double(0), "HcalNoise");
   cmstree->column("nHFPMT", NumberOfHFPMTHits, int(0), "HcalNoise");
   cmstree->column("sumEHFPMT",SumEOfHFPMTHits, double(0), "HcalNoise");
   cmstree->column("failHPDHits", FailHPDHits, false, "HcalNoise");
   cmstree->column("failHPDNoOtherHits", FailHPDNoOtherHits, false, "HcalNoise");
   cmstree->column("failMaxZeros", FailMaxZeros, false, "HcalNoise");
   cmstree->column("failE2E10", FailE2E10, false, "HcalNoise");
   cmstree->column("failICHEPFilter", FailICHEPFilter, false, "HcalNoise");
}










