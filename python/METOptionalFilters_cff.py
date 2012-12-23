import MyAnalysis.METFlags.EcalDeadCellEventFlagProducer_cfi
EcalDeadCellEventFilterSkim = MyAnalysis.METFlags.EcalDeadCellEventFlagProducer_cfi.EcalDeadCellEventFlagProducer.clone()
EcalDeadCellEventFilterSkim.taggingMode = False

import MyAnalysis.METFlags.simpleDRFlagProducer_cfi
simpleDRFilterSkim = MyAnalysis.METFlags.simpleDRFlagProducer_cfi.simpleDRFlagProducer.clone()
simpleDRFilterSkim.taggingMode = False

# the CSC Halo filter only works in the RECO
#from MyAnalysis.METFlags.CSCHaloFlagProducer_cfi import *

import RecoMET.METFilters.trackingFailureFilter_cfi
trackingFailureFilterSkim = RecoMET.METFilters.trackingFailureFilter_cfi.trackingFailureFilter.clone()
trackingFailureFilterSkim.taggingMode = False
trackingFailureFilterSkim.VertexSource = "goodPrimaryVertices"

import PhysicsTools.EcalAnomalousEventFilter.ecalanomalouseventfilter_cfi
EcalAnomalousEventFilterSkim = PhysicsTools.EcalAnomalousEventFilter.ecalanomalouseventfilter_cfi.EcalAnomalousEventFilter.clone()
EcalAnomalousEventFilterSkim.recHitsEB = 'reducedEcalRecHitsEB'
EcalAnomalousEventFilterSkim.recHitsEE = 'reducedEcalRecHitsEE'

from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *


# in 42X, it has to be customized to use run/lumi/event filter, otherwise it needs AOD
import RecoMET.METFilters.hcalLaserEventFilter_cfi
hcalLaserEventFilterSkim = RecoMET.METFilters.hcalLaserEventFilter_cfi.hcalLaserEventFilter.clone()
hcalLaserEventFilterSkim.taggingMode = False
hcalLaserEventFilterSkim.vetoByRunEventNumber = False
hcalLaserEventFilterSkim.vetoByHBHEOccupancy= True

import RecoMET.METFilters.eeBadScFilter_cfi
eeBadScFilterSkim = RecoMET.METFilters.eeBadScFilter_cfi.eeBadScFilter.clone()
eeBadScFilterSkim.taggingMode = False

import MyAnalysis.METFlags.logErrorAnalysisProducer_cff
tooManySeedsFilterSkim = MyAnalysis.METFlags.logErrorAnalysisProducer_cff.tooManySeeds.clone()
tooManyClustersFilterSkim = MyAnalysis.METFlags.logErrorAnalysisProducer_cff.tooManyClusters.clone()
tooManyTripletsPairsFilterSkim = MyAnalysis.METFlags.logErrorAnalysisProducer_cff.tooManyTripletsPairs.clone()
tooManyTripletsPairsMainIterationsFilterSkim = MyAnalysis.METFlags.logErrorAnalysisProducer_cff.tooManyTripletsPairsMainIterations.clone()
tooManySeedsMainIterationsFilterSkim = MyAnalysis.METFlags.logErrorAnalysisProducer_cff.tooManySeedsMainIterations.clone()


metOptionalFiltersSkim = cms.Sequence( EcalDeadCellEventFilterSkim
                                       * HBHENoiseFilter
                                       * hcalLaserEventFilterSkim
                                       * eeBadScFilterSkim
                                       #* simpleDRFilterSkim  # this is super-efficient!
                                       * trackingFailureFilterSkim
                                       * EcalAnomalousEventFilterSkim
                                       * tooManySeedsFilterSkim
                                       * tooManyClustersFilterSkim
                                       * tooManyTripletsPairsFilterSkim
                                       * tooManyTripletsPairsMainIterationsFilterSkim
                                       * tooManySeedsMainIterationsFilterSkim
                                       )

