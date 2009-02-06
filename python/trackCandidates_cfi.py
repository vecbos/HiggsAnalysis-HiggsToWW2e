import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

trackCandidates = cms.EDProducer("ConcreteChargedCandidateProducer",
                                 src = cms.InputTag("generalTracks"),
                                 particleType = cms.string('mu-')
                             )

