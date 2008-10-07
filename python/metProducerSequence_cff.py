import FWCore.ParameterSet.Config as cms

from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.Configuration.RecoGenMET_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import *

pfMET = cms.EDProducer("PFMET",
                       PFCandidates = cms.InputTag("particleFlow"),
                       verbose = cms.untracked.bool(False)
                       )

metSequence = cms.Sequence ( pfMET )
