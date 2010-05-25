import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.DefaultJEC_cff import *

from HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff import *

jetVertexAlpha1 = cms.EDProducer("JetVertexAssociation",
                                 JV_deltaZ = cms.double(0.3),
                                 JV_sigmaZ = cms.double(9.5),
                                 JV_alpha_threshold = cms.double(0.2),
                                 JV_cone_size = cms.double(0.5),
                                 JV_type_Algo = cms.int32(1),
                                 JET_ALGO = cms.string('ak5CaloJetsL2L3'),
                                 TRACK_ALGO = cms.string('generalTracks'),
                                 VERTEX_ALGO = cms.string('offlinePrimaryVertices'),
                                 JV_cutType = cms.string('delta') )

jetVertexAlpha2 = cms.EDProducer("JetVertexAssociation",
                                 JV_deltaZ = cms.double(0.3),
                                 JV_sigmaZ = cms.double(9.5),
                                 JV_alpha_threshold = cms.double(0.2),
                                 JV_cone_size = cms.double(0.5),
                                 JV_type_Algo = cms.int32(1),
                                 JET_ALGO = cms.string('ak5CaloJets'),
                                 TRACK_ALGO = cms.string('generalTracks'),
                                 VERTEX_ALGO = cms.string('offlinePrimaryVertices'),
                                 JV_cutType = cms.string('delta') )

jetSequence = cms.Sequence( ak5CaloJetsL2L3 + jetVertexAlpha1 + jetVertexAlpha2 )
pfjetAK5Sequence = cms.Sequence( ak5PFJetsL2L3 )

ourJetSequence = cms.Sequence( jetSequence * pfjetAK5Sequence )
