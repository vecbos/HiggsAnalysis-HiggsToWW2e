import FWCore.ParameterSet.Config as cms

# produces  a collection with all the charged candidates to be used for chMET
reducedPFCands = cms.EDProducer("ReducedCandidatesProducer",
                                srcCands = cms.InputTag("particleFlow"),
                                srcVertices = cms.InputTag("offlinePrimaryVertices"),
                                dz = cms.double(0.1)
                                )

# produces the reduced collection, subset of the above, of all the candidates in a cone 0.1 from reco leptons (to be removed offline from the chMET)
reducedPFCandsToSave = cms.EDFilter("leptonCandidateFilter",
                                    src = cms.InputTag("reducedPFCands"),
                                    ElectronLabel = cms.InputTag("gsfElectrons"),
                                    MuonLabel = cms.InputTag("muons")
                                    )

ourChPFMet = cms.EDProducer("ChargedPFMetProducer",
                            collectionTag = cms.InputTag("reducedPFCands")
                            )

metSequence = cms.Sequence ( reducedPFCands * ourChPFMet * reducedPFCandsToSave )
