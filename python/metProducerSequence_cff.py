import FWCore.ParameterSet.Config as cms

# produces  a collection with all the charged candidates to be used for chMET
reducedPFCands = cms.EDProducer("HWWReducedCandidatesProducer",
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

from JetMETCorrections.Type1MET.pfMETCorrections_cff import *
# this does not exist in the release
pfType0CorrectedMet = pfType1CorrectedMet.clone()
pfType0CorrectedMet.applyType1Corrections = False

# merge the PF met and the corrected ones into a unique collection of candidates
pfmets = cms.EDProducer("PfMetCombiner",
                        labels = cms.VInputTag(cms.InputTag("pfMet"),
                                               cms.InputTag("pfType0CorrectedMet"),
                                               cms.InputTag("pfType1CorrectedMet"),
                                               cms.InputTag("pfType1p2CorrectedMet")
                                               )
                        )

metSequence = cms.Sequence ( reducedPFCands * ourChPFMet * reducedPFCandsToSave * pfCandsNotInJet * pfJetMETcorr * pfCandMETcorr * pfchsMETcorr * pfType0CorrectedMet * pfType1CorrectedMet * pfType1p2CorrectedMet * pfmets)
