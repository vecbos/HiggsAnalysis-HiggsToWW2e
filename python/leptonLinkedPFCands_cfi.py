import FWCore.ParameterSet.Config as cms

leptonLinkedChargedPFCands = cms.EDFilter("LeptonChargedPFCandidateFilter",
                                          src = cms.InputTag("particleFlow"),
                                          ElectronLabel = cms.InputTag("gsfElectrons"),
                                          MuonLabel = cms.InputTag("muons")
                                          )
