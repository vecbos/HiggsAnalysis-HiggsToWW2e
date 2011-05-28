import FWCore.ParameterSet.Config as cms

leptonLinkedTracks = cms.EDFilter("leptonTrackFilter",
                                  src = cms.InputTag("generalTracks"),
                                  ElectronLabel = cms.InputTag("gsfElectrons"),
                                  MuonLabel = cms.InputTag("muons"),
                                  JetLabel = cms.InputTag("ak5CaloJets")
                                  )
