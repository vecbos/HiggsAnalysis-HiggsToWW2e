import FWCore.ParameterSet.Config as cms

# jet id mva 
from CMGTools.External.puJetIDAlgo_cff import PhilV1
mvaJetIDMapProd = cms.EDProducer("mvaJetIDMapProd",
                                 puJetIDAlgo = PhilV1,
                                 vtxLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
                                 uncorrJetLabel = cms.untracked.InputTag("ak5PFJets"),
                                 corrJetLabel = cms.untracked.InputTag("ak5PFJetsL1FastL2L3")
                                 )                                      

