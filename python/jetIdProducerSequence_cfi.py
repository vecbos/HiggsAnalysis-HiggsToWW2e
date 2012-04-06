import FWCore.ParameterSet.Config as cms

# jet id mva 
from CMGTools.External.puJetIDAlgo_cff import PhilV0
mvaJetIDMapProd = cms.EDProducer("mvaJetIDMapProd",
                                 puJetIDAlgo = PhilV0,
                                 vtxLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
                                 uncorrJetLabel = cms.untracked.InputTag("ak5PFJets"),
                                 corrJetLabel = cms.untracked.InputTag("ak5PFJetsL1FastL2L3")
                                 )                                      

