import FWCore.ParameterSet.Config as cms
#from JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_cff import *

from JetMETCorrections.Configuration.DefaultJEC_cff import *

from HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff import *

CaloJetSequence = cms.Sequence( ak5CaloJetsL2L3 )
PFJetAK5Sequence = cms.Sequence( ak5PFJetsL2L3 )
JPTjetsAK5Sequence = cms.Sequence( ZSPJetCorrectionsAntiKt5 * ZSPrecoJetAssociationsAntiKt5 * ak5JPTJets * ak5JPTJetsL2L3 )

ourJetSequence = cms.Sequence( CaloJetSequence * PFJetAK5Sequence * JPTjetsAK5Sequence )
