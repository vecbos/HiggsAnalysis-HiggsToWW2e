import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *
#ak5JPTJetsL2L3   = cms.EDProducer('JPTJetCorrectionProducer',
#    src         = cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
#    correctors  = cms.vstring('ak5JPTL2L3')
#    )

from HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff import *

CaloJetSequence = cms.Sequence( ak5CaloJetsL2L3 )
PFJetAK5Sequence = cms.Sequence( ak5PFJetsL2L3 )
#JPTjetsAK5Sequence = cms.Sequence( ak5JPTJetsL2L3 )

#ourJetSequence = cms.Sequence( CaloJetSequence * PFJetAK5Sequence * JPTjetsAK5Sequence )
ourJetSequence = cms.Sequence( CaloJetSequence * PFJetAK5Sequence )
