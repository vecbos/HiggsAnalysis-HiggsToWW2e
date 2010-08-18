import FWCore.ParameterSet.Config as cms

# ESSource from file (until JPT corrections are not in DB)
JetCorrectionEra = cms.PSet( era = cms.string('Spring10') )
ak5JPTL2Relative = cms.ESSource(
    'LXXXCorrectionService',
    JetCorrectionEra,
    section   = cms.string(''),
    level     = cms.string('L2Relative'),
    algorithm = cms.string('AK5JPT')
    )

ak5JPTL3Absolute = cms.ESSource(
    'LXXXCorrectionService',
    JetCorrectionEra,
    section   = cms.string(''),
    level     = cms.string('L3Absolute'),
    algorithm = cms.string('AK5JPT')
    )

ak5JPTL2L3 = cms.ESSource(
    'JetCorrectionServiceChain',
    correctors = cms.vstring('ak5JPTL2Relative','ak5JPTL3Absolute')
    )

ak5JPTJetsL2L3   = cms.EDProducer('JPTJetCorrectionProducer',
                                  src         = cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
                                  correctors  = cms.vstring('ak5JPTL2L3')
                                  )
