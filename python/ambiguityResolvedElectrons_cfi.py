import FWCore.ParameterSet.Config as cms

ambiguityResolvedElectrons = cms.EDFilter("AmbResolver",
                                          src = cms.InputTag("pixelMatchGsfElectrons"),
                                          filter = cms.bool(False),
                                          reducedElectronsRefCollectionLabel = cms.InputTag("isolatedElectronsRef"),
                                          doRefCheck = cms.bool(True)
                                          )
