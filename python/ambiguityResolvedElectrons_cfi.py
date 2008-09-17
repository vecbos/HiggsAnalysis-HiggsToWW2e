import FWCore.ParameterSet.Config as cms

ambiguityResolvedElectrons = cms.EDFilter("AmbResolver",
                                          src = cms.InputTag("pixelMatchGsfElectrons"),
                                          filter = cms.bool(False),
                                          reducedElectronsRefCollectionLabel = cms.InputTag("selectedElectronsRef"),
                                          doRefCheck = cms.bool(True)
                                          )
