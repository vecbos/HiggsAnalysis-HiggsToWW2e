import FWCore.ParameterSet.Config as cms

mergedSuperClusters = cms.EDProducer("SuperClusterCombiner",
                                     labels = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                                                            cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"))
                                     )

electronAndPhotonSuperClusters = cms.EDProducer("ElectronAndPhotonSuperClusterProducer",
                                                ElectronLabel = cms.InputTag("gsfElectrons"),
                                                PhotonLabel = cms.InputTag("photons"),
                                                SuperClusterLabel = cms.InputTag("mergedSuperClusters")
                                                )

electronAndPhotonSuperClustersSequence = cms.Sequence( mergedSuperClusters
                                                       * electronAndPhotonSuperClusters )
