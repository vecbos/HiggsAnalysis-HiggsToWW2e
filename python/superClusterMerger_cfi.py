import FWCore.ParameterSet.Config as cms

mergedSuperClusters = cms.EDProducer("SuperClusterCombiner",
                                     labels = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                                                            cms.InputTag("multi5x5SuperClustersWithPreshower"))
                                     )
