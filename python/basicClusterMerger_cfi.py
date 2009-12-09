import FWCore.ParameterSet.Config as cms

mergedBasicClusters = cms.EDProducer("BasicClusterCombiner",
                                     labels = cms.VInputTag(cms.InputTag("islandBasicClusters","islandBarrelBasicClusters"),
                                                            cms.InputTag("islandBasicClusters","islandEndcapBasicClusters"))
                                     )
