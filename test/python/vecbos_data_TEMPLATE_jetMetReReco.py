import FWCore.ParameterSet.Config as cms

process = cms.Process("VecBosAnalysis")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR10_P_V7::All'
process.GlobalTag.globaltag = 'GR_R_36X_V12B::All'
#process.GlobalTag.globaltag = 'GR_R_36X_V12A::All'

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'vecbos_CHANGERUNNUMBER.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = True
process.treeDumper.dumpTree = True

process.options = cms.untracked.PSet(
    fileMode =  cms.untracked.string('NOMERGE'),
    wantSummary = cms.untracked.bool(True) 
    )

    #####################################################################################################
    ####
    ####  Top level replaces for handling strange scenarios of early collisions
    ####

    ## TRACKING:
    ## Skip events with HV off
process.newSeedFromTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.newSeedFromPairs.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
process.secTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.fifthSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 20000
process.fourthPLSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters= 20000
process.thTripletsA.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000
process.thTripletsB.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000

## local tracker strip reconstruction
process.OutOfTime.TOBlateBP=0.071
process.OutOfTime.TIBlateBP=0.036

###### FIXES TRIPLETS FOR LARGE BS DISPLACEMENT ######

### prevent bias in pixel vertex
process.pixelVertices.useBeamConstraint = False

### pixelTracks
#---- new parameters ----
process.pixelTracks.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.06)
process.pixelTracks.RegionFactoryPSet.RegionPSet.originHalfLength = cms.double(40.6)

### 0th step of iterative tracking
#---- new parameters ----
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ   = cms.double(4.06)
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 40.6

### 2nd step of iterative tracking
#---- new parameters ----
process.secTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.47)
process.secTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 44.7

## Primary Vertex
process.offlinePrimaryVerticesWithBS.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minPixelLayersWithHits = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minSiliconLayersWithHits = 5
process.offlinePrimaryVerticesWithBS.TkClusParameters.TkGapClusParameters.zSeparation = 1
process.offlinePrimaryVertices.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVertices.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVertices.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVertices.TkFilterParameters.minPixelLayersWithHits = 2
process.offlinePrimaryVertices.TkFilterParameters.minSiliconLayersWithHits = 5
process.offlinePrimaryVertices.TkClusParameters.TkGapClusParameters.zSeparation = 1

## ECAL
process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]

##Preshower
process.ecalPreshowerRecHit.ESBaseline = 0

## HCAL temporary fixes
process.hfreco.firstSample  = 3
process.hfreco.samplesToAdd = 4

## EGAMMA
process.photons.minSCEtBarrel = 5.
process.photons.minSCEtEndcap =5.
process.photonCore.minSCEt = 5.
process.conversionTrackCandidates.minSCEt =5.
process.conversions.minSCEt =5.
process.trackerOnlyConversions.rCut = 2.
process.trackerOnlyConversions.vtxChi2 = 0.0005

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.L1BSC=process.hltLevel1GTSeed.clone()
process.L1BSC.L1TechTriggerSeeding = cms.bool(True)
#process.L1BSC.L1SeedsLogicalExpression = cms.string('0 AND (34) AND NOT (36,37,38,39)')
#Only BPTX gating
process.L1BSC.L1SeedsLogicalExpression = cms.string('0')

process.HFCoincidence =  cms.EDFilter("HFCoincidence",
                                      threshold = cms.untracked.double(2.0),
                                      HFRecHitCollection = cms.InputTag("hfreco"),
                                      #                                      maskedChannels = cms.vint32( 8137, 8141, 8146, 8149, 8150, 8153 )
                                      )

process.L1HF=process.hltLevel1GTSeed.clone()
process.L1HF.L1TechTriggerSeeding = cms.bool(True)
process.L1HF.L1SeedsLogicalExpression = cms.string('9')

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15),	
                                           maxd0 = cms.double(2)	
                                           )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
#                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = cms.untracked.vstring('/store/data/Run2010A/EG/RECO/v4/000/140/059/50AEDE30-738E-DF11-8EEA-003048F11C5C.root')
                            )

process.goodElectrons = cms.EDFilter("GsfElectronSelector",
                                         filter = cms.bool(True),
                                         src = cms.InputTag("gsfElectrons::RECO"),
                                         cut = cms.string('et > 10.')
                                     )


process.goodElectronsCounter = cms.EDFilter("CandViewCountFilter",
                                                                                        src = cms.InputTag("goodElectrons"),
                                                                                        minNumber = cms.uint32(1)
                                                                                        )


process.p = cms.Path (
    process.L1BSC *
#    process.HFCoincidence * process.noscraping * process.primaryVertexFilter *
#    process.goodElectrons * process.goodElectronsCounter *
    process.siPixelRecHits * process.siStripMatchedRecHits *
    process.ckftracks * process.vertexreco * process.electronGsfTracking *
    process.particleFlowCluster * process.caloTowersRec * process.jetGlobalReco * process.particleFlowReco *
    process.egammarecoFull *process.jetHighLevelReco *
    process.tautagging * process.metreco* process.btagging * process.recoPFMET *
    process.lowThrCaloTowers * process.mergedSuperClusters * process.mergedBasicClusters *
    process.ourJetSequence * process.newBtaggingSequence *
#    process.gsfElectrons *
    process.eIdSequence *
    process.eleIsolationSequence *
    process.ambiguityResolvedElectrons *
    process.treeDumper
    )
