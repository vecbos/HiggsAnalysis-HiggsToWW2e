import FWCore.ParameterSet.Config as cms

process = cms.Process("VecBosAnalysis")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.TFileService = cms.Service("TFileService", fileName = cms.string("lumi.root") )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR10_P_V7::All'
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
#process.treeDumper.saveJet1BTag = False
process.treeDumper.dumpJets = False
process.treeDumper.dumpTree = True
process.treeDumper.dumpHLTObjects = cms.untracked.bool(True)
process.treeDumper.dumpHLTObjectsInfo = cms.untracked.PSet(
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    processName = cms.string("HLT"),
    triggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD","","HLT")
    )

process.options = cms.untracked.PSet(
    fileMode =  cms.untracked.string('NOMERGE'),
    #wantSummary = cms.untracked.bool(True) 
    )

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
#                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = cms.untracked.vstring('file:46040E22-37D1-DF11-8F23-0030487C6062.root')
                            )

process.lumiAna = cms.EDAnalyzer("LumiAnalyzer")

process.p = cms.Path (
     process.lumiAna *
    process.L1BSC *
#    process.HFCoincidence * process.noscraping * process.primaryVertexFilter *
    process.lowThrCaloTowers * process.mergedSuperClusters * process.mergedBasicClusters *
    process.ourJetSequence * #process.newBtaggingSequence *
    process.gsfElectrons *
    process.eIdSequence *
    process.eleIsolationSequence *
    process.ambiguityResolvedElectrons *
    process.treeDumper
    )
