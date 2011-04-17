import FWCore.ParameterSet.Config as cms

runOnAOD = 1
useL1Offset = 0 # 1=L1Offset with vtx correction 0=FastJet

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_311_V2::All'

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFPUcorrJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.dinamicAnnealingVertexing_cff")
#process.load("HiggsAnalysis.HiggsToWW2e.btagJPTJetsProducerSequence_cff")

if (useL1Offset == 1) :
    process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
    process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL2L3Residual'
    process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL2L3Residual'
    process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL2L3Residual'
    process.newPFJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL2L3Residual'
    process.newPFJetsSoftElectronTagInfos.jets = 'ak5PFJetsL2L3Residual'
    process.newPFJetsSoftMuonTagInfos.jets = 'ak5PFJetsL2L3Residual'
    process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1L2L3Residual'
    process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1L2L3Residual'
    process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1L2L3Residual'
else:
    process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")
    process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL2L3Residual'
    process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL2L3Residual'
    process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL2L3Residual'
    process.newPFJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL2L3Residual'
    process.newPFJetsSoftElectronTagInfos.jets = 'ak5PFJetsL2L3Residual'
    process.newPFJetsSoftMuonTagInfos.jets = 'ak5PFJetsL2L3Residual'
    process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1FastL2L3Residual'
    process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1FastL2L3Residual'
    process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1FastL2L3Residual'

# to correct calo met ---
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
#process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
#process.metMuonJESCorAK5.inputUncorMetLabel = "met"
#process.metMuonJESCorAK5.useTypeII = True
#process.metMuonJESCorAK5.hasMuonsCorr = False

# --- electron sequences ---
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
process.treeDumper.nameFile = 'default_data.root'
if (useL1Offset == 1) :
    process.treeDumper.PFpuCorrJetCollection1 = 'ak5PFJetsL1L2L3Residual'
    process.treeDumper.PFpuCorrJetCollection2 = 'ak5PFJetsL1L2L3Residual'
else:
    process.treeDumper.PFpuCorrJetCollection1 = 'ak5PFJetsL1FastL2L3Residual'
    process.treeDumper.PFpuCorrJetCollection2 = 'ak5PFJetsL1FastL2L3Residual'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = False
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.dumpTree = True
if (runOnAOD == 1) :
    process.treeDumper.saveFatTrk = False
    process.treeDumper.saveTrackDeDx = False
    process.treeDumper.dumpPFlowElectrons = False
    process.treeDumper.dumpHcalNoiseFlags = True
    process.treeDumper.AODHcalNoiseFlags = True
else :
    process.treeDumper.saveFatTrk = True
    process.treeDumper.saveTrackDeDx = True
    process.treeDumper.dumpPFlowElectrons = True
    process.treeDumper.dumpHcalNoiseFlags = True
    process.treeDumper.AODHcalNoiseFlags = False

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.lumiAna = cms.EDAnalyzer("LumiAnalyzer")
process.TFileService = cms.Service("TFileService", fileName = cms.string("lumi.root") )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23_2/emanuele/data/AOD_Run2011A.root')
                            )

process.p = cms.Path ( process.mergedSuperClusters * process.mergedBasicClusters *
                       process.ourJetSequenceData *
                       process.offlinePrimaryVertices *
                       process.newBtaggingSequence * process.newPFJetBtaggingSequence * process.newPFPUcorrJetBtaggingSequence *
                       process.eIdSequence * process.FastjetForIsolation  *
                       process.ambiguityResolvedElectrons *
                       process.lumiAna # save lumi info by LS
                       )

process.q = cms.EndPath ( process.treeDumper )
