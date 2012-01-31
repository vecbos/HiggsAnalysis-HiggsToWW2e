import FWCore.ParameterSet.Config as cms

runOnAOD = 1

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V21::All'

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFPUcorrJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFNoPUJetsProducerSequence_cff")
#process.load("HiggsAnalysis.HiggsToWW2e.btagJPTJetsProducerSequence_cff")
process.load("WWAnalysis.Tools.chargedMetProducer_cfi")
process.chargedMetProducer.collectionTag = "particleFlow"
process.chargedMetProducer.vertexTag = "offlinePrimaryVertices"

# --- noise filters ---
process.load("HiggsAnalysis.HiggsToWW2e.METOptionalFilterFlags_cff")

# --- tracker failures ---
process.load("MyAnalysis.METFlags.logErrorAnalysisProducer_cff")

process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")
process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL2L3Residual'
process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL2L3Residual'
process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL2L3Residual'
process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1FastL2L3Residual'
process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1FastL2L3Residual'
process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1FastL2L3Residual'
#process.newPFJetNoPUTracksAssociatorAtVertex.jets = 'ak5PFJetsNoPUL1FastL2L3Residual'
#process.newPFJetsNoPUSoftElectronTagInfos.jets = 'ak5PFJetsNoPUL1FastL2L3Residual'
#process.newPFJetsNoPUSoftMuonTagInfos.jets = 'ak5PFJetsNoPUL1FastL2L3Residual'
#process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL1FastL2L3'
#process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL1FastL2L3'
#process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL1FastL2L3'

# to correct calo met ---
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
#process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
#process.metMuonJESCorAK5.inputUncorMetLabel = "met"
#process.metMuonJESCorAK5.useTypeII = True
#process.metMuonJESCorAK5.hasMuonsCorr = False

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonLinkedTracks_cfi")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- pf isolation sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonPFIsoSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_data.root'
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
process.treeDumper.dumpPFCandidates = False
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

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc24_2/emanuele/data/reRecoMay10File.root')
                            )

process.p = cms.Path ( process.leptonLinkedTracks
                       * process.mergedSuperClusters 
                       * process.chargedMetProducer
                       * process.metOptionalFilterSequence
                       * process.metSequence
                       * process.pfIsolationAllSequence
                       * process.ourJetSequenceDataReduced
                       * process.newBtaggingSequence * process.newPFPUcorrJetBtaggingSequence
                       * process.eIdSequence
                       * process.FastjetForIsolation
                       * process.logErrorAnalysis
                       * process.treeDumper
                       )

