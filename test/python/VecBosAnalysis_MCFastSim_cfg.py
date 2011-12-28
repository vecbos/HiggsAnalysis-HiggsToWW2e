import FWCore.ParameterSet.Config as cms

runOnAOD = 1

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START44_V10::All'
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

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

# --- tracker failures ---
process.load("MyAnalysis.METFlags.logErrorAnalysisProducer_cff")

# do not use residual corrections in MC
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")
process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL1FastL2L3'
process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL1FastL2L3'
process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL1FastL2L3'
process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1FastL2L3'
process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1FastL2L3'
process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1FastL2L3'

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
process.treeDumper.nameFile = 'default_MC.root'
process.treeDumper.jetCollection1 = 'ak5CaloJetsL1FastL2L3'
process.treeDumper.jetCollection2 = 'ak5CaloJets::HLT'
process.treeDumper.jetCollection3 = 'ak5CaloJetsL2L3'
process.treeDumper.JPTjetCollection1 = 'ak5JPTJetsL2L3'
process.treeDumper.PFjetCollection1 = 'ak5PFJetsNoPUL1FastL2L3'
process.treeDumper.PFpuCorrJetCollection1 = 'ak5PFJetsL1FastL2L3'
process.treeDumper.PFpuCorrJetCollection3 = 'ak5PFJetsL2L3'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpLHE = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = False
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpGenJets = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.dumpPFCandidates = False
process.treeDumper.dumpTree = True
if (runOnAOD == 1) :
    process.treeDumper.saveFatTrk = False
    process.treeDumper.saveTrackDeDx = False
    process.treeDumper.dumpPFlowElectrons = False
    process.treeDumper.dumpHcalNoiseFlags = False
    process.treeDumper.AODHcalNoiseFlags = False
else :
    process.treeDumper.saveFatTrk = True
    process.treeDumper.saveTrackDeDx = True
    process.treeDumper.dumpPFlowElectrons = True
    process.treeDumper.dumpHcalNoiseFlags = False
    process.treeDumper.AODHcalNoiseFlags = False

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23/emanuele/data/Pool/jpsiEE_Fall10.root') # RECO
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23_2/emanuele/Pool/AODSIM_Winter10_FlatPU.root')
                            fileNames = cms.untracked.vstring('file:batch5_msugra_600_600_10_0_1_razor.root')
                            )

process.p = cms.Path ( process.leptonLinkedTracks
                       * process.mergedSuperClusters
                       * process.chargedMetProducer
                       * process.metSequence
                       * process.pfIsolationAllSequence
                       * process.kt6CaloJets
                       * process.ourJetSequenceMCReduced
                       * process.newBtaggingSequence * process.newPFPUcorrJetBtaggingSequence
                       * process.eIdSequence
                       * process.FastjetForIsolation
                       * process.treeDumper
                       )
