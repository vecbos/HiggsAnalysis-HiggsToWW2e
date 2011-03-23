import FWCore.ParameterSet.Config as cms

runOnAOD = 1
useL1Offset = 0 # 1=L1Offset with vtx correction 0=FastJet

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START311_V2::All'
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFPUcorrJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.dinamicAnnealingVertexing_cff")
#process.load("HiggsAnalysis.HiggsToWW2e.btagJPTJetsProducerSequence_cff")

# do not use residual corrections in MC
if (useL1Offset == 1) :
    process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
    process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL2L3'
    process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL2L3'
    process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL2L3'
    process.newPFJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL2L3'
    process.newPFJetsSoftElectronTagInfos.jets = 'ak5PFJetsL2L3'
    process.newPFJetsSoftMuonTagInfos.jets = 'ak5PFJetsL2L3'
    process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1L2L3'
    process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1L2L3'
    process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1L2L3'
else:
    process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")
    process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL2L3'
    process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL2L3'
    process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL2L3'
    process.newPFJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL2L3'
    process.newPFJetsSoftElectronTagInfos.jets = 'ak5PFJetsL2L3'
    process.newPFJetsSoftMuonTagInfos.jets = 'ak5PFJetsL2L3'
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

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- to recover the ak5 GenJets that are not re-recoed in 33X samples ---
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load("RecoJets.Configuration.RecoGenJets_cff")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_MC.root'
process.treeDumper.jetCollection1 = 'ak5CaloJetsL2L3'
process.treeDumper.JPTjetCollection1 = 'ak5JPTJetsL2L3'
process.treeDumper.PFjetCollection1 = 'ak5PFJetsL2L3'
if (useL1Offset == 1) :
    process.treeDumper.PFpuCorrJetCollection1 = 'ak5PFJetsL1L2L3'
    process.treeDumper.PFpuCorrJetCollection2 = 'ak5PFJetsL1L2L3'
else:
    process.treeDumper.PFpuCorrJetCollection1 = 'ak5PFJetsL1FastL2L3'
    process.treeDumper.PFpuCorrJetCollection2 = 'ak5PFJetsL1FastL2L3'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpGenInfo = True
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
process.treeDumper.dumpTree = True
if (runOnAOD == 1) :
    process.treeDumper.saveFatTrk = False
    process.treeDumper.saveTrackDeDx = False
    process.treeDumper.dumpPFlowElectrons = False
else :
    process.treeDumper.saveFatTrk = True
    process.treeDumper.saveTrackDeDx = True
    process.treeDumper.dumpPFlowElectrons = True
    process.treeDumper.dumpHcalNoiseFlags = True

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23/emanuele/data/Pool/jpsiEE_Fall10.root') # RECO
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23_2/emanuele/Pool/AODSIM_Winter10_FlatPU.root')
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23_2/emanuele/data/AOD_HWW_Spring11.root')
                            )

process.p = cms.Path ( process.mergedBasicClusters * process.mergedSuperClusters *
                       process.genParticlesForJets * process.ak5GenJets * # added for re-recoed V9 Summer09 samples where the ak5GenJet collection was dropped
                       process.offlinePrimaryVertices *
                       process.ourJetSequenceMC *
                       process.newBtaggingSequence * process.newPFJetBtaggingSequence * process.newPFPUcorrJetBtaggingSequence *
                       process.eIdSequence * process.FastjetForIsolation *
                       process.ambiguityResolvedElectrons )

process.q = cms.EndPath ( process.treeDumper )
