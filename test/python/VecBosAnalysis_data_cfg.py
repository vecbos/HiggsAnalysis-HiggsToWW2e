import FWCore.ParameterSet.Config as cms

is42X = True
is52X = False
runOnAOD = 1

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
if(is42X): process.load("Configuration.StandardSequences.Reconstruction_cff")
else:      process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(is42X):
    process.GlobalTag.globaltag = 'GR_R_42_V25::All'
elif(is52X):
    process.GlobalTag.globaltag = 'GR_R_52_V9D::All'
else:
    process.GlobalTag.globaltag = 'GR_P_V41_AN1::All'

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
if(is42X):
    process.hcalLaserEventFilter.vetoByRunEventNumber = True
    process.hcalLaserEventFilter.vetoByHBHEOccupancy= False

# --- tracker failures ---
process.load("MyAnalysis.METFlags.logErrorAnalysisProducer_cff")

process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonLinkedTracks_cfi")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- pf isolation sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonPFIsoSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.electronAndPhotonSuperClusters_cff")
process.load("HiggsAnalysis.HiggsToWW2e.seedBasicClusters_cff")

# --- good vertex filter ---
process.load("HiggsAnalysis.HiggsToWW2e.vertexFiltering_cff")

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
process.treeDumper.dumpTracks = False
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = True
process.treeDumper.dumpConversions = True
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

# this replaces the collections with existing ones (even if they are saved already)
if(is42X==True):
    process.treeDumper.hpsTauDiscrByVLooseCombinedIsolationDBSumPtCorrTag = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation")
    process.treeDumper.hpsTauDiscrByLooseCombinedIsolationDBSumPtCorrTag = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")
    process.treeDumper.hpsTauDiscrByMediumCombinedIsolationDBSumPtCorrTag = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation")
    process.treeDumper.hpsTauDiscrByTightCombinedIsolationDBSumPtCorrTag = cms.InputTag("hpsPFTauDiscriminationByTightIsolation")   

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

if(is42X):
    inputfile = cms.untracked.vstring('file:/cmsrm/pc24_2/emanuele/data/reRecoMay10File.root')
elif(is52X):
    inputfile = cms.untracked.vstring('/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/191/700/00327508-AF8B-E111-8151-BCAEC53296F4.root')
else:
    inputfile = cms.untracked.vstring('/store/data/Run2012C/DoubleElectron/AOD/PromptReco-v2/000/201/278/F8F44E0A-ACED-E111-8621-E0CB4E553651.root')

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = inputfile
                            )

process.prejets = cms.Sequence( process.leptonLinkedTracks
                                * process.electronAndPhotonSuperClustersSequence
                                * process.chargedMetProducer
                                * process.metOptionalFilterSequence
                                * process.pfIsolationAllSequence )

if(process.treeDumper.dumpBCs==True):
    process.prejets *= process.seedBasicClustersSequence

process.jets = cms.Sequence( process.ourJetSequenceDataReduced
                             * process.newBtaggingSequence 
                             * process.newPFPUcorrJetBtaggingSequence
                             * process.newPFNoPUJetBtaggingSequence
                             * process.metSequence )

process.postjets = cms.Sequence( process.eIdSequence
                                 * process.FastjetForIsolation
                                 * process.logErrorAnalysis
                                 * process.treeDumper )

# In order to use the good primary vertices everywhere (It would be nicer to set the proper inputTags in the first place)
from PhysicsTools.PatAlgos.tools.helpers import *
massSearchReplaceAnyInputTag(process.prejets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.jets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.postjets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)

process.p = cms.Path ( process.goodPrimaryVertices
                       * process.prejets
                       * process.jets
                       * process.postjets
                       )
