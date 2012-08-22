import FWCore.ParameterSet.Config as cms

is42X = False
runOnAOD = 1

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if(is42X):
    process.GlobalTag.globaltag = 'START42_V14B::All'
else:
    process.GlobalTag.globaltag = 'START52_V9::All'
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

# process.IgProfService = cms.Service("IgProfService",
#                                     reportFirstEvent            = cms.untracked.int32(0),
#                                     reportEventInterval         = cms.untracked.int32(50),
#                                     reportToFileAtPostEvent     = cms.untracked.string("| gzip -c > YYYY.%I.gz")
#                                     )

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonLinkedTracks_cfi")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- pf isolation sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonPFIsoSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection and seeds filtering ---
process.load("HiggsAnalysis.HiggsToWW2e.electronAndPhotonSuperClusters_cff")
process.load("HiggsAnalysis.HiggsToWW2e.seedBasicClusters_cff")

# --- good vertex filter ---
process.load("HiggsAnalysis.HiggsToWW2e.vertexFiltering_cff")

#PDF systematics
# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
                                    # Fix POWHEG if buggy (this PDF set will also appear on output,
                                    # so only two more PDF sets can be added in PdfSetNames if not "")
                                    #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
                                    GenTag = cms.untracked.InputTag("genParticles"),
                                    PdfInfoTag = cms.untracked.InputTag("generator"),
                                    PdfSetNames = cms.untracked.vstring("cteq66.LHgrid", "MRST2006nnlo.LHgrid", "NNPDF10_100.LHgrid")
                                    )

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_MC.root'
process.treeDumper.PFJetCorrectionService = 'ak5PFL1FastL2L3'
process.treeDumper.JetCorrectionService = 'ak5CaloL1FastL2L3'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpLHE = False
process.treeDumper.dumpPdfWeight = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpTracks = False
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = True
process.treeDumper.dumpConversions = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpGenJets = True
#process.treeDumper.dumpParticleFlowObjects = False
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

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc21_2/emanuele/data/AOD_WWSummer12.root')
                            )

process.prejets = cms.Sequence( process.leptonLinkedTracks 
                                * process.electronAndPhotonSuperClustersSequence
                                * process.chargedMetProducer
                                * process.pfIsolationAllSequence )

if(process.treeDumper.dumpBCs==True):
    process.prejets *= process.seedBasicClustersSequence 

process.jets = cms.Sequence( process.ourJetSequenceMCReduced
                             * process.newBtaggingSequence 
                             * process.newPFPUcorrJetBtaggingSequence
                             * process.newPFNoPUJetBtaggingSequence
                             * process.metSequence )

process.postjets = cms.Sequence( process.eIdSequence
                                 * process.FastjetForIsolation
                                 * process.treeDumper )

# In order to use the good primary vertices everywhere (It would be nicer to set the proper inputTags in the first place)
from PhysicsTools.PatAlgos.tools.helpers import *
massSearchReplaceAnyInputTag(process.prejets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.jets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.postjets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)

if(process.treeDumper.dumpPdfWeight == False) :
    process.p = cms.Path ( process.goodPrimaryVertices
                           * process.prejets
                           * process.jets
                           * process.postjets
                           )
else :
    process.p = cms.Path ( process.pdfWeights
                           * process.goodPrimaryVertices
                           * process.prejets
                           * process.jets
                           * process.postjets
                           )
