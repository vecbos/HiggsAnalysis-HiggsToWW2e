import FWCore.ParameterSet.Config as cms

runOnAOD = 1

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
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
process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL1FastL2L3'
process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL1FastL2L3'
process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL1FastL2L3'
process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1FastL2L3'
process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1FastL2L3'
process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1FastL2L3'
process.newPFNoPUJetTracksAssociatorAtVertex.jets = 'ak5PFNoPUJetsL1FastL2L3'
process.newPFNoPUJetsSoftElectronTagInfos.jets = 'ak5PFNoPUJetsL1FastL2L3'
process.newPFNoPUJetsSoftMuonTagInfos.jets = 'ak5PFNoPUJetsL1FastL2L3'

# to correct calo met ---
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
#process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
#process.metMuonJESCorAK5.inputUncorMetLabel = "met"
#process.metMuonJESCorAK5.useTypeII = True
#process.metMuonJESCorAK5.hasMuonsCorr = False

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

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

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
process.treeDumper.jetCollection1 = 'ak5CaloJetsL1FastL2L3'
process.treeDumper.JPTjetCollection1 = 'ak5JPTJetsL2L3'
process.treeDumper.PFjetCollection1 = 'ak5PFNoPUJetsL1FastL2L3'
process.treeDumper.PFpuCorrJetCollection1 = 'ak5PFJetsL1FastL2L3'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpLHE = False
process.treeDumper.dumpPdfWeight = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = False
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

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23/emanuele/data/Pool/jpsiEE_Fall10.root') # RECO
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23_2/emanuele/Pool/AODSIM_Winter10_FlatPU.root')
                           fileNames = cms.untracked.vstring('file:/cmsrm/pc24_2/emanuele/data/DYeeSummer11.root')
                            #fileNames = cms.untracked.vstring('file:/cmsrm/pc25/emanuele/data/DYToEE_Fall11_44X.root')
                            )

if(process.treeDumper.dumpPdfWeight == False) :
    process.p = cms.Path ( process.leptonLinkedTracks
                           * process.mergedSuperClusters
                           * process.chargedMetProducer
                           * process.pfIsolationAllSequence
                           * process.ourJetSequenceMCReduced
                           * process.newBtaggingSequence 
                           * process.newPFPUcorrJetBtaggingSequence
                           * process.newPFNoPUJetBtaggingSequence
                           * process.metSequence
                           * process.eIdSequence
                           * process.FastjetForIsolation
                           * process.treeDumper
                           )
else :
    process.p = cms.Path ( process.pdfWeights
                           * process.leptonLinkedTracks
                           * process.mergedSuperClusters
                           * process.chargedMetProducer
                           * process.pfIsolationAllSequence
                           * process.ourJetSequenceMCReduced
                           * process.newBtaggingSequence 
                           * process.newPFPUcorrJetBtaggingSequence
                           * process.newPFNoPUJetBtaggingSequence
                           * process.metSequence
                           * process.eIdSequence
                           * process.FastjetForIsolation
                           * process.treeDumper
                           )
