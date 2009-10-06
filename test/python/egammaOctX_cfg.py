import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V5::All'

# --- common code to comput gg fusion signal k factor
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWKFactorProducer_cfi")

# --- apply the muon correction from the common code
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWMetCorrector_cfi")

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.trackCandidates_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default.root'
process.treeDumper.metCollection = 'tcMet'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpMuons = False
process.treeDumper.dumpTracks = False
process.treeDumper.dumpVertices = False
process.treeDumper.dumpParticleFlowObjects = False
process.treeDumper.saveFatTrk = False
process.treeDumper.saveJet1BTag = True
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpTree = True

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
#                            fileNames = cms.untracked.vstring('file:/cmsrm/pc17/emanuele/data/Pool/EGM2ELE15Skim_2.root')
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/e/emanuele/EGMSEL_Oct09X/EGM2ELE15Skim_2.root')
                            )

process.p = cms.Path ( process.jetSequence * process.pfjetSCSequence * process.newBtaggingSequence *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons *
                       process.trackCandidates)
                       
process.q = cms.EndPath ( process.treeDumper )
