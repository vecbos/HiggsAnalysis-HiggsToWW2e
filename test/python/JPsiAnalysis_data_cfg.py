import FWCore.ParameterSet.Config as cms

process = cms.Process("JPsiAnalysis")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_P_V6::All'

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- track sequences ---
process.load("RecoTracker.DeDx.dedxEstimatorsFromRefitter_cff")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_data.root'
process.treeDumper.dumpTriggerResults = False
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpSCs = True
process.treeDumper.dumpTracks = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = False
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpTree = True

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring(
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_1.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_11.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_12.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_13.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_16.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_17.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_2.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_23.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_24.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_25.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_26.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_27.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_29.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_3.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_30.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_46.root',
'file:/cmsrm/pc21/emanuele/data/900GeVBscFilterAndReReco_Paolo/reRecoOutput_123151_6.root'
                            )
                            )

process.dedx = cms.Sequence (process.RefitterForDeDx * process.dedxTruncated40)

process.p = cms.Path ( process.jetSequence * process.pfjetSCSequence * process.newBtaggingSequence *
#                      process.doAlldEdXEstimators *
#                      process.dedx *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons )
                       
process.q = cms.EndPath ( process.treeDumper )
