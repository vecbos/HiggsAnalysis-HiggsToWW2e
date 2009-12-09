import FWCore.ParameterSet.Config as cms

process = cms.Process("KsPi0Pi0Analysis")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_P_V6::All'

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_pi0pi0_data1.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpSCs = False
process.treeDumper.dumpBCs = True
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = False
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = False
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = False
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpElectrons = False
process.treeDumper.dumpPFlowElectrons = False
process.treeDumper.dumpMuons = False
process.treeDumper.dumpJets = False
process.treeDumper.dumpMet = False
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
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123591_1.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_2.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_36.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123615_2.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_19.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_6.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_22.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_39.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_9.root',
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123591_2.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_20.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_37.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123615_3.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_2.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_7.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_23.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_4.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123778_3.root',
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123592_1.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_21.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_38.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123615_4.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_20.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_8.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_24.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_40.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123791_2.root',
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123592_2.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_22.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_39.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123615_5.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_21.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_9.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_25.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_41.root',
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123592_3.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_23.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_4.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123615_6.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_22.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_1.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_26.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_42.root',
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123592_4.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_24.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_40.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123615_7.root',   'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123732_23.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_10.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_27.root',  'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123734_43.root'
)
                            )

process.p = cms.Path ( process.mergedBasicClusters )                       
process.q = cms.EndPath ( process.treeDumper )
