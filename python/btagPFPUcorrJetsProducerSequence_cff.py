import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newPFPUcorrJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newPFPUcorrJetTracksAssociatorAtVertex.jets = "ak5PFJetsL1FastL2L3Residual"
newPFPUcorrJetTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
import RecoBTag.Configuration.RecoBTag_cff
newPFPUcorrJetsImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
newPFPUcorrJetsImpactParameterTagInfos.jetTracks = "newPFPUcorrJetTracksAssociatorAtVertex"
newTrackCountingHighEffBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos") )
newTrackCountingHighPurBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos") )
newJetProbabilityBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
newJetProbabilityBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos") )
newJetBProbabilityBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
newJetBProbabilityBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos") )
newImpactParameterMVABPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.impactParameterMVABJetTags.clone()
newImpactParameterMVABPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos") )

# impact parameter done with the first track instead of the second
from RecoBTag.ImpactParameter.trackCounting3D1stComputer_cfi import *
import RecoBTag.ImpactParameter.trackCountingVeryHighEffBJetTags_cfi
newTrackCountingVeryHighEffBPFPUcorrJetTags = RecoBTag.ImpactParameter.trackCountingVeryHighEffBJetTags_cfi.trackCountingVeryHighEffBJetTags.clone()
newTrackCountingVeryHighEffBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos") )

# secondary vertex b-tag
newPFPUcorrJetsSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
newPFPUcorrJetsSecondaryVertexTagInfos.trackIPTagInfos = "newPFPUcorrJetsImpactParameterTagInfos"
newSimpleSecondaryVertexHighEffBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
newSimpleSecondaryVertexHighEffBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSecondaryVertexTagInfos") )
newSimpleSecondaryVertexHighPurBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
newSimpleSecondaryVertexHighPurBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos"), cms.InputTag("newPFPUcorrJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsImpactParameterTagInfos"), cms.InputTag("newPFPUcorrJetsSecondaryVertexTagInfos") )

# soft electron b-tag
newPFPUcorrJetsSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newPFPUcorrJetsSoftElectronTagInfos.jets = "ak5PFJetsL1FastL2L3Residual"
newSoftElectronBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSoftElectronTagInfos") )
newSoftElectronByIP3dBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByIP3dBJetTags.clone()
newSoftElectronByIP3dBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSoftElectronTagInfos") )
newSoftElectronByPtBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByPtBJetTags.clone()
newSoftElectronByPtBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSoftElectronTagInfos") )

# soft muon b-tag
newPFPUcorrJetsSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newPFPUcorrJetsSoftMuonTagInfos.jets = "ak5PFJetsL1FastL2L3Residual"
newSoftMuonBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSoftMuonTagInfos") )
newSoftMuonByIP3dBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByIP3dBJetTags.clone()
newSoftMuonByIP3dBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSoftMuonTagInfos") )
newSoftMuonByPtBPFPUcorrJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByPtBJetTags.clone()
newSoftMuonByPtBPFPUcorrJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFPUcorrJetsSoftMuonTagInfos") )

# prepare a path running the new modules
newPFPUcorrJetTracksAssociator = cms.Sequence(
    newPFPUcorrJetTracksAssociatorAtVertex
    )

newPFPUcorrJetBtaggingIP = cms.Sequence(
    newPFPUcorrJetsImpactParameterTagInfos * (
       newTrackCountingVeryHighEffBPFPUcorrJetTags +
       newTrackCountingHighEffBPFPUcorrJetTags +
       newTrackCountingHighPurBPFPUcorrJetTags +
       newJetProbabilityBPFPUcorrJetTags +
       newJetBProbabilityBPFPUcorrJetTags )
    )

newPFPUcorrJetBtaggingSV = cms.Sequence(
    newPFPUcorrJetsImpactParameterTagInfos *
    newPFPUcorrJetsSecondaryVertexTagInfos * (
       newSimpleSecondaryVertexHighEffBPFPUcorrJetTags +
       newSimpleSecondaryVertexHighPurBPFPUcorrJetTags +
       newCombinedSecondaryVertexBPFPUcorrJetTags +
       newCombinedSecondaryVertexMVABPFPUcorrJetTags )
    )

newPFPUcorrJetBtaggingEle = cms.Sequence(
    softElectronCands * # already run for calojets
    newPFPUcorrJetsSoftElectronTagInfos *
       newSoftElectronBPFPUcorrJetTags +
       newSoftElectronByIP3dBPFPUcorrJetTags +
       newSoftElectronByPtBPFPUcorrJetTags
    )

newPFPUcorrJetBtaggingMu = cms.Sequence(
    newPFPUcorrJetsSoftMuonTagInfos * (
       newSoftMuonBPFPUcorrJetTags +
       newSoftMuonByIP3dBPFPUcorrJetTags +
       newSoftMuonByPtBPFPUcorrJetTags )
    )

newPFPUcorrJetBtagging = cms.Sequence(
    newPFPUcorrJetBtaggingIP +
    newPFPUcorrJetBtaggingSV +
    newPFPUcorrJetBtaggingEle +
    newPFPUcorrJetBtaggingMu )

newPFPUcorrJetBtaggingSequence = cms.Sequence(
    newPFPUcorrJetTracksAssociator *
       newPFPUcorrJetBtagging )
