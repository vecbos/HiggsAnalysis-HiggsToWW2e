import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newPFNoPUJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newPFNoPUJetTracksAssociatorAtVertex.jets = "ak5PFNoPUJetsL1FastL2L3Residual"
newPFNoPUJetTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
import RecoBTag.Configuration.RecoBTag_cff
newPFNoPUJetsImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
newPFNoPUJetsImpactParameterTagInfos.jetTracks = "newPFNoPUJetTracksAssociatorAtVertex"
newTrackCountingHighEffBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos") )
newTrackCountingHighPurBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos") )
newJetProbabilityBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
newJetProbabilityBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos") )
newJetBProbabilityBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
newJetBProbabilityBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos") )
newImpactParameterMVABPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.impactParameterMVABJetTags.clone()
newImpactParameterMVABPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos") )

# secondary vertex b-tag
newPFNoPUJetsSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
newPFNoPUJetsSecondaryVertexTagInfos.trackIPTagInfos = "newPFNoPUJetsImpactParameterTagInfos"
newSimpleSecondaryVertexHighEffBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
newSimpleSecondaryVertexHighEffBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSecondaryVertexTagInfos") )
newSimpleSecondaryVertexHighPurBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
newSimpleSecondaryVertexHighPurBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos"), cms.InputTag("newPFNoPUJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsImpactParameterTagInfos"), cms.InputTag("newPFNoPUJetsSecondaryVertexTagInfos") )

# soft electron b-tag
newPFNoPUJetsSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newPFNoPUJetsSoftElectronTagInfos.jets = "ak5PFJetsL2L3Residual"
newSoftElectronBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSoftElectronTagInfos") )
newSoftElectronByIP3dBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByIP3dBJetTags.clone()
newSoftElectronByIP3dBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSoftElectronTagInfos") )
newSoftElectronByPtBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByPtBJetTags.clone()
newSoftElectronByPtBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSoftElectronTagInfos") )

# soft muon b-tag
newPFNoPUJetsSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newPFNoPUJetsSoftMuonTagInfos.jets = "ak5PFJetsL2L3Residual"
newSoftMuonBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSoftMuonTagInfos") )
newSoftMuonByIP3dBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByIP3dBJetTags.clone()
newSoftMuonByIP3dBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSoftMuonTagInfos") )
newSoftMuonByPtBPFNoPUJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByPtBJetTags.clone()
newSoftMuonByPtBPFNoPUJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFNoPUJetsSoftMuonTagInfos") )

# prepare a path running the new modules
newPFNoPUJetTracksAssociator = cms.Sequence(
    newPFNoPUJetTracksAssociatorAtVertex
    )

newPFNoPUJetBtaggingIP = cms.Sequence(
    newPFNoPUJetsImpactParameterTagInfos * (
       newTrackCountingHighEffBPFNoPUJetTags +
       newTrackCountingHighPurBPFNoPUJetTags +
       newJetProbabilityBPFNoPUJetTags +
       newJetBProbabilityBPFNoPUJetTags )
    )

newPFNoPUJetBtaggingSV = cms.Sequence(
    newPFNoPUJetsImpactParameterTagInfos *
    newPFNoPUJetsSecondaryVertexTagInfos * (
       newSimpleSecondaryVertexHighEffBPFNoPUJetTags +
       newSimpleSecondaryVertexHighPurBPFNoPUJetTags +
       newCombinedSecondaryVertexBPFNoPUJetTags +
       newCombinedSecondaryVertexMVABPFNoPUJetTags )
    )

newPFNoPUJetBtaggingEle = cms.Sequence(
    softElectronCands * # already run for calojets
    newPFNoPUJetsSoftElectronTagInfos *
       newSoftElectronBPFNoPUJetTags +
       newSoftElectronByIP3dBPFNoPUJetTags +
       newSoftElectronByPtBPFNoPUJetTags
    )

newPFNoPUJetBtaggingMu = cms.Sequence(
    newPFNoPUJetsSoftMuonTagInfos * (
       newSoftMuonBPFNoPUJetTags +
       newSoftMuonByIP3dBPFNoPUJetTags +
       newSoftMuonByPtBPFNoPUJetTags )
    )

newPFNoPUJetBtagging = cms.Sequence(
    newPFNoPUJetBtaggingIP +
    newPFNoPUJetBtaggingSV +
    newPFNoPUJetBtaggingEle +
    newPFNoPUJetBtaggingMu )

newPFNoPUJetBtaggingSequence = cms.Sequence(
    newPFNoPUJetTracksAssociator *
       newPFNoPUJetBtagging )
