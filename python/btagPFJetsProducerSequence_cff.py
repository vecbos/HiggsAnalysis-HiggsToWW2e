import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newPFJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newPFJetTracksAssociatorAtVertex.jets = "ak5PFJetsL2L3Residual"
newPFJetTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
import RecoBTag.Configuration.RecoBTag_cff
newPFJetsImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
newPFJetsImpactParameterTagInfos.jetTracks = "newPFJetTracksAssociatorAtVertex"
newTrackCountingHighEffBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos") )
newTrackCountingHighPurBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos") )
newJetProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
newJetProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos") )
newJetBProbabilityBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
newJetBProbabilityBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos") )
newImpactParameterMVABPFJetTags = RecoBTag.Configuration.RecoBTag_cff.impactParameterMVABJetTags.clone()
newImpactParameterMVABPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos") )

# secondary vertex b-tag
newPFJetsSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
newPFJetsSecondaryVertexTagInfos.trackIPTagInfos = "newPFJetsImpactParameterTagInfos"
newSimpleSecondaryVertexHighEffBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
newSimpleSecondaryVertexHighEffBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSecondaryVertexTagInfos") )
newSimpleSecondaryVertexHighPurBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
newSimpleSecondaryVertexHighPurBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos"), cms.InputTag("newPFJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABPFJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsImpactParameterTagInfos"), cms.InputTag("newPFJetsSecondaryVertexTagInfos") )

# soft electron b-tag
newPFJetsSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newPFJetsSoftElectronTagInfos.jets = "ak5PFJetsL2L3Residual"
newSoftElectronBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSoftElectronTagInfos") )
newSoftElectronByIP3dBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByIP3dBJetTags.clone()
newSoftElectronByIP3dBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSoftElectronTagInfos") )
newSoftElectronByPtBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByPtBJetTags.clone()
newSoftElectronByPtBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSoftElectronTagInfos") )

# soft muon b-tag
newPFJetsSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newPFJetsSoftMuonTagInfos.jets = "ak5PFJetsL2L3Residual"
newSoftMuonBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSoftMuonTagInfos") )
newSoftMuonByIP3dBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByIP3dBJetTags.clone()
newSoftMuonByIP3dBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSoftMuonTagInfos") )
newSoftMuonByPtBPFJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByPtBJetTags.clone()
newSoftMuonByPtBPFJetTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsSoftMuonTagInfos") )

# prepare a path running the new modules
newPFJetTracksAssociator = cms.Sequence(
    newPFJetTracksAssociatorAtVertex
    )

newPFJetBtaggingIP = cms.Sequence(
    newPFJetsImpactParameterTagInfos * (
       newTrackCountingHighEffBPFJetTags +
       newTrackCountingHighPurBPFJetTags +
       newJetProbabilityBPFJetTags +
       newJetBProbabilityBPFJetTags )
    )

newPFJetBtaggingSV = cms.Sequence(
    newPFJetsImpactParameterTagInfos *
    newPFJetsSecondaryVertexTagInfos * (
       newSimpleSecondaryVertexHighEffBPFJetTags +
       newSimpleSecondaryVertexHighPurBPFJetTags +
       newCombinedSecondaryVertexBPFJetTags +
       newCombinedSecondaryVertexMVABPFJetTags )
    )

newPFJetBtaggingEle = cms.Sequence(
    softElectronCands * # already run for calojets
    newPFJetsSoftElectronTagInfos *
       newSoftElectronBPFJetTags +
       newSoftElectronByIP3dBPFJetTags +
       newSoftElectronByPtBPFJetTags
    )

newPFJetBtaggingMu = cms.Sequence(
    newPFJetsSoftMuonTagInfos * (
       newSoftMuonBPFJetTags +
       newSoftMuonByIP3dBPFJetTags +
       newSoftMuonByPtBPFJetTags )
    )

newPFJetBtagging = cms.Sequence(
    newPFJetBtaggingIP +
    newPFJetBtaggingSV +
    newPFJetBtaggingEle +
    newPFJetBtaggingMu )

newPFJetBtaggingSequence = cms.Sequence(
    newPFJetTracksAssociator *
       newPFJetBtagging )
