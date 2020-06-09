import FWCore.ParameterSet.Config as cms

from RecoBTag.TensorFlow.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags

from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel
from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import _pfDeepFlavourJetTags
run2_miniAOD_devel.toModify(pfDeepFlavourJetTags, _pfDeepFlavourJetTags)

pfNegativeDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
        src = cms.InputTag('pfNegativeDeepFlavourTagInfos')
        )
