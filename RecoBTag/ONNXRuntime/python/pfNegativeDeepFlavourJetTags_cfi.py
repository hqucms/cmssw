import FWCore.ParameterSet.Config as cms

from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags

pfNegativeDeepFlavourJetTags = pfDeepFlavourJetTags.deep_clone(
    src = 'pfNegativeDeepFlavourTagInfos'
    )
