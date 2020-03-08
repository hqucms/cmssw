import FWCore.ParameterSet.Config as cms

from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags

pfNegativeDeepFlavourJetTags = pfDeepFlavourJetTags.clone(
    native = {'src': 'pfNegativeDeepFlavourTagInfos'},
    onnx = {'src': 'pfNegativeDeepFlavourTagInfos'},
    )
