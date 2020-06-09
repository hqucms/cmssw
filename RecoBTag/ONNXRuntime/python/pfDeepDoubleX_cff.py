from RecoBTag.TensorFlow.pfDeepDoubleBvLJetTags_cfi import pfDeepDoubleBvLJetTags
from RecoBTag.TensorFlow.pfDeepDoubleCvBJetTags_cfi import pfDeepDoubleCvBJetTags
from RecoBTag.TensorFlow.pfDeepDoubleCvLJetTags_cfi import pfDeepDoubleCvLJetTags

pfMassIndependentDeepDoubleBvLJetTags = pfDeepDoubleBvLJetTags.clone(
    graph_path = 'RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDB_mass_independent.pb')
pfMassIndependentDeepDoubleCvLJetTags = pfDeepDoubleCvLJetTags.clone(
    graph_path = 'RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDC_mass_independent.pb')
pfMassIndependentDeepDoubleCvBJetTags = pfDeepDoubleCvBJetTags.clone(
    graph_path = 'RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDCvB_mass_independent.pb')

# ONNXRuntime-based implementation
from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel
from RecoBTag.ONNXRuntime.pfDeepDoubleBvLJetTags_cfi import pfDeepDoubleBvLJetTags as _pfDeepDoubleBvLJetTags
from RecoBTag.ONNXRuntime.pfDeepDoubleCvBJetTags_cfi import pfDeepDoubleCvBJetTags as _pfDeepDoubleCvBJetTags
from RecoBTag.ONNXRuntime.pfDeepDoubleCvLJetTags_cfi import pfDeepDoubleCvLJetTags as _pfDeepDoubleCvLJetTags

run2_miniAOD_devel.toModify(pfDeepDoubleBvLJetTags, _pfDeepDoubleBvLJetTags)
run2_miniAOD_devel.toModify(pfDeepDoubleCvBJetTags, _pfDeepDoubleCvBJetTags)
run2_miniAOD_devel.toModify(pfDeepDoubleCvLJetTags, _pfDeepDoubleCvLJetTags)

run2_miniAOD_devel.toModify(
    pfMassIndependentDeepDoubleBvLJetTags,
    _pfDeepDoubleBvLJetTags.clone(
        model_path = 'RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDB_mass_independent.onnx')
)
run2_miniAOD_devel.toModify(
    pfMassIndependentDeepDoubleCvLJetTags,
    _pfDeepDoubleCvLJetTags.clone(
        model_path='RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDC_mass_independent.onnx')
)
run2_miniAOD_devel.toModify(
    pfMassIndependentDeepDoubleCvBJetTags,
    _pfDeepDoubleCvBJetTags.clone(
        model_path='RecoBTag/Combined/data/DeepDoubleX/94X/V01/DDCvB_mass_independent.onnx')
)
