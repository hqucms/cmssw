from RecoBTag.ONNXRuntime.SwitchProducerONNX import SwitchProducerONNX
from RecoBTag.ONNXRuntime.deepFlavourONNXJetTagsProducer_cfi import deepFlavourONNXJetTagsProducer
from RecoBTag.TensorFlow.deepFlavourTFJetTagsProducer_cfi import deepFlavourTFJetTagsProducer

pfDeepFlavourJetTags = SwitchProducerONNX(
    native = deepFlavourTFJetTagsProducer.clone(),
    onnx = deepFlavourONNXJetTagsProducer.clone(),
)
