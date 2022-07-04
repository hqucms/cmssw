import FWCore.ParameterSet.Config as cms

from RecoBTag.FeatureTools.pfDeepBoostedTauJetTagInfos_cfi import pfDeepBoostedTauJetTagInfos
from RecoBTag.ONNXRuntime.boostedJetTauONNXJetTagsProducer_cfi import boostedJetTauONNXJetTagsProducer
#from RecoBTag.ONNXRuntime.pfParticleNetDiscriminatorsJetTags_cfi import pfParticleNetDiscriminatorsJetTags
from RecoBTag.ONNXRuntime.pfMassDecorrelatedParticleNetTauDiscriminatorsJetTags_cfi import pfMassDecorrelatedParticleNetTauDiscriminatorsJetTags

pfParticleNetTauTagInfos = pfDeepBoostedTauJetTagInfos.clone(
    use_puppiP4 = False
)

pfMassDecorrelatedParticleNetTauJetTags = boostedJetTauONNXJetTagsProducer.clone(
    src = 'pfParticleNetTauTagInfos',
    preprocess_json = 'RecoBTag/ONNXRuntime/data/pnet_ak8_tmp/preprocess.json',
    model_path = 'RecoBTag/ONNXRuntime/data/pnet_ak8_tmp/particle-net.onnx',
    flav_names = ["probXtt", "probXtm", "probXte", "probXmm", "probXee",
                  "probXbb", "probXcc", "probXqq", "probXgg", "probQCD",
                  "mass"],
)



from CommonTools.PileupAlgos.Puppi_cff import puppi
from PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi import primaryVertexAssociation

# This task is not used, useful only if we run it from RECO jets (RECO/AOD)
pfParticleNetTauTask = cms.Task(puppi, primaryVertexAssociation, pfParticleNetTauTagInfos,
                                pfMassDecorrelatedParticleNetTauJetTags, pfMassDecorrelatedParticleNetTauDiscriminatorsJetTags)

# declare all the discriminators
# mass-decorrelated: probs
_pfMassDecorrelatedParticleNetTauJetTagsProbs = ['pfMassDecorrelatedParticleNetTauJetTags:' + flav_name
                                                 for flav_name in pfMassDecorrelatedParticleNetTauJetTags.flav_names]
# mass-decorrelated: meta-taggers
_pfMassDecorrelatedParticleNetTauJetTagsMetaDiscrs = ['pfMassDecorrelatedParticleNetTauDiscriminatorsJetTags:' + disc.name.value()
                                                      for disc in pfMassDecorrelatedParticleNetTauDiscriminatorsJetTags.discriminators]

_pfParticleNetTauJetTagsAll = _pfMassDecorrelatedParticleNetTauJetTagsProbs + _pfMassDecorrelatedParticleNetTauJetTagsMetaDiscrs
