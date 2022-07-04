import FWCore.ParameterSet.Config as cms

from RecoBTag.FeatureTools.pfDeepBoostedTauJetTagInfos_cfi import pfDeepBoostedTauJetTagInfos
from RecoBTag.ONNXRuntime.boostedJetTauONNXJetTagsProducer_cfi import boostedJetTauONNXJetTagsProducer
from RecoBTag.ONNXRuntime.pfParticleNetAK4TauDiscriminatorsJetTags_cfi import pfParticleNetAK4TauDiscriminatorsJetTags

pfParticleNetAK4TauTagInfos = pfDeepBoostedTauJetTagInfos.clone(
    jet_radius = 0.4,
    min_jet_pt = 15,
    min_puppi_wgt = -1,
    use_puppiP4 = False,
    jets = "ak4PFJetsCHS"
)

pfParticleNetAK4TauJetTags = boostedJetTauONNXJetTagsProducer.clone(
    src = 'pfParticleNetAK4TauTagInfos',
    preprocess_json = 'RecoBTag/ONNXRuntime/data/pnet_ak4_tmp/preprocess.json',
    model_path = 'RecoBTag/ONNXRuntime/data/pnet_ak4_tmp/particle-net.onnx',
    flav_names = ["probmu", "probele", "probtau1h0p", "probtau1h1or2p", "probtau3h0p", "probtau3h1p",
                  "probb", "probc", "probuds", "probg",
                  "ptvis", "ptnu"],
)

from CommonTools.PileupAlgos.Puppi_cff import puppi
from PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi import primaryVertexAssociation

# This task is not used, useful only if we run it from RECO jets (RECO/AOD)
pfParticleNetAK4TauTask = cms.Task(puppi, primaryVertexAssociation, pfParticleNetAK4TauTagInfos,
                                pfParticleNetAK4TauJetTags, pfParticleNetAK4TauDiscriminatorsJetTags)

# declare all the discriminators
# probs
_pfParticleNetAK4TauJetTagsProbs = ['pfParticleNetAK4TauJetTags:' + flav_name
                                 for flav_name in pfParticleNetAK4TauJetTags.flav_names]
# meta-taggers
_pfParticleNetAK4TauJetTagsMetaDiscrs = ['pfParticleNetAK4TauDiscriminatorsJetTags:' + disc.name.value()
                                      for disc in pfParticleNetAK4TauDiscriminatorsJetTags.discriminators]

_pfParticleNetAK4TauJetTagsAll = _pfParticleNetAK4TauJetTagsProbs + _pfParticleNetAK4TauJetTagsMetaDiscrs

