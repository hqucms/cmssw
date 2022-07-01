import FWCore.ParameterSet.Config as cms

pfParticleNetAK4TauDiscriminatorsJetTags = cms.EDProducer(
   'BTagProbabilityToDiscriminator',
   discriminators = cms.VPSet(
      cms.PSet(
         name = cms.string('BvsAll'),
         numerator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probb'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probbb'),
            ),
         denominator=cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probb'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probbb'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probcc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probuds'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probg'),
         ),
      ),
      cms.PSet(
         name = cms.string('CvsL'),
         numerator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probcc'),
            ),
         denominator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probcc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probuds'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probg'),
            ),
         ),
      cms.PSet(
         name = cms.string('CvsB'),
         numerator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probcc'),
            ),
         denominator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probcc'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probb'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probbb'),
            ),
         ),
      cms.PSet(
         name = cms.string('QvsG'),
         numerator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probuds'),
            ),
         denominator = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probuds'),
            cms.InputTag('pfParticleNetAK4TauJetTags', 'probg'),
            ),
         ),

      )
   )
