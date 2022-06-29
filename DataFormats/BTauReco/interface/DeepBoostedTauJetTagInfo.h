#ifndef DataFormats_BTauReco_DeepBoostedTauJetTagInfo_h
#define DataFormats_BTauReco_DeepBoostedTauJetTagInfo_h

#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"
#include "DataFormats/BTauReco/interface/DeepBoostedTauJetFeatures.h"

namespace reco {

typedef  FeaturesTagInfo<btagbtvdeep::DeepBoostedTauJetFeatures> DeepBoostedTauJetTagInfo;

DECLARE_EDM_REFS( DeepBoostedTauJetTagInfo )

}

#endif // DataFormats_BTauReco_DeepBoostedTauJetTagInfo_h
