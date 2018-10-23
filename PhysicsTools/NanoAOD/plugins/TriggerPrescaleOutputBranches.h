#ifndef PhysicsTools_NanoAOD_TriggerPrescaleOutputBranches_h
#define PhysicsTools_NanoAOD_TriggerPrescaleOutputBranches_h

#include <string>
#include <vector>
#include <TTree.h>
#include "FWCore/Framework/interface/EventForOutput.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

class TriggerPrescaleOutputBranches {
public:
  TriggerPrescaleOutputBranches(const edm::BranchDescription *desc, const edm::EDGetToken & token ) :
    m_token(token), m_lastRun(-1),m_fills(0)
{
    if (desc->className() != "pat::PackedTriggerPrescales") throw cms::Exception("Configuration", "NanoAODOutputModule/TriggerPrescaleOutputBranches can only write out pat::PackedTriggerPrescales objects");
}

  void updateTriggerNames(TTree &tree, const edm::TriggerNames & names, const edm::TriggerResults & ta);
  void fill(const edm::EventForOutput &iEvent, TTree & tree) ;

private:
  edm::TriggerNames triggerNames(const edm::TriggerResults& triggerResults); //FIXME: if we have to keep it local we may use PsetID check per event instead of run boundary

  edm::EDGetToken m_token;
  std::string  m_baseName;
  bool         m_singleton;
  UInt_t       m_counter;
  struct NamedBranchPtr {
    std::string name, title;
    int idx;
    TBranch * branch;
    uint16_t buffer;
    NamedBranchPtr(const std::string & aname, const std::string & atitle, TBranch *branchptr = nullptr) :
      name(aname), title(atitle), idx(-1), branch(branchptr), buffer(-1) {}
  };
  std::vector<NamedBranchPtr> m_triggerPrescaleBranches;
  long m_lastRun;
  unsigned long m_fills;

  template<typename T>
  void fillColumn(NamedBranchPtr & nb, const pat::PackedTriggerPrescales & triggerPrescales) {
    if(nb.idx>=0) nb.buffer=triggerPrescales.getPrescaleForIndex(nb.idx);
    nb.branch->SetAddress(&(nb.buffer)); // Can be improved: this is not reallt needed at each event
    //but we should be sure that resize of vectors of TriggerPrescaleOutputBranches do not mess up things
  }

};

#endif

