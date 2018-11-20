#include "PhysicsTools/NanoAOD/plugins/TriggerPrescaleOutputBranches.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include <iostream>

void
TriggerPrescaleOutputBranches::updateTriggerNames(TTree & tree, const edm::TriggerNames & names, const edm::TriggerResults & triggers)
{
   std::vector<std::string> newNames(triggers.getTriggerNames());
   if(newNames.empty()) {
       for(unsigned int j=0;j<triggers.size();j++) {
           newNames.push_back(names.triggerName(j));
       }
   }

   for(auto & existing : m_triggerPrescaleBranches)
   {
       existing.idx=-1;// reset all triggers as not found
       for(unsigned int j=0;j<newNames.size();j++) {
	 std::string name=newNames[j]; // no const & as it will be modified below!
	 std::size_t vfound = name.rfind("_v");
	 if (vfound!=std::string::npos){
           name.replace(vfound,name.size()-vfound,"");
	 }
	 if(name==existing.name) existing.idx=j;
       }
   }
   // Find new ones
   for(unsigned int j=0;j<newNames.size();j++) {
       std::string name=newNames[j]; // no const & as it will be modified below!
       std::size_t vfound = name.rfind("_v");
       if (vfound!=std::string::npos){
           name.replace(vfound,name.size()-vfound,"");
       }
       bool found=false;
       if(name.compare(0,3,"HLT")==0 || name.compare(0,2,"L1")==0 ){
           for(auto & existing : m_triggerPrescaleBranches) {if(name==existing.name) found=true;}
           if(!found){
                std::string branchName = "Prescale_"+name;
                NamedBranchPtr nb(branchName,"Trigger prescale"); //FIXME: If the title can be updated we can use it to list the versions _v* that were seen in this file
                uint16_t backFillValue=0;
                nb.branch= tree.Branch(branchName.c_str(), &backFillValue, (branchName + "/s").c_str());
                nb.branch->SetTitle(nb.title.c_str());
                nb.idx=j;
                m_triggerPrescaleBranches.push_back(nb);
                for(size_t i=0;i<m_fills;i++) nb.branch->Fill(); // Back fill
           }
       }
   }
}

edm::TriggerNames TriggerPrescaleOutputBranches::triggerNames(const edm::TriggerResults& triggerResults){
    edm::pset::Registry* psetRegistry = edm::pset::Registry::instance();
    edm::ParameterSet const* pset=nullptr;
    if (nullptr!=(pset=psetRegistry->getMapped(triggerResults.parameterSetID()))) {

         if (pset->existsAs<std::vector<std::string> >("@trigger_paths", true)) {
            edm::TriggerNames triggerNames(*pset);

            // This should never happen
            if (triggerNames.size() != triggerResults.size()) {
               throw cms::Exception("LogicError")
                  << "edm::EventBase::triggerNames_ Encountered vector\n"
                     "of trigger names and a TriggerResults object with\n"
                     "different sizes.  This should be impossible.\n"
                     "Please send information to reproduce this problem to\n"
                     "the edm developers.\n";
            }
            return triggerNames;
         }
      }
    return edm::TriggerNames();
}

void TriggerPrescaleOutputBranches::fill(const edm::EventForOutput &iEvent, TTree & tree)
{
    edm::Handle<pat::PackedTriggerPrescales> handle;
    iEvent.getByToken(m_token, handle);
    const pat::PackedTriggerPrescales & triggerPrescales = *handle;
    const edm::TriggerNames &names = triggerNames(triggerPrescales.triggerResults());

    if(m_lastRun!=iEvent.id().run()) {
        m_lastRun=iEvent.id().run();
        updateTriggerNames(tree, names, triggerPrescales.triggerResults());
    }
    for (auto & pair : m_triggerPrescaleBranches) fillColumn<uint16_t>(pair, triggerPrescales);
    m_fills++;
}



