////////////////////////////////////////////////////////////////////////
// Class:       SimpleAna
// Plugin Type: analyzer (art v2_05_00)
// File:        SimpleAna_module.cc
//
// Generated at Tue Oct 31 16:36:23 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "uboone/UBXSec/DataTypes/SelectionResult.h"

class SimpleAna;


class SimpleAna : public art::EDAnalyzer {
public:
  explicit SimpleAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleAna(SimpleAna const &) = delete;
  SimpleAna(SimpleAna &&) = delete;
  SimpleAna & operator = (SimpleAna const &) = delete;
  SimpleAna & operator = (SimpleAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

};


SimpleAna::SimpleAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void SimpleAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  art::Handle<std::vector<ubana::SelectionResult>> selection_h;
  e.getByLabel("UBXSec",selection_h);
  if (!selection_h.isValid() || selection_h->empty()) {
    std::cout << "[SimpleAna] SelectionResult handle is not valid or empty." << std::endl;
  }

  std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
  art::fill_ptr_vector(selection_v, selection_h);

  if (selection_v.at(0)->GetSelectionStatus()) {
    std::cout << "[SimpleAna] Event is selected" << std::endl;
  } else {
    std::cout << "[SimpleAna] Event is not selected" << std::endl;
    std::cout << "[SimpleAna] Failure reason " << selection_v.at(0)->GetFailureReason()  << std::endl;
  }

}

DEFINE_ART_MODULE(SimpleAna)
