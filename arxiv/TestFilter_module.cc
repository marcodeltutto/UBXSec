////////////////////////////////////////////////////////////////////////
// Class:       TestFilter
// Plugin Type: filter (art v2_05_00)
// File:        TestFilter_module.cc
//
// Generated at Tue Jan 23 08:27:38 2018 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

class TestFilter;


class TestFilter : public art::EDFilter {
public:
  explicit TestFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestFilter(TestFilter const &) = delete;
  TestFilter(TestFilter &&) = delete;
  TestFilter & operator = (TestFilter const &) = delete;
  TestFilter & operator = (TestFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.

};


TestFilter::TestFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  std::cout << __PRETTY_FUNCTION__ << " Called" << std::endl;
}

bool TestFilter::filter(art::Event & e)
{
  // Implementation of required member function here.
  std::cout << __PRETTY_FUNCTION__ << " Called >>>>>>>>>>>>>>>>>>>> FILTER" << std::endl;
  std::cout << "\t Filtering run: " << e.id().run() << ", subrun: " << e.id().subRun() << ", event: " << e.event() << std::endl;
  return false;
}


DEFINE_ART_MODULE(TestFilter)
