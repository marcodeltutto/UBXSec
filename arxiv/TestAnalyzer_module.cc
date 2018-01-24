////////////////////////////////////////////////////////////////////////
// Class:       TestAnalyzer
// Plugin Type: analyzer (art v2_05_00)
// File:        TestAnalyzer_module.cc
//
// Generated at Tue Jan 23 08:32:47 2018 by Marco Del Tutto using cetskelgen
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

class TestAnalyzer;


class TestAnalyzer : public art::EDAnalyzer {
public:
  explicit TestAnalyzer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestAnalyzer(TestAnalyzer const &) = delete;
  TestAnalyzer(TestAnalyzer &&) = delete;
  TestAnalyzer & operator = (TestAnalyzer const &) = delete;
  TestAnalyzer & operator = (TestAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  void beginRun(const art::Run&) override;
  void endSubRun(const art::SubRun &sr) override;
  void beginSubRun(const art::SubRun &sr) override;

private:

  // Declare member data here.

};


TestAnalyzer::TestAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void TestAnalyzer::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  std::cout << __PRETTY_FUNCTION__ << " Called" << std::endl;
  std::cout << "\t Analyzing run: " << e.id().run() << ", subrun: " << e.id().subRun() << ", event: " << e.event() << std::endl;
}


void TestAnalyzer::beginRun(const art::Run&)
{
    std::cout << __PRETTY_FUNCTION__ << " Called" << std::endl;
}

void TestAnalyzer::beginSubRun(const art::SubRun& sr)
{
    std::cout << __PRETTY_FUNCTION__ << " Called, run " << sr.run() << ", subrun " << sr.subRun() << std::endl;
}

void TestAnalyzer::endSubRun(const art::SubRun& sr)
{
    std::cout << __PRETTY_FUNCTION__ << " Called, run " << sr.run() << ", subrun " << sr.subRun() << std::endl;
}
DEFINE_ART_MODULE(TestAnalyzer)
