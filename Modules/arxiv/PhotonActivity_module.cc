////////////////////////////////////////////////////////////////////////
// Class:       PhotonActivity
// Plugin Type: producer (art v2_05_00)
// File:        PhotonActivity_module.cc
//
// Generated at Fri Apr 21 12:50:07 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "TTree.h"

#include <memory>

class PhotonActivity;


class PhotonActivity : public art::EDProducer {
public:
  explicit PhotonActivity(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonActivity(PhotonActivity const &) = delete;
  PhotonActivity(PhotonActivity &&) = delete;
  PhotonActivity & operator = (PhotonActivity const &) = delete;
  PhotonActivity & operator = (PhotonActivity &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string _mctruth_label;
  std::string _trigger_label;
  std::string _simphot_label;
  bool        _debug;

  double _beam_window_start, _beam_window_end, _before_spill_window;

  int _run, _subrun, _event;
  int _n_ph_beam_spill = 0;
  int _n_ph_before_beam_spill = 0;

  void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);

  TTree* _tree1;
};


PhotonActivity::PhotonActivity(fhicl::ParameterSet const & p)
{
  _mctruth_label = p.get<std::string>("MCTruthProduct", "generator");
  _trigger_label = p.get<std::string>("TriggerProduct", "triggersim");
  _simphot_label = p.get<std::string>("SimPhotProduct", "largeant");    
  _debug         = p.get<bool>       ("DebugMode",      false);

  _beam_window_start   = p.get<double> ("BeamWindowStart",   3.2);
  _beam_window_end     = p.get<double> ("BeamWindowEnd",     4.8);
  _before_spill_window = p.get<double> ("BeforeSpillWindow", 2.0);

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",                &_run,                "run/I");
  _tree1->Branch("subrun",             &_subrun,             "subrun/I");
  _tree1->Branch("event",              &_event,              "event/I");
  _tree1->Branch("n_ph_beam_spill",    &_n_ph_beam_spill,    "n_ph_beam_spill/I");
  _tree1->Branch("n_ph_before_beam_spill", &_n_ph_before_beam_spill, "n_ph_before_beam_spill/I");

  //produces< std::vector<recob::OpFlash> >();
}

void PhotonActivity::produce(art::Event & e)
{
  if (_debug) std::cout << "***** PhotonActivity starts." << std::endl;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // produce OpFlash data-product to be filled within module
  //std::unique_ptr< std::vector<recob::OpFlash> > opflashes(new std::vector<recob::OpFlash>);

  art::Handle<std::vector<raw::Trigger> > evt_trigger_h;
  e.getByLabel(_trigger_label,evt_trigger_h);

  if( !evt_trigger_h.isValid() || evt_trigger_h->empty() ) {
    std::cerr << "Trigger product is not valid or empty." << std::endl;
    return;
  }

  art::Handle<std::vector<simb::MCTruth> > evt_mctruth_h;
  e.getByLabel(_mctruth_label,evt_mctruth_h);

  if( !evt_mctruth_h.isValid() || evt_mctruth_h->empty() ) {
    std::cerr << "MCTruth product is not valid or empty." << std::endl;
    return;
  }

  art::Handle<std::vector<sim::SimPhotons> > evt_simphot_h;
  e.getByLabel(_simphot_label,evt_simphot_h);

  if( !evt_simphot_h.isValid() || evt_simphot_h->empty() ) {
    std::cerr << "SimPhotons product is not valid or empty." << std::endl;
    return;
  }

  ::art::ServiceHandle<geo::Geometry> geo; 

  if(evt_simphot_h->size() != geo->NOpDets()) {
    std::cerr << "Unexpected # of channels in simphotons!" << std::endl;
    return;
  }

  // opdet=>opchannel mapping
  std::vector<size_t> opdet2opch(geo->NOpDets(),0);
  for(size_t opch=0; opch<opdet2opch.size(); ++opch){
    opdet2opch[geo->OpDetFromOpChannel(opch)] = opch;
  }

  
  auto const & evt_trigger = (*evt_trigger_h)[0];
  auto const trig_time = evt_trigger.TriggerTime();
  auto const * ts = lar::providerFrom<detinfo::DetectorClocksService>();

  /*
  double nuTime = 0.;
  if (_debug) std::cout << "We have " << evt_mctruth_h->size() << " mctruth events." << std::endl;
  for (size_t n = 0; n < evt_mctruth_h->size(); n++) {

    simb::MCTruth const& evt_mctruth = (*evt_mctruth_h)[n];
    if (_debug) std::cout << "Origin: " << evt_mctruth.Origin() << std::endl;
    if (evt_mctruth.Origin() != 1 ) continue;
    if (_debug) std::cout << "We have " << evt_mctruth.NParticles() << " particles." << std::endl;
    for (int p = 0; p < evt_mctruth.NParticles(); p++) {
   
      simb::MCParticle const& par = evt_mctruth.GetParticle(p);
      //if (par.PdgCode() != 14) continue;
      if (_debug){
        std::cout << "Particle pdg: " << par.PdgCode() << std::endl;
        std::cout << "Particle time: " << par.Trajectory().T(0) << std::endl;
        std::cout << "    converted: " << ts->G4ToElecTime(par.Trajectory().T(0)) - trig_time << std::endl;
        std::cout << "new Particle time: " << par.T() << std::endl;
        std::cout << "new    converted: " << ts->G4ToElecTime(par.T()) - trig_time << std::endl;
        std::cout << std::endl;
      }
      if (par.PdgCode() == 14) nuTime = par.T();//ts->G4ToElecTime(par.T()) - trig_time;
    }
  }

  std::cout << "Neutrino G4 interaction time: "  << nuTime << std::endl; 

  std::vector<std::vector<double> > pmt_v(1,std::vector<double>(geo->NOpDets(),0));

  */

  _n_ph_beam_spill = 0;

  std::cout << "JJJJ " << _beam_window_start << "  " << _beam_window_end << "  " << _before_spill_window << std::endl; 

  for(size_t opdet=0; opdet<32; ++opdet) {

    // Getting simulated photons from opdet
    sim::SimPhotons const& simph = (*evt_simphot_h)[opdet];

    if (_debug) std::cout << "Opdet " << opdet << std::endl;

    for(auto const& oneph : simph) {

      double t = ts->G4ToElecTime(oneph.Time) - trig_time;
      if (_debug) std::cout << "Photon G4 time is " << oneph.Time << "  while the TPC time is " << t << std::endl;

      if (t > _beam_window_start && t < _beam_window_end) {
        _n_ph_beam_spill += 1;
      }  

      if (t > _beam_window_start - _before_spill_window && t < _beam_window_start) {
        _n_ph_before_beam_spill += 1;
      }    

    }
  }

  _tree1->Fill();

  if (_debug) std::cout << "***** PhotonActivity ends." << std::endl;

}



//_______________________________________________________________________
void PhotonActivity::GetFlashLocation(std::vector<double> pePerOpChannel,
                                       double& Ycenter,
                                       double& Zcenter,
                                       double& Ywidth,
                                       double& Zwidth)
{

  // Reset variables
  Ycenter = Zcenter = 0.;
  Ywidth  = Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {

    // Get physical detector location for this opChannel
    double PMTxyz[3];
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);

    // Add up the position, weighting with PEs
    sumy    += pePerOpChannel[opch]*PMTxyz[1];
    sumy2   += pePerOpChannel[opch]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpChannel[opch]*PMTxyz[2];
    sumz2   += pePerOpChannel[opch]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpChannel[opch];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. )
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

  if ( (sumz2*totalPE - sumz*sumz) > 0. )
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}


DEFINE_ART_MODULE(PhotonActivity)
