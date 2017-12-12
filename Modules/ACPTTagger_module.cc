////////////////////////////////////////////////////////////////////////
// Class:       ACPTTagger
// Plugin Type: producer (art v2_06_03)
// File:        ACPTTagger_module.cc
//
// Generated at Thu Jun  8 06:01:57 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v2_03_00.
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

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h" 

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "TVector3.h"
#include "TTree.h"
#include "TH1D.h"

#include <memory>
#include <fstream>


class ACPTTagger;


class ACPTTagger : public art::EDProducer {
public:
  explicit ACPTTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ACPTTagger(ACPTTagger const &) = delete;
  ACPTTagger(ACPTTagger &&) = delete;
  ACPTTagger & operator = (ACPTTagger const &) = delete;
  ACPTTagger & operator = (ACPTTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  /// Claculates the (track reco time) - (flash time) from the _end specified point, that is the closest to the _value given
  bool GetClosestDtDz(TVector3 _end, double _value, double trk_z_center, std::vector<double> &_dt, std::vector<double> &_dz, bool fill_histo);
  /// Returns true if sign is positive, negative otherwise
  bool GetSign(std::vector<TVector3> sorted_points);
  void SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_points);
  void SortSpacePoints(std::vector<art::Ptr<recob::SpacePoint>> sp_v, std::vector<TVector3>& sorted_points);
  void SortHitPoints(std::vector<art::Ptr<recob::Hit>> hit_v, std::vector<TVector3>& sorted_points, TVector3 highest_point);
  TVector3 ContainPoint(TVector3 p);
 
  std::string _flash_producer;
  std::string _pfp_producer;
  std::string _track_producer;
  std::string _spacepoint_producer;
  std::string _cluster_producer;
  std::string _swtrigger_producer;

  ::detinfo::DetectorProperties const* fDetectorProperties;
  ::art::ServiceHandle<geo::Geometry> geo;

  bool _use_tracks; ///< If true, uses tracks to get start and end points of a PFP
  bool _use_spacepoints; ///< If true, uses spacepoints to get start and end points of a PFP
  bool _use_hits; ///< If true, uses hits on the collection plane to get start x and end x points of a PFP 

  double _min_track_length;

  double _anodeTime;
  double _cathodeTime;

  double _dt_resolution_a, _dz_resolution_a;
  double _dt_resolution_c, _dz_resolution_c;

  // vector to hold flash-times for the event
  std::vector<double> _flash_times;
  std::vector<size_t> _flash_idx_v;
  std::vector<double> _flash_zcenter;
  std::vector<double> _flash_zwidth;
  double _pe_min;

  bool _debug, _create_histo;

  const std::vector<float> endPt1 = {-9999., -9999., -9999.};
  const std::vector<float> endPt2 = {-9999., -9999., -9999.};

  TH1D *_h_diff, *_h_diff_a, *_h_diff_c;
  TH1D *_h_dt_u_anode, *_h_dt_d_anode, *_h_dt_u_cathode, *_h_dt_d_cathode;
  TH1D *_h_dz_u_anode, *_h_dz_d_anode, *_h_dz_u_cathode, *_h_dz_d_cathode;

  //std::ofstream _csvfile;
  TTree* _tree1;
  int _run, _subrun, _event;
  bool _sw_trigger;
  double _drift_vel;
  std::vector<double> _trk_x_up, _trk_x_down, _trk_len, _trk_z_center;
  std::vector<double> _dt_u_anode, _dz_u_anode, _dt_d_anode, _dz_d_anode;
  std::vector<double> _dt_u_cathode, _dz_u_cathode, _dt_d_cathode, _dz_d_cathode;

  TTree* _tree2;
  double _tree2_dt;
  double _tree2_flstime;
  double _tree2_tracktime;
};


ACPTTagger::ACPTTagger(fhicl::ParameterSet const & p) {

  _flash_producer      = p.get<std::string>("FlashProducer", "simpleFlashCosmic");
  _pfp_producer        = p.get<std::string>("PFPartProducer", "pandoraCosmic");
  _track_producer      = p.get<std::string>("TrackProducer", "pandoraCosmic");
  _spacepoint_producer = p.get<std::string>("SpacePointProducer", "pandoraCosmic");
  _cluster_producer    = p.get<std::string>("ClusterProducer", "pandoraCosmic");
  _swtrigger_producer  = p.get<std::string>("SWTriggerProducer", "swtrigger");

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

  _use_tracks          = p.get<bool>("UseSpaceTracks", true);
  _use_spacepoints     = p.get<bool>("UseSpacePoints", false);
  _use_hits            = p.get<bool>("UseHits", false);

  _min_track_length    = p.get<double>("MinTrackLength", 0.0);

  _anodeTime           = p.get<double>("AnodeTime",   0.53);
  _cathodeTime         = p.get<double>("CathodeTime", 2291);

  _dt_resolution_a     = p.get<double>("DtResolutionAnode", 5); 
  _dz_resolution_a     = p.get<double>("DzResolutionAnode", 80);
  _dt_resolution_c     = p.get<double>("DtResolutionCathode", 5);
  _dz_resolution_c     = p.get<double>("DzResolutionCathode", 80);

  _pe_min              = p.get<double> ("PEMin", 0);
  _debug               = p.get<bool>("Debug", false);
  _create_histo        = p.get<bool>("CreateHisto", false); 

  //_csvfile.open ("acpt.csv", std::ofstream::out | std::ofstream::trunc);
  //_csvfile << "trk_x_up,trk_x_down,fls_time" << std::endl;

  art::ServiceHandle<art::TFileService> fs;

  if (_create_histo) {
    _h_diff   = fs->make<TH1D>("h_diff",  ";Reco time - Flash time;Events", 200, -6000, 6000);
    _h_diff_a = fs->make<TH1D>("h_diff_a",";Reco time - Flash time;Events", 100, -20, 20);
    _h_diff_c = fs->make<TH1D>("h_diff_c",";Reco time - Flash time;Events", 300, 2200, 3000);
 
    _h_dt_u_anode   = fs->make<TH1D>("h_dt_u_anode",   ";Reco time - Flash time [#mus];Events", 200, -20, 20);
    _h_dt_d_anode   = fs->make<TH1D>("h_dt_d_anode",   ";Reco time - Flash time [#mus];Events", 200, -20, 20);
    _h_dt_u_cathode = fs->make<TH1D>("h_dt_u_cathode", ";Reco time - Flash time [#mus];Events", 300, 2200, 3000);
    _h_dt_d_cathode = fs->make<TH1D>("h_dt_d_cathode", ";Reco time - Flash time [#mus];Events", 300, 2200, 3000);
    _h_dz_u_anode   = fs->make<TH1D>("h_dz_u_anode",   ";Track Z - Flash Z [cm];Events", 200, -200, 200);
    _h_dz_d_anode   = fs->make<TH1D>("h_dz_d_anode",   ";Track Z - Flash Z [cm];Events", 200, -200, 200);
    _h_dz_u_cathode = fs->make<TH1D>("h_dz_u_cathode", ";Track Z - Flash Z [cm];Events", 200, -200, 200);
    _h_dz_d_cathode = fs->make<TH1D>("h_dz_d_cathode", ";Track Z - Flash Z [cm];Events", 200, -200, 200);
  }

  if (_debug) {
    _tree1 = fs->make<TTree>("tree","");
    _tree1->Branch("run",           &_run,                 "run/I");
    _tree1->Branch("subrun",        &_subrun,              "subrun/I");
    _tree1->Branch("event",         &_event,               "event/I");
    _tree1->Branch("sw_trigger",    &_sw_trigger,          "sw_trigger/O");
    _tree1->Branch("drift_vel",     &_drift_vel,           "drift_vel/D");
    _tree1->Branch("trk_x_up",      "std::vector<double>", &_trk_x_up);
    _tree1->Branch("trk_x_down",    "std::vector<double>", &_trk_x_down);
    _tree1->Branch("trk_len",       "std::vector<double>", &_trk_len);
    _tree1->Branch("trk_z_center",  "std::vector<double>", &_trk_z_center);
    _tree1->Branch("flash_times",   "std::vector<double>", &_flash_times);
    _tree1->Branch("flash_zcenter", "std::vector<double>", &_flash_zcenter);
    _tree1->Branch("flash_zwidth",  "std::vector<double>", &_flash_zwidth);
  
    _tree1->Branch("dt_u_anode",    "std::vector<double>", &_dt_u_anode);
    _tree1->Branch("dz_u_anode",    "std::vector<double>", &_dz_u_anode);
    _tree1->Branch("dt_d_anode",    "std::vector<double>", &_dt_d_anode);
    _tree1->Branch("dz_d_anode",    "std::vector<double>", &_dz_d_anode);
    _tree1->Branch("dt_u_cathode",  "std::vector<double>", &_dt_u_cathode);
    _tree1->Branch("dz_u_cathode",  "std::vector<double>", &_dz_u_cathode);
    _tree1->Branch("dt_d_cathode",  "std::vector<double>", &_dt_d_cathode);
    _tree1->Branch("dz_d_cathode",  "std::vector<double>", &_dz_d_cathode);
 
    _tree2 = fs->make<TTree>("tree2","");
    _tree1->Branch("run",             &_run,                 "run/I");
    _tree1->Branch("subrun",          &_subrun,              "subrun/I");
    _tree1->Branch("event",           &_event,               "event/I");
    _tree1->Branch("tree2_dt",        &_tree2_dt,            "tree2_dt/D");
    _tree1->Branch("tree2_flstime",   &_tree2_flstime,       "tree2_flstime/D");
    _tree1->Branch("tree2_tracktime", &_tree2_tracktime,     "tree2_tracktime/D");
  }

  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<anab::CosmicTag,   recob::Track>>();
  produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

void ACPTTagger::produce(art::Event & e)
{

  if (_debug) {
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
   
    std::cout << "AnodeTime set to " << _anodeTime << std::endl;
    std::cout << "CathodeTime set to " << _cathodeTime << std::endl;
  }

  // Instantiate the output
  std::unique_ptr< std::vector< anab::CosmicTag>>                  cosmicTagTrackVector      (new std::vector<anab::CosmicTag>);
  std::unique_ptr< art::Assns<anab::CosmicTag, recob::Track>>      assnOutCosmicTagTrack     (new art::Assns<anab::CosmicTag, recob::Track>);
  std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag>> assnOutCosmicTagPFParticle(new art::Assns<recob::PFParticle, anab::CosmicTag>);

  // sw trigger
  art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
  e.getByLabel(_swtrigger_producer, softwareTriggerHandle);

  if (!softwareTriggerHandle.isValid() || softwareTriggerHandle.failedToGet()){
    //std::cerr << "Failed to get software trigget data product with label " << _swtrigger_producer << std::endl;
  } else {
    std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
    size_t trigger = 0;
    _sw_trigger = softwareTriggerHandle->passedAlgo(algoNames[trigger]);
  }

  // load Flash
  if (_debug) { std::cout << "Loading flashes from producer " << _flash_producer << std::endl; }
  art::Handle<std::vector<recob::OpFlash> > flash_h;
  e.getByLabel(_flash_producer,flash_h);
  if (_debug) { std::cout << "Initially we have " << flash_h->size() << " flashes." << std::endl; }

  // make sure flash look good
  if(!flash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    e.put(std::move(cosmicTagTrackVector));
    e.put(std::move(assnOutCosmicTagTrack));
    e.put(std::move(assnOutCosmicTagPFParticle));
    return; //throw std::exception();
  }

  // load PFParticles for which T0 reconstruction should occur
  if (_debug) { std::cout << "Loading PFParticles from producer " << _pfp_producer << std::endl; }
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);

  // make sure pfparticles look good
  if(!pfp_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate PFParticles!"<<std::endl;
    e.put(std::move(cosmicTagTrackVector));
    e.put(std::move(assnOutCosmicTagTrack));
    e.put(std::move(assnOutCosmicTagPFParticle));
    return; //throw std::exception();
  }

  std::vector<art::Ptr<recob::PFParticle> > PFPVec;
  art::fill_ptr_vector(PFPVec, pfp_h);

  // grab tracks associated with PFParticles
  art::FindManyP<recob::Track> pfp_track_assn_v(pfp_h, e, _track_producer);
  if (_debug) {
    std::cout << "There are " << pfp_track_assn_v.size() << " pfpart -> track associations" << std::endl;
  }

  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, _spacepoint_producer);
  if (_debug) {
    std::cout << "There are " << pfp_spacepoint_assn_v.size() << " pfpart -> spacepoint associations" << std::endl;
  }

  // grab hits associated with PFParticles
  lar_pandora::PFParticleVector particleVector;
  lar_pandora::PFParticlesToClusters particlesToClusters;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, particleVector, particlesToClusters);
  
  lar_pandora::ClusterVector clusterVector;
  lar_pandora::ClustersToHits clustersToHits;
  lar_pandora::LArPandoraHelper::CollectClusters(e, _cluster_producer, clusterVector, clustersToHits);


  // prepare a vector of optical flash times, if flash above some PE cut value

  _flash_times.clear();
  _flash_idx_v.clear();
  _flash_zcenter.clear();
  _flash_zwidth.clear();

  size_t flash_ctr = 0;
  for (auto const& flash : *flash_h){
    if (flash.TotalPE() > _pe_min){
      _flash_times.push_back( flash.Time() );
      _flash_idx_v.push_back(flash_ctr);
      _flash_zcenter.push_back(flash.ZCenter());
      _flash_zwidth.push_back(flash.ZWidth());
      if (_debug) { 
        std::cout << "\t flash time : " << flash.Time() << ", PE : " << flash.TotalPE() << ", ZCenter : " << flash.ZCenter() << " +- " << flash.ZWidth() << std::endl; 
      }
    }
    flash_ctr += 1;
  } // for all flashes

  if (_debug) { 
    std::cout << __PRETTY_FUNCTION__ << " Selected a total of " << _flash_times.size() << " OpFlashes" << std::endl; 
  }

  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield     = _detp->Efield();
  double temp       = _detp->Temperature();
  _drift_vel        = _detp->DriftVelocity(efield,temp);

  _trk_len.clear();
  _trk_x_up.clear();
  _trk_x_down.clear();

  _dt_u_anode.clear();
  _dz_u_anode.clear();
  _dt_d_anode.clear();
  _dz_d_anode.clear();
  _dt_u_cathode.clear();
  _dz_u_cathode.clear();
  _dt_d_cathode.clear();
  _dz_d_cathode.clear();

  for (size_t i = 0; i < PFPVec.size(); i++) {

    if (_debug) {
      std::cout << "[ACPTTagger] Looping through pfpart number " << i << std::endl;
    }

    bool isCosmic = false;
    auto pfp = PFPVec.at(i);

    // grab associated tracks
    std::vector<art::Ptr<recob::Track>> track_v = pfp_track_assn_v.at(i);

    // grab associated spacepoints
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_v = pfp_spacepoint_assn_v.at(i);

    // grab associated hits
    std::vector<art::Ptr<recob::Hit>> hit_v;
    auto iter = particlesToClusters.find(PFPVec.at(i));
    if (iter != particlesToClusters.end()) {
      for (auto c : iter->second) {
        auto iter2 = clustersToHits.find(c);
        if (iter2 != clustersToHits.end()) {
          hit_v.reserve(hit_v.size() + iter2->second.size()); 
          hit_v.insert(hit_v.end(), iter2->second.begin(), iter2->second.end());
        }
      }
    }

    // Will store sorted points for the object [assuming downwards going]
    // The first vector is two consider end points estimated via different methods
    // (spacepoints, hits, tracks). The second vector has length==2, and
    // contains the start and end point of the track
    std::vector<std::vector<TVector3>> sorted_points_v;
    sorted_points_v.clear();

    if (_use_spacepoints) {
      if(_debug) std::cout << "Using Spacepoints" << std::endl;
      if(spacepoint_v.size() >= 2) {
        std::vector<TVector3> pts;
        this->SortSpacePoints(spacepoint_v, pts);
        sorted_points_v.emplace_back(pts);
      }
    }

    if (_use_hits) {
      if(_debug) std::cout << "Using Hits" << std::endl;
      std::vector<TVector3> pts;
      this->SortSpacePoints(spacepoint_v, pts);
      if (pts.size() > 0) {
        TVector3 hp = this->ContainPoint(pts[0]);
        this->SortHitPoints(hit_v, pts, hp);
        if (pts.size() >= 2) {
          sorted_points_v.emplace_back(pts);
        }
      }
    }

    if (_use_tracks) {
      if(_debug) std::cout << "Using Tracks" << std::endl;
      for (auto t : track_v) {
        if (t->Length() < _min_track_length) continue;   
        std::vector<TVector3> pts;
        this->SortTrackPoints(*t, pts);
        if (pts.size() >= 2) {
          _trk_len.emplace_back(t->Length());
          _trk_x_up.emplace_back(pts[0].X());
          _trk_x_down.emplace_back(pts[pts.size()-1].X());
          double z_center = pts[0].Z();
          z_center += pts[pts.size()-1].Z();
          z_center /= 2.;
          _trk_z_center.emplace_back(z_center);

          TVector3 start = pts[0];

          TVector3 end = pts[pts.size()-1];
          pts.resize(2);
          pts.at(0) = start;
          pts.at(1) = end;
          sorted_points_v.emplace_back(pts);
        } 
      }
    }
      

    size_t loop_max = sorted_points_v.size();

    if(_debug) std::cout << "The loop will be " << loop_max << std::endl;

    for (size_t i = 0; i < loop_max; i ++) {

      if(_debug) std::cout << "At loop stage " << i << std::endl;

      std::vector<TVector3> sorted_points = sorted_points_v.at(i);
     
      double z_center = sorted_points[0].Z();
      z_center += sorted_points[sorted_points.size()-1].Z();
      z_center /= 2.;

      this->GetClosestDtDz(sorted_points[0],                      _anodeTime,   z_center, _dt_u_anode,   _dz_u_anode,   true);
      this->GetClosestDtDz(sorted_points[sorted_points.size()-1], _anodeTime,   z_center, _dt_d_anode,   _dz_d_anode,   true);
      this->GetClosestDtDz(sorted_points[0],                      _cathodeTime, z_center, _dt_u_cathode, _dz_u_cathode, false);
      this->GetClosestDtDz(sorted_points[sorted_points.size()-1], _cathodeTime, z_center, _dt_d_cathode, _dz_d_cathode, false);

      bool sign = this->GetSign(sorted_points);

      // A
      if (_dt_u_anode.back() > _anodeTime - _dt_resolution_a && _dt_u_anode.back() < _anodeTime + _dt_resolution_a 
       && _dz_u_anode.back() > -_dz_resolution_a && _dz_u_anode.back() < _dz_resolution_a
       && sign) isCosmic = true;

      // B
      if (_dt_d_anode.back() > _anodeTime - _dt_resolution_a && _dt_d_anode.back() < _anodeTime + _dt_resolution_a 
       && _dz_d_anode.back() > -_dz_resolution_a && _dz_d_anode.back() < _dz_resolution_a
       && !sign) isCosmic = true;

      // C
      if (_dt_u_cathode.back() > _cathodeTime - _dt_resolution_c  && _dt_u_cathode.back() < _cathodeTime + _dt_resolution_c 
       && _dz_u_cathode.back() > -_dz_resolution_a && _dz_u_cathode.back() < _dz_resolution_a
       && !sign) isCosmic = true;

      // D
      if (_dt_d_cathode.back() > _cathodeTime - _dt_resolution_c && _dt_d_cathode.back() < _cathodeTime + _dt_resolution_c 
       && _dz_d_cathode.back() > -_dz_resolution_a && _dz_d_cathode.back() < _dz_resolution_a
       && sign) isCosmic = true;

      if (_create_histo) {
        _h_dt_u_anode->Fill(_dt_u_anode.back());
        _h_dt_d_anode->Fill(_dt_d_anode.back());
        _h_dt_u_cathode->Fill(_dt_u_cathode.back());
        _h_dt_d_cathode->Fill(_dt_d_cathode.back());
        
        _h_dz_u_anode->Fill(_dz_u_anode.back());
        _h_dz_d_anode->Fill(_dz_d_anode.back());
        _h_dz_u_cathode->Fill(_dz_u_cathode.back());
        _h_dz_d_cathode->Fill(_dz_d_cathode.back());
      }

      if (_debug) {

        std::cout << "_dt_u_anode " << _dt_u_anode.back() << std::endl;
        std::cout << "_dt_d_anode " << _dt_d_anode.back() << std::endl;
        std::cout << "_dt_u_cathode " << _dt_u_cathode.back() << std::endl;
        std::cout << "_dt_d_cathode " << _dt_d_cathode.back() << std::endl;
        std::cout << "_dz_u_anode " << _dz_u_anode.back() << std::endl;
        std::cout << "_dz_d_anode " << _dz_d_anode.back() << std::endl;
        std::cout << "_dz_u_cathode " << _dz_u_cathode.back() << std::endl;
        std::cout << "_dz_d_cathode " << _dz_d_cathode.back() << std::endl << std::endl;
        if (isCosmic) std::cout << "Tagged!" << std::endl;
      }
      
    } // Points loop

    float cosmicScore = 0;
    if (isCosmic) {
      cosmicScore = 1;
    }
     
    cosmicTagTrackVector->emplace_back(endPt1, endPt2, cosmicScore, anab::CosmicTagID_t::kGeometry_XY);
    util::CreateAssn(*this, e, *cosmicTagTrackVector, track_v, *assnOutCosmicTagTrack );
    util::CreateAssn(*this, e, *cosmicTagTrackVector, pfp, *assnOutCosmicTagPFParticle); 
    
  } // PFP loop

  e.put(std::move(cosmicTagTrackVector));
  e.put(std::move(assnOutCosmicTagTrack));
  e.put(std::move(assnOutCosmicTagPFParticle));

  if (_debug) _tree1->Fill();
}


bool ACPTTagger::GetClosestDtDz(TVector3 _end, double _value, double trk_z_center, std::vector<double> &_dt, std::vector<double> &_dz, bool fill_histo) {

  double dist = 1e9;
  double min_dist = 1e9;
  double min_diff = 1e9;
  double theflash = -1;

  for (size_t f = 0; f < _flash_times.size(); f++) {

    double diff = _end.X() / _drift_vel - _flash_times[f];

    dist = abs(diff - _value);

    if (dist < min_dist) {
      min_dist = dist;
      min_diff = diff;
      theflash = f;
    }

    if (fill_histo && _create_histo){
      _h_diff->Fill(diff);
      _h_diff_a->Fill(diff);
      _h_diff_c->Fill(diff);
    }
    if (_debug && fill_histo) {
      _tree2_dt = diff;
      _tree2_flstime = _flash_times[theflash];
      _tree2_tracktime = _end.X() / _drift_vel;
      _tree2->Fill();
    }

  } // flash loop

  if (theflash == -1) return false;

  std::cout << "flash index " << theflash << ", flashtime " << _flash_times.at(theflash) << ", flashcenter " << _flash_zcenter.at(theflash) << ", trackcenter " << trk_z_center << std::endl;
 
  _dt.emplace_back(min_diff);
  _dz.emplace_back(trk_z_center - _flash_zcenter[theflash]); 

  return true;

}

void ACPTTagger::SortSpacePoints(std::vector<art::Ptr<recob::SpacePoint>> sp_v, std::vector<TVector3>& sorted_points) {

  sorted_points.clear();

  if (sp_v.size() < 2) 
    return;

  // Sort SpacePoints by y position
  std::sort(sp_v.begin(), sp_v.end(),
            [](art::Ptr<recob::SpacePoint> a, art::Ptr<recob::SpacePoint> b) -> bool
            {
              return a->XYZ()[1] > b->XYZ()[1];
            });

  // Just save start and end point
  sorted_points.resize(2);
  TVector3 pt1(sp_v.at(0)->XYZ()[0], 
               sp_v.at(0)->XYZ()[1], 
               sp_v.at(0)->XYZ()[2]);
  TVector3 pt2(sp_v.at(sp_v.size()-1)->XYZ()[0], 
               sp_v.at(sp_v.size()-1)->XYZ()[1], 
               sp_v.at(sp_v.size()-1)->XYZ()[2]);
  sorted_points.at(0) = std::move(pt1);
  sorted_points.at(1) = std::move(pt2);
  

}



void ACPTTagger::SortHitPoints(std::vector<art::Ptr<recob::Hit>> hit_v, std::vector<TVector3>& sorted_points, TVector3 highest_point) {

  sorted_points.clear();
  if (hit_v.size() < 2)
    return;
  // Only collection plane hits
  std::vector<art::Ptr<recob::Hit>> temp_v;
  temp_v.clear();
  for (size_t i = 0; i < hit_v.size(); i++) {
    if (hit_v.at(i)->View() == 2) {
      temp_v.push_back(hit_v.at(i));
    }
  }
  std::swap(temp_v, hit_v);
  if (hit_v.size() < 2)
    return;
  for (size_t i = 0; i < hit_v.size(); i++) {

  }

  // Sort Hit by x position
  std::sort(hit_v.begin(), hit_v.end(),
            [](art::Ptr<recob::Hit> a, art::Ptr<recob::Hit> b) -> bool
            {
              return a->PeakTime() > b->PeakTime();
            });

  double time_highest = fDetectorProperties->ConvertXToTicks(highest_point.X(), geo::PlaneID(0,0,2));
  int wire_highest = geo->NearestWire(highest_point, 2);

  std::cout << "[ACPTTagger] wire_highest " << wire_highest << ", time_highest " << time_highest << std::endl;

  TVector3 pt0 (wire_highest,                             time_highest,                         0);
  TVector3 pt_1 (hit_v.at(0)->WireID().Wire,              hit_v.at(0)->PeakTime(),              0);
  TVector3 pt_2 (hit_v.at(hit_v.size()-1)->WireID().Wire, hit_v.at(hit_v.size()-1)->PeakTime(), 0);

  if ( (pt0-pt_1).Mag() > (pt0-pt_2).Mag()) {
    auto temp = hit_v.at(0);
    hit_v.at(0) = hit_v.at(hit_v.size()-1);
    hit_v.at(hit_v.size()-1) = temp;
  }

  std::cout << "first pt wire " << hit_v.at(0)->WireID().Wire << ", time " << hit_v.at(0)->PeakTime() << ", which is x " << fDetectorProperties->ConvertTicksToX(hit_v.at(0)->PeakTime(), geo::PlaneID(0,0,2)) << std::endl;   
  std::cout << "second pt wire " << hit_v.at(hit_v.size()-1)->WireID().Wire << ", time " << hit_v.at(hit_v.size()-1)->PeakTime() << ", which is x " << fDetectorProperties->ConvertTicksToX(hit_v.at(hit_v.size()-1)->PeakTime(), geo::PlaneID(0,0,2)) << std::endl; 
     
  // Just save start and end point
  sorted_points.resize(2);
  TVector3 pt1(fDetectorProperties->ConvertTicksToX(hit_v.at(0)->PeakTime(), geo::PlaneID(0,0,2)),
               0.,
               hit_v.at(0)->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,2)));
  TVector3 pt2(fDetectorProperties->ConvertTicksToX(hit_v.at(hit_v.size()-1)->PeakTime(), geo::PlaneID(0,0,2)),
               0.,
               hit_v.at(hit_v.size()-1)->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,2)));
  sorted_points.at(0) = std::move(pt1);
  sorted_points.at(1) = std::move(pt2);


}



void ACPTTagger::SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_points) {

  // vector to store 3D coordinates of                                                                                                                                           
  // ordered track                              
  sorted_points.clear();

  // take the reconstructed 3D track                                                                                                                                           
  // and assuming it is downwards                                                                                                    
  // going, sort points so that                                                                                                              
  // the track starts at the top                                                                                                     
  // which point is further up in Y coord?                                                                                                                  
  // start or end?                                                                                                                 
  auto const&N = track.NumberTrajectoryPoints();
  auto const&start = track.LocationAtPoint(0);
  auto const&end   = track.LocationAtPoint( N - 1 );

  if (_debug) {
    std::cout << "[ACPTTagger] Track start " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
    std::cout << "[ACPTTagger] Track start " << end.X() << " " << end.Y() << " " << end.Z() << std::endl;
  }

  // if points are ordered correctly                                                                                                                                       
  if (start.Y() > end.Y()){
    for (size_t i=0; i < N; i++)
      sorted_points.push_back( track.LocationAtPoint(i) );
  }

  // otherwise flip order                                                                                                                                                 
  else {
    for (size_t i=0; i < N; i++)
      sorted_points.push_back( track.LocationAtPoint( N - i - 1) );
  }
}

bool ACPTTagger::GetSign(std::vector<TVector3> sorted_points)
{

  double t_down = sorted_points[sorted_points.size()-1].X();
  double t_up = sorted_points[0].X();

  bool is_positive = (t_down - t_up) > 0.;

  return is_positive;

}

TVector3 ACPTTagger::ContainPoint(TVector3 p) {

  TVector3 p_out = p;

  double x = p.X();
  double y = p.Y();
  double z = p.Z();

  double e = std::numeric_limits<double>::epsilon();

  if (x < 0. + e)
    p_out.SetX(0. + e);
  if (x > 2.*geo->DetHalfWidth() - e)
    p_out.SetX(2.*geo->DetHalfWidth() - e);

  if (y < -geo->DetHalfWidth() + e)
    p_out.SetY(-geo->DetHalfWidth() + e);
  if (y > geo->DetHalfWidth() - e)
    p_out.SetY(geo->DetHalfWidth() - e);

  if (z < 0. + e)
    p_out.SetZ(0.+ e);
  if (z > geo->DetLength() - e)
    p_out.SetZ(geo->DetLength() - e);

  return p_out;
}

DEFINE_ART_MODULE(ACPTTagger)
