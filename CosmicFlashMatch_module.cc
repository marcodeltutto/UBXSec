////////////////////////////////////////////////////////////////////////
// Class:       CosmicFlashMatch
// Plugin Type: producer (art v2_05_00)
// File:        CosmicFlashMatch_module.cc
//
// Generated at Wed Jan 25 10:00:40 2017 by Marco Del Tutto using cetskelgen
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
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "uboone/UBXSec/FlashMatch.h" // new!
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"

#include "uboone/UBXSec/UBXSecHelper.h"

#include "TTree.h"

#include <memory>

class CosmicFlashMatch;


class CosmicFlashMatch : public art::EDProducer {
public:
  explicit CosmicFlashMatch(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicFlashMatch(CosmicFlashMatch const &) = delete;
  CosmicFlashMatch(CosmicFlashMatch &&) = delete;
  CosmicFlashMatch & operator = (CosmicFlashMatch const &) = delete;
  CosmicFlashMatch & operator = (CosmicFlashMatch &&) = delete;

  void GetTPCObjects(lar_pandora::PFParticleVector, lar_pandora::PFParticlesToTracks, lar_pandora::PFParticlesToVertices, std::vector<lar_pandora::PFParticleVector> &, std::vector<lar_pandora::TrackVector> &);
  flashana::QCluster_t GetQCluster(std::vector<art::Ptr<recob::Track>>);
  void CollectTracksAndPFP(lar_pandora::PFParticlesToTracks, lar_pandora::PFParticleVector, art::Ptr<recob::PFParticle>,  lar_pandora::PFParticleVector &, lar_pandora::TrackVector &);
  int  GetTrajectory(std::vector<art::Ptr<recob::Track>>, ::geoalgo::Trajectory &);
  flashana::Flash_t Trial(std::vector<art::Ptr<recob::Track>> track_v, flashana::Flash_t flashBeam, double & _chi2, double & _ll); 

  // Required functions.
  void produce(art::Event & e) override;

private:
  std::string _particleLabel;          ///<
  std::string _opflash_producer_beam;
  std::string _opflash_producer_cosmic;
  double _flash_trange_start;
  double _flash_trange_end;
  bool _debug;

  std::vector<::flashana::Flash_t>    beam_flashes;
  std::vector<::flashana::Flash_t>    cosmic_flashes;
  ::flashana::FlashMatchManager       _mgr;
  std::vector<flashana::FlashMatch_t> _result;

  std::vector<double>    _xfixed_hypo_spec;
  double _xfixed_chi2, _xfixed_ll;

  TTree* _tree1;
  int _run, _subrun, _event, _matchid, _flashid;
  std::vector<double>               _score, _t0;
  std::vector<double>               _qll_xmin, _tpc_xmin;
  std::vector<double>              _beam_flash_spec;
  std::vector<std::vector<double>> _hypo_flash_spec;
  std::vector<double>              _numc_flash_spec;
  int _fv, _ccnc, _nupdg;
};


CosmicFlashMatch::CosmicFlashMatch(fhicl::ParameterSet const & p)
{
  _particleLabel           = p.get<std::string>("PFParticleModule",      "pandoraNu");
  _debug                   = p.get<bool>       ("DebugMode",             true);
  _opflash_producer_beam   = p.get<std::string>("BeamOpFlashProducer",   "simpleFlashBeam");
  _opflash_producer_cosmic = p.get<std::string>("CosmicOpFlashProducer", "simpleFlashCosmic");
  _flash_trange_start      = p.get<double>     ("FlashVetoTimeStart",    -1000000);
  _flash_trange_end        = p.get<double>     ("FlashVetoTimeEnd",      1000000);
    
  _mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));

  if (_debug) {
    art::ServiceHandle<art::TFileService> fs;
    _tree1 = fs->make<TTree>("flashmatchtree","");
    _tree1->Branch("run",             &_run,                             "run/I");
    _tree1->Branch("subrun",          &_subrun,                          "subrun/I");
    _tree1->Branch("event",           &_event,                           "event/I");
    _tree1->Branch("beam_flash_spec", "std::vector<double>",             &_beam_flash_spec);
    _tree1->Branch("hypo_flash_spec", "std::vector<std::vector<double>>",&_hypo_flash_spec);
    _tree1->Branch("numc_flash_spec", "std::vector<double>",             &_numc_flash_spec);
    _tree1->Branch("score",           "std::vector<double>",             &_score);
    _tree1->Branch("t0",              "std::vector<double>",             &_t0);
    _tree1->Branch("qll_xmin",        "std::vector<double>",             &_qll_xmin);
    _tree1->Branch("tpc_xmin",        "std::vector<double>",             &_tpc_xmin);
    _tree1->Branch("xfixed_hypo_spec","std::vector<double>",             &_xfixed_hypo_spec);
    _tree1->Branch("fv",              &_fv,                              "fv/I");
    _tree1->Branch("ccnc",            &_ccnc,                            "ccnc/I");
    _tree1->Branch("nupdg",           &_nupdg,                           "nupdg/I");
  }


  produces< std::vector<ubana::FlashMatch>>();
  produces< art::Assns<ubana::FlashMatch,   recob::Track>>();
  produces< art::Assns<ubana::FlashMatch,   recob::PFParticle>>();
}

void CosmicFlashMatch::produce(art::Event & e)
{

  std::cout << "CosmicFlashMatch starts." << std::endl;

  // Instantiate the output
  std::unique_ptr< std::vector<ubana::FlashMatch>>                   flashMatchTrackVector      (new std::vector<ubana::FlashMatch>);
  std::unique_ptr< art::Assns<ubana::FlashMatch, recob::Track>>      assnOutFlashMatchTrack     (new art::Assns<ubana::FlashMatch, recob::Track>     );
  std::unique_ptr< art::Assns<ubana::FlashMatch, recob::PFParticle>> assnOutFlashMatchPFParticle(new art::Assns<ubana::FlashMatch, recob::PFParticle>);

  ::art::ServiceHandle<geo::Geometry> geo;

  _mgr.Reset();
  _result.clear();
  if(_debug) _mgr.PrintConfig();

  // Get Beam Flashes from the ART event
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cerr << "Don't have good beam flashes." << std::endl;
    e.put(std::move(flashMatchTrackVector));
    e.put(std::move(assnOutFlashMatchTrack));
    e.put(std::move(assnOutFlashMatchPFParticle));
    return;
  }
  int nBeamFlashes = 0;

  for (size_t n = 0; n < beamflash_h->size(); n++) {

    auto const& flash = (*beamflash_h)[n];

    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      continue;
    }
    nBeamFlashes++;

    // Construct a Flash_t
    ::flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }
    f.time = flash.Time();
    f.idx = nBeamFlashes-1;
    beam_flashes.resize(nBeamFlashes);
    beam_flashes[nBeamFlashes-1] = f;

    _mgr.Emplace(std::move(f));
  } // flash loop

  // Get Cosmic Flashes from the ART event
  ::art::Handle<std::vector<recob::OpFlash>> cosmicflash_h;
  e.getByLabel(_opflash_producer_cosmic,cosmicflash_h);
  if( !cosmicflash_h.isValid() || cosmicflash_h->empty() ) {
    std::cerr << "Don't have good cosmic flashes." << std::endl;
    e.put(std::move(flashMatchTrackVector));
    e.put(std::move(assnOutFlashMatchTrack));
    e.put(std::move(assnOutFlashMatchPFParticle));
    return;
  }
  int nCosmicFlashes = 0;

  for (size_t n = 0; n < cosmicflash_h->size(); n++) {

    auto const& flash = (*cosmicflash_h)[n];

    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      continue;
    }
    nCosmicFlashes++;

    // Construct a Flash_t
    ::flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }
    f.time = flash.Time();
    f.idx = nBeamFlashes + nCosmicFlashes - 1;
    cosmic_flashes.resize(nCosmicFlashes);
    cosmic_flashes[nCosmicFlashes-1] = f;

    _mgr.Emplace(std::move(f));

  } // flash loop


  if (nBeamFlashes == 0 && nCosmicFlashes == 0) {
    std::cout << "Zero beam and cosmic flashes in this event." << std::endl;
    e.put(std::move(flashMatchTrackVector));
    e.put(std::move(assnOutFlashMatchTrack));
    e.put(std::move(assnOutFlashMatchPFParticle));
    return;
  }


  //Vectors and maps we will use to store Pandora information
  lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::PFParticlesToClusters pfParticleToClusterMap; //PFParticle-to-cluster map

  //Use LArPandoraHelper functions to collect Pandora information
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _particleLabel, pfParticleList, pfParticleToClusterMap); //collect PFParticles and build map PFParticles to Clusters

  // Collect vertices and tracks
  lar_pandora::VertexVector           allPfParticleVertices;
  lar_pandora::PFParticlesToVertices  pfParticleToVertexMap;
  lar_pandora::LArPandoraHelper::CollectVertices(e, _particleLabel, allPfParticleVertices, pfParticleToVertexMap);
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _particleLabel, allPfParticleTracks, pfParticleToTrackMap);


  // ********************
  // Construct TPC Objects
  // ********************

  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;

  UBXSecHelper::GetTPCObjects(pfParticleList, pfParticleToTrackMap, pfParticleToVertexMap, pfp_v_v, track_v_v);

  if(_debug) std::cout << " For this event we have " << track_v_v.size() << " pandora slices." << std::endl;

  for (unsigned int tpcObj = 0; tpcObj < track_v_v.size(); tpcObj++) {

    lar_pandora::TrackVector tpcObjTrk_v = track_v_v[tpcObj];

    // Get QCluster for this TPC Object
    flashana::QCluster_t qcluster = this->GetQCluster(tpcObjTrk_v);

    qcluster.idx = tpcObj;

    // Emplace the QCluster to the FlashMatching Manager
    _mgr.Emplace(std::move(qcluster));
  }

  if(_debug) std::cout << "Finished emplacing beam flash and tpc objects" << std::endl;

  // ********************
  // Run Flash Matching
  // ********************

  _result = _mgr.Match();


  // ********************
  // Save the results
  // ********************
  
  if(_debug) {
    std::cout << "Number of matches: " << _result.size() << std::endl;
    _hypo_flash_spec.resize(_result.size());
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
  }

  _score.resize(_result.size());
  _t0.resize(_result.size()); 
  _qll_xmin.resize(_result.size());
  _tpc_xmin.resize(_result.size());

  for(_matchid=0; _matchid < (int)(_result.size()); ++_matchid) {

    auto const& match = _result[_matchid];

    _flashid         = match.flash_id;
    _score[_matchid] = match.score;

    auto const& flash = _mgr.FlashArray()[_flashid];
    _t0[_matchid] = flash.time;

    if(_debug) std::cout << "For this match, the score is " << match.score << ", the tpc id is " << match.tpc_id << std::endl;

    // Get the TPC obj 
    lar_pandora::TrackVector      track_v = track_v_v[match.tpc_id];
    lar_pandora::PFParticleVector pfp_v   = pfp_v_v[match.tpc_id];

    // Get hypo spec
    if(_debug){
      _hypo_flash_spec[_matchid].resize(geo->NOpDets());
      for(size_t pmt=0; pmt<_hypo_flash_spec[_matchid].size(); ++pmt) _hypo_flash_spec[_matchid][pmt] = match.hypothesis[pmt];
    }

    // Save x position
    _qll_xmin[_matchid] = match.tpc_point.x;

    _tpc_xmin[_matchid] = 1.e4;
    for(auto const& pt : _mgr.QClusterArray()[match.tpc_id]) {
      if(pt.x < _tpc_xmin[_matchid]) _tpc_xmin[_matchid] = pt.x;
    }

    ubana::FlashMatch fm;
    fm.SetScore               ( _score[_matchid] );
    fm.SetTPCX                ( _tpc_xmin[_matchid] );
    fm.SetEstimatedX          ( _qll_xmin[_matchid] );
    fm.SetT0                  ( _t0[_matchid] );
    fm.SetHypoFlashSpec       ( _hypo_flash_spec[_matchid] );
    fm.SetRecoFlashSpec       ( _beam_flash_spec );
    fm.SetMCFlashSpec         ( _numc_flash_spec );
    fm.SetXFixedHypoFlashSpec ( _xfixed_hypo_spec );
    fm.SetXFixedChi2          ( _xfixed_chi2 );
    fm.SetXFixedLl            ( _xfixed_ll );

    flashMatchTrackVector->emplace_back(std::move(fm));
    util::CreateAssn(*this, e, *flashMatchTrackVector, track_v, *assnOutFlashMatchTrack);
    util::CreateAssn(*this, e, *flashMatchTrackVector, pfp_v,   *assnOutFlashMatchPFParticle);
  }


  e.put(std::move(flashMatchTrackVector));
  e.put(std::move(assnOutFlashMatchTrack));
  e.put(std::move(assnOutFlashMatchPFParticle));


  if (_debug) _tree1->Fill();

  std::cout << "CosmicFlashMatch ends." << std::endl;
}






//______________________________________________________________________________________________________________________________________
flashana::QCluster_t CosmicFlashMatch::GetQCluster(std::vector<art::Ptr<recob::Track>> track_v) {

  flashana::QCluster_t summed_qcluster;
  summed_qcluster.clear();

  for (unsigned int trk = 0; trk < track_v.size(); trk++) {

    art::Ptr<recob::Track> trk_ptr = track_v.at(trk);

    ::geoalgo::Trajectory track_geotrj;
    track_geotrj.resize(trk_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));

    for (size_t pt_idx=0; pt_idx < trk_ptr->NumberTrajectoryPoints(); ++pt_idx) {
      auto const& pt = trk_ptr->LocationAtPoint(pt_idx);
      track_geotrj[pt_idx][0] = pt[0];
      track_geotrj[pt_idx][1] = pt[1];
      track_geotrj[pt_idx][2] = pt[2];
    }

    auto qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(track_geotrj);
    summed_qcluster += qcluster;

  } // track loop

  return summed_qcluster;
}



DEFINE_ART_MODULE(CosmicFlashMatch)
