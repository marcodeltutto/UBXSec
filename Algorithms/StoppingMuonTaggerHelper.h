/**
 * \file StoppingMuonTaggerHelper.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class StoppingMuonTaggerHelper
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef STOPPINGMUONTAGGERHELPER_H
#define STOPPINGMUONTAGGERHELPER_H

#include <iostream>
#include <fstream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

namespace ubana {

  /// Struct to represent a recob Hit
  struct SimpleHit {

    public: 

    double time;     ///< in cm, basically x
    int wire;        ///< in cm, given the wire pitch
    int plane;       ///< plane number
    double integral; ///< integral (~charge)

    double t;        ///< in time ticks
    int w;           ///< in wire number

    SimpleHit () {
      time = t = -1;
      wire = w = -1;
      plane = -1;
    }

    bool operator==(const SimpleHit& x) const {
      return (time == x.time) && 
             (wire == x.wire) &&
             (plane == x.plane);
    }
  };

  typedef std::vector<ubana::SimpleHit> SimpleHitVector;

  enum StopMuAlgoType{
    kAlgoUnknown = -1,
    kAlgoMichel = 0,
    kAlgoBragg = 1,
    kAlgoCurvature = 2,
    kAlgoSimpleMIP = 3,
  };

}


namespace ubana{
  
  /**
   \class StoppingMuonTaggerHelper
   User defined class StoppingMuonTaggerHelper ... these comments are used to generate
   doxygen documentation!
 */

  class StoppingMuonTaggerHelper {
    
  public:
    
    /// Default constructor
    StoppingMuonTaggerHelper();

    /// Default destructor
    ~StoppingMuonTaggerHelper(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Printd the current configuration
    void PrintConfig();

    ///
    void PrintOnFile(int index);

    /// Sets the Simple Hit Vector
    void Emplace(ubana::SimpleHitVector v) {_s_hit_v = v;}
    
    /// Filter hits keeping only in selected plane
    size_t FilterByPlane(int);

    ///
    size_t FilterOnSingleWire();

    /// Sets the start hit given time, wire and plane
    void SetStartHit(double t, double w, int p);

    /// Orde hits from start hit (has to be set)
    size_t OrderHits();

    /// Calulates the dQds given ordered hits
    void CalculatedQds();

    /// Calculated the dqds averaging neibouring hits
    bool PerformdQdsSlider();

    /// Calculates the local linearity along the hits
    void CalculateLocalLinearity();

    ///
    double cov (const std::vector<double>& data1,
                const std::vector<double>& data2) const;

    ///
    double stdev(const std::vector<double>& data) const;

    ///
    double mean (const std::vector<double>& data) const;

    /// Returns true if this object has been id'ed as a Stopping Muon given an algo
    bool MakeDecision(ubana::StopMuAlgoType algo);

    /// Algo, returns true if is a stopping muon decaying to Michel
    bool IsStopMuMichel();

    /// Algo, returns true if is a stopping muon, looks at the Bragg peak
    bool IsStopMuBragg();

    /// Algo, returns true if is is a crossing simple MIP
    bool IsSimpleMIP();

    /// Algo, returns true if is a stopping muon, looks at the curvature
    bool IsStopMuCurvature(); 

    /// Restores flags
    void Clear();

    /// Set the conversion from wire number to cm
    void SetW2Cm(double x) {_w2cm = x;}

    ///Set the conversion from time to cm
    void SetT2Cm(double x) {_t2cm = x;}

    /// Sets the max distance between consecutive hits
    void SetMaxAllowedHitDistance(double x) {_max_allowed_hit_distance = x;}

    /// Sets the vertex, or highest point
    void SetVertexPoint(double *p) {_vertex.SetX(p[0]); _vertex.SetY(p[1]); _vertex.SetZ(p[2]);}

    /// Sets the FV object
    void SetFiducialVolume(ubana::FiducialVolume fv) {_fv = fv;}

    /// Removes max and min value and returns the median
    double GetTruncMedian(std::vector<double> v);

    /// Returns the curvature of a circle passing trough 3 points (2D, z value of points has to be zero)
    double GetCurvature(TVector3 pt1, TVector3 pt2, TVector3 pt3);
      
    /// Returns selected elements in window
    template<typename T>
    std::vector<std::vector<T> > get_windows(const std::vector<T>& the_thing,
                                             const size_t window_size) const;
  protected:

    ubana::SimpleHitVector _s_hit_v;
    int _start_index = -1;
    std::vector<double> _dqds_v;
    std::vector<double> _ds_v;
    std::vector<double> _dqds_slider;
    std::vector<double> _linearity_v;
    bool _linearity_is_set = false;

    ubana::FiducialVolume _fv;
    TVector3 _vertex;

    bool _hits_ordered = false;

    double _w2cm = 0.3;    // to cm 
    double _t2cm = 0.0557; // to cm
    double _dqds_calib = 198.;

    size_t _slider_window = 10;
    double _max_allowed_hit_distance = 15.;
    double _slope_threshold = 0.15;
    int _hits_to_remove = 3;
    int _pre_post_window = 5;
    double _perc_diff_cut = 20;
    double _local_linearity_threshold = 0.85;
    int _min_muon_hits = 20;
    int _max_muon_hits = 2000;
    int _min_michel_hits = 2;
    int _max_michel_hits = 50;
    int _max_end_hits = 120;

    bool _debug = false;

    std::ofstream _csvfile; 
  };
}

#endif
/** @} */ // end of doxygen group 

