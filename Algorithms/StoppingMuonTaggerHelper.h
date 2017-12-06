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
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


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

    /// Sets the Simple Hit Vector
    void Emplace(ubana::SimpleHitVector v) {_s_hit_v = v;}
    
    /// Filter hits keeping only in selected plane
    size_t FilterByPlane(int);

    /// Sets the start hit given time, wire and plane
    void SetStartHit(double t, double w, int p);

    /// Orde hits from start hit (has to be set)
    void OrderHits();

    /// Calulates the dQds given ordered hits
    void CalculatedQds();

    /// Calculated the dqds averaging neibouring hits
    void PerformdQdsSlider();

    /// Returns true if this object has been id'ed as a Stopping Muon
    bool MakeDecision();

    /// Restores flags
    void Clear();

    /// Set the conversion from wire number to cm
    void SetW2Cm(double x) {_w2cm = x;}

    ///Set the conversion from time to cm
    void SetT2Cm(double x) {_t2cm = x;}

    /// Sets the max distance between consecutive hits
    void SetMaxAllowedHitDistance(double x) {_max_allowed_hit_distance = x;}

    /// Removes max and min value and returns the median
    double GetTruncMedian(std::vector<double> v);

  protected:

    ubana::SimpleHitVector _s_hit_v;
    int _start_index = -1;
    std::vector<double> _dqds_v;
    std::vector<double> _ds_v;
    std::vector<double> _dqds_slider;

    double _w2cm = 0.3;    // to cm 
    double _t2cm = 0.0557; // to cm
    double _dqds_calib = 198.;

    double _max_allowed_hit_distance = 15.;

    bool _hits_ordered = false;
  };
}

#endif
/** @} */ // end of doxygen group 

