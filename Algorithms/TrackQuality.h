/**
 * \file TrackQuality.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class TrackQuality
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef TRACKQUALITY_H
#define TRACKQUALITY_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace ubana{
  
  /**
   \class TrackQuality
   User defined class TrackQuality ... these comments are used to generate
   doxygen documentation!
 */

  class TrackQuality {
    
  public:
    
    /// Default constructor
    TrackQuality();

    /// Default destructor
    ~TrackQuality(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Printd the current configuration
    void PrintConfig();

    /// Set track points (projected to collection plane)
    void SetTrackPoints(std::vector<TVector3>);
    
    /// Set hits in collection plane assns to track
    void SetHitCollection(std::vector<TVector3>);

    /// Restores flags
    void Clear();

    /// Returns the residuals value (mean, std)
    std::pair<double, double> GetResiduals();

    ///
    std::pair<double,int> GetTrackGap(); 

    ///
    double GetR();

  
  protected:

    double mean(const std::vector<double>& data);
    double stdev(const std::vector<double>& data);
    double cov (const std::vector<double>& data1, const std::vector<double>& data2);

    std::vector<TVector3> _track_v; ///< The track
    std::vector<TVector3> _hit_v;   ///< The hits

    bool _track_is_set = false;
    bool _hits_are_set = false;

    bool _debug = false;

  };
}

#endif
/** @} */ // end of doxygen group 

