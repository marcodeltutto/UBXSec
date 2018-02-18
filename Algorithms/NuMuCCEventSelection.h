/**
 * \file NuMuCCEventSelection.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class NuMuCCEventSelection
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef NUMUCCEVENTSELECTION_H
#define NUMUCCEVENTSELECTION_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"

namespace ubana{
  
  /**
   \class NuMuCCEventSelection
   User defined class NuMuCCEventSelection ... these comments are used to generate
   doxygen documentation!
 */

  class NuMuCCEventSelection {
    
  public:
    
    /// Default constructor
    NuMuCCEventSelection();

    /// Default destructor
    ~NuMuCCEventSelection(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Print the current configuration
    void PrintConfig();

    /// Set the UBXSecEvent
    void SetEvent(UBXSecEvent*);

    /// Returns true if this event is selected
    bool IsSelected(size_t & slice_index, std::map<std::string,bool> & failure_map);

  protected:

    UBXSecEvent* _ubxsec_event;
    bool _event_is_set = false;
    bool _configured = false;

    bool _verbose;

    double _deltax_cut_down;
    double _deltax_cut_up;
    double _deltaz_cut_down;
    double _deltaz_cut_up;
    double _flsmatch_score_cut;
    double _vtxcheck_angle_cut_down;
    double _vtxcheck_angle_cut_up;
    double _mcs_length_cut;
    double _ntrack_cut;
    double _residuals_std_down_cut;
    double _residuals_std_up_cut;
    double _residuals_mean_down_cut;
    double _residuals_mean_up_cut;
    double _perc_used_hits_in_cluster_cut;

    double _pe_cut;
    double _beamSpillStarts;
    double _beamSpillEnds;
  };
}

#endif
/** @} */ // end of doxygen group 

