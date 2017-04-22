/**
 * \file FindDeadRegions.h
 *
 * \ingroup 
 * 
 * \brief Class def header for a class FindDeadRegions
 *
 * @author Mike Mooney, Marco Del Tutto
 */

/** \addtogroup 

    @{*/

#ifndef FINDDEADREGIONS_H
#define FINDDEADREGIONS_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

struct BoundaryWire {
  unsigned int wire_num;
  float y_start;
  float z_start;
  float y_end;
  float z_end;
  bool isLowWire;
};

class FindDeadRegions;

class FindDeadRegions {
public:
  explicit FindDeadRegions();//fhicl::ParameterSet const & p, art::ActivityRegistry & areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  /// Configure function parameters
  void Configure(fhicl::ParameterSet const& p);

  /// Returns true if the passed point is close to a dead region given a tolerance considering two planes only
  bool NearDeadReg2P(float yVal, float zVal, float tolerance);

  /// Returns true if the passed point is close to a dead region given a tolerance considering all three planes
  bool NearDeadReg3P(float yVal, float zVal, float tolerance);

  /// Returns a root 2D histogram (y v.s. z) containing the detector dead regions considering two planes only
  TH2F* GetDeadRegionHisto2P();

  /// Returns a root 2D histogram (y v.s. z) containing the detector dead regions considering two planes only
  void GetDeadRegionHisto2P(TH2F *);

  /// Returns a root 2D histogram (y v.s. z) containing the detector dead regions considering all three planes
  TH2F* GetDeadRegionHisto3P();

  /// Returns a root 2D histogram (y v.s. z) containing the detector dead regions considering all three planes
  void GetDeadRegionHisto3P(TH2F *);

private:

  void LoadBWires();

  std::vector<BoundaryWire> BWires_U; ///< Contains list of wires marking the boundaries of dead regions (U plane)
  std::vector<BoundaryWire> BWires_V; ///< Contains list of wires marking the boundaries of dead regions (V plane)
  std::vector<BoundaryWire> BWires_Y; ///< Contains list of wires marking the boundaries of dead regions (Y plane)

  bool _use_file = true;  ///< If true, uses input files instad of geometry and database
  double _tolerance = 0.6; ///< Tolerance in cm to claim a point is in a dead region
  int _ch_thres = 4;       ///< Channels with status _less_ than threshold are considered as bad (only if using database)
};

#endif
/** @} */ // end of doxygen group 
