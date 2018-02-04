#ifndef TRACKQUALITY_CXX
#define TRACKQUALITY_CXX

#include "TrackQuality.h"
#include <iostream>

namespace ubana {

  TrackQuality::TrackQuality()
  {
    _track_is_set = false;
    _hits_are_set = false;
  }


  void TrackQuality::Configure(fhicl::ParameterSet const& pset)
  {
    //_max_distance   = pset.get< double > ( "MaxDistance" );
  }

  void TrackQuality::PrintConfig() 
  {

    std::cout << "--- TrackQuality configuration:" << std::endl;
    //std::cout << "---   _max_distance  = " << _max_distance << std::endl;

  }

  void TrackQuality::SetTrackPoints(std::vector<TVector3> v) 
  {
    _track_v = v;
    _track_is_set = true;
  }

  void TrackQuality::SetHitCollection(std::vector<TVector3> v) 
  {
    _hit_v = v;
    _hits_are_set = true;
  }

  void TrackQuality::Clear()
  {
    _track_is_set = _hits_are_set = false;
  }

  double TrackQuality::GetR()
  {

    std::vector<double> X;
    std::vector<double> Y;
    X.reserve(_hit_v.size());
    Y.reserve(_hit_v.size());

    for(auto hit : _hit_v) {
      X.push_back(hit.X()); 
      Y.push_back(hit.Y());
    }

    auto c  = cov(X,Y);
    auto sX = stdev(X);
    auto sY = stdev(Y);
    auto r  = std::abs(c/(sX * sY));

    return r;

  }

  std::pair<double,int> TrackQuality::GetTrackGap() 
  {

    std::vector<std::pair<double,int>> distance_to_previous_wire;

    for (size_t i = 1; i < _track_v.size(); i++) {

      double distance = _track_v.at(i).X() - _track_v.at(i-1).X();
      distance_to_previous_wire.emplace_back(std::make_pair(distance, _track_v.at(i).X()));

    }

    if (distance_to_previous_wire.size() == 0) {
      if (_debug) std::cout << "[TrackQuality] is zero." << std::endl;
      return std::make_pair(-9999, -9999);
    }

    for (auto p : distance_to_previous_wire) {
      if (_debug) std::cout << "[TrackQuality] Distance: " << p.first << ", wire at end: " << p.second << std::endl;
    }

    // Order based on wire distance
    std::sort(distance_to_previous_wire.begin(), distance_to_previous_wire.end(),
      [](std::pair<double,int> a, std::pair<double,int> b) -> bool
      {
        return std::abs(a.first) > std::abs(b.first);
      });

    return distance_to_previous_wire.at(0);

  }

  std::pair<double, double> TrackQuality::GetResiduals() 
  {

    if (!_track_is_set || !_hits_are_set){
      std::cerr << "[TrackQuality] The Track or the hits were not set. Exiting now." << std::endl;
      throw std::exception();
    }

    std::vector<double> residuals;
    std::vector<double> angles;

    std::vector<bool> used_track_pts;
    used_track_pts.resize(_track_v.size());
    for (size_t i = 0; i < used_track_pts.size(); i++) used_track_pts.at(i) = false;

    //
    // Loop over hit points
    //
      for (size_t i = 0; i < _hit_v.size(); i++) {
        TVector3 hit = _hit_v.at(i);

        if (_debug) std::cout << "[TrackQuality]This is hit " << hit.X() << ", " << hit.Y() << std::endl;

      //
      // Loop over track points and find the closest one
      //
        double min_distance = 1e9;
        int min_sign = 1;
        int min_track_index = -1;

        for (size_t j = 0; j < _track_v.size(); j++) {

          TVector3 track_pt = _track_v.at(j);

          double distance = (hit - track_pt).Mag();
          double sign = (hit.Y() < track_pt.Y() ? 1 : -1);

          if (distance < min_distance) {
            min_distance = distance;
            min_sign = sign;
            min_track_index = j;
          }

      } // track loop

      if (_debug) std::cout << "[TrackQuality]\t matched to track " << _track_v.at(min_track_index).X() << ", " << _track_v.at(min_track_index).Y() << std::endl;

      double angle = -9999;
      if (min_track_index > 1 && min_track_index < (int)_track_v.size()-2) {
        TVector3 one = _track_v.at(min_track_index+2) - _track_v.at(min_track_index-2);
        TVector3 two = hit - _track_v.at(min_track_index);
        angle = one.Angle(two);
      } 

      if (min_track_index != -1 && min_track_index != 0 && min_track_index != (int)_track_v.size()-1) {
        residuals.emplace_back(min_distance * min_sign);
        angles.emplace_back(angle);

        used_track_pts.at(min_track_index) = true;
      }
    }

    for (size_t i = 0; i < residuals.size(); i++) {
      auto r = residuals.at(i);
      auto a = angles.at(i);
      if (_debug) std::cout << "[TrackQuality] Residual at " << i <<": " << r << ", angle " << a << std::endl;
    } 

    //for (size_t i = 0; i < used_track_pts.size(); i++) {
      //std::cout << "Track point " << i << " is " << (used_track_pts.at(i) ? "shadowed" : "not shadowed") << std::endl; 
    //}



    return std::make_pair(this->mean(residuals), this->stdev(residuals));

  }


  double TrackQuality::mean(const std::vector<double>& data)
  {
    if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to mean" << std::endl;

    double result = 0.0;

    for(const auto& d : data) 
      result += d;

    return (result / ((double)data.size()));
  }



  double TrackQuality::stdev(const std::vector<double>& data)
  {
    if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to stdev" << std::endl;

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return std::sqrt(result/((double)data.size()));
  }

  double TrackQuality::cov (const std::vector<double>& data1,
    const std::vector<double>& data2)
  {
    if(data1.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;
    if(data2.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(size_t i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());

  }

}


#endif
