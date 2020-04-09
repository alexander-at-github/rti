#pragma once

// a logger for the rays we are using
//
// NOT THREAD SAVE: use with a single thread only!
// CAUSES SEGMENTATION FAULT WHEN USED IN MULTIPLE THREADS.
//

#include "rti/util/utils.hpp"

namespace rti { namespace util {
//
//#define RAYLOG_ON
//#define RAYSRCLOG_ON

  // vector of 3D line segments
  static auto sRayLogVec = std::vector<rti::util::pair<rti::util::triple<float> > > {};
  // vector of 3D points
  static auto sPointLogVec = std::vector<rti::util::triple<float> > {};
  static constexpr auto sLineMaxLength = 64.f;


#ifdef RAYLOG_ON
    //auto tfar = std::min((rh).ray.tfar, rti::util::sLineMaxLength);
    //
    // This macro expects a semicolon at the end
#define RAYLOG(rh, tfar)                             \
  { \
    auto p1 = rti::util::triple<float> {(rh).ray.org_x, (rh).ray.org_y, (rh).ray.org_z}; \
    tfar = ((tfar) > 10 ? 10 : (tfar));                                 \
    auto p2 = rti::util::triple<float> {p1[0] + tfar * (rh).ray.dir_x,  \
                                        p1[1] + tfar * (rh).ray.dir_y,  \
                                        p1[2] + tfar * (rh).ray.dir_z}; \
    rti::util::sRayLogVec.push_back({p1, p2}); \
  } \
  do {} while(false) // allows us to put a semicolon after the macro (without warnings)
#else
#define RAYLOG(rh, tfar)                                                     \
  do {} while(false) // allows us to put a semicolon after the macro (without warnings)
#endif


#ifdef RAYLOG_ON
    // This macro expects a semicolon at the end
#define RAYLOG_GET_PTR() \
    [](){return &rti::util::sRayLogVec;}()
#else
#define RAYLOG_GET_PTR() \
    [](){return (std::vector<rti::util::pair<rti::util::triple<float> > >*) nullptr;}()
#endif


#ifdef RAYSRCLOG_ON
#define RAYSRCLOG(rh) \
  { \
  rti::util::sPointLogVec.push_back({rh.ray.org_x, rh.ray.org_y, rh.ray.org_z}); \
  } \
  do {} while(false) // allows to put a semicolon after macro (without warnings)
#else
#define RAYSRCLOG(rh) \
  do {} while(false) // allows to put a semicolon after macro (without warnings)
#endif


#ifdef RAYSRCLOG_ON
#define RAYSRCLOG_GET_PTR(rh) \
  [](){return &rti::util::sPointLogVec;}()
#else
#define RAYSRCLOG_GET_PTR(rh) \
  [](){return (std::vector<rti::util::triple<float> >*)  nullptr;}()
#endif


}}
