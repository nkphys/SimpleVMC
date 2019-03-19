/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <stdexcept>
#include <chrono>
#include <array>
#include <map>
#include <random>
#include <time.h> 
//#include "../lattice/lattice.h" 

class RandomGenerator : public std::mt19937_64
{
public:
  RandomGenerator();
  RandomGenerator(const unsigned& seed_type);
  ~RandomGenerator() {};
  void set_site_generator(const unsigned& min, const unsigned& max);
  void set_upspin_generator(const unsigned& min, const unsigned& max);
  void set_dnspin_generator(const unsigned& min, const unsigned& max);
  void set_uphole_generator(const unsigned& min, const unsigned& max);
  void set_dnhole_generator(const unsigned& min, const unsigned& max);
  void seed(const int& seed_type);
  void time_seed(void);
  //unsigned random_idx(const unsigned& site_type) {return state_dist_map[site_type](*this); }
  //unsigned random_idx(const unsigned& site_type) { return state_generators[site_type](*this); }
  //unsigned random_site(void) { return site_dist[0](*this); }
  unsigned random_site(void) { return site_generator(*this); }
  unsigned random_upspin(void) { return upspin_generator(*this); }
  unsigned random_dnspin(void) { return dnspin_generator(*this); }
  unsigned random_uphole(void) { return uphole_generator(*this); }
  unsigned random_dnhole(void) { return dnhole_generator(*this); }
  double random_real(void) { return real_generator(*this); }
private:
  using int_generator = std::uniform_int_distribution<unsigned>;
  using myclock = std::chrono::high_resolution_clock;
  int seed_type_;
  int_generator site_generator; 
  int_generator upspin_generator; 
  int_generator dnspin_generator; 
  int_generator uphole_generator; 
  int_generator dnhole_generator; 
  std::uniform_real_distribution<double> real_generator;
};


#endif
