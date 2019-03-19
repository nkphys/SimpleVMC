/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-15 22:23:46
*----------------------------------------------------------------------------*/
#include "random.h"

RandomGenerator::RandomGenerator() : seed_type_(0), real_generator(0.0, 1.0) 
{
  //for (auto& g : state_generators) g = int_dist(-1,-1);
}
  
RandomGenerator::RandomGenerator(const unsigned& seed_type) 
  : seed_type_(seed_type), real_generator(0.0, 1.0)
{
  if (seed_type_==1) time_seed();
  //for (auto& g : state_generators) g = int_dist(-1,-1);
}

void RandomGenerator::seed(const int& seed_type) {
  seed_type_ = seed_type;
  if (seed_type_==1) time_seed();
} 

void RandomGenerator::time_seed(void) 
{
  myclock::time_point now = myclock::now();
  myclock::duration till_now = now.time_since_epoch();
  unsigned itc = till_now.count();
  this->std::mt19937_64::seed(itc);
}

void RandomGenerator::set_site_generator(const unsigned& min, const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomGenerator::set_site_generator: invalid input");
  site_generator = int_generator(min, max);
}

void RandomGenerator::set_upspin_generator(const unsigned& min, const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomGenerator::set_upspin_generator: invalid input");
  upspin_generator = int_generator(min, max);
}

void RandomGenerator::set_dnspin_generator(const unsigned& min, const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomGenerator::set_dnspin_generator: invalid input");
  dnspin_generator = int_generator(min, max);
}

void RandomGenerator::set_uphole_generator(const unsigned& min, const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomGenerator::set_uphole_generator: invalid input");
  uphole_generator = int_generator(min, max);
}

void RandomGenerator::set_dnhole_generator(const unsigned& min, const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomGenerator::set_dnhole_generator: invalid input");
  dnhole_generator = int_generator(min, max);
}



