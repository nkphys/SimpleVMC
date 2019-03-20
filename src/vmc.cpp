/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-20 13:07:50
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-20 16:44:06
*----------------------------------------------------------------------------*/
// File: vmc.cpp

#include "vmc.h"

int VMC::init(void) 
{
  lattice.construct(lattice_id::SQUARE, lattice_size(4,4));
  config.init(lattice);
  num_vparams = config.num_vparams();
  vparams.resize(num_vparams);

  //config.build(lattice);
  // run parameters
  num_samples = 100;
  warmup_steps = 100;
  interval = 3;

  return 0;
}

int VMC::run_simulation(void) 
{
  // set variational parameters
  vparams.setOnes();
  config.build(lattice, vparams);

  // warmup run
  config.reset();
  for (int n=0; n<warmup_steps; ++n) {
    config.update_state();
  } 
  std::cout << " warmup done\n";
  // measuring run
  int sample = 0;
  int skip_count = interval;
  // Initialize observables
  while (sample < num_samples) {
    if (skip_count == interval) {
      skip_count = 0;
      ++sample;
      int progress = int((100.0*sample)/num_samples);
      if (progress%10==0) std::cout<<" done = "<<progress<<"% \n";
      // Make measurements
    }
    config.update_state();
    skip_count++;
  }
  // Finalize observables
  std::cout << " simulation done\n";

  return 0;
}