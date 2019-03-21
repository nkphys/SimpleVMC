/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-20 13:07:50
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-20 13:08:22
*----------------------------------------------------------------------------*/
// File: vmc.h
#ifndef VMC_H
#define VMC_H

#include <iostream>
#include "sysconfig.h"
#include "mcdata/mc_observable.h"

class VMC
{
public:
	VMC() {}
	~VMC() {}
	int init(void);
	int run_simulation(void);
private:
	SysConfig config;
	RealVector vparams;
	int num_vparams;
	int num_samples;
	int warmup_steps;
	int interval;

	// observables
	mcdata::MC_Observable energy;
};



#endif