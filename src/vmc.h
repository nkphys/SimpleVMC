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
#include "lattice.h"
#include "sysconfig.h"

class VMC
{
public:
	VMC() {}
	~VMC() {}
	int init(void);
	int run_simulation(void);
private:
	Lattice lattice;
	SysConfig config;
	RealVector vparams;
	int num_vparams;
	int num_samples;
	int warmup_steps;
	int interval;
};



#endif