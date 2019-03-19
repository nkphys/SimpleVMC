/*---------------------------------------------------------------------------
* @Author: amedhi
* @Date:   2019-03-19 13:12:20
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-19 16:52:32
*----------------------------------------------------------------------------*/

#include <iostream>
//#include "vmc.h"
#include "lattice.h"
#include "basis.h"
#include "wavefunction.h"

int main(int argc, const char *argv[])
{
  Lattice lattice(lattice_id::SQUARE, lattice_size(4,4));
  /*for (int n=0; n<lattice.num_bonds(); ++n) {
  	int i = lattice.bond(n).src();
  	int j = lattice.bond(n).tgt();
  	// calculate <c^dag_i c_j + hc>
  	std::cout << "i,j=" << i << "  " << j << "\n";
  }*/

  FockBasis basis(lattice.num_sites());
  basis.init_spins(4,4);
  std::cout << basis.state().transpose() << "\n";
  for (int n=0; n<100; ++n) {
  	//basis.gen_upspin_hop();
  	// calculate the W
  	//basis.commit_last_move();
  	//basis.undo_last_move();
  	//std::cout << basis.state().transpose() << "\n";
  	//basis.gen_upspin_hop();
  }

  Wavefunction wf(wf_id::BCS, lattice, 0.0);

  //SysConfig config

  // VMC vmc

  //vmc.run_simulation()


}
