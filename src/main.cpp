/*---------------------------------------------------------------------------
* @Author: amedhi
* @Date:   2019-03-19 13:12:20
* @Last Modified by:   amedhi
* @Last Modified time: 2019-03-19 14:06:02
*----------------------------------------------------------------------------*/

#include <iostream>
//#include "vmc.h"
#include "lattice.h"
#include "basis.h"

int main(int argc, const char *argv[])
{
  Lattice lattice(lattice_id::SQUARE, lattice_size(4,4));
  FockBasis basis(lattice.num_sites());
  basis.init_spins(4,4);
  std::cout << basis << "\n";
}
