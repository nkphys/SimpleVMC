/*---------------------------------------------------------------------------
* @Author: amedhi
* @Date:   2019-03-19 13:12:20
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-22 00:07:49
*----------------------------------------------------------------------------*/
// File: lattice.cpp

#include <Eigen/Dense>
#include "constants.h"
#include "lattice.h"

void Lattice::construct(const lattice_id& id, const lattice_size& size)
{
	id_ = id;
	size_ = size;
	switch (id_) {
		case lattice_id::SQUARE: 
			lattice_dim_ = 2;
      bc_ = {bc_t::PERIODIC, bc_t::ANTIPERIODIC, bc_t::OPEN};
			construct_square(size);
			break;
		default: 
			throw std::range_error("This lattice not implemented\n");
			break;
	}
}

void Lattice::construct_square(const lattice_size& size)
{
	num_basis_sites_ = 1;
	num_sites_ = num_basis_sites_ * size_.L1() * size_.L2();
	a1_ = Vector3d(1,0,0);
	a2_ = Vector3d(0,1,0);
	a3_ = Vector3d(0,0,0);
	construct_kpoints();

  //------- Nearest Neighbout Table
  // Numbering scheme for the SQUARE lattice:
  /*
  *   12   13   14   15      
  *    8    9   10   11      
  *    4    5    6    7     
  *    0    1    2    3    
  *-----------------------------------------*/
  num_neighbs_ = 4; 
  nn_table_.resize(num_sites_);
  for (int i=0; i<num_sites_; ++i) nn_table_[i].resize(num_neighbs_);

  int L1 = size_.L1();
  int L2 = size_.L2();
  enum nn_dir {right_nn, top_nn, left_nn, bottom_nn};

  // Right NNs
  for (int i=0; i<num_sites_; ++i) nn_table_[i][right_nn] = i+1;
  // Top NNs
  for (int i=0; i<num_sites_; ++i) nn_table_[i][top_nn] = i+L1;
  // Left NNs
  for (int i=0; i<num_sites_; ++i) nn_table_[i][left_nn] = i-1;
  // Bottom NNs
  for (int i=0; i<num_sites_; ++i) nn_table_[i][bottom_nn] = i-L1;

  // Take care of the boundary points
  int m1 = L1-1; 
  int m2 = L1 * (L2-1);

  // bottom line
  for (int i=0; i<L1; ++i) nn_table_[i][bottom_nn] = i + m2;
  // right line
  for (int i=m1; i<num_sites_; i += L1) nn_table_[i][right_nn] = i - m1;
  // top line
  for (int i=m2; i<num_sites_; ++i) nn_table_[i][top_nn] = i - m2;
  // left line
  for (int i=0; i<num_sites_; i += L1) nn_table_[i][left_nn] = i + m1;
  // print
  /*
  for (int i=0; i<num_sites_; ++i) {
    std::cout << "nn_table[i][nn] = " << i << " " << nn_table_[i][0] 
      << " " << nn_table_[i][1] << " "
      << nn_table_[i][2] << " " << nn_table_[i][3] << std::endl;
  }*/

  // Position Coordinates
  /*Vector3i n = {0,0,0};
  rpoints_.clear();
  for (int i=0; i<num_sites_; ++i) {
    Vector3d R = n(0) * a1_ + n(1) * a2_ + n(2) * a3_;
    rpoints_.push_back(R);
    n = get_next_bravindex(n);
  }*/
  // check
  /*for (int i=0; i<num_sites_; ++i) {
    std::cout << i << ": " << rpoints_[i].transpose() << "\n";
  }*/

  // Sites in the lattice
  sites_.clear();
  Vector3i n = {0,0,0};
  for (int i=0; i<num_sites_; ++i) {
    Vector3d R = n(0)*a1_ + n(1)*a2_ + n(2)*a3_;
    sites_.push_back(Site(i,0,R));
    n = get_next_bravindex(n);
  }

  // Bonds in the lattice
  bonds_.clear();
  int id = 0;
  int nn, phase;
  Vector3d R;
  for (int i=0; i<num_sites_; ++i) {
    nn = nn_table_[i][right_nn];
    phase = 1;
    if (bc_.L1_bc()==bc_t::ANTIPERIODIC && nn<i) phase = -1;
    R = sites_[nn].cell_coord()-sites_[i].cell_coord();
    bonds_.push_back(Bond(id,i,nn,phase,R));
    id++;

    nn = nn_table_[i][top_nn];
    phase = 1;
    if (bc_.L2_bc()==bc_t::ANTIPERIODIC && nn<i) phase = -1;
    R = sites_[nn].cell_coord()-sites_[i].cell_coord();
    bonds_.push_back(Bond(id,i,nn,phase,R));
    id++;
  }
  num_bonds_ = bonds_.size();
}

void Lattice::construct_kpoints(void)
{
	// reciprocal lattice vectors
  b1_ = Vector3d(0,0,0);
  b2_ = Vector3d(0,0,0);
  b3_ = Vector3d(0,0,0);
  double v;
  Vector3d n3;
  switch (lattice_dim_) {
  	case 1:
    	b1_ = TWO_PI * a1_/a1_.dot(a1_); 
    	break;
  	case 2:
      n3 = a1_.cross(a2_);
      v = a1_.dot(a2_.cross(n3));
      b1_ = TWO_PI * a2_.cross(n3)/v;
      b2_ = TWO_PI * n3.cross(a1_)/v;
    	break;
  	case 3:
      v = a1_.dot(a2_.cross(a3_));
      b1_ = TWO_PI * a2_.cross(a3_) / v;
      b2_ = TWO_PI * a3_.cross(a1_) / v;
      b3_ = TWO_PI * a1_.cross(a2_) / v;
      break;
  	default: break;
  }

  // set 'antiperiodic boundary' along 'y'
  Vector3d antipb_shift(0.0,0.0,0.0);
  if (bc_.L1_bc()==bc_t::ANTIPERIODIC) {
    antipb_shift(0) = 0.5/size_.L1();
  }
  if (bc_.L2_bc()==bc_t::ANTIPERIODIC) {
    antipb_shift(1) = 0.5/size_.L2();
  }
  if (bc_.L3_bc()==bc_t::ANTIPERIODIC) {
    antipb_shift(2) = 0.5/size_.L3();
  }

  // kpoints
  num_kpoints_ = num_sites_/num_basis_sites_;
  Vector3i n = {0,0,0};
  Vector3i m = {-size_.L1()/2, -size_.L2()/2, -size_.L3()/2};
  kpoints_.clear();
  double x1, x2, x3;
  for (int i=0; i<num_kpoints_; ++i) {
    x1 = static_cast<double>(m(0)+n(0))/size_.L1() + antipb_shift(0);
    x2 = static_cast<double>(m(1)+n(1))/size_.L2() + antipb_shift(1);
    x3 = static_cast<double>(m(2)+n(2))/size_.L3() + antipb_shift(2);
    kpoints_.push_back(x1*b1_ + x2*b2_ + x3*b3_);
    n = get_next_bravindex(n);
  }
  // check
  /*
  for (int i=0; i<num_kpoints_; ++i) {
    std::cout << i << ": " << kpoints_[i].transpose() << "\n";
  }*/

}

Vector3i Lattice::get_next_bravindex(const Vector3i& current_index) const
{
  /* Returns the next Bravais lattice index. 
  ! Index for first unit cell = (0,0,0)
  ! Index for last unit cell = (N1-1, N2-1, N3-1)
   */
  Vector3i next_index = current_index;
	next_index[0] += 1;
  if (next_index[0] >= size_.L1()) {
    next_index[0] = 0;
		next_index[1] += 1;
    if (next_index[1] >= size_.L2()) {
      next_index[1] = 0;
			next_index[2] += 1;
      if (next_index[2] >= size_.L3()) {
        next_index[2] = 0;
      }
    }
  }
  return next_index;
}



