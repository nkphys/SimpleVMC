/*---------------------------------------------------------------------------
* @Author: amedhi
* @Date:   2019-03-19 13:12:20
* @Last Modified by:   amedhi
* @Last Modified time: 2019-03-19 13:42:09
*----------------------------------------------------------------------------*/
// File: lattcie.h
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include "matrix.h"

enum class lattice_id {
  CHAIN, SQUARE, HONEYCOMB, SIMPLECUBIC
};

class lattice_size
{
public:
	lattice_size() : L1_{1}, L2_{1}, L3_{1}
	{
	}
	lattice_size(const int& L1, const int& L2=1, const int& L3=1) 
		: L1_{L1}, L2_{L2}, L3_{L3}
	{
		if (L1_<1 || L2_<1 || L3_<1) {
			throw std::invalid_argument("lattice_size:: Invalid lattice size\n");
		}
	}
	~lattice_size() {}
	const int& L1(void) const { return L1_; }
	const int& L2(void) const { return L2_; }
	const int& L3(void) const { return L3_; }
private:
	int L1_;
	int L2_;
	int L3_;
};

class Bond
{
public:
	Bond() {}
	Bond(const int& src, const int& tgt, const Vector3d& vec)
		: src_{src}, tgt_{tgt}, vector_{vec} {}
	~Bond() {}
	const int& src(void) const { return src_; }
	const int& tgt(void) const { return tgt_; }
	const Vector3d& vector(void) const { return vector_; }
private:
	int src_{0};
	int tgt_{0};
	Vector3d vector_{0,0,0};
};

class Lattice
{
public:
	Lattice() : id_{lattice_id::CHAIN} { num_sites_=1; }
	Lattice(const lattice_id& id, const lattice_size& size) { construct(id, size); }
	void construct(const lattice_id& id, const lattice_size& size);
	~Lattice() {}
	const int& size_L1(void) const { return size_.L1(); }
	const int& size_L2(void) const { return size_.L2(); }
	const int& size_L3(void) const { return size_.L3(); }
	const int& num_sites(void) const { return num_sites_; }
	const int& num_basis_sites(void) const { return num_basis_sites_; }
	const int& num_kpoints(void) const { return num_kpoints_; }
	const int& num_neighbs(void) const { return num_neighbs_; }
	const Bond& bond(const int& i) const { return bonds_[i]; }
	const std::vector<int>& site_nn(const int& site) const { return nn_table_[site]; }
	const Vector3d& kpoint(const int& i) { return kpoints_[i]; }
	const std::vector<Vector3d>& kpoints(void) { return kpoints_; }
	const Vector3d& site_coord(const int& i) const { return rpoints_[i]; }
private:
	lattice_id id_;
	lattice_size size_;
	int lattice_dim_;
	int num_basis_sites_; // number of sites per unit cell
	int num_sites_; // total number of sites
	int num_bonds_; // total number of bonds
	int num_kpoints_; 
	int num_neighbs_;
	Vector3d a1_;
	Vector3d a2_;
	Vector3d a3_;
	Vector3d b1_;
	Vector3d b2_;
	Vector3d b3_;
	std::vector<Vector3d> kpoints_;
	std::vector<Vector3d> rpoints_; // position coordinates
	std::vector<std::vector<int> > nn_table_;
	std::vector<Bond> bonds_;
	void construct_square(const lattice_size& size);
	void construct_kpoints(void);
	Vector3i get_next_bravindex(const Vector3i& current_index) const;
};


#endif