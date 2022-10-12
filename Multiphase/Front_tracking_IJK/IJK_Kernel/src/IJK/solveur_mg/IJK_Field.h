//TRUST_NO_INDENT
/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Field.h
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file IJK_Field.h.P
//
#ifndef IJK_Field_H
#define IJK_Field_H

#include <Vect.h>
#include <Static_Int_Lists.h>
#include <TRUSTLists.h>
#include <TRUSTVect.h>
#include <TRUSTArray.h>
#include <IJK_Splitting.h>  
// #include <ArrOfFloat.h>

#define DIRECTION_I 0
#define DIRECTION_J 1
#define DIRECTION_K 2

// .Description
//  Class to store lists of boundary condition cells
//  In an IJK field, the list of all cells on the k layer is in K_Layer_Cells_Lists(k, ..)
//  It contains "linear_index" of cells in the IJ_Layout
//  See prepare_cond_lim
class K_Layer_Cells_Lists : public Objet_U
{
  Declare_instanciable(K_Layer_Cells_Lists);
public:
  // list is a k index, with -ghost_k_ < k < nk + ghost_k_
  int operator()(int list, int index) const { return data_(list + ghost_k_, index); }
  int get_list_size(int list) const { return data_.get_list_size(list + ghost_k_); }
  
  int ghost_k_;
  Static_Int_Lists data_;
};

Declare_vect(K_Layer_Cells_Lists);

// .Description
//  This is an array with [] operator.
//  Allocate array with resize(n, ghost).
//  tab[i] is valid for "-ghost <= i < n + ghost"
class ArrOfFloat_with_ghost : public Objet_U
{
  Declare_instanciable_sans_constructeur(ArrOfFloat_with_ghost);
public:
  ArrOfFloat_with_ghost() : ghost_(0) {};
  void resize(int n, int new_ghost) {
    tab_.resize_array(n + 2 * new_ghost);
    ghost_ = new_ghost;
  }
  float operator[](int i) const { return tab_[i + ghost_]; }
  float & operator[](int i) { return tab_[i + ghost_]; }
  const ArrOfFloat & data() const { return tab_; }
  int size() const { return tab_.size_array() - 2 * ghost_; }
  int ghost() const { return ghost_; }
  void echange_espace_virtuel(int pe_min, int pe_max);
protected:
  ArrOfFloat tab_;
  int ghost_;
};
Declare_vect(ArrOfFloat_with_ghost);
// .Description
//  This is an array with [] operator.
//  Allocate array with resize(n, ghost).
//  tab[i] is valid for "-ghost <= i < n + ghost"
class ArrOfDouble_with_ghost : public Objet_U
{
  Declare_instanciable_sans_constructeur(ArrOfDouble_with_ghost);
public:
  ArrOfDouble_with_ghost() : ghost_(0) {};
  void resize(int n, int new_ghost) {
    tab_.resize_array(n + 2 * new_ghost);
    ghost_ = new_ghost;
  }
  double operator[](int i) const { return tab_[i + ghost_]; }
  double & operator[](int i) { return tab_[i + ghost_]; }
  const ArrOfDouble & data() const { return tab_; }
  int size() const { return tab_.size_array() - 2 * ghost_; }
  int ghost() const { return ghost_; }
  void echange_espace_virtuel(int pe_min, int pe_max);
protected:
  ArrOfDouble tab_;
  int ghost_;
};
Declare_vect(ArrOfDouble_with_ghost);



// .Description
// This class describes a scalar field in an ijk box without any parallel information.
// The scalar field can be accessed by
//  - field(i,j,k) with "-ghost() <= i < ni() + ghost()", same for j and k
//  - field.data()[linear_index(i,j,k)]
class IJK_Field_local_float : public Objet_U
{
  Declare_instanciable(IJK_Field_local_float);
public:
  void allocate(int ni, int nj, int nk, int ghosts, int additional_k_layers = 0, int nb_compo = 1);
  void shift_k_origin(int n);
  void ref_ij(IJK_Field_local_float & src, int k_layer);

  int linear_index(int i, int j, int k) const {
    assert(nb_compo_ == 1); // otherwise, must specify component
    assert(i >= -ghost_size_ && i < ni_ + ghost_size_
	   && j >= -ghost_size_ && j < nj_ + ghost_size_
	   && k >= -ghost_size_ && k < nk_ + ghost_size_);
    return offset_ + k * compo_stride_ + j * j_stride_ + i;
  }
  int linear_index(int i, int j, int k, int compo) const {
    assert(i >= -ghost_size_ && i < ni_ + ghost_size_
	   && j >= -ghost_size_ && j < nj_ + ghost_size_
	   && k >= -ghost_size_ && k < nk_ + ghost_size_
	   && compo >= 0 && compo < nb_compo_);
    return offset_ + (k * nb_compo_ + compo) * compo_stride_ + j * j_stride_ + i;
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  float & operator()(int i, int j, int k) {
    int idx = linear_index(i, j, k);
    return data_[idx];
  }
  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  const float & operator()(int i, int j, int k) const {
    int idx = linear_index(i, j, k);
    return data_[idx];
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  float & operator()(int i, int j, int k, int compo) {
    int idx = linear_index(i, j, k, compo);
    return data_[idx];
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  const float & operator()(int i, int j, int k, int compo) const {
    int idx = linear_index(i, j, k, compo);
    return data_[idx];
  }

  int linear_index_relaxed_test(int i, int j, int k) const {
    assert(nb_compo_ == 1); // otherwise, must specify component
    assert(k >= -ghost_size_-k_layer_shift_ && k < nk_ + ghost_size_ + additional_k_layers_ - k_layer_shift_);
    int x = offset_ + k * compo_stride_ + j * j_stride_ + i;
    assert(x >= 0 && x < data_.size_array());
    return x;
  }
  int linear_index_relaxed_test(int i, int j, int k, int compo) const {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_-k_layer_shift_ && k < nk_ + ghost_size_ + additional_k_layers_ - k_layer_shift_);
    int x = offset_ + (k * nb_compo_ + compo) * compo_stride_ + j * j_stride_ + i;
    assert(x >= 0 && x < data_.size_array());
    return x;
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  float & get_in_allocated_area(int i, int j, int k) {
    int idx = linear_index_relaxed_test(i, j, k);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  const float & get_in_allocated_area(int i, int j, int k) const {
    int idx = linear_index_relaxed_test(i, j, k);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  float & get_in_allocated_area(int i, int j, int k, int compo) {
    int idx = linear_index_relaxed_test(i, j, k, compo);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  const float & get_in_allocated_area(int i, int j, int k, int compo) const {
    int idx = linear_index_relaxed_test(i, j, k, compo);
    return data_[idx];
  }
  
  float *k_layer(int k) {
    assert(nb_compo_ == 1);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + k * compo_stride_;
  }
  const float *k_layer(int k) const {
    assert(nb_compo_ == 1);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + k * compo_stride_;
  }
  float *k_layer(int k, int compo) {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + (k * nb_compo_ + compo) * compo_stride_;
  }
  const float *k_layer(int k, int compo) const {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + (k * nb_compo_ + compo) * compo_stride_;
  }

  int ni() const { return ni_; }
  int nj() const { return nj_; }
  int nk() const { return nk_; }
  int nb_elem_local(int dir) const { return (dir==0)?ni_:((dir==1)?nj_:nk_); }
  int nb_compo() const { return nb_compo_; }
  int j_stride() const { return j_stride_; }
  int compo_stride() const { return compo_stride_; }
  int k_stride() const { return compo_stride_ * nb_compo_; }
  int ghost() const { return ghost_size_; }
  int k_shift() const { return k_layer_shift_; }
  int k_shift_max() const { return additional_k_layers_; }
  ArrOfFloat & data() { return data_; }
  const ArrOfFloat & data() const { return data_; }
protected:
  // local size on this proc: (real items)
  // ni_ nj_ nk_ do not include the ghost size
  int ni_, nj_, nk_, ghost_size_, nb_compo_;
  int j_stride_; // how to jump to next j
  int compo_stride_; // how to jump to next component, k_stride is compo_stride_ * nb_compo_
  int offset_; // offset to first non ghost cell
  int k_layer_shift_; // current shift value of the origin in the k direction
  int additional_k_layers_;
  ArrOfFloat data_;
};


// .Description
// This class is an IJK_Field_local with parallel informations.
// Each processor has a sub_box of the global box,
// and echange_espace_virtuel(n) exchanges n layers of ghost cells,
// echange_espace_virtuel handles periodicity by copying the first layer
//  into the ghost layer on the opposite side.
class IJK_Field_float : public IJK_Field_local_float
{
  Declare_instanciable(IJK_Field_float);
public:
  void allocate(const IJK_Splitting &, IJK_Splitting::Localisation, int ghost_size, int additional_k_layers = 0, int nb_compo = 1);

  const IJK_Splitting & get_splitting() const { return splitting_ref_.valeur(); }
  IJK_Splitting::Localisation get_localisation() const { return localisation_; }
  void echange_espace_virtuel(int ghost);
protected:
  REF(IJK_Splitting) splitting_ref_;
  IJK_Splitting::Localisation localisation_;

  void exchange_data(int pe_imin_, /* processor to send to */
		     int is, int js, int ks, /* ijk coordinates of first data to send */
		     int pe_imax_, /* processor to recv from */
		     int ir, int jr, int kr, /* ijk coordinates of first data to recv */
		     int isz, int jsz, int ksz); /* size of block data to send/recv */
};

Declare_vect(IJK_Field_float);

double norme_ijk(const IJK_Field_float & x);
float max_ijk(const IJK_Field_float & x);
float prod_scal_ijk(const IJK_Field_float & x, const IJK_Field_float & y);
double somme_ijk(const IJK_Field_float & residu);

// .Description
// This class describes a scalar field in an ijk box without any parallel information.
// The scalar field can be accessed by
//  - field(i,j,k) with "-ghost() <= i < ni() + ghost()", same for j and k
//  - field.data()[linear_index(i,j,k)]
class IJK_Field_local_double : public Objet_U
{
  Declare_instanciable(IJK_Field_local_double);
public:
  void allocate(int ni, int nj, int nk, int ghosts, int additional_k_layers = 0, int nb_compo = 1);
  void shift_k_origin(int n);
  void ref_ij(IJK_Field_local_double & src, int k_layer);

  int linear_index(int i, int j, int k) const {
    assert(nb_compo_ == 1); // otherwise, must specify component
    assert(i >= -ghost_size_ && i < ni_ + ghost_size_
	   && j >= -ghost_size_ && j < nj_ + ghost_size_
	   && k >= -ghost_size_ && k < nk_ + ghost_size_);
    return offset_ + k * compo_stride_ + j * j_stride_ + i;
  }
  int linear_index(int i, int j, int k, int compo) const {
    assert(i >= -ghost_size_ && i < ni_ + ghost_size_
	   && j >= -ghost_size_ && j < nj_ + ghost_size_
	   && k >= -ghost_size_ && k < nk_ + ghost_size_
	   && compo >= 0 && compo < nb_compo_);
    return offset_ + (k * nb_compo_ + compo) * compo_stride_ + j * j_stride_ + i;
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  double & operator()(int i, int j, int k) {
    int idx = linear_index(i, j, k);
    return data_[idx];
  }
  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  const double & operator()(int i, int j, int k) const {
    int idx = linear_index(i, j, k);
    return data_[idx];
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  double & operator()(int i, int j, int k, int compo) {
    int idx = linear_index(i, j, k, compo);
    return data_[idx];
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  const double & operator()(int i, int j, int k, int compo) const {
    int idx = linear_index(i, j, k, compo);
    return data_[idx];
  }

  int linear_index_relaxed_test(int i, int j, int k) const {
    assert(nb_compo_ == 1); // otherwise, must specify component
    assert(k >= -ghost_size_-k_layer_shift_ && k < nk_ + ghost_size_ + additional_k_layers_ - k_layer_shift_);
    int x = offset_ + k * compo_stride_ + j * j_stride_ + i;
    assert(x >= 0 && x < data_.size_array());
    return x;
  }
  int linear_index_relaxed_test(int i, int j, int k, int compo) const {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_-k_layer_shift_ && k < nk_ + ghost_size_ + additional_k_layers_ - k_layer_shift_);
    int x = offset_ + (k * nb_compo_ + compo) * compo_stride_ + j * j_stride_ + i;
    assert(x >= 0 && x < data_.size_array());
    return x;
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  double & get_in_allocated_area(int i, int j, int k) {
    int idx = linear_index_relaxed_test(i, j, k);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  const double & get_in_allocated_area(int i, int j, int k) const {
    int idx = linear_index_relaxed_test(i, j, k);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  double & get_in_allocated_area(int i, int j, int k, int compo) {
    int idx = linear_index_relaxed_test(i, j, k, compo);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  const double & get_in_allocated_area(int i, int j, int k, int compo) const {
    int idx = linear_index_relaxed_test(i, j, k, compo);
    return data_[idx];
  }
  
  double *k_layer(int k) {
    assert(nb_compo_ == 1);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + k * compo_stride_;
  }
  const double *k_layer(int k) const {
    assert(nb_compo_ == 1);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + k * compo_stride_;
  }
  double *k_layer(int k, int compo) {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + (k * nb_compo_ + compo) * compo_stride_;
  }
  const double *k_layer(int k, int compo) const {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + (k * nb_compo_ + compo) * compo_stride_;
  }

  int ni() const { return ni_; }
  int nj() const { return nj_; }
  int nk() const { return nk_; }
  int nb_elem_local(int dir) const { return (dir==0)?ni_:((dir==1)?nj_:nk_); }
  int nb_compo() const { return nb_compo_; }
  int j_stride() const { return j_stride_; }
  int compo_stride() const { return compo_stride_; }
  int k_stride() const { return compo_stride_ * nb_compo_; }
  int ghost() const { return ghost_size_; }
  int k_shift() const { return k_layer_shift_; }
  int k_shift_max() const { return additional_k_layers_; }
  ArrOfDouble & data() { return data_; }
  const ArrOfDouble & data() const { return data_; }
protected:
  // local size on this proc: (real items)
  // ni_ nj_ nk_ do not include the ghost size
  int ni_, nj_, nk_, ghost_size_, nb_compo_;
  int j_stride_; // how to jump to next j
  int compo_stride_; // how to jump to next component, k_stride is compo_stride_ * nb_compo_
  int offset_; // offset to first non ghost cell
  int k_layer_shift_; // current shift value of the origin in the k direction
  int additional_k_layers_;
  ArrOfDouble data_;
};


// .Description
// This class is an IJK_Field_local with parallel informations.
// Each processor has a sub_box of the global box,
// and echange_espace_virtuel(n) exchanges n layers of ghost cells,
// echange_espace_virtuel handles periodicity by copying the first layer
//  into the ghost layer on the opposite side.
class IJK_Field_double : public IJK_Field_local_double
{
  Declare_instanciable(IJK_Field_double);
public:
  void allocate(const IJK_Splitting &, IJK_Splitting::Localisation, int ghost_size, int additional_k_layers = 0, int nb_compo = 1);

  const IJK_Splitting & get_splitting() const { return splitting_ref_.valeur(); }
  IJK_Splitting::Localisation get_localisation() const { return localisation_; }
  void echange_espace_virtuel(int ghost);
protected:
  REF(IJK_Splitting) splitting_ref_;
  IJK_Splitting::Localisation localisation_;

  void exchange_data(int pe_imin_, /* processor to send to */
		     int is, int js, int ks, /* ijk coordinates of first data to send */
		     int pe_imax_, /* processor to recv from */
		     int ir, int jr, int kr, /* ijk coordinates of first data to recv */
		     int isz, int jsz, int ksz); /* size of block data to send/recv */
};

Declare_vect(IJK_Field_double);

double norme_ijk(const IJK_Field_double & x);
double max_ijk(const IJK_Field_double & x);
double prod_scal_ijk(const IJK_Field_double & x, const IJK_Field_double & y);
double somme_ijk(const IJK_Field_double & residu);

// .Description
// This class describes a scalar field in an ijk box without any parallel information.
// The scalar field can be accessed by
//  - field(i,j,k) with "-ghost() <= i < ni() + ghost()", same for j and k
//  - field.data()[linear_index(i,j,k)]
class IJK_Field_local_int : public Objet_U
{
  Declare_instanciable(IJK_Field_local_int);
public:
  void allocate(int ni, int nj, int nk, int ghosts, int additional_k_layers = 0, int nb_compo = 1);
  void shift_k_origin(int n);
  void ref_ij(IJK_Field_local_int & src, int k_layer);

  int linear_index(int i, int j, int k) const {
    assert(nb_compo_ == 1); // otherwise, must specify component
    assert(i >= -ghost_size_ && i < ni_ + ghost_size_
	   && j >= -ghost_size_ && j < nj_ + ghost_size_
	   && k >= -ghost_size_ && k < nk_ + ghost_size_);
    return offset_ + k * compo_stride_ + j * j_stride_ + i;
  }
  int linear_index(int i, int j, int k, int compo) const {
    assert(i >= -ghost_size_ && i < ni_ + ghost_size_
	   && j >= -ghost_size_ && j < nj_ + ghost_size_
	   && k >= -ghost_size_ && k < nk_ + ghost_size_
	   && compo >= 0 && compo < nb_compo_);
    return offset_ + (k * nb_compo_ + compo) * compo_stride_ + j * j_stride_ + i;
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  int & operator()(int i, int j, int k) {
    int idx = linear_index(i, j, k);
    return data_[idx];
  }
  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  const int & operator()(int i, int j, int k) const {
    int idx = linear_index(i, j, k);
    return data_[idx];
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  int & operator()(int i, int j, int k, int compo) {
    int idx = linear_index(i, j, k, compo);
    return data_[idx];
  }

  // Operator() checks if the requested i,j,k index lies within the valid
  // range [-ghost,n+ghost-1]
  const int & operator()(int i, int j, int k, int compo) const {
    int idx = linear_index(i, j, k, compo);
    return data_[idx];
  }

  int linear_index_relaxed_test(int i, int j, int k) const {
    assert(nb_compo_ == 1); // otherwise, must specify component
    assert(k >= -ghost_size_-k_layer_shift_ && k < nk_ + ghost_size_ + additional_k_layers_ - k_layer_shift_);
    int x = offset_ + k * compo_stride_ + j * j_stride_ + i;
    assert(x >= 0 && x < data_.size_array());
    return x;
  }
  int linear_index_relaxed_test(int i, int j, int k, int compo) const {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_-k_layer_shift_ && k < nk_ + ghost_size_ + additional_k_layers_ - k_layer_shift_);
    int x = offset_ + (k * nb_compo_ + compo) * compo_stride_ + j * j_stride_ + i;
    assert(x >= 0 && x < data_.size_array());
    return x;
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  int & get_in_allocated_area(int i, int j, int k) {
    int idx = linear_index_relaxed_test(i, j, k);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  const int & get_in_allocated_area(int i, int j, int k) const {
    int idx = linear_index_relaxed_test(i, j, k);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  int & get_in_allocated_area(int i, int j, int k, int compo) {
    int idx = linear_index_relaxed_test(i, j, k, compo);
    return data_[idx];
  }
  // This method allows to access padding cells outside of the valid data range
  // but inside the allocated data block.
  const int & get_in_allocated_area(int i, int j, int k, int compo) const {
    int idx = linear_index_relaxed_test(i, j, k, compo);
    return data_[idx];
  }
  
  int *k_layer(int k) {
    assert(nb_compo_ == 1);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + k * compo_stride_;
  }
  const int *k_layer(int k) const {
    assert(nb_compo_ == 1);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + k * compo_stride_;
  }
  int *k_layer(int k, int compo) {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + (k * nb_compo_ + compo) * compo_stride_;
  }
  const int *k_layer(int k, int compo) const {
    assert(compo >= 0 && compo < nb_compo_);
    assert(k >= -ghost_size_ && k < nk_ + ghost_size_);
    return data_.addr() + offset_ + (k * nb_compo_ + compo) * compo_stride_;
  }

  int ni() const { return ni_; }
  int nj() const { return nj_; }
  int nk() const { return nk_; }
  int nb_elem_local(int dir) const { return (dir==0)?ni_:((dir==1)?nj_:nk_); }
  int nb_compo() const { return nb_compo_; }
  int j_stride() const { return j_stride_; }
  int compo_stride() const { return compo_stride_; }
  int k_stride() const { return compo_stride_ * nb_compo_; }
  int ghost() const { return ghost_size_; }
  int k_shift() const { return k_layer_shift_; }
  int k_shift_max() const { return additional_k_layers_; }
  ArrOfInt & data() { return data_; }
  const ArrOfInt & data() const { return data_; }
protected:
  // local size on this proc: (real items)
  // ni_ nj_ nk_ do not include the ghost size
  int ni_, nj_, nk_, ghost_size_, nb_compo_;
  int j_stride_; // how to jump to next j
  int compo_stride_; // how to jump to next component, k_stride is compo_stride_ * nb_compo_
  int offset_; // offset to first non ghost cell
  int k_layer_shift_; // current shift value of the origin in the k direction
  int additional_k_layers_;
  ArrOfInt data_;
};


// .Description
// This class is an IJK_Field_local with parallel informations.
// Each processor has a sub_box of the global box,
// and echange_espace_virtuel(n) exchanges n layers of ghost cells,
// echange_espace_virtuel handles periodicity by copying the first layer
//  into the ghost layer on the opposite side.
class IJK_Field_int : public IJK_Field_local_int
{
  Declare_instanciable(IJK_Field_int);
public:
  void allocate(const IJK_Splitting &, IJK_Splitting::Localisation, int ghost_size, int additional_k_layers = 0, int nb_compo = 1);

  const IJK_Splitting & get_splitting() const { return splitting_ref_.valeur(); }
  IJK_Splitting::Localisation get_localisation() const { return localisation_; }
  void echange_espace_virtuel(int ghost);
protected:
  REF(IJK_Splitting) splitting_ref_;
  IJK_Splitting::Localisation localisation_;

  void exchange_data(int pe_imin_, /* processor to send to */
		     int is, int js, int ks, /* ijk coordinates of first data to send */
		     int pe_imax_, /* processor to recv from */
		     int ir, int jr, int kr, /* ijk coordinates of first data to recv */
		     int isz, int jsz, int ksz); /* size of block data to send/recv */
};

Declare_vect(IJK_Field_int);

double norme_ijk(const IJK_Field_int & x);
int max_ijk(const IJK_Field_int & x);
int prod_scal_ijk(const IJK_Field_int & x, const IJK_Field_int & y);
double somme_ijk(const IJK_Field_int & residu);


typedef IJK_Field_local_double IJK_Field_local;
typedef IJK_Field_double IJK_Field;
typedef VECT(IJK_Field_double) VECT(IJK_Field);
#ifdef INT_is_64_
using IJK_Field_int = IJK_Field_long;
#endif

// .Description
//  This class describes an IJ plane of an IJK_Field_local.
//  It just encapsulates the offset computation:
//    j * j_stride_ + i
//  with bound checking in debug...
//  This class is usefull to safely optimize the code !
class IJ_layout
{
public:
  int ni() const { return ni_; }
  int nj() const { return nj_; }
  int ghost() const { return ghost_; }
  int j_stride() const { return j_stride_; }

  IJ_layout(const IJK_Field_local_double & f) :
    ni_(f.ni()), 
    nj_(f.nj()), 
    ghost_(f.ghost()), 
    j_stride_(f.j_stride()) {}

  // origin points to x(0,0,k) for an IJK_Field
  double & operator()(double * origin, int i, int j) const {
    assert(i >= -ghost_ && i < ni_ + ghost_);
    assert(j >= -ghost_ && j < nj_ + ghost_);
    return origin[j * j_stride_ + i];
  }

  // origin points to x(0,0,k) for an IJK_Field
  const double & operator()(const double * origin, int i, int j) const {
    assert(i >= -ghost_ && i < ni_ + ghost_);
    assert(j >= -ghost_ && j < nj_ + ghost_);
    return origin[j * j_stride_ + i];
  }
  // linear index of cell (i=0,j=0) is zero
  // origin points to x(0,0,k) for an IJK_Field
  double & linear(double * origin, int linear_index) const {
    assert(linear_index >= -ghost_*(j_stride_+1));
    assert(linear_index < (ni_ + ghost_) + j_stride_ * (nj_ + ghost_ - 1));
    return origin[linear_index];
  }

  const double & linear(const double * origin, int linear_index) const {
    assert(linear_index >= -ghost_*(j_stride_+1));
    assert(linear_index < (ni_ + ghost_) + j_stride_ * (nj_ + ghost_ - 1));
    return origin[linear_index];
  }
  IJ_layout(const IJK_Field_local_float & f) :
    ni_(f.ni()), 
    nj_(f.nj()), 
    ghost_(f.ghost()), 
    j_stride_(f.j_stride()) {}

  // origin points to x(0,0,k) for an IJK_Field
  float & operator()(float * origin, int i, int j) const {
    assert(i >= -ghost_ && i < ni_ + ghost_);
    assert(j >= -ghost_ && j < nj_ + ghost_);
    return origin[j * j_stride_ + i];
  }

  // origin points to x(0,0,k) for an IJK_Field
  const float & operator()(const float * origin, int i, int j) const {
    assert(i >= -ghost_ && i < ni_ + ghost_);
    assert(j >= -ghost_ && j < nj_ + ghost_);
    return origin[j * j_stride_ + i];
  }
  // linear index of cell (i=0,j=0) is zero
  // origin points to x(0,0,k) for an IJK_Field
  float & linear(float * origin, int linear_index) const {
    assert(linear_index >= -ghost_*(j_stride_+1));
    assert(linear_index < (ni_ + ghost_) + j_stride_ * (nj_ + ghost_ - 1));
    return origin[linear_index];
  }

  const float & linear(const float * origin, int linear_index) const {
    assert(linear_index >= -ghost_*(j_stride_+1));
    assert(linear_index < (ni_ + ghost_) + j_stride_ * (nj_ + ghost_ - 1));
    return origin[linear_index];
  }
  void linear_to_ij(int linear_index, int &i, int &j) const {
    assert(linear_index >= -ghost_*(j_stride_+1));
    assert(linear_index < (ni_ + ghost_) + j_stride_ * (nj_ + ghost_ - 1));
    int l = linear_index + ghost_*(j_stride_+1);
    j = (l / j_stride_) - ghost_;
    i = (l % j_stride_) - ghost_;
  }

  // ghost_ is the number of ghost cells allocated in memory in I and J directions
  int ni_, nj_, ghost_, j_stride_;
};


#endif
