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
// File      : Multigrille_base.h
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file Multigrille_base.h.P
//
#ifndef Multigrille_base_H
#define Multigrille_base_H

#include <SolveurSys.h>
#include <TRUSTLists.h>
#include <Matrice_Grossiere.h>

class Param;

class Multigrille_base : public SolveurSys_base
{
  Declare_base_sans_constructeur(Multigrille_base);
 public:
  Multigrille_base();
  int get_flag_updated_input() const override;
  int resoudre_systeme(const Matrice_Base & a, const DoubleVect & b, DoubleVect & x, int ) override;
  int resoudre_systeme(const Matrice_Base & a, const DoubleVect & b, DoubleVect & x) override;
  int solveur_direct() const override { return 0 ; };
  void resoudre_systeme_IJK(const IJK_Field_double & rhs, IJK_Field_double & x);

  void convert_from_ijk(const IJK_Field_double & ijk_x, DoubleVect & x);
  void convert_to_ijk(const DoubleVect & x, IJK_Field_double & ijk_x);

  double multigrille(IJK_Field_double & x, const IJK_Field_double & b, IJK_Field_double & residu);

  virtual void prepare_secmem(IJK_Field_double & x) const = 0;
  virtual void alloc_field(IJK_Field_double & x, int level, bool with_additional_k_layers = false) const = 0;

  virtual void dump_lata(const Nom & field, const IJK_Field_double & data, int tstep) const = 0;
  void resoudre_systeme_IJK(const IJK_Field_float & rhs, IJK_Field_float & x);

  void convert_from_ijk(const IJK_Field_float & ijk_x, DoubleVect & x);
  void convert_to_ijk(const DoubleVect & x, IJK_Field_float & ijk_x);

  double multigrille(IJK_Field_float & x, const IJK_Field_float & b, IJK_Field_float & residu);

  virtual void prepare_secmem(IJK_Field_float & x) const = 0;
  virtual void alloc_field(IJK_Field_float & x, int level, bool with_additional_k_layers = false) const = 0;

  virtual void dump_lata(const Nom & field, const IJK_Field_float & data, int tstep) const = 0;
  virtual int nb_grid_levels() const = 0;

 protected:
  virtual void ajouter_param(Param &);

  Matrice_Grossiere & set_coarse_matrix() { 
    solveur_grossier_.valeur().reinit(); 
    return coarse_matrix_; 
  }
  enum StorageId { STORAGE_RHS, STORAGE_X, STORAGE_RESIDUE };
  virtual int needed_kshift_for_jacobi(int level) const = 0;
  int nsweeps_jacobi_residu(int level) const;
  void coarse_solver(IJK_Field_double & x, const IJK_Field_double & b);

  void resoudre_avec_gcp(IJK_Field_double & x, IJK_Field_double & b, IJK_Field_double & residu);
  void resoudre_avec_gmres(IJK_Field_double & x, IJK_Field_double & b, IJK_Field_double & residu);
  void resoudre_systeme_double(const Matrice_Base & a, const DoubleVect & b, DoubleVect & x); 
  void resoudre_systeme_IJK_double(const IJK_Field_double & rhs, IJK_Field_double & x);
  void resoudre_systeme_IJK_double(const IJK_Field_float & rhs, IJK_Field_float & x);
  void solve_ijk_in_storage_double();

  virtual void jacobi_residu(IJK_Field_double & x,
			     const IJK_Field_double *secmem, /* if null pointer, take secmem = 0 (to compute A*x) */
			     const int grid_level,
			     const int n_jacobi,
			     IJK_Field_double *residu) const = 0;

  virtual void coarsen(const IJK_Field_double & fine, IJK_Field_double & coarse, int fine_level) const = 0;
  virtual void interpolate_sub_shiftk(const IJK_Field_double & coarse, IJK_Field_double & fine, int fine_level) const = 0;
  
  double multigrille_(IJK_Field_double & x, const IJK_Field_double & b, IJK_Field_double & residu, 
		      int grid_level, int iter_number);

  virtual IJK_Field_double & get_storage_double(StorageId, int level) = 0;
  void coarse_solver(IJK_Field_float & x, const IJK_Field_float & b);

  void resoudre_avec_gcp(IJK_Field_float & x, IJK_Field_float & b, IJK_Field_float & residu);
  void resoudre_avec_gmres(IJK_Field_float & x, IJK_Field_float & b, IJK_Field_float & residu);
  void resoudre_systeme_float(const Matrice_Base & a, const DoubleVect & b, DoubleVect & x); 
  void resoudre_systeme_IJK_float(const IJK_Field_double & rhs, IJK_Field_double & x);
  void resoudre_systeme_IJK_float(const IJK_Field_float & rhs, IJK_Field_float & x);
  void solve_ijk_in_storage_float();

  virtual void jacobi_residu(IJK_Field_float & x,
			     const IJK_Field_float *secmem, /* if null pointer, take secmem = 0 (to compute A*x) */
			     const int grid_level,
			     const int n_jacobi,
			     IJK_Field_float *residu) const = 0;

  virtual void coarsen(const IJK_Field_float & fine, IJK_Field_float & coarse, int fine_level) const = 0;
  virtual void interpolate_sub_shiftk(const IJK_Field_float & coarse, IJK_Field_float & fine, int fine_level) const = 0;
  
  double multigrille_(IJK_Field_float & x, const IJK_Field_float & b, IJK_Field_float & residu, 
		      int grid_level, int iter_number);

  virtual IJK_Field_float & get_storage_float(StorageId, int level) = 0;
  double relax_jacobi(int level) const { 
    int i;
    if (level < relax_jacobi_.size_array())
      i = level;
    else
      i = relax_jacobi_.size_array()-1;
    return relax_jacobi_[i]; 
  }
  
  const int precision_double_;
  const int precision_float_;
  const int precision_mix_;
  int solver_precision_;
 private:


  // Solveur systeme pour la grille grossiere
  SolveurSys solveur_grossier_;
  Matrice_Grossiere coarse_matrix_;

  // PARAMETRES MULTIGRILLE
  // Coefficient de relaxation du lisseur Jacobi
  //  d'apres les publis, l'optimum est 0.71 en 2D et 0.65 en 3D
  ArrOfDouble relax_jacobi_;
  // Pour chaque niveau, nombre d'iterations lissage pre/post et multigrille
  ArrOfInt pre_smooth_steps_;
  ArrOfInt smooth_steps_;
  ArrOfInt nb_full_mg_steps_; // number of multigrid iterations at each level

  double seuil_; // exit from main multigrid iterations at fine level if below seuil_
  int max_iter_gcp_; // if > 0, use gcp solver and fix max iterations
  int n_krilov_; // number of krilov vectors in gmres
  int max_iter_gmres_;
  int solv_jacobi_;
  int impr_;
  int impr_gmres_;
  int max_iter_mixed_solver_; // maximum number of iterations in mixed precision solver
  int check_residu_;
};

#endif