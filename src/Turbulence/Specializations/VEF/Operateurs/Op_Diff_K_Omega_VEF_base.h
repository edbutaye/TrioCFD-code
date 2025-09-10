/****************************************************************************
* Copyright (c) 2023, CEA
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
/*! @brief class Op_Diff_K_Omega_VEF_Face Cette classe represente l'operateur de diffusion turbulente
 *
 *    La discretisation est VEF
 *   Les methodes pour l'implicite sont codees.
 *
 */
#ifndef Op_Diff_K_Omega_VEF_base_included
#define Op_Diff_K_Omega_VEF_base_included

#include <Op_Dift_VEF_base.h>
#include <Op_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <K_Omega_Eps_constants.h>

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_K_Omega_VEF_base
//
//////////////////////////////////////////////////////////////////////////////

class Op_Diff_K_Omega_VEF_base : public Op_Dift_VEF_base
{

  Declare_instanciable_sans_constructeur(Op_Diff_K_Omega_VEF_base);

public:

  Op_Diff_K_Omega_VEF_base(double Sigma_K = SIGMA_K,
                           double Sigma_Omega = SIGMA_OMEGA,
                           double Sigma_K1 = SIGMA_K1,
                           double Sigma_K2 = SIGMA_K2,
                           double Sigma_OMEGA1 = SIGMA_OMEGA1,
                           double Sigma_OMEGA2 = SIGMA_OMEGA2
                          ) : Sigma_K_(Sigma_K), Sigma_Omega_(Sigma_Omega), Sigma_K1_(Sigma_K1), Sigma_K2_(Sigma_K2), Sigma_OMEGA1_(Sigma_OMEGA1), Sigma_OMEGA2_(Sigma_OMEGA2)
  { }

  void completer() override;
  inline void associer_tab_F1(const DoubleTab& tabF1)
  {
    tab_F1_face_.ref(tabF1);
  }
  inline void initialize_tab_F1(const int nb_faces_tot)
  {
    tab_F1_face_.resize(nb_faces_tot);
    tab_F1_face_=0;
  }

protected:
  double Sigma_K_;
  double Sigma_Omega_;
  double Sigma_K1_,Sigma_K2_,Sigma_OMEGA1_,Sigma_OMEGA2_;
  int is_SST_or_BSL_=0;
  DoubleTab tab_F1_face_;
  DoubleTab nu_turb_m_;
  OBS_PTR(Modele_turbulence_hyd_K_Omega) turbulence_model;
};

#endif
