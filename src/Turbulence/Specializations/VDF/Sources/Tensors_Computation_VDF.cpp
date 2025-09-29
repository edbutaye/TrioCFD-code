/****************************************************************************
 * Copyright (c) 2025, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Tensors_Computation_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Tensors_Computation_VDF.h>
#include <TRUSTTab.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Champ_Face_VDF.h>
#include <Probleme_base.h>

using TCV = Tensors_Computation_VDF;

void TCV::compute_enstrophy(const Domaine_VDF& dom_VDF,
                            const Domaine_Cl_VDF& dom_BC_VDF,
                            const DoubleTab& velocity,
                            Champ_Face_VDF& ch_vit,
                            Navier_Stokes_Turbulent& ns_turb,
                            DoubleTab& enstrophy) const
{
  // Get gradient of velocity
  const int nb_elem_tot = dom_VDF.nb_elem_tot();
  const int dimension = Objet_U::dimension;
  DoubleTab gradient_elem(nb_elem_tot, dimension*dimension);
  gradient_elem = 0.;

  if (ns_turb.has_champ("gradient_vitesse")) //if the gradient is already computed (post_processing)
    gradient_elem.ref(ns_turb.get_champ("gradient_vitesse").valeurs());
  else
    ch_vit.calcul_grad_u(velocity, gradient_elem, dom_BC_VDF);

  antisym_enstrophy(dom_VDF, gradient_elem, enstrophy);
}

void TCV::compute_strain_invariant(const Domaine_VDF& dom_VDF,
                                   const Domaine_Cl_VDF& dom_BC_VDF,
                                   const DoubleTab& velocity,
                                   Champ_Face_VDF& ch_vit,
                                   Navier_Stokes_Turbulent& ns_turb,
                                   DoubleTab& strain_invariant) const
{
  // Get gradient of velocity
  const int nb_elem_tot = dom_VDF.nb_elem_tot();
  const int dimension = Objet_U::dimension;
  DoubleTab gradient_elem(nb_elem_tot, dimension*dimension);
  gradient_elem = 0.;

  if (ns_turb.has_champ("gradient_vitesse")) //if the gradient is already computed (post_processing)
    gradient_elem.ref(ns_turb.get_champ("gradient_vitesse").valeurs());
  else
    ch_vit.calcul_grad_u(velocity, gradient_elem, dom_BC_VDF);

  antisym_strain_invariant(dom_VDF, gradient_elem, strain_invariant);
}

void TCV::antisym_strain_invariant(const Domaine_VDF& dom_VDF,
                                   const DoubleTab& gradient_elem,
                                   DoubleTab& strain_invariant) const
{
  const int nb_elem = dom_VDF.nb_elem();
  for (int elem=0; elem<nb_elem; elem++)
    {
      double SijSij=0;
      if (Objet_U::dimension==2)
        {
          const double du_dx= gradient_elem(elem,0);
          const double du_dy= gradient_elem(elem,1);
          const double dv_dx= gradient_elem(elem,2);
          const double dv_dy= gradient_elem(elem,3);
          SijSij = du_dx*du_dx
                   +0.5*(du_dy+dv_dx)*(du_dy+dv_dx)
                   +dv_dy*dv_dy ;
        }
      if(Objet_U::dimension==3)
        {
          const double du_dz= gradient_elem(elem,2);
          const double dv_dz= gradient_elem(elem,5);
          const double dw_dx= gradient_elem(elem,6);
          const double dw_dy= gradient_elem(elem,7);
          const double dw_dz= gradient_elem(elem,8);

          SijSij += 0.5*(du_dz+dw_dx)*(du_dz+dw_dx)
                    + 0.5*(dv_dz+dw_dy)*(dv_dz+dw_dy)
                    + dw_dz*dw_dz ;
        }
      strain_invariant(elem) = sqrt(2*SijSij);
    }
}

void TCV::antisym_enstrophy(const Domaine_VDF& dom_VDF,
                            const DoubleTab& gradient_elem,
                            DoubleTab& enstrophy) const
{
  const int nb_elem = dom_VDF.nb_elem();
  for (int elem=0; elem<nb_elem; elem++)
    {
      const double du_dy= gradient_elem(elem,1);
      const double dv_dx= gradient_elem(elem,2);

      double tmp =  0.5*((du_dy - dv_dx)*(du_dy - dv_dx) +
                         (dv_dx - du_dy)*(dv_dx - du_dy));

      if(Objet_U::dimension==3)
        {
          const double du_dz= gradient_elem(elem,2);
          const double dv_dz= gradient_elem(elem,5);
          const double dw_dx= gradient_elem(elem,6);
          const double dw_dy= gradient_elem(elem,7);

          tmp += 0.5*((du_dz - dw_dx)*(du_dz - dw_dx) +
                      (dv_dz - dw_dy)*(dv_dz - dw_dy) +
                      (dw_dx - du_dz)*(dw_dx - du_dz) +
                      (dw_dy - dv_dz)*(dw_dy - dv_dz));


        }
      enstrophy(elem) = sqrt(tmp);
    }
}

