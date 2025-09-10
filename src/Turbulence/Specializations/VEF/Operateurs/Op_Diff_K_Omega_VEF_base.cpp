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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_K_Omega_VEF_base.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Omega_VEF_base.h>
#include <Modele_turbulence_hyd_K_Omega.h>

Implemente_instanciable_sans_constructeur(Op_Diff_K_Omega_VEF_base,
                                          "Op_Diff_K_Omega_VEF_base",
                                          Op_Dift_VEF_base);

Sortie& Op_Diff_K_Omega_VEF_base::printOn(Sortie& s) const
{
  return s << que_suis_je();
}

Entree& Op_Diff_K_Omega_VEF_base::readOn(Entree& s)
{
  return s;
}

void Op_Diff_K_Omega_VEF_base::completer()
{
  Operateur_base::completer();

  Op_Dift_VEF_base::completer();

  const RefObjU& modele_turbulence = equation().get_modele(TURBULENCE);
  const Modele_turbulence_hyd_K_Omega& mod_turb = ref_cast(Modele_turbulence_hyd_K_Omega, modele_turbulence.valeur());

  turbulence_model = mod_turb;

  Sigma_K1_ = mod_turb.get_Sigma_K1();
  Sigma_K2_ = mod_turb.get_Sigma_K2();
  Sigma_OMEGA1_ = mod_turb.get_Sigma_Omega1();
  Sigma_OMEGA2_ = mod_turb.get_Sigma_Omega2();
  if (mod_turb.is_SST_or_BSL())
    {
      is_SST_or_BSL_=1;
      associer_tab_F1(mod_turb.get_tabF1());
    }
  else
    {
      const Domaine_VF& domain_VF = ref_cast(Domaine_VF,mod_turb.get_eq_transport().domaine_dis());
      const int nb_faces_tot = domain_VF.nb_faces_tot();
      initialize_tab_F1(nb_faces_tot);
    }

  const Domaine_VEF& domain_VEF = le_dom_vef.valeur();
  nu_turb_m_.resize(domain_VEF.nb_elem_tot(), 2);
  domain_VEF.domaine().creer_tableau_elements(nu_turb_m_,RESIZE_OPTIONS::NOCOPY_NOINIT);
  nu_turb_m_=0.;
}
