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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_K_Omega_VDF_Elem_base.h
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_Omega_VDF_Elem_base_included
#define Source_Transport_K_Omega_VDF_Elem_base_included

#include <Modele_turbulence_hyd_K_Omega.h>
#include <Calcul_Production_K_VDF.h>
#include <Source_Transport_VDF_Elem_base.h>
#include <Domaine_Cl_VDF.h>
#include <Domaine_VDF.h>
#include <TRUST_Ref.h>

class Source_Transport_K_Omega_VDF_Elem_base: public Source_Transport_VDF_Elem_base
{
  Declare_base_sans_constructeur(Source_Transport_K_Omega_VDF_Elem_base);
public:
  Source_Transport_K_Omega_VDF_Elem_base() { }

  DoubleTab& calculer(DoubleTab&) const override;
  void associer_domaines(const Domaine_dis_base&, const Domaine_Cl_dis_base&) override;
  void associer_pb(const Probleme_base&) override;

  inline int has_interface_blocs() const override { return 1; }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const override {}
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override = 0;

  inline void mettre_a_jour(double t) override { Calcul_Production_K_VDF::mettre_a_jour(t); }

protected:
  DoubleTab& ajouter_komega(DoubleTab&) const;

  // Constants for the classic k-omega model Wilcox 1988
  static constexpr double BETA_K = 0.09; // Cmu or BETA_STAR, but clearer with _K
  static constexpr double BETA_OMEGA = 3./40.; // BETA
  static constexpr double SIGMA_K = 0.5; // SIGMA_STAR
  static constexpr double SIGMA_OMEGA = 0.5; // SIGMA
  static constexpr double ALPHA_OMEGA = 5./9.; // ALPHA
};

#endif /* Source_Transport_K_Omega_VDF_Elem_base_included */
