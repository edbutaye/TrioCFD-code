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
// File:        Transport_K_Omega_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Transport_K_Omega_base_included
#define Transport_K_Omega_base_included

#include <Modele_turbulence_hyd_RANS_K_Omega_base.h>
#include <Transport_2eq_base.h>
#include <TRUST_Ref.h>


class Modele_turbulence_hyd_RANS_K_Omega_base;
class Champ_Inc_base;
class Milieu_base;

/*! @brief Classe Transport_K_Omega_base Classe de base pour les equations
 *
 *     de transport des modeles k_Omega.
 *
 */
class Transport_K_Omega_base: public Transport_2eq_base
{
  Declare_base(Transport_K_Omega_base);

public:

  virtual void associer_modele_turbulence(const Modele_turbulence_hyd_RANS_K_Omega_base& )=0;
  void discretiser() override;

  int controler_K_Omega();
  void valider_iteration() override;
  inline const Champ_Inc_base& inconnue() const override;
  inline Champ_Inc_base& inconnue() override;
  inline const Modele_turbulence_hyd_RANS_K_Omega_base& modele_turbulence() const;
  inline Modele_turbulence_hyd_RANS_K_Omega_base& modele_turbulence();

  void get_position_cells(Nom&, int&);
  void get_position_faces(Nom&, int&);


protected:

  OWN_PTR(Champ_Inc_base) le_champ_K_Omega;
  REF(Modele_turbulence_hyd_RANS_K_Omega_base) mon_modele;
};

/*! @brief Renvoie le champ inconnue de l'equation.
 *
 * Un champ vecteur contenant K et epsilon.
 *
 * @return (Champ_Inc_base&) le champ inconnue de l'equation
 */
inline Champ_Inc_base& Transport_K_Omega_base::inconnue() { return le_champ_K_Omega; }


/*! @brief Renvoie le champ inconnue de l'equation.
 *
 * Un champ vecteur contenant K et epsilon.
 *     (version const)
 *
 * @return (Champ_Inc_base&) le champ inconnue de l'equation
 */
inline const Champ_Inc_base& Transport_K_Omega_base::inconnue() const { return le_champ_K_Omega; }

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * (version const)
 *
 * @return (Modele_turbulence_hyd_K_Omega&) le modele de turbulence associe a l'equation
 */
inline const Modele_turbulence_hyd_RANS_K_Omega_base& Transport_K_Omega_base::modele_turbulence() const
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * @return (Modele_turbulence_hyd_K_Omega&) le modele de turbulence associe a l'equation
 */
inline Modele_turbulence_hyd_RANS_K_Omega_base& Transport_K_Omega_base::modele_turbulence()
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

#endif
