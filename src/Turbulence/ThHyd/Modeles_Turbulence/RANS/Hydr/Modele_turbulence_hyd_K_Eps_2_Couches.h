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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_turbulence_hyd_K_Eps_2_Couches.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_K_Eps_2_Couches_included
#define Modele_turbulence_hyd_K_Eps_2_Couches_included

#include <Transport_K_KEps.h>

/*! @brief Classe Modele_turbulence_hyd_K_Eps_2_Couches Cette classe represente le modele de turbulence (k,eps) pour les
 *
 *     equations de Navier-Stokes.
 *
 * @sa Modele_turbulence_hyd_base Modele_turbulence_hyd_LES_base
 */
class Modele_turbulence_hyd_K_Eps_2_Couches: public Modele_turbulence_hyd_RANS_K_Eps_base, public Modele_turbulence_hyd_RANS_Gen<Modele_turbulence_hyd_K_Eps_2_Couches>
{
  Declare_instanciable(Modele_turbulence_hyd_K_Eps_2_Couches);
public:

  void completer() override;
  int preparer_calcul() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double) override;
  inline int get_nbcouches() const;
  inline int get_yswitch() const;
  inline int get_nutswitch() const;
  inline int get_switch() const;
  inline int get_impr() const;
  void controler() override;
  inline Transport_K_KEps& get_set_eq_transport() override;
  inline const Transport_K_KEps& get_eq_transport() const override;

  Champ_Fonc_base& calculer_viscosite_turbulente(double temps);

private:
  void fill_turbulent_viscosity_tab(const int , const DoubleTab&, DoubleTab& );
};

inline Transport_K_KEps& Modele_turbulence_hyd_K_Eps_2_Couches::get_set_eq_transport()
{
  Transport_K_KEps& eq_transport = ref_cast(Transport_K_KEps,ptr_eq_transport_.valeur());
  return eq_transport;
}

inline const Transport_K_KEps& Modele_turbulence_hyd_K_Eps_2_Couches::get_eq_transport() const
{
  const Transport_K_KEps& eq_transport = ref_cast(Transport_K_KEps,ptr_eq_transport_.valeur());
  return eq_transport;
}

/*! @brief Renvoie le nombre couches utilises pour les lois de paroi Ce nombre est porte par l'equation de transport K_K-eps.
 *
 * @return (int) le nombre de couches
 */
inline int Modele_turbulence_hyd_K_Eps_2_Couches::get_nbcouches() const
{
  return get_eq_transport().get_nbcouches();
}

/*! @brief Renvoie le y* de switch entre les deux couches pour le modele a deux couches.
 *
 * (version const)
 *
 * @return (int) le nombre de couches
 */
inline int Modele_turbulence_hyd_K_Eps_2_Couches::get_yswitch() const
{
  return get_eq_transport().get_yswitch();
}

/*! @brief Renvoie 0 si on choisit le switch par y*, 1 par nu_t.
 *
 * (version const)
 *
 * @return (int)
 */
inline int Modele_turbulence_hyd_K_Eps_2_Couches::get_switch() const
{
  return get_eq_transport().get_switch();
}

/*! @brief Renvoie la valeur de nut/nu qui delimite les deux couches.
 *
 * (version const)
 *
 * @return (int) valeur limite de nu_t/nu
 */
inline int Modele_turbulence_hyd_K_Eps_2_Couches::get_nutswitch() const
{
  return get_eq_transport().get_nutswitch();
}

/*! @brief indique si on doit ecrire le domaine des 2 couches dans le .
 *
 * out. (version const)
 *
 * @return (int) 1 si on imprime , 0 sinon.
 */
inline int Modele_turbulence_hyd_K_Eps_2_Couches::get_impr() const
{
  return get_eq_transport().get_impr();
}

#endif /* Modele_turbulence_hyd_K_Eps_2_Couches_included */
