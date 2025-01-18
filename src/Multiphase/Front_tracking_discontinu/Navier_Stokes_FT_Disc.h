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
// File:        Navier_Stokes_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_FT_Disc_included
#define Navier_Stokes_FT_Disc_included

#include <Navier_Stokes_Turbulent.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Champ_Don.h>
#include <TRUST_Ref.h>

class Probleme_FT_Disc_gen;
class Navier_Stokes_FT_Disc_interne;
class Maillage_FT_Disc;
class Fluide_Diphasique;

class Navier_Stokes_FT_Disc : public Navier_Stokes_Turbulent
{
  Declare_instanciable_sans_constructeur(Navier_Stokes_FT_Disc);
public:

  Navier_Stokes_FT_Disc();
  //
  // Methodes surchargees
  //
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  const Milieu_base& milieu() const override;
  Milieu_base&        milieu() override;
  void                associer_pb_base(const Probleme_base& probleme) override;
  void                discretiser() override;
  int              preparer_calcul() override;
  void                preparer_pas_de_temps();
  void mettre_a_jour(double temps) override;
  void                calculer_la_pression_en_pa() override;
  DoubleTab&          derivee_en_temps_inco(DoubleTab& vpoint) override;
  void                projeter() override;
  virtual const Champ_base& calculer_div_normale_interface();
  void correct_at_exit_bad_gradient(DoubleTab& u0) const;
  void calculer_delta_u_interface(Champ_base& u0, int phase_pilote, int ordre);
  const Champ_Don& diffusivite_pour_transport() const override;

  virtual const Champ_base * get_delta_vitesse_interface() const;
  virtual const Fluide_Diphasique&     fluide_diphasique() const;

  void compute_boussinesq_additional_gravity(
    const Convection_Diffusion_Temperature_FT_Disc& eq,
    const Fluide_Diphasique& fluide_diphasique,
    const IntTab& face_voisins,
    const DoubleVect& volumes_entrelaces,
    const IntVect& orientation,
    const DoubleTab& indicatrice,
    const ArrOfDouble& g,
    DoubleTab& gravite_face) const;

  int is_terme_gravite_rhog() const;
  const Champ_Fonc& champ_rho_faces() const;
  const Champ_Fonc& get_num_compo() const;
  Champ_Fonc& get_num_compo();
  void reprendre_num_compo(Entree& is) ;

  virtual void calculer_dI_dt(DoubleVect& dI_dt); // const;
  const int& get_is_penalized() const;
  const int& get_new_mass_source() const;
  const DoubleTab& get_interfacial_area() const;
  DoubleTab& get_set_interfacial_area();  // Open access  in write-mode..
  const DoubleTab& get_mpoint() const;
  DoubleTab& get_set_mpoint(); // Open access to mpoint in write-mode...
  //void corriger_mpoint(); // Apply correction based on TCL model

  const SolveurSys& get_solveur_pression() const;

  const DoubleTab& get_force_pression_interf() const; // EB
  const DoubleTab& get_force_frottements_interf() const; // EB
  const DoubleTab& get_pression_interf() const; // EB
  const DoubleTab& get_force_tot_pression_interf() const;  // EB
  const DoubleTab& get_force_tot_frottements_interf() const;  // EB
  const DoubleVect& get_surface_tot_interf() const; // EB
  const DoubleTab& get_force_pression_tot_interf_stokes_th() const; // EB
  const DoubleTab& get_force_frottements_tot_interf_stokes_th() const; // EB
  const DoubleTab& get_force_pression_tot_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_force_frottements_tot_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_xx_interf() const; // EB
  const DoubleTab& get_sigma_xy_interf() const; // EB
  const DoubleTab& get_sigma_xz_interf() const; // EB
  const DoubleTab& get_sigma_yx_interf() const; // EB
  const DoubleTab& get_sigma_yy_interf() const; // EB
  const DoubleTab& get_sigma_yz_interf() const; // EB
  const DoubleTab& get_sigma_zx_interf() const; // EB
  const DoubleTab& get_sigma_zy_interf() const; // EB
  const DoubleTab& get_sigma_zz_interf() const; // EB

  DoubleTab& get_force_pression_interf();  // EB
  DoubleTab& get_force_frottements_interf();  // EB
  DoubleTab& get_pression_interf(); // EB
  DoubleTab& get_force_tot_pression_interf();  // EB
  DoubleTab& get_force_tot_frottements_interf();  // EB
  DoubleVect& get_surface_tot_interf(); // EB
  DoubleTab& get_force_pression_tot_interf_stokes_th(); // EB
  DoubleTab& get_force_frottements_tot_interf_stokes_th(); // EB
  DoubleTab& get_force_pression_tot_interf_stokes_th_dis(); // EB
  DoubleTab& get_force_frottements_tot_interf_stokes_th_dis(); // EB

  DoubleTab& get_sigma_xx_interf(); // EB
  DoubleTab& get_sigma_xy_interf(); // EB
  DoubleTab& get_sigma_xz_interf(); // EB
  DoubleTab& get_sigma_yx_interf(); // EB
  DoubleTab& get_sigma_yy_interf(); // EB
  DoubleTab& get_sigma_yz_interf(); // EB
  DoubleTab& get_sigma_zx_interf(); // EB
  DoubleTab& get_sigma_zy_interf(); // EB
  DoubleTab& get_sigma_zz_interf(); // EB

  const DoubleTab& get_force_pression_stokes_th() const; // EB
  const DoubleTab& get_force_frottements_stokes_th() const; // EB
  const DoubleTab& get_force_pression_stokes_th_dis() const; // EB
  const DoubleTab& get_force_frottements_stokes_th_dis() const; // EB
  const DoubleTab& get_pression_interf_stokes_th_dis() const; // EB

  const DoubleTab& get_sigma_xx_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_xy_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_xz_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_yx_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_yy_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_yz_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_zx_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_zy_interf_stokes_th_dis() const; // EB
  const DoubleTab& get_sigma_zz_interf_stokes_th_dis() const; // EB

  DoubleTab& get_force_pression_stokes_th(); // EB
  DoubleTab& get_force_frottements_stokes_th(); // EB
  DoubleTab& get_force_pression_stokes_th_dis(); // EB
  DoubleTab& get_force_frottements_stokes_th_dis(); // EB
  DoubleTab& get_pression_interf_stokes_th_dis(); // EB


  DoubleTab& get_sigma_xx_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_xy_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_xz_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_yx_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_yy_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_yz_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_zx_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_zy_interf_stokes_th_dis(); // EB
  DoubleTab& get_sigma_zz_interf_stokes_th_dis(); // EB

  const DoubleTab& get_sigma_xx_interf_stokes_th() const; // EB
  const DoubleTab& get_sigma_xy_interf_stokes_th() const; // EB
  const DoubleTab& get_sigma_xz_interf_stokes_th() const; // EB
  const DoubleTab& get_sigma_yy_interf_stokes_th() const; // EB
  const DoubleTab& get_sigma_yz_interf_stokes_th() const; // EB
  const DoubleTab& get_sigma_zz_interf_stokes_th() const; // EB

  DoubleTab& get_sigma_xx_interf_stokes_th(); // EB
  DoubleTab& get_sigma_xy_interf_stokes_th(); // EB
  DoubleTab& get_sigma_xz_interf_stokes_th(); // EB
  DoubleTab& get_sigma_yy_interf_stokes_th(); // EB
  DoubleTab& get_sigma_yz_interf_stokes_th(); // EB
  DoubleTab& get_sigma_zz_interf_stokes_th(); // EB


  const DoubleTab& get_dUdx_P1() const; // EB
  const DoubleTab& get_dUdy_P1() const; // EB
  const DoubleTab& get_dUdz_P1() const; // EB
  const DoubleTab& get_dVdx_P1() const; // EB
  const DoubleTab& get_dVdy_P1() const; // EB
  const DoubleTab& get_dVdz_P1() const; // EB
  const DoubleTab& get_dWdx_P1() const; // EB
  const DoubleTab& get_dWdy_P1() const; // EB
  const DoubleTab& get_dWdz_P1() const; // EB

  DoubleTab& get_dUdx_P1(); // EB
  DoubleTab& get_dUdy_P1(); // EB
  DoubleTab& get_dUdz_P1(); // EB
  DoubleTab& get_dVdx_P1(); // EB
  DoubleTab& get_dVdy_P1(); // EB
  DoubleTab& get_dVdz_P1(); // EB
  DoubleTab& get_dWdx_P1(); // EB
  DoubleTab& get_dWdy_P1(); // EB
  DoubleTab& get_dWdz_P1(); // EB


  const DoubleTab& get_dUdx_P2() const; // EB
  const DoubleTab& get_dUdy_P2() const; // EB
  const DoubleTab& get_dUdz_P2() const; // EB
  const DoubleTab& get_dVdx_P2() const; // EB
  const DoubleTab& get_dVdy_P2() const; // EB
  const DoubleTab& get_dVdz_P2() const; // EB
  const DoubleTab& get_dWdx_P2() const; // EB
  const DoubleTab& get_dWdy_P2() const; // EB
  const DoubleTab& get_dWdz_P2() const; // EB

  DoubleTab& get_dUdx_P2(); // EB
  DoubleTab& get_dUdy_P2(); // EB
  DoubleTab& get_dUdz_P2(); // EB
  DoubleTab& get_dVdx_P2(); // EB
  DoubleTab& get_dVdy_P2(); // EB
  DoubleTab& get_dVdz_P2(); // EB
  DoubleTab& get_dWdx_P2(); // EB
  DoubleTab& get_dWdy_P2(); // EB
  DoubleTab& get_dWdz_P2(); // EB

  const DoubleTab& get_dUdx_P1_th_dis() const; // EB
  const DoubleTab& get_dUdy_P1_th_dis() const; // EB
  const DoubleTab& get_dUdz_P1_th_dis() const; // EB
  const DoubleTab& get_dVdx_P1_th_dis() const; // EB
  const DoubleTab& get_dVdy_P1_th_dis() const; // EB
  const DoubleTab& get_dVdz_P1_th_dis() const; // EB
  const DoubleTab& get_dWdx_P1_th_dis() const; // EB
  const DoubleTab& get_dWdy_P1_th_dis() const; // EB
  const DoubleTab& get_dWdz_P1_th_dis() const; // EB

  DoubleTab& get_dUdx_P1_th_dis(); // EB
  DoubleTab& get_dUdy_P1_th_dis(); // EB
  DoubleTab& get_dUdz_P1_th_dis(); // EB
  DoubleTab& get_dVdx_P1_th_dis(); // EB
  DoubleTab& get_dVdy_P1_th_dis(); // EB
  DoubleTab& get_dVdz_P1_th_dis(); // EB
  DoubleTab& get_dWdx_P1_th_dis(); // EB
  DoubleTab& get_dWdy_P1_th_dis(); // EB
  DoubleTab& get_dWdz_P1_th_dis(); // EB

  const DoubleTab& get_dUdx_P2_th_dis() const; // EB
  const DoubleTab& get_dUdy_P2_th_dis() const; // EB
  const DoubleTab& get_dUdz_P2_th_dis() const; // EB
  const DoubleTab& get_dVdx_P2_th_dis() const; // EB
  const DoubleTab& get_dVdy_P2_th_dis() const; // EB
  const DoubleTab& get_dVdz_P2_th_dis() const; // EB
  const DoubleTab& get_dWdx_P2_th_dis() const; // EB
  const DoubleTab& get_dWdy_P2_th_dis() const; // EB
  const DoubleTab& get_dWdz_P2_th_dis() const; // EB

  DoubleTab& get_dUdx_P2_th_dis(); // EB
  DoubleTab& get_dUdy_P2_th_dis(); // EB
  DoubleTab& get_dUdz_P2_th_dis(); // EB
  DoubleTab& get_dVdx_P2_th_dis(); // EB
  DoubleTab& get_dVdy_P2_th_dis(); // EB
  DoubleTab& get_dVdz_P2_th_dis(); // EB
  DoubleTab& get_dWdx_P2_th_dis(); // EB
  DoubleTab& get_dWdy_P2_th_dis(); // EB
  DoubleTab& get_dWdz_P2_th_dis(); // EB

  const DoubleTab& get_dUdx_P1_th() const; // EB
  const DoubleTab& get_dUdy_P1_th() const; // EB
  const DoubleTab& get_dUdz_P1_th() const; // EB
  const DoubleTab& get_dVdx_P1_th() const; // EB
  const DoubleTab& get_dVdy_P1_th() const; // EB
  const DoubleTab& get_dVdz_P1_th() const; // EB
  const DoubleTab& get_dWdx_P1_th() const; // EB
  const DoubleTab& get_dWdy_P1_th() const; // EB
  const DoubleTab& get_dWdz_P1_th() const; // EB

  DoubleTab& get_dUdx_P1_th(); // EB
  DoubleTab& get_dUdy_P1_th(); // EB
  DoubleTab& get_dUdz_P1_th(); // EB
  DoubleTab& get_dVdx_P1_th(); // EB
  DoubleTab& get_dVdy_P1_th(); // EB
  DoubleTab& get_dVdz_P1_th(); // EB
  DoubleTab& get_dWdx_P1_th(); // EB
  DoubleTab& get_dWdy_P1_th(); // EB
  DoubleTab& get_dWdz_P1_th(); // EB

  const DoubleTab& get_dUdx_P2_th() const; // EB
  const DoubleTab& get_dUdy_P2_th() const; // EB
  const DoubleTab& get_dUdz_P2_th() const; // EB
  const DoubleTab& get_dVdx_P2_th() const; // EB
  const DoubleTab& get_dVdy_P2_th() const; // EB
  const DoubleTab& get_dVdz_P2_th() const; // EB
  const DoubleTab& get_dWdx_P2_th() const; // EB
  const DoubleTab& get_dWdy_P2_th() const; // EB
  const DoubleTab& get_dWdz_P2_th() const; // EB

  DoubleTab& get_dUdx_P2_th(); // EB
  DoubleTab& get_dUdy_P2_th(); // EB
  DoubleTab& get_dUdz_P2_th(); // EB
  DoubleTab& get_dVdx_P2_th(); // EB
  DoubleTab& get_dVdy_P2_th(); // EB
  DoubleTab& get_dVdz_P2_th(); // EB
  DoubleTab& get_dWdx_P2_th(); // EB
  DoubleTab& get_dWdy_P2_th(); // EB
  DoubleTab& get_dWdz_P2_th(); // EB


  const DoubleTab& get_U_P1() const; // EB
  const DoubleTab& get_U_P2() const; // EB
  const DoubleTab& get_U_P1_th() const; // EB
  const DoubleTab& get_U_P2_th() const; // EB
  const DoubleTab& get_U_P1_th_dis() const; // EB
  const DoubleTab& get_U_P2_th_dis() const; // EB

  DoubleTab& get_U_P1(); // EB
  DoubleTab& get_U_P2(); // EB
  DoubleTab& get_U_P1_th(); // EB
  DoubleTab& get_U_P2_th(); // EB
  DoubleTab& get_U_P1_th_dis(); // EB
  DoubleTab& get_U_P2_th_dis(); // EB
  const IntTab& get_list_elem_P1() const; // EB
  IntTab& get_list_elem_P1(); // EB
  const IntTab& get_list_elem_diph() const; // EB
  IntTab& get_list_elem_diph(); // EB
  const IntTab& get_list_elem_P1_all() const; // EB
  IntTab& get_list_elem_P1_all(); // EB
  void init_champs_forces_interf(); // EB
  void compute_num_compo(DoubleTab& num_compo,const DoubleTab& indicatrice) const; // EB

  void calcul_forces_interface();
  void calcul_forces_interface_stokes_th();
  void calcul_forces_interface_taylor_lagrange(DoubleVect& surface_tot_interf, DoubleTab& force_pression_tot_interf, DoubleTab& force_frottements_tot_interf);

  // Interpolation trilineaire de valeurs_champs aux facettes du maillage lagrangien. Valeurs_champs contient les infos aux faces du maillage eulerien.
  int trilinear_interpolation_face(const DoubleTab& indicatrice_faces, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu); // EB // on le declare public car on en a besoin dans Transport_Interfaces_FT_Disc:calculer_vitesse_transport_interpolee
  int trilinear_interpolation_elem(const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu);
  int trilinear_interpolation_elem(const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu, const int is_P2, const int discr);
  // Interpolation trilineaire de valeurs_champs aux sommets du maillage lagrangien. Valeurs_champs contient les infos aux faces du maillage eulerien.
  int trilinear_interpolation_face_sommets(const DoubleTab& indicatrice_faces, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu); // EB // on le declare public car on en a besoin dans Transport_Interfaces_FT_Disc:calculer_vitesse_transport_interpolee

protected:
  // Methode surchargee de Navier_Stokes_std :
  void discretiser_assembleur_pression() override;
  void associer_milieu_base(const Milieu_base& fluide) override;

  // Nouvelles methodes
  virtual const Probleme_FT_Disc_gen& probleme_ft() const;
  virtual Probleme_FT_Disc_gen&        probleme_ft() ;
  virtual void calculer_champ_forces_superficielles(const Maillage_FT_Disc& maillage,
                                                    const Champ_base& gradient_indicatrice,
                                                    Champ_base& potentiel_elements,
                                                    Champ_base& potentiel_faces,
                                                    Champ_base& champ);
  virtual void calculer_gradient_indicatrice(const Champ_base& indicatrice,
                                             const DoubleTab& distance_interface_sommets,
                                             Champ_base& gradient_i);

  REF(Probleme_FT_Disc_gen)  probleme_ft_;

  // Masse volumique calculee aux elements
  Champ_Fonc champ_rho_elem_;
  // Masse volumique calculee pour les volumes de controle de la vitesse
  // (pour division   v = (rho.v) / rho et pour matrice de pression)
  Champ_Fonc champ_rho_faces_;
  // Viscosite dynamique (calcul dans preparer_pas_de_temps)
  // champ du type requis pour l'operateur diffusion.
  Champ_Don champ_mu_;
  // Viscosite cinematique pour le calcul du pas de temps de diffusion
  Champ_Don champ_nu_;

protected:


private:
  const Navier_Stokes_FT_Disc_interne& variables_internes() const;
  Navier_Stokes_FT_Disc_interne& variables_internes();

  // Ne pas utiliser ce pointeur : utiliser variables_internes() a la place !
  Navier_Stokes_FT_Disc_interne *variables_internes_;

  double minx,maxx,pente;
  int is_repulsion;
};


#endif
