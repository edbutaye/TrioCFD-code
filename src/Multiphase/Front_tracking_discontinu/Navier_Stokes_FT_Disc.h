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
#include <Modele_Collision_FT.h>
//EB
#include <Domaine_VDF.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Transport_Interfaces_FT_Disc.h>

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
  //void ouvrir_fichier(SFichier& ,const Nom& , const int& , const Navier_Stokes_FT_Disc&);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int sauvegarder(Sortie&) const override; //EB
  int reprendre(Entree&) override; // EB
  void imprimer(Sortie& os) const override;
  virtual int impr_fpi(Sortie& os) const override; // EB
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

  // debut EB
  int trilinear_interpolation_face(const DoubleTab& indicatrice_faces, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu) const;
  int trilinear_interpolation_gradU_face(const DoubleTab& indicatrice_face, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu);
  int trilinear_interpolation_gradU_elem_P1(const DoubleTab& indicatrice_face, const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu);
  int trilinear_interpolation_gradU_elem(const DoubleTab& indicatrice_face, const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu);
  double calculer_viscosite_arete(int face1, int face2, int compo);
  inline double chercher_elem_voisins(const DoubleTab& indicatrice, DoubleVect& coord_elem_interp, IntVect& elem_voisins, const int sauv_list_P1=0, const int num_fa7=-1); // sauv_list_P1 : on sauvegarde la liste des elements auxquels appartiennent les points P1
  inline void chercher_faces_voisines (DoubleVect& coord_elem_interp, IntVect& faces_voisines, int orientation);
  inline void chercher_faces_voisines_xyz (DoubleVect& coord_elem_interp, IntTab& faces_voisines);
  // fin EB

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

  void calculer_champ_forces_collisions(const DoubleTab& indicatrice, DoubleTab& valeurs_champ,  const Transport_Interfaces_FT_Disc& eq_transport,Transport_Interfaces_FT_Disc& eq_transport_non_const, REF(Transport_Interfaces_FT_Disc)& refeq_transport, const Maillage_FT_Disc& maillage); // HMS
  void calculer_correction_trainee(DoubleTab& valeurs_champ, const Transport_Interfaces_FT_Disc& eq_transport,Transport_Interfaces_FT_Disc& eq_transport_non_const, REF(Transport_Interfaces_FT_Disc)& refeq_transport, const Maillage_FT_Disc& maillage);// EB
  void init_positions_vitesses_FT();

};

// Description
/*! @brief
*  Dans cette fonction, les conditions if qui retournent 0 permettent de ne pas passer dans le assert qui ferait planter le programme.
*  On teste donc si on a acces aux elements distants, si ce n'est pas le cas, la fonction renvoie 0.
*  Ainsi, si chercher_elem_voisins renvoie 0, on ne calculera pas la force de pression.
*  En entree : l'indicatrice, les coodonnees xyz coord_elem_interp
*  Le tableau elem_voisins est rempli avec les 8 elements euleriens les plus proches du point coord_elem_interp. Ces elements
*  seront utilises pour les interpolations trilineaires.
*/
inline double Navier_Stokes_FT_Disc::chercher_elem_voisins(const DoubleTab& indicatrice, DoubleVect& coord_elem_interp, IntVect& elem_voisins, const int sauv_list_P1, const int num_fa7)
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  const Domaine& domaine = domaine_vdf.domaine();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());

  DoubleVect coord_elem_eulerien(dimension);
  IntTab& list_elem_P1=get_list_elem_P1();
  int elem_eulerien=domaine.chercher_elements(coord_elem_interp(0), coord_elem_interp(1),coord_elem_interp(2));
  if (sauv_list_P1)  list_elem_P1(num_fa7,0)=elem_eulerien;
  if (elem_eulerien<0)
    {
      elem_voisins=-1;
      return -1;
    }
  //if (indicatrice(elem_eulerien)<1) Cerr << "Elem eulerien " <<elem_eulerien << " d'indicatrice " << indicatrice(elem_eulerien) << finl;

  for (int dim=0; dim<dimension; dim++) coord_elem_eulerien(dim)=domaine_vdf.xp(elem_eulerien,dim);
  IntVect direction_interp(dimension); // Pour chaque direction, on regarde de quel cote de l'element le point se trouve. 0 : a gauche (pos_i_point<=pos_i_elem), 1 : a droite(pos_i_point>pos_i_elem)
  IntVect faces_elem_interp(2*dimension);
  for (int dim=0; dim<dimension; dim++)
    {
      faces_elem_interp(dim)=domaine_vdf.elem_faces_pour_interp(elem_eulerien,dim);
      faces_elem_interp(dimension+dim)=domaine_vdf.elem_faces_pour_interp(elem_eulerien,dimension+dim);
      if (coord_elem_interp(dim)<=coord_elem_eulerien(dim))
        {
          direction_interp(dim)=0;
        }
      else
        {
          direction_interp(dim)=1;
        }
    }

  // Les 8 elements forment un grand cube avec 2 cubes dans chaque direction (1 cube = 1 elem).
  // Soit un cube de cote 1. Chaque sommet represente un element voisin. On remplit la liste elem_voisins de commencant
  // par les coordonnees z=0, y=0 puis z=0, y=1, puis z=1, y=0 puis z=1, y=1
  // On a ainsi (0,0,0) ; (1,0,0) ; (0,1,0) ; (1,1,0) ; (0,0,+-1) ; (1,0,+-1) ; (0,1,+-1) ; (1,1,+-1)


  if (direction_interp(0)==1)
    {
      if (direction_interp(1)==1)
        {
          elem_voisins(0)=elem_eulerien;
          elem_voisins(1)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(0+dimension),1);
          elem_voisins(2)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(1+dimension),1);
          elem_voisins(3)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(1),1+dimension),1);
        }
      else
        {
          elem_voisins(0)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(1),0);
          elem_voisins(1)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(0),0+dimension),1);
          elem_voisins(2)=elem_eulerien;
          elem_voisins(3)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(0+dimension),1);

        }
    }
  else
    {
      if (direction_interp(1)==1)
        {
          elem_voisins(0)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(0),0);
          elem_voisins(1)=elem_eulerien;
          elem_voisins(2)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(0),1+dimension),1);
          elem_voisins(3)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(1+dimension),1);

        }
      else
        {
          elem_voisins(1)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(1),0);
          elem_voisins(2)=domaine_vf.face_voisins_pour_interp(faces_elem_interp(0),0);
          elem_voisins(3)=elem_eulerien;
          elem_voisins(0)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(1),0),0);

        }
    }

  if (direction_interp(2)==1)
    {
      elem_voisins(4)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(0),2+dimension),1);
      elem_voisins(5)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(1),2+dimension),1);
      elem_voisins(6)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(2),2+dimension),1);
      elem_voisins(7)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(3),2+dimension),1);
    }
  else
    {
      elem_voisins(4)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(0),2),0);
      elem_voisins(5)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(1),2),0);
      elem_voisins(6)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(2),2),0);
      elem_voisins(7)=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_voisins(3),2),0);
      for (int i=0; i<4; i++)
        {
          int tmp=elem_voisins(i);
          elem_voisins(i)=elem_voisins(i+4);
          elem_voisins(i+4)=tmp;
        }
    }
  return (indicatrice(elem_eulerien));
}
/*! @brief meme principe que Navier_Stokes_FT_Disc::chercher_elem_voisins mais pour les faces
 * identifie 8 faces voisines d'orientation "orientation"
 */
inline void Navier_Stokes_FT_Disc::chercher_faces_voisines (DoubleVect& coord_elem_interp, IntVect& faces_voisines, int orientation)
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  const Domaine& domaine = domaine_vdf.domaine();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  DoubleVect coord_elem_eulerien(dimension);

  int elem_eulerien=domaine.chercher_elements(coord_elem_interp(0), coord_elem_interp(1),coord_elem_interp(2));
  if (elem_eulerien<0)
    {
      faces_voisines=-1;
      return;
    }
  for (int dim=0; dim<dimension; dim++) coord_elem_eulerien(dim)=domaine_vdf.xp(elem_eulerien,dim);
  IntVect direction_interp(dimension); // Pour chaque direction, on regarde de quel cote de l'element le point se trouve. 0 : a gauche (pos_i_point<=pos_i_elem), 1 : a droite(pos_i_point>pos_i_elem)
  IntVect faces_elem_interp(2*dimension);

  for (int dim=0; dim<dimension; dim++)
    {
      faces_elem_interp(dim)=domaine_vdf.elem_faces_pour_interp(elem_eulerien,dim);
      faces_elem_interp(dimension+dim)=domaine_vdf.elem_faces_pour_interp(elem_eulerien,dimension+dim);
      if (coord_elem_interp(dim)<=coord_elem_eulerien(dim))
        {
          direction_interp(dim)=0;
        }
      else
        {
          direction_interp(dim)=1;
        }
    }

  if (orientation==0)
    {
      if (direction_interp(1)==1)
        {
          faces_voisines(0)=faces_elem_interp(orientation);
          faces_voisines(1)=faces_elem_interp(orientation+dimension);
          faces_voisines(2)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(1+dimension),1),orientation);
          faces_voisines(3)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(1+dimension),1),orientation+dimension);
        }
      else
        {
          faces_voisines(0)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(1),0),orientation);
          faces_voisines(1)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(1),0),orientation+dimension);
          faces_voisines(2)=faces_elem_interp(orientation);
          faces_voisines(3)=faces_elem_interp(orientation+dimension);
        }
    }
  if (orientation==1)
    {
      if (direction_interp(0)==1)
        {
          faces_voisines(0)=faces_elem_interp(orientation);
          faces_voisines(1)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0+dimension),1),orientation);
          faces_voisines(2)=faces_elem_interp(orientation+dimension);
          faces_voisines(3)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0+dimension),1),orientation+dimension);
        }
      else
        {
          faces_voisines(0)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0),0),orientation);
          faces_voisines(1)=faces_elem_interp(orientation);
          faces_voisines(2)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0),0),orientation+dimension);
          faces_voisines(3)=faces_elem_interp(orientation+dimension);
        }
    }
  if (direction_interp(2)==1 && orientation!=2)
    {
      faces_voisines(4)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(0),1),2+dimension),1),orientation);
      faces_voisines(5)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(1),1),2+dimension),1),orientation);
      faces_voisines(6)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(2),1),2+dimension),1),orientation);
      faces_voisines(7)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(3),1),2+dimension),1),orientation);
    }
  else if (direction_interp(2)==0 && orientation!=2)
    {
      faces_voisines(4)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(0),1),2),0),orientation);
      faces_voisines(5)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(1),1),2),0),orientation);
      faces_voisines(6)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(2),1),2),0),orientation);
      faces_voisines(7)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(3),1),2),0),orientation);
      for (int i=0; i<4; i++)
        {
          int tmp=faces_voisines(i);
          faces_voisines(i)=faces_voisines(i+4);
          faces_voisines(i+4)=tmp;
        }
    }

  if (orientation==2)
    {
      if (direction_interp(0)==1)
        {
          faces_voisines(0)=faces_elem_interp(orientation);
          faces_voisines(1)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0+dimension),1),orientation);
          faces_voisines(4)=faces_elem_interp(orientation+dimension);
          faces_voisines(5)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0+dimension),1),orientation+dimension);
        }
      else
        {
          faces_voisines(0)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0),0),orientation);
          faces_voisines(1)=faces_elem_interp(orientation);
          faces_voisines(4)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_elem_interp(0),0),orientation+dimension);
          faces_voisines(5)=faces_elem_interp(orientation+dimension);

        }

      if (direction_interp(1)==1)
        {
          faces_voisines(2)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(0),1),1+dimension),1),orientation);
          faces_voisines(3)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(1),1),1+dimension),1),orientation);
          faces_voisines(6)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(4),1),1+dimension),1),orientation);
          faces_voisines(7)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(5),1),1+dimension),1),orientation);
        }
      else
        {
          faces_voisines(2)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(0),1),1),0),orientation);
          faces_voisines(3)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(1),1),1),0),orientation);
          faces_voisines(6)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(4),1),1),0),orientation);
          faces_voisines(7)=domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(domaine_vf.face_voisins_pour_interp(faces_voisines(5),1),1),0),orientation);

          int tmp0=faces_voisines(0);
          int tmp1=faces_voisines(1);
          int tmp4=faces_voisines(4);
          int tmp5=faces_voisines(5);

          faces_voisines(0)=faces_voisines(2);
          faces_voisines(1)=faces_voisines(3);
          faces_voisines(4)=faces_voisines(6);
          faces_voisines(5)=faces_voisines(7);

          faces_voisines(2)=tmp0;
          faces_voisines(3)=tmp1;
          faces_voisines(6)=tmp4;
          faces_voisines(7)=tmp5;

        }
    }
}

/*! @brief Pour un point de coordonnees coord_elem_interp, identifie les 8 faces euleriennes les plus proches,
 * pour chaque orientation (8 faces de normale x, 8 de normale y et 8 de normale z).
 * voir aussi Navier_Stokes_FT_Disc::chercher_faces_voisines
 */
inline void Navier_Stokes_FT_Disc::chercher_faces_voisines_xyz (DoubleVect& coord_elem_interp, IntTab& faces_voisines)
{
  int nb_faces_voisines=8;
  IntVect faces_voisines_x(nb_faces_voisines),faces_voisines_y(nb_faces_voisines),faces_voisines_z(nb_faces_voisines);
  chercher_faces_voisines(coord_elem_interp,faces_voisines_x,0);
  chercher_faces_voisines(coord_elem_interp,faces_voisines_y,1);
  chercher_faces_voisines(coord_elem_interp,faces_voisines_z,2);

  for (int i=0; i<8; i++)
    {
      faces_voisines(0,i)=faces_voisines_x(i);
      faces_voisines(1,i)=faces_voisines_y(i);
      faces_voisines(2,i)=faces_voisines_z(i);
    }

}
#endif
