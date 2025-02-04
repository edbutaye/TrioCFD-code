/****************************************************************************
* Copyright (c) 2024, CEA
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
// File:        Navier_Stokes_FT_Disc_interne.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_FT_Disc_interne_included
#define Navier_Stokes_FT_Disc_interne_included

class Navier_Stokes_FT_Disc_interne
{
public:
  Navier_Stokes_FT_Disc_interne() :
    correction_courbure_ordre_(0), // Par defaut, pas de prise en compte de la courbure pour corriger le champ etendu delta_vitesse
    mpoint_inactif(0),   // Par defaut, mpoint cree un saut de vitesse
    mpointv_inactif(0),   // Par defaut, mpointv  cree un saut de vitesse
    matrice_pression_invariante(0),   // Par defaut, recalculer la matrice pression
    clipping_courbure_interface(1e40),   // Par defaut, pas de clipping
    flag_correction_trainee_(0), // EB : Par defaut, pas de correction de trainee
    alpha_correction_trainee_(0.), // EB
    beta_correction_trainee_(0.), // EB
    faces_diphasiques_(1), // EB
    extension_reynolds_(0), // EB
    proportionnel_(0), // EB
    terme_gravite_(GRAVITE_GRAD_I),   // Par defaut terme gravite ft sans courants parasites
    is_explicite(1),                  // Par defaut, calcul explicite de vpoint etape predicition
    is_boussinesq_(0),                // Par defaut, l'hypothese de Boussinesq n'est pas utilisee pour la flottabilite dans les phases.
    new_mass_source_(0),              // Par defaut, on utilise la methode historique pour imposer le saut de vitesse du changement de phase.
    type_interpol_indic_pour_dI_dt_(INTERP_STANDARD), // Default is the historical interpolation
    OutletCorrection_pour_dI_dt_(NO_CORRECTION),   // Default is the historical
    is_penalized(0),                  // Par defaut, pas de penalisation L2 du forcage
    eta(1.0),                         // Par defaut, coefficient de penalisation L2 = 1.
    p_ref_pena(-1.e40),               // Par defaut, pas de penalisation L2 de la pression sinon valeur reference
    is_pfl_flottant(0),               // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
    x_pfl_imp(-1.e40),                // Par defaut, x, y, z du point de modification de la pression fluide
    y_pfl_imp(-1.e40), z_pfl_imp(-1.e40)
  { }

  int correction_courbure_ordre_;
  int mpoint_inactif;
  int mpointv_inactif;

  OWN_PTR(Champ_Fonc_base)  second_membre_projection;
  OWN_PTR(Champ_Fonc_base)  second_membre_projection_jump_;
  OWN_PTR(Champ_Fonc_base)  derivee_u_etoile;
  OWN_PTR(Champ_Fonc_base)  gradient_pression;
  OWN_PTR(Champ_Fonc_base)  terme_diffusion;
  OWN_PTR(Champ_Fonc_base)  terme_convection;
  OWN_PTR(Champ_Fonc_base)  terme_source;
  OWN_PTR(Champ_Fonc_base)  terme_source_interfaces;
  OWN_PTR(Champ_Fonc_base)  indicatrice_p1b;
  OWN_PTR(Champ_Fonc_base)  gradient_indicatrice;
  OWN_PTR(Champ_Fonc_base)  potentiel_faces;
  OWN_PTR(Champ_Fonc_base)  potentiel_elements;
  // delta_u_interface = la partie "saut de vitesse" du champ de vitesse a l'interface
  OWN_PTR(Champ_Inc_base) delta_u_interface;
  OWN_PTR(Champ_Fonc_base)  laplacien_d;
  OWN_PTR(Champ_Fonc_base)  mpoint;
  OWN_PTR(Champ_Fonc_base)  mpoint_vap;
  // Variation temporelle indicatrice de phase
  OWN_PTR(Champ_Fonc_base)  derivee_temporelle_indicatrice;
  OWN_PTR(Champ_Fonc_base)  ai; // Eulerian interfacial area.
  OWN_PTR(Champ_Inc_base) vitesse_jump0_; // Extended Velocity of phase 0.

  LIST(OBS_PTR(Champ_base)) liste_champs_compris;
  OWN_PTR(Champ_Fonc_base) terme_source_collisions; // HMS
  OWN_PTR(Champ_Fonc_base)  num_compo; // HMS
  OWN_PTR(Champ_Fonc_base) vitesse_stokes_th_; // EB
  OWN_PTR(Champ_Fonc_base) pression_stokes_th_; // EB
  OWN_PTR(Champ_Fonc_base) terme_correction_trainee; // EB

  // Si matrice_pression_invariante != 0,
  //   on ne recalcule pas la matrice de pression a chaque pas de temps.
  int matrice_pression_invariante;
  // Si on veut ajouter une interface a vitesse imposee :
  //  reference a l'equation de transport correspondante :
  VECT(OBS_PTR(Transport_Interfaces_FT_Disc)) ref_eq_interf_vitesse_imposee;
  // Si le fluide est diphasique, c'est l'indicatrice de l'equation suivante
  // qui est utilisee pour determiner les proprietes du fluide:
  // (masse volumique, viscosite, tension superficielle, ...)
  OBS_PTR(Transport_Interfaces_FT_Disc) ref_eq_interf_proprietes_fluide;
  // Si le fluide est diphasique, la reference au fluide:
  OBS_PTR(Fluide_Diphasique) ref_fluide_diphasique;

  OBS_PTR(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_;
  OBS_PTR(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_vap_;

  // Valeur maximale de courbure autorisee pour calculer le
  // terme source de tension de surface (clipping si valeur superieur)
  double clipping_courbure_interface;

  int flag_correction_trainee_; // EB
  double alpha_correction_trainee_; // EB
  double beta_correction_trainee_; // EB
  int faces_diphasiques_; // EB
  int extension_reynolds_; // EB
  int proportionnel_; // EB

  enum Terme_Gravite
  {
    GRAVITE_RHO_G, GRAVITE_GRAD_I
  };
  Terme_Gravite terme_gravite_;
  Noms equations_concentration_source_fluide_;
  // Si is_explicite != 0,
  //   on calcul vpoint de facon explicite dans l etape de prediction des vitesses.
  int is_explicite;
  // Si is_boussinesq_ != 0, on calcul une force par le modele de boussinesq
  int is_boussinesq_;
  // Flag pour la method du saut de vitesse :
  int new_mass_source_;

  enum Type_interpol_indic_pour_dI_dt
  {
    INTERP_STANDARD, INTERP_MODIFIEE, INTERP_AI_BASED, INTERP_STANDARD_UVEXT, INTERP_MODIFIEE_UVEXT, INTERP_AI_BASED_UVEXT, INTERP_STANDARD_UIEXT, INTERP_MODIFIEE_UIEXT, INTERP_AI_BASED_UIEXT
  };
  Type_interpol_indic_pour_dI_dt type_interpol_indic_pour_dI_dt_;

  enum OutletCorrection_pour_dI_dt
  {
    NO_CORRECTION, CORRECTION_GHOST_INDIC, ZERO_NET_FLUX_ON_MIXED_CELLS, ZERO_OUT_FLUX_ON_MIXED_CELLS
  };
  OutletCorrection_pour_dI_dt OutletCorrection_pour_dI_dt_;

  // Si is_penalized != 0,
  //   on penalise L2 le terme de forcage.
  int is_penalized;
  // Valeur de l'inverse du coefficient de penalisation L2 du terme de forcage.
  double eta;
  // Valeur pour la penalisation L2 de la pression.
  double p_ref_pena;
  // Point de penalisation L2 de la pression du fluide
  int is_pfl_flottant; // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
  double x_pfl_imp;
  double y_pfl_imp;
  double z_pfl_imp;

  DoubleTab force_pression_interf_; // EB pression locale pour chaque fa7
  DoubleTab force_frottements_interf_; // EB force de frottement locale pour chaque fa7
  DoubleTab pression_interf_; // EB pression a l'interface locale pour chaque fa7
  DoubleVect surface_tot_interf_; // EB
  DoubleTab force_pression_tot_interf_; // EB
  DoubleTab force_frottements_tot_interf_; // EB
  DoubleTab force_pression_tot_interf_stokes_th_; // EB
  DoubleTab force_frottements_tot_interf_stokes_th_; // EB
  DoubleTab force_pression_tot_interf_stokes_th_dis_; // EB
  DoubleTab force_frottements_tot_interf_stokes_th_dis_; // EB


  DoubleTab sigma_xx_interf_, sigma_xy_interf_, sigma_xz_interf_, sigma_yx_interf_, sigma_yy_interf_, sigma_yz_interf_, sigma_zx_interf_, sigma_zy_interf_, sigma_zz_interf_; // EB tenseur des contraintes local pour c
  DoubleTab sigma_xx_interf_stokes_th_dis_, sigma_xy_interf_stokes_th_dis_, sigma_xz_interf_stokes_th_dis_, sigma_yx_interf_stokes_th_dis_, sigma_yy_interf_stokes_th_dis_, sigma_yz_interf_stokes_th_dis_, sigma_zx_interf_stokes_th_dis_, sigma_zy_interf_stokes_th_dis_, sigma_zz_interf_stokes_th_dis_; // EB tenseur des contraintes local pour c

  DoubleTab sigma_xx_interf_stokes_th_, sigma_xy_interf_stokes_th_, sigma_xz_interf_stokes_th_, sigma_yy_interf_stokes_th_, sigma_yz_interf_stokes_th_, sigma_zz_interf_stokes_th_;
  DoubleTab dUdx_P1_, dUdy_P1_, dUdz_P1_, dVdx_P1_, dVdy_P1_, dVdz_P1_, dWdx_P1_, dWdy_P1_, dWdz_P1_;
  DoubleTab dUdx_P2_, dUdy_P2_, dUdz_P2_, dVdx_P2_, dVdy_P2_, dVdz_P2_, dWdx_P2_, dWdy_P2_, dWdz_P2_;
  DoubleTab dUdx_P1_th_, dUdy_P1_th_, dUdz_P1_th_, dVdx_P1_th_, dVdy_P1_th_, dVdz_P1_th_, dWdx_P1_th_, dWdy_P1_th_, dWdz_P1_th_;
  DoubleTab dUdx_P2_th_, dUdy_P2_th_, dUdz_P2_th_, dVdx_P2_th_, dVdy_P2_th_, dVdz_P2_th_, dWdx_P2_th_, dWdy_P2_th_, dWdz_P2_th_;
  DoubleTab dUdx_P1_th_dis_, dUdy_P1_th_dis_, dUdz_P1_th_dis_, dVdx_P1_th_dis_, dVdy_P1_th_dis_, dVdz_P1_th_dis_, dWdx_P1_th_dis_, dWdy_P1_th_dis_, dWdz_P1_th_dis_;
  DoubleTab dUdx_P2_th_dis_, dUdy_P2_th_dis_, dUdz_P2_th_dis_, dVdx_P2_th_dis_, dVdy_P2_th_dis_, dVdz_P2_th_dis_, dWdx_P2_th_dis_, dWdy_P2_th_dis_, dWdz_P2_th_dis_;

  DoubleTab force_pression_stokes_th_, force_pression_stokes_th_dis_;
  DoubleTab force_frottements_stokes_th_, force_frottements_stokes_th_dis_;
  DoubleTab pression_interf_stokes_th_dis_;

  DoubleTab U_P1_, U_P2_, U_P1_th_, U_P2_th_, U_P1_th_dis_, U_P2_th_dis_;
  DoubleTab U_P2_moy_; // vitesse moyenne en P2 pour les points de calcul le permettant (fluide ou solide)
  DoubleTab Indic_elem_P2_, Prop_P2_fluide_compo_; // Indic de l'element dans lequel se trouve P2, proportion de points P2 fluide par compo
  DoubleTab Proportion_fa7_ok_UP2_;  // Nb_fa7_ok_prop_ : pourcentage de fa7 pour lesquelles on a pu calculer la vitesse moyenne
  IntTab list_elem_P1_, list_elem_diph_, list_elem_P1_all_; // EB list_elem_P1_ : liste des elements dans lesquels se trouvent les points P1, list_elem_diph_ : liste des elements traverses par l'interface,
  // list_elem_P1_all_ : liste de tous les elements - PUREMENT FLUIDE UNIQUEMENT - ayant servis a l'interpolation des champs en P1
};

#endif /* Navier_Stokes_FT_Disc_interne_included */
