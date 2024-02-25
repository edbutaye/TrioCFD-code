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
// File:        Navier_Stokes_FT_Disc.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/55
//
//////////////////////////////////////////////////////////////////////////////
#include <Parametre_implicite.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Probleme_FT_Disc_gen.h>
#include <Modele_turbulence_hyd_null.h>
#include <Discret_Thyd.h>
#include <Operateur_Diff_base.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Assembleur_base.h>
#include <TRUST_Vector.h>
#include <Schema_Temps_base.h>
#include <Domaine_VDF.h>
#include <Debog.h>
#include <Op_Dift_VDF_Face_leaves.h>
// Ces includes seront a retirer quand on aura clairement separe les operations
// specifiques au VDF et VEF
#include <Domaine_VF.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Terme_Source_Constituant_Vortex_VEF_Face.h>
#include <TRUSTTrav.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice_Bloc.h>
#include <Param.h>
#include <Vecteur3.h>
#include <Domaine_Cl_VDF.h>
#if TCL_MODEL
//#include <Dirichlet_paroi_fixe.h>
//#include <Dirichlet_paroi_defilante.h>
//#include <Neumann_sortie_libre.h>
//#include <Periodique.h>
//#include <Symetrie.h>
//#include <Dirichlet.h>
//#include <Dirichlet_homogene.h>
#endif
#define NS_VERBOSE 0 // To activate verbose mode on err ...
// EB
#include <Modele_Collision_FT.h>
#include <Particule_Solide.h>
#include <Matrice_Dense.h>
#include <Matrice_Diagonale.h>
#include <Maillage_FT_Disc.h>
#include <Postraitement_Forces_Interfaces_FT.h>
#include <MD_Vector_tools.h>
#include <EcritureLectureSpecial.h>
#include <Avanc.h>
#include <iostream>
#include <vector>
#include <Statistiques.h>
#include <Connex_components_FT.h>
#include <Connex_components.h>

#include <fstream>
#include <iomanip>
#include <SFichier.h>
#include <EFichier.h>
#include <sys/stat.h>
// fin EB
// #include <chrono>
// using namespace std::chrono;

Implemente_instanciable_sans_constructeur_ni_destructeur(Navier_Stokes_FT_Disc,"Navier_Stokes_FT_Disc",Navier_Stokes_Turbulent);

class Navier_Stokes_FT_Disc_interne
{
public:
  Navier_Stokes_FT_Disc_interne() :
    correction_courbure_ordre_(0), // Par defaut, pas de prise en compte de la courbure pour corriger le champ etendu delta_vitesse
    mpoint_inactif(0),   // Par defaut, mpoint cree un saut de vitesse
    mpointv_inactif(0),   // Par defaut, mpointv  cree un saut de vitesse
    matrice_pression_invariante(0),   // Par defaut, recalculer la matrice pression
    clipping_courbure_interface(1e40),// Par defaut, pas de clipping
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
    y_pfl_imp(-1.e40),
    z_pfl_imp(-1.e40)
  {};
  int correction_courbure_ordre_;
  int mpoint_inactif;
  int mpointv_inactif;

  Champ_Fonc second_membre_projection;
  Champ_Fonc second_membre_projection_jump_;
  Champ_Fonc derivee_u_etoile;
  Champ_Fonc gradient_pression;
  Champ_Fonc terme_diffusion;
  Champ_Fonc terme_convection;
  Champ_Fonc terme_source;
  Champ_Fonc terme_source_interfaces;
  Champ_Fonc indicatrice_p1b;
  Champ_Fonc gradient_indicatrice;
  Champ_Fonc potentiel_faces;
  Champ_Fonc potentiel_elements;
  // delta_u_interface = la partie "saut de vitesse" du champ de vitesse a l'interface
  Champ_Inc delta_u_interface;
  Champ_Fonc laplacien_d;
  Champ_Fonc mpoint;
  Champ_Fonc mpoint_vap;
  // Variation temporelle indicatrice de phase
  Champ_Fonc derivee_temporelle_indicatrice;
  Champ_Fonc ai; // Eulerian interfacial area.
  Champ_Inc vitesse_jump0_; // Extended Velocity of phase 0.
  Champ_Fonc terme_source_collisions; // HMS
  Champ_Fonc  num_compo; // HMS
  Champ_Fonc vitesse_stokes_th_; // EB
  Champ_Fonc pression_stokes_th_; // EB
  Champ_Fonc terme_correction_trainee; // EB
  LIST(REF(Champ_base)) liste_champs_compris;

  // Si matrice_pression_invariante != 0,
  //   on ne recalcule pas la matrice de pression a chaque pas de temps.
  int matrice_pression_invariante;
  // Si on veut ajouter une interface a vitesse imposee :
  //  reference a l'equation de transport correspondante :
  VECT(REF(Transport_Interfaces_FT_Disc)) ref_eq_interf_vitesse_imposee;
  // Si le fluide est diphasique, c'est l'indicatrice de l'equation suivante
  // qui est utilisee pour determiner les proprietes du fluide:
  // (masse volumique, viscosite, tension superficielle, ...)
  REF(Transport_Interfaces_FT_Disc) ref_eq_interf_proprietes_fluide;
  // Si le fluide est diphasique, la reference au fluide:
  REF(Fluide_Diphasique) ref_fluide_diphasique;

  REF(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_;
  REF(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_vap_;

  // Valeur maximale de courbure autorisee pour calculer le
  // terme source de tension de surface (clipping si valeur superieur)
  double clipping_courbure_interface;

  int flag_correction_trainee_; // EB
  double alpha_correction_trainee_; // EB
  double beta_correction_trainee_; // EB
  int faces_diphasiques_; // EB
  int extension_reynolds_; // EB
  int proportionnel_; // EB

  //-----------
  //enum Modele_collisions { RESSORT_AMORTI_VVA, RESSORT_AMORTI_ESI, RESSORT_AMORTI_EE, MOHAGHEG,HYBRID_EE, HYBRID_ESI, HYBRID_RK3, HYBRID_VVA, BREUGEM };
  //Modele_collisions modele_collisions;
  //-----------
  enum Terme_Gravite { GRAVITE_RHO_G, GRAVITE_GRAD_I };
  Terme_Gravite terme_gravite_;
  Noms equations_concentration_source_fluide_;
  // Si is_explicite != 0,
  //   on calcul vpoint de facon explicite dans l etape de prediction des vitesses.
  int is_explicite;
  // Si is_boussinesq_ != 0, on calcul une force par le modele de boussinesq
  int is_boussinesq_;
  // Flag pour la method du saut de vitesse :
  int new_mass_source_;

  enum Type_interpol_indic_pour_dI_dt { INTERP_STANDARD, INTERP_MODIFIEE, INTERP_AI_BASED,
                                        INTERP_STANDARD_UVEXT, INTERP_MODIFIEE_UVEXT, INTERP_AI_BASED_UVEXT,
                                        INTERP_STANDARD_UIEXT, INTERP_MODIFIEE_UIEXT, INTERP_AI_BASED_UIEXT
                                      };
  Type_interpol_indic_pour_dI_dt type_interpol_indic_pour_dI_dt_;

  enum OutletCorrection_pour_dI_dt { NO_CORRECTION, CORRECTION_GHOST_INDIC, ZERO_NET_FLUX_ON_MIXED_CELLS, ZERO_OUT_FLUX_ON_MIXED_CELLS };
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

/*! @brief Calcul de champ_rho_faces_ et champ_rho_elem_ en fonction de l'indicatrice: rho_elem_ = indicatrice * ( rho(phase_1) - rho(phase_0) ) + rho(phase_0)
 *
 *    rho_faces_ = 0.5 * (rho_elem(voisin_0) + rho_elem(voisin_1))
 *   Le calcul des viscosites cinematique et dynamique est pour l'instant le suivant:
 *    nu_elem = indicatrice * ( nu(phase_1) - nu(phase_0) ) + nu(phase_0)
 *    mu_elem = nu_elem * rho_elem
 *
 */
static void FT_disc_calculer_champs_rho_mu_nu_dipha(const Domaine_dis_base&      domaine_dis_base,
                                                    const Fluide_Diphasique& fluide,
                                                    const DoubleVect&         indicatrice_elem,
                                                    Transport_Interfaces_FT_Disc&         eq_transport,
                                                    DoubleVect& rho_elem,
                                                    DoubleVect& nu_elem,
                                                    DoubleVect& mu_elem,
                                                    DoubleVect& rho_faces, const double temps_courant, const int calcul_precis)
{
  const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
  const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
  const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
  const double rho_phase_0 = tab_rho_phase_0(0,0);
  const double rho_phase_1 = tab_rho_phase_1(0,0);
  const double delta_rho = rho_phase_1 - rho_phase_0;
  const DoubleTab& tab_nu_phase_0 = phase_0.viscosite_cinematique().valeur().valeurs();
  const DoubleTab& tab_nu_phase_1 = phase_1.viscosite_cinematique().valeur().valeurs();
  const double nu_phase_0 = tab_nu_phase_0(0,0);
  const double nu_phase_1 = tab_nu_phase_1(0,0);
  const double delta_nu = nu_phase_1 - nu_phase_0;
  const double mu_phase_0 = nu_phase_0 * rho_phase_0;
  const double mu_phase_1 = nu_phase_1 * rho_phase_1;
  const double delta_mu = mu_phase_1 - mu_phase_0;
  double mu = 0.;
  const int formule_mu = fluide.formule_mu();

  // Calcul de rho, nu, mu aux elements
  {
    const int n = indicatrice_elem.size();
    assert(n == rho_elem.size());
    assert(n == nu_elem.size());
    assert(n == mu_elem.size());
    for (int i = 0; i < n; i++)
      {
        const double indic = indicatrice_elem[i];
        const double rho = indic * delta_rho + rho_phase_0;
        rho_elem[i] = rho;
        double nu  = indic * delta_nu  + nu_phase_0;
        switch(formule_mu)
          {
          case 0: // standard default method (will be removed)
            {
              mu  = nu*rho ;
            }
            break;
          case 1: // Arithmetic average
            {
              mu  = indic * delta_mu  + mu_phase_0;
            }
            break;
          case 2: // Harmonic average
            {
              mu  = (mu_phase_0 * mu_phase_1)/(mu_phase_1 - indic * delta_mu);
            }
            break;
          case 3: // Staircase Average
            {
              mu =  indic==0 ? mu_phase_0 : mu_phase_1;
            }
            break;
          default:
            {
              Cerr << "The method specified for formule_mu in not recognized. \n" << finl;
              Cerr << "you can choose : standard, arithmetic or harmonic. \n" << finl;
              //Cerr << "We should not be here Navier_Stokes_FT_Disc::FT_disc_calculer_champs_rho_mu_nu_dipha" << finl;
              Process::exit();
            }
          }
        mu_elem[i]  = mu;
        nu  = mu / rho;
        nu_elem[i]  = nu;
      }
    rho_elem.echange_espace_virtuel();
    mu_elem.echange_espace_virtuel();
    nu_elem.echange_espace_virtuel();
  }

// Calcul de rho aux faces (on suppose que la vitesse est aux faces)
  if (calcul_precis && temps_courant>0)
    {
      const DoubleTab& indicatrice_face_elem=eq_transport.get_compute_indicatrice_faces().valeurs();
      assert(rho_elem.size() == domaine_dis_base.nb_elem());
      const IntTab& face_voisins = domaine_dis_base.face_voisins();
      const int nfaces = face_voisins.dimension(0);
      assert(rho_faces.size() == nfaces);
      for (int i = 0; i < nfaces; i++)
        {
          const double indic =  indicatrice_face_elem[i];
          rho_faces[i] = indic * delta_rho  + rho_phase_0;
        }
    }
  else
    {
      assert(rho_elem.size() == domaine_dis_base.nb_elem());
      const IntTab& face_voisins = domaine_dis_base.face_voisins();
      const int nfaces = face_voisins.dimension(0);
      assert(rho_faces.size() == nfaces);
      for (int i = 0; i < nfaces; i++)
        {
          const int elem0 = face_voisins(i, 0);
          const int elem1 = face_voisins(i, 1);
          double rho = 0.;
          if (elem0 >= 0)
            rho = rho_elem[elem0];
          if (elem1 >= 0)
            rho += rho_elem[elem1];
          if (elem0 >= 0 && elem1 >= 0)
            rho *= 0.5;
          rho_faces[i] = rho;
        }
    }
  rho_faces.echange_espace_virtuel();
}

static void FT_disc_calculer_champs_rho_mu_nu_mono(const Domaine_dis_base& zdis,
                                                   const Fluide_Incompressible& fluide,
                                                   Champ_Fonc& champ_rho_elem_,
                                                   Champ_Don& champ_mu_,
                                                   Champ_Don& champ_nu_,
                                                   Champ_Fonc& champ_rho_faces_)
{

  if (sub_type(Champ_Uniforme,champ_rho_elem_.valeur()) && (sub_type(Champ_Uniforme,champ_nu_.valeur())))
    {
      const DoubleTab& tab_rho_phase_0 = fluide.masse_volumique().valeurs();
      const double rho = tab_rho_phase_0(0,0);
      const DoubleTab& tab_nu_phase_0 = fluide.viscosite_cinematique().valeur().valeurs();
      const double nu = tab_nu_phase_0(0,0);
      const double mu = nu * rho;
      champ_rho_elem_. valeur().valeurs() = rho;
      champ_nu_.       valeur().valeurs() = nu;
      champ_mu_.       valeur().valeurs() = mu;
      champ_rho_faces_.valeur().valeurs() = rho;
    }
  else
    {
      const int nb_el = zdis.domaine().nb_elem();
      const Domaine_VF& zvf = ref_cast(Domaine_VF,zdis);
      const IntTab& face_vois = zvf.face_voisins();
      const DoubleVect& vol = zvf.volumes();
      const int nb_faces = zvf.nb_faces();
      int elem1,elem2;
      double volume;
      //


      IntVect les_polys(nb_el);
      for (int i=0; i<nb_el; i++)
        les_polys(i) = i;

      const DoubleTab& cg = zvf.xp();
      DoubleTab& val_rho = champ_rho_elem_.valeur().valeurs();
      DoubleTab& val_nu = champ_nu_.valeur().valeurs();
      DoubleTab& val_mu = champ_mu_.valeur().valeurs();
      DoubleTab& val_rho_faces = champ_rho_faces_.valeur().valeurs();

      fluide.masse_volumique()->valeur_aux_elems(cg,les_polys,val_rho);
      fluide.viscosite_dynamique()->valeur_aux_elems(cg,les_polys,val_mu);
      val_rho.echange_espace_virtuel();
      val_mu.echange_espace_virtuel();

      for (int el=0; el<nb_el; el++)
        val_nu(el) = val_mu(el)/val_rho(el);
      val_nu.echange_espace_virtuel();

      for (int fac=0; fac<nb_faces; fac++)
        {
          elem1 = face_vois(fac,0);
          elem2 = face_vois(fac,1);
          val_rho_faces(fac) = 0.;
          volume = 0.;
          if (elem1!=-1)
            {
              val_rho_faces(fac) += val_rho(elem1)*vol(elem1);
              volume += vol(elem1);
            }
          if (elem2!=-1)
            {
              val_rho_faces(fac) += val_rho(elem2)*vol(elem2);
              volume += vol(elem2);
            }
          val_rho_faces(fac) /=volume;
        }

      val_rho_faces.echange_espace_virtuel();
    }
}
/*! @brief constructeur par defaut
 *
 */
Navier_Stokes_FT_Disc::Navier_Stokes_FT_Disc()
{
  variables_internes_ = new Navier_Stokes_FT_Disc_interne;
  // terme_diffusif.associer_diffusivite(champ_mu_);
  is_repulsion=0;
}

/*! @brief le destructeur qui va avec
 *
 */
Navier_Stokes_FT_Disc::~Navier_Stokes_FT_Disc()
{
  delete variables_internes_;
}

Sortie& Navier_Stokes_FT_Disc::printOn(Sortie& os) const
{
  return Navier_Stokes_std::printOn(os);
}

Entree& Navier_Stokes_FT_Disc::readOn(Entree& is)
{
  terme_diffusif.associer_diffusivite(champ_mu_);
  return Navier_Stokes_Turbulent::readOn(is);
}

void Navier_Stokes_FT_Disc::set_param(Param& param)
{
  Navier_Stokes_Turbulent::set_param(param);
  param.ajouter_flag("boussinesq_approximation",&variables_internes().is_boussinesq_);
  param.ajouter_flag("new_mass_source",&variables_internes().new_mass_source_);
  param.ajouter_flag("matrice_pression_invariante",&variables_internes().matrice_pression_invariante);
  param.ajouter_non_std("equation_interfaces_vitesse_imposee",(this) );
  param.ajouter_non_std("equations_interfaces_vitesse_imposee",(this) );
  param.ajouter_non_std("equation_interfaces_proprietes_fluide",(this));
  param.ajouter("clipping_courbure_interface",&variables_internes().clipping_courbure_interface);
  param.ajouter_non_std("equation_temperature_mpoint",(this));
  param.ajouter_non_std("equation_temperature_mpoint_vapeur",(this));
  param.ajouter_non_std("terme_gravite",(this));
  param.ajouter("equations_concentration_source_vortex",&variables_internes().equations_concentration_source_fluide_);
  param.ajouter_non_std("repulsion_aux_bords",(this));
  param.ajouter_non_std("penalisation_forcage",(this));
  param.ajouter_flag("mpoint_inactif_sur_qdm",&variables_internes().mpoint_inactif);
  param.ajouter_flag("mpoint_vapeur_inactif_sur_qdm",&variables_internes().mpointv_inactif);
  param.ajouter("correction_courbure_ordre", &variables_internes().correction_courbure_ordre_);
  param.ajouter_non_std("interpol_indic_pour_dI_dt", (this));
  param.ajouter_non_std("OutletCorrection_pour_dI_dt", (this));
  param.ajouter_non_std("correction_trainee",(this)); // EB
}

int Navier_Stokes_FT_Disc::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      // Si on a lu le modele de turbulence et qu'il est nul,
      // alors on utilise l'operateur de diffusion standard.
      if (le_modele_turbulence.non_nul() // L'operateur a ete type (donc lu)
          && sub_type(Modele_turbulence_hyd_null, le_modele_turbulence.valeur()))
        {
          is >> terme_diffusif; // Operateur de diffusion standard (non turbulent)
          if (Process::je_suis_maitre())
            {
              Cerr << " WARNING : standard diffusion operator used for front_tracking\n";
              Cerr << "           the transposed term grad(v) is missing !!!" << finl;
            }
        }
      else
        {
          // debut EB
          const Nom dis=discretisation().que_suis_je();

          if (dis=="VDF+")
            {
              Nom type="Op_Dift_VDF_Face_FT";
              terme_diffusif.typer(type);
              terme_diffusif.l_op_base().associer_eqn(*this);
              Cerr << terme_diffusif.valeur().que_suis_je() << finl;
              terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
              Motcle motbidon;
              is >>  motbidon; // on passe les 2 accolades
              is >>  motbidon;
              // fin EB
            }
          else
            lire_op_diff_turbulent(is);
        }
      // Le coefficient de diffusion est une viscosite dynamique.
      // Il faut le diviser par rho pour calculer le pas de temps de stabilite.
      terme_diffusif.valeur().
      associer_champ_masse_volumique(champ_rho_elem_.valeur());

      // Pareil avec le modele de turbulence
      if (le_modele_turbulence.non_nul())
        le_modele_turbulence.valeur().
        associer_champ_masse_volumique(champ_rho_elem_.valeur());
      return 1;
    }
  else if ((mot=="equation_interfaces_vitesse_imposee") || (mot=="equation_interfaces_proprietes_fluide"))
    {
      Motcle nom_equation;
      is >> nom_equation;
      Cerr << mot << " equation : " << nom_equation << finl;
      const Transport_Interfaces_FT_Disc& eq =
        probleme_ft().equation_interfaces(nom_equation);

      if (mot=="equation_interfaces_vitesse_imposee")
        {
          if (variables_internes().ref_eq_interf_vitesse_imposee.size()>0)
            {
              Cerr<<"Error: You have already a equation_interfaces_vitesse_imposee defined."<<finl;
              Cerr<<"Since 1.6.4 TRUST version, you need to use the following syntax when"<<finl;
              Cerr<<"using several equation_interfaces_vitesse_imposee equations:"<<finl;
              Cerr<<"equations_interfaces_vitesse_imposee number_of_equations equation_name_one equation_name_two ..."<<finl;
              exit();
            }
          else
            {
              Cerr << "===============================================================================" << finl;
              Cerr << "Warning: You are using a future obsolete syntax for defining a solid interface:" << finl;
              Cerr << "equation_interfaces_vitesse_imposee " << nom_equation << finl;
              Cerr << "Should be written:" << finl;
              Cerr << "equations_interfaces_vitesse_imposee 1 " << nom_equation << finl;
              Cerr << "===============================================================================" << finl;
            }
          variables_internes().ref_eq_interf_vitesse_imposee.add(eq);
        }
      else if (mot=="equation_interfaces_proprietes_fluide")
        variables_internes().ref_eq_interf_proprietes_fluide = eq;
      return 1;
    }
  else if (mot=="equations_interfaces_vitesse_imposee")
    {
      Noms na;
      is >> na;
      variables_internes().ref_eq_interf_vitesse_imposee.reset();
      for (int i=0; i<na.size(); i++)
        {
          const Transport_Interfaces_FT_Disc& eq =
            probleme_ft().equation_interfaces(na[i]);
          variables_internes().ref_eq_interf_vitesse_imposee.add(eq);
        }
      return 1;
    }
  else if (mot=="equation_temperature_mpoint")
    {
      Nom nom_equation;
      is >> nom_equation;
      Cerr << " equation : " << nom_equation << finl;
      const Equation_base& eq =
        probleme_ft().get_equation_by_name(nom_equation);
      if (!sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
        {
          Cerr << " Error : equation is not of type Convection_Diffusion_Temperature_FT_Disc"
               << finl;
          exit();
        }
      variables_internes().ref_equation_mpoint_ =
        ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
      return 1;
    }
  else if (mot=="equation_temperature_mpoint_vapeur")
    {
      Nom nom_equation;
      is >> nom_equation;
      Cerr << " equation : " << nom_equation << finl;
      const Equation_base& eq =
        probleme_ft().get_equation_by_name(nom_equation);
      if (!sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
        {
          Cerr << " Error : equation is not of type Convection_Diffusion_Temperature_FT_Disc"
               << finl;
          Process::exit();
        }
      variables_internes().ref_equation_mpoint_vap_ =
        ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
      return 1;
    }
  else if (mot=="terme_gravite")
    {
      Motcles mots;
      mots.add("rho_g"); // Terme source volumique traditionnel avec courants parasites
      mots.add("grad_i"); // Terme source "front-tracking" sans courants parasites mais pbs en VEF.
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading terme_gravite : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          variables_internes().terme_gravite_ = Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G;
          break;
        case 1:
          variables_internes().terme_gravite_ = Navier_Stokes_FT_Disc_interne::GRAVITE_GRAD_I;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot=="repulsion_aux_bords")
    {
      is_repulsion=1;
      is >> minx;
      is >> maxx;
      is >> pente;
      Cerr <<"Interfaces repulsion on the boundaries for : "<<minx<<" "<<maxx<<finl;
      return 1;
    }
  else if (mot=="penalisation_forcage")
    {
      variables_internes().is_penalized = 1;
      variables_internes().is_explicite = 0;
      variables_internes().eta = 1.e-12;
      Cerr << "Navier_Stokes_FT_Disc : option penalisation_forcage" << finl;
      Motcles mots;
      mots.add("pression_reference"); // Valeur reference pour penalisation L2 pression
      mots.add("domaine_flottant_fluide"); // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
      Motcle motbis;
      is >> motbis;
      Motcle accouverte = "{" , accfermee = "}" ;
      if (motbis == accouverte)
        {
          is >> motbis;
          while (motbis != accfermee)
            {
              int rang = mots.search(motbis);
              switch(rang)
                {
                case 0:
                  {
                    is >> variables_internes().p_ref_pena;
                    Cerr << "Lecture de la pression de reference : " << variables_internes().p_ref_pena << finl;
                    break;
                  }
                case 1:
                  {
                    variables_internes().is_pfl_flottant = 1;
                    is >> variables_internes().x_pfl_imp ;
                    is >> variables_internes().y_pfl_imp ;
                    if( Objet_U::dimension == 3 )
                      is >> variables_internes().z_pfl_imp ;
                    Cerr <<"Domaine flottant fluide"<< finl;
                    Cerr <<"Lecture de la position du point de reference pression fluide : "<<variables_internes().x_pfl_imp <<" "<<variables_internes().y_pfl_imp<<" "<<variables_internes().z_pfl_imp<< finl;
                    break;
                  }
                default:
                  Cerr << "Erreur, on attendait " << mots << "On a trouve : " << motbis << finl;
                  barrier();
                  exit();
                }
              is >> motbis;
            }
        }
      else
        {
          Cerr << "Erreur a la lecture des parametres de la penalisation du forcage " << finl;
          Cerr << "On attendait : " << accouverte << finl;
          Process::exit();
        }
    }
  else if (mot=="correction_trainee")
    {
      Cerr << "Lecture des parametres de la correction de la trainee : ";
      Motcles mots;
      mots.add("alpha"); // constante correlation
      mots.add("beta"); // puissance correlation
      mots.add("faces_diphasiques"); // 1 : correction appliquee sur les faces d'indicatrice strictement inferieure a 1,
      //0: correction appliquee sur les faces d'indicatrice nulle -> initialisee a 1
      mots.add("extension_reynolds"); // 1 : on multiplie la correction par la correlation d'Abraham, 0 : on n'en tient pas compte -> intialise a 0
      mots.add("proportionnel"); // 1 : la correction est definie comme proportionnelle a la force de trainee calculee, 0 : la correction est egale
      // a la force de stokes a un facteur multiplicatif pres. Dans les deux cas, les expressions sont obtenues par une fonction d'interpolation python
      // -> initialisee a 0
      Motcle motbis;
      Motcle accouverte = "{" , accfermee = "}" ;
      is >> motbis;
      if (motbis == accouverte)
        {
          is >> variables_internes().flag_correction_trainee_;
          is >> motbis;
          while (motbis != accfermee)
            {
              int rang = mots.search(motbis);
              switch(rang)
                {
                case 0:
                  {
                    is >> variables_internes().alpha_correction_trainee_;
                    Cerr << "\talpha = " << variables_internes().alpha_correction_trainee_;
                    break;
                  }
                case 1:
                  {
                    is >> variables_internes().beta_correction_trainee_;
                    Cerr << "\tbeta = " << variables_internes().beta_correction_trainee_ << finl;
                    break;
                  }
                case 2:
                  {
                    is >> variables_internes().faces_diphasiques_;
                    if (variables_internes().faces_diphasiques_) Cerr << "La correction de la trainee sera appliquee sur les faces"
                                                                        "solides et diphasiques." << finl;
                    else Cerr << "La correction de la trainee sera appliquee sur les faces solides uniquement." << finl;
                    break;
                  }
                case 3:
                  {
                    is >> variables_internes().extension_reynolds_;
                    if (variables_internes().extension_reynolds_) Cerr << "On utilise la correlation d'Abraham pour etendre"
                                                                         "la correction de la force de trainee a plus haut reynolds." << finl;
                    break;
                  }
                case 4:
                  {
                    is >> variables_internes().proportionnel_;
                    if (variables_internes().proportionnel_) Cerr << "On utilise la forme F_PR-DNS = F_PR-SCS x g(D/delta_x)" << finl;
                    else  Cerr << "On utilise la forme F_PR-DNS = F_PR-SCS + F_stokes x (1-h(D/delta_x))." << finl;
                    break;
                  }
                default:
                  Cerr << "Erreur, on attendait " << mots << " On a trouve : " << motbis << finl;
                  barrier();
                  exit();
                }
              is >> motbis;
            }
        }
      else
        {
          Cerr << "Erreur, on attendait " << accouverte << "On a trouve : " << motbis << finl;
          barrier();
          exit();
        }
    }
  else if (mot =="interpol_indic_pour_dI_dt")
    {
      Motcles motcles2(9);
      motcles2[0] = "interp_standard";
      motcles2[1] = "interp_modifiee";
      motcles2[2] = "interp_ai_based";
      motcles2[3] = "interp_standard_uvext";
      motcles2[4] = "interp_modifiee_uvext";
      motcles2[5] = "interp_ai_based_uvext";
      motcles2[6] = "interp_standard_uiext";
      motcles2[7] = "interp_modifiee_uiext";
      motcles2[8] = "interp_ai_based_uiext";
      Motcle motlu;
      is >> motlu;
      Cerr << mot << " " << motlu << finl;
      Cout << "Setting the type of interpolation for calculer_dI_dt to " << motlu << "." << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_STANDARD;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the historical way"
                   << " where a mean value + upwind is used." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the field indicatrice_faces"
                   << " as defined by the interfacial transport option." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area"
                   << " and on the normal to the interface." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the historical way"
                   << " where a mean value + upwind is used."
                   << " Additionally, uv_ext is used." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the field indicatrice_faces"
                   << " as defined by the interfacial transport option."
                   << " Additionally, uv_ext is used."  << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area"
                   << " and on the normal to the interface."
                   << " Additionally, uv_ext is used."  << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the historical way"
                   << " where a mean value + upwind is used."
                   << " Additionally, ui_ext is used." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the field indicatrice_faces"
                   << " as defined by the interfacial transport option."
                   << " Additionally, ui_ext is used."  << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT:
          {
            variables_internes_->type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area"
                   << " and on the normal to the interface."
                   << " Additionally, ui_ext is used."  << finl;
            return 1;
          }

        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << "The options for methode_transport are :\n"
               << motcles2;
          Process::exit();
        }
    }
  else if (mot =="OutletCorrection_pour_dI_dt")
    {
      Motcles motcles2(4);
      motcles2[0] = "no_correction";
      motcles2[1] = "CORRECTION_GHOST_INDIC";
      motcles2[2] = "zero_net_flux_on_mixed_cells";
      motcles2[3] = "zero_out_flux_on_mixed_cells";
      Motcle motlu;
      is >> motlu;
      Cerr << mot << " " << motlu << finl;
      Cout << "Setting the type of correction at outlet BC for calculer_dI_dt to " << motlu << "." << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case Navier_Stokes_FT_Disc_interne::NO_CORRECTION:
          {
            variables_internes_->OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::NO_CORRECTION;
            if (Process::je_suis_maitre())
              Cerr << " No correction of div(chi u) at exit (historical way)" << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC:
          {
            variables_internes_->OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC;
            if (Process::je_suis_maitre())
              Cerr << " Correction of chi in ghost cells (virtually)" << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS:
          {
            variables_internes_->OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS;
            if (Process::je_suis_maitre())
              Cerr << " correction of div(chi u) at exit : zero divergence on cells touching outlet." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS:
          {
            variables_internes_->OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS;
            if (Process::je_suis_maitre())
              {
                Cerr << " correction of div(chi u) at exit : zero vapour mass flux on cells touching outlet." << finl;
                Cerr << " This is a bad option because it does not let vapour get out explicitly (prevents interface contact)" << finl;
                Cerr << " Should not be used or with great care." << finl;
                // Process::exit();
              }
            return 1;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << "The options for methode_transport are :\n"
               << motcles2;
          Process::exit();
        }
    }

  else
    return Navier_Stokes_Turbulent::lire_motcle_non_standard(mot,is);
  return 1;
}

// debut EB

// Description:
//    Sauvegarde num_compo
//    sur un flot de sortie.
// Precondition:
// Parametre:Sortie& os
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: int
//    Signification: renvoie toujours 1
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
int Navier_Stokes_FT_Disc::sauvegarder(Sortie& os) const
{
  int bytes=0;
  bytes += Navier_Stokes_Turbulent::sauvegarder(os);
  const REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  if (refeq_transport.non_nul())
    {
      const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
      if (eq_transport.is_solid_particle())
        {
          bytes+=variables_internes().num_compo.sauvegarder(os);
          assert(bytes % 4 == 0);
        }
    }
  return bytes;
}

int Navier_Stokes_FT_Disc::reprendre(Entree& is)
{
  Navier_Stokes_Turbulent::reprendre(is);
  const REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  if (eq_transport.is_solid_particle())
    {
      double temps = schema_temps().temps_courant();
      Nom ident_num_compo(variables_internes().num_compo.le_nom());
      ident_num_compo += variables_internes().num_compo.valeur().que_suis_je();
      ident_num_compo+= probleme().domaine().le_nom();
      ident_num_compo+=Nom(temps,probleme().reprise_format_temps());

      avancer_fichier(is,ident_num_compo);
      variables_internes().num_compo.reprendre(is);
    }
  return 1;
}
void ouvrir_fichier(SFichier& os,const Nom& type, const int& flag, const Navier_Stokes_FT_Disc& equation)
{
  // flag nul on n'ouvre pas le fichier
  if (flag==0)
    return ;
  int rang = -1;
  Nom fichier=Objet_U::nom_du_cas();
  if (type=="forces_particule")
    {
      fichier+="_Forces_Particule_";
      rang = 0;
    }
  else if( type=="forces_particule_th" )
    {
      fichier+= "_Forces_Particule_theoriques_";
      rang = 1;
    }
  else if (type=="forces_particule_lit")
    {
      fichier+= "_Forces_Particules_Lit_";
      rang = 2;
    }
  else if (type=="liste_collision")
    {
      fichier+= "_Liste_Collision_Lit_";
      rang = 3;
    }


  else
    {
      Cerr << "Le fichier " << type << " n est pas compris par Navier_Stokes_FT_Disc::ouvrir_fichier. "
           << "Cela semble du a une erreur d'implementation au sein de votre BALTIK." << finl;
    }

  fichier+=equation.le_nom();
  fichier+=".out";
  const Schema_Temps_base& sch=equation.probleme().schema_temps();
  const int& precision=sch.precision_impr();
  // On cree le fichier a la premiere impression avec l'en tete ou si le fichier n'existe pas

  struct stat f;
  //const int rang=fichier.search(type);
  if (stat(fichier,&f) && (sch.nb_impr_fpi()==1 && !equation.probleme().reprise_effectuee()))
    {
      os.ouvrir(fichier,ios::app);
      SFichier& file=os;
      Nom espace="\t";

      if (rang==0)
        {
          file << "###################################" << finl;
          file << "# Hydrodynamic force computation #"  << finl;
          file << "###################################" << finl;
          file << finl;
          file << "# Time [s]"<< espace << "Particle surface [m^2]" << espace << "Pressure force [N] (fpx fpy fpz)" << espace << "Friction force [N] (ffx ffy ffz)" << finl;
          file << finl;
        }
      if (rang==1)
        {
          file << "#####################################################################" << finl;
          file << "# Hydrodynamic force computation - Stokes theoretical configuration #"  << finl;
          file << "#####################################################################" << finl;
          file << "# Time [s]"<< finl;
          file << "# Stokes theoretical PRESSURE FORCE computed from the integration, on the lagrangian mesh, of the discretized analytical solution [N] (fpx_th fpy_th fpz_th)" << finl;
          file << "# Stokes PRESSURE FORCE computed with the developed method on the theoretical pressure field discretized on the eulerian mesh [N] (fpx_th_interp fpy_th_interp fpz_th_interp)" << finl;
          file << "# Stokes theoretical FRICTION FORCE computed from the integration, on the lagrangian mesh, of the discretized analytical solution [N] (ffx_th ffy_th ffz_th)" << finl;
          file << "# Stokes FRICTION FORCE computed with the developed method on the theoretical velocity field discretized on the eulerian mesh [N] (ffx_th_interp ffy_th_interp ffz_th_interp)" << finl;
          file << finl;
          file << "# Time" << espace << "fpx_th fpy_th fpz_th" << espace << "fpx_th_interp fpy_th_interp fpz_th_interp" << espace << "ffx_th ffy_th ffz_th" << espace << "ffx_th_interp ffy_th_interp ffz_th_interp" << finl;
          file << finl;
        }
      if (rang==2)
        {
          file << "#########################################################" << finl;
          file << "# Hydrodynamic force computation in a particle assembly #"  << finl;
          file << "#########################################################" << finl;
          file << finl;
          file << "# Time [s]" << espace << "num_compo" << "Pressure force [N] (fpx fpy fpz)" << espace << "Friction force [N] (ffx ffy ffz)" << espace << "Percentage of facets for which forces were computable" << espace << "Average fluid velocity in P2" << "Percentage of purely fluid cells in P2"<<  finl;
          file << finl;
        }
      if (rang==3)
        {
          file << "# List of colliding particle pairs" << finl;
          file << finl;
        }
    }
  else
    {
      os.ouvrir(fichier,ios::app);
    }
  os.precision(precision);
  os.setf(ios::scientific);
}

void Navier_Stokes_FT_Disc::imprimer(Sortie& os) const
{
  Navier_Stokes_Turbulent::imprimer(os);
}

// EB
/*! @brief imprile les forces de pression et de frottements, la surface totale de l'interface
 *  Pour le multi-particule (>5), imprime egalement le pourcentage des fa7, par particule, pour lesquels le calcul des efforts est possible,
 *  la vitesse moyenne aux points P1 et P2 (en ne prenant en compte que les facettes pour lesquelles l'interpolation est possible),
 *  ainsi que la proportion des elements P2 (elements contenant les points P2) qui sont purement fluide.
 */
int Navier_Stokes_FT_Disc::impr_fpi(Sortie& os) const
{
  const REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  if (refeq_transport.non_nul())
    {
      const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
      if (eq_transport.is_solid_particle())
        {
          const Postraitement_Forces_Interfaces_FT& les_post_interf=eq_transport.postraitement_forces_interf();

          if (les_post_interf.postraiter_forces())
            {
              if (Process::je_suis_maitre())
                {
                  const DoubleTab& force_pression_tot_interf=get_force_tot_pression_interf(); // EB
                  const DoubleTab& force_frottements_tot_interf=get_force_tot_frottements_interf(); // EB
                  const DoubleVect& surface_tot_interf=get_surface_tot_interf();

                  Nom espace= " ";
                  int dim_max_impr=5; // on imprime pas les valeurs si il y a plus de 5 particules dans le domaine

                  int nb_compo = force_pression_tot_interf.dimension(0); //
                  // On imprime les forces exercees par le fluide sur les particules
                  if (nb_compo<dim_max_impr)
                    {
                      SFichier Force_Particule;
                      const Navier_Stokes_FT_Disc& mon_eq = *this;
                      ouvrir_fichier(Force_Particule,"forces_particule",1,mon_eq);
                      schema_temps().imprimer_temps_courant(Force_Particule);

                      for (int compo=0; compo<nb_compo; compo++)
                        {
                          Force_Particule << espace << surface_tot_interf(compo) << espace;
                          for (int dim=0; dim<dimension; dim++) Force_Particule << espace << force_pression_tot_interf(compo,dim);
                          Force_Particule << espace;
                          for (int dim=0; dim<dimension; dim++) Force_Particule << espace << force_frottements_tot_interf(compo,dim);
                          Force_Particule << espace;
                        }

                      Force_Particule << finl;
                    }
                  else
                    {
                      SFichier Forces_Particule_Lit;
                      const Navier_Stokes_FT_Disc& mon_eq = *this;
                      ouvrir_fichier(Forces_Particule_Lit,"forces_particule_lit",1,mon_eq);
                      Forces_Particule_Lit << "TIME ";
                      schema_temps().imprimer_temps_courant(Forces_Particule_Lit);
                      Forces_Particule_Lit << finl;
                      const DoubleTab& U_P2_moy = variables_internes().U_P2_moy_;
                      const DoubleTab& Nb_fa7_ok_prop = variables_internes().Proportion_fa7_ok_UP2_;
                      const DoubleTab& Prop_indic_fluide_P2 = variables_internes().Prop_P2_fluide_compo_;
                      for (int compo = 0; compo < nb_compo; compo++)
                        {
                          int dim;
                          Forces_Particule_Lit << compo ;
                          Forces_Particule_Lit << espace << surface_tot_interf(compo) << espace;
                          for (dim=0; dim<dimension; dim++) Forces_Particule_Lit << espace << force_pression_tot_interf(compo,dim);
                          Forces_Particule_Lit << espace;
                          for (dim=0; dim<dimension; dim++) Forces_Particule_Lit << espace << force_frottements_tot_interf(compo,dim);
                          Forces_Particule_Lit << espace;
                          for (dim=0; dim<dimension; dim++) Forces_Particule_Lit << espace << Nb_fa7_ok_prop(compo,dim);
                          for (dim=0; dim<dimension; dim++) Forces_Particule_Lit << espace << U_P2_moy(compo,dim);
                          Forces_Particule_Lit << espace << Prop_indic_fluide_P2(compo);
                          Forces_Particule_Lit << finl;
                        }
                    }

                  // On imprime les forces exercees par le fluide sur les particules - donnees theoriques issues de la resoluton de
                  // l'ecoulement de Stokes
                  if (les_post_interf.calcul_forces_theoriques_stokes_ && schema_temps().nb_pas_dt()==1)
                    {
                      const DoubleTab& force_pression_tot_interf_stokes_th=get_force_pression_tot_interf_stokes_th(); // EB
                      const DoubleTab& force_frottements_tot_interf_stokes_th=get_force_frottements_tot_interf_stokes_th(); // EB
                      const DoubleTab& force_pression_tot_interf_stokes_th_dis=get_force_pression_tot_interf_stokes_th_dis(); // EB
                      const DoubleTab& force_frottements_tot_interf_stokes_th_dis=get_force_frottements_tot_interf_stokes_th_dis(); // EB

                      if (nb_compo<dim_max_impr)
                        {
                          SFichier Force_Particule_th;
                          const Navier_Stokes_FT_Disc& mon_eq = *this;
                          ouvrir_fichier(Force_Particule_th,"forces_particule_th",1,mon_eq);
                          schema_temps().imprimer_temps_courant(Force_Particule_th);
                          for (int compo=0; compo<nb_compo; compo++)
                            {
                              Force_Particule_th << espace;
                              for (int dim=0; dim<dimension; dim++) Force_Particule_th << espace << force_pression_tot_interf_stokes_th(compo,dim);
                              Force_Particule_th << espace;
                              for (int dim=0; dim<dimension; dim++) Force_Particule_th << espace << force_pression_tot_interf_stokes_th_dis(compo,dim);
                              Force_Particule_th << espace;
                              for (int dim=0; dim<dimension; dim++) Force_Particule_th << espace << force_frottements_tot_interf_stokes_th(compo,dim);
                              Force_Particule_th << espace;
                              for (int dim=0; dim<dimension; dim++) Force_Particule_th << espace << force_frottements_tot_interf_stokes_th_dis(compo,dim);
                              Force_Particule_th << finl;
                            }
                        }
                    }
                }
            }
        }
    }
  return 1;
}


// fin EB


const Champ_Don& Navier_Stokes_FT_Disc::diffusivite_pour_transport() const
{
  return champ_mu_;
}

void Navier_Stokes_FT_Disc::associer_pb_base(const Probleme_base& un_probleme)
{
  if (! sub_type(Probleme_FT_Disc_gen, un_probleme))
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::associer_pb_base\n";
      Cerr << " Navier_Stokes_FT_Disc equation must be associated to\n";
      Cerr << " a Probleme_FT_Disc_gen problem type\n";
      exit();
    }
  probleme_ft_ = ref_cast(Probleme_FT_Disc_gen, un_probleme);
  Navier_Stokes_std::associer_pb_base(un_probleme);
}

const Fluide_Diphasique& Navier_Stokes_FT_Disc::fluide_diphasique() const
{
  const REF(Fluide_Diphasique) & fluide_dipha = variables_internes().ref_fluide_diphasique;
  if (! fluide_dipha.non_nul())
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::fluide_diphasique()\n";
      Cerr << " The fluid has not been associated\n";
      Cerr << " (check that the fluid is of Fluide_Diphasique type)" << finl;
      assert(0);
      exit();
    }
  return fluide_dipha.valeur();
}

void Navier_Stokes_FT_Disc::associer_milieu_base(const Milieu_base& un_fluide)
{
  if (sub_type(Fluide_Diphasique, un_fluide))
    variables_internes().ref_fluide_diphasique = ref_cast(Fluide_Diphasique, un_fluide);
  else
    Navier_Stokes_Turbulent::associer_milieu_base(un_fluide);
}

const Milieu_base& Navier_Stokes_FT_Disc::milieu() const
{
  if (variables_internes().ref_fluide_diphasique.non_nul())
    return variables_internes().ref_fluide_diphasique.valeur();
  else
    return Navier_Stokes_Turbulent::milieu();
}

Milieu_base& Navier_Stokes_FT_Disc::milieu()
{
  if (variables_internes().ref_fluide_diphasique.non_nul())
    return variables_internes().ref_fluide_diphasique.valeur();
  else
    return Navier_Stokes_Turbulent::milieu();
}

/*! @brief Discretisation des champs utilises dans cette equation.
 *
 * Fonction appelee par Probleme_base::discretiser.
 *   B. Mathieu : a titre experimental, au lieu de dupliquer les noms
 *                des champs ici et dans "a_pour_champ_Fonc", on stocke
 *                les champs dans une liste. (voir a_pour_champ_fonc).
 *
 */
void Navier_Stokes_FT_Disc::discretiser()
{
  Navier_Stokes_Turbulent::discretiser();
  const Discretisation_base& dis = discretisation();
  const double temps = schema_temps().temps_courant();
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  LIST(REF(Champ_base)) & champs_compris = variables_internes().liste_champs_compris;

  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        "diffusivite", "m^2/s",
                        1 /* composantes */, temps,
                        champ_nu_);
  champs_compris.add(champ_nu_.valeur());
  //Nouvelle formulation
  champs_compris_.ajoute_champ(champ_nu_);

  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        "viscosite_dynamique", "kg/m.s",
                        1 /* composantes */, temps,
                        champ_mu_);
  champs_compris.add(champ_mu_.valeur());
  champs_compris_.ajoute_champ(champ_mu_);

  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        "masse_volumique", "kg/m3",
                        1 /* composantes */, temps,
                        champ_rho_elem_);
  champs_compris.add(champ_rho_elem_.valeur());
  champs_compris_.ajoute_champ(champ_rho_elem_);

  // La masse volumique discretisee sur les volumes de controle de la vitesse
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "masse_volumique_vnodes", "kg/m3",
                        1 /* composantes */, temps,
                        champ_rho_faces_);
  champs_compris.add(champ_rho_faces_.valeur());
  champs_compris_.ajoute_champ(champ_rho_faces_);

  // Variables internes
  dis.discretiser_champ("pression", mon_dom_dis,
                        "second_membre_projection", "",
                        1 /* composantes */, temps,
                        variables_internes().second_membre_projection);
  champs_compris.add(variables_internes().second_membre_projection.valeur());
  champs_compris_.ajoute_champ(variables_internes().second_membre_projection);
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        "second_membre_projection_jump", "",
                        1 /* composantes */, temps,
                        variables_internes().second_membre_projection_jump_);
  champs_compris.add(variables_internes().second_membre_projection_jump_.valeur());
  champs_compris_.ajoute_champ(variables_internes().second_membre_projection_jump_);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "gradient_pression", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().gradient_pression);
  champs_compris.add(variables_internes().gradient_pression.valeur());
  champs_compris_.ajoute_champ(variables_internes().gradient_pression);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "derivee_u_etoile", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().derivee_u_etoile);
  champs_compris.add(variables_internes().derivee_u_etoile.valeur());
  champs_compris_.ajoute_champ(variables_internes().derivee_u_etoile);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "terme_diffusion_vitesse", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_diffusion);
  champs_compris.add(variables_internes().terme_diffusion.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_diffusion);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "terme_convection_vitesse", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_convection);
  champs_compris.add(variables_internes().terme_convection.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_convection);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "terme_source_vitesse", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_source);
  champs_compris.add(variables_internes().terme_source.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "terme_source_interfaces", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_source_interfaces);
  champs_compris.add(variables_internes().terme_source_interfaces.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source_interfaces);
  if (dis.que_suis_je() == "VEFPreP1B")
    {
      dis.discretiser_champ("pression", mon_dom_dis,
                            "indicatrice_p1b", "",
                            1 /* composantes */, temps,
                            variables_internes().indicatrice_p1b);
      champs_compris.add(variables_internes().indicatrice_p1b.valeur());
      champs_compris_.ajoute_champ(variables_internes().indicatrice_p1b);
    }

  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "gradient_indicatrice", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().gradient_indicatrice);
  champs_compris.add(variables_internes().gradient_indicatrice.valeur());
  champs_compris_.ajoute_champ(variables_internes().gradient_indicatrice);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "potentiel_faces", "",
                        1 /* composante */, temps,
                        variables_internes().potentiel_faces);
  champs_compris.add(variables_internes().potentiel_faces.valeur());
  champs_compris_.ajoute_champ(variables_internes().potentiel_faces);
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        "potentiel_elements", "",
                        1 /* composante */, temps,
                        variables_internes().potentiel_elements);
  champs_compris.add(variables_internes().potentiel_elements.valeur());
  champs_compris_.ajoute_champ(variables_internes().potentiel_elements);
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "vitesse_delta_interface", "m/s",
                        -1 /* nb composantes par defaut */, 1, temps,
                        variables_internes().delta_u_interface);
  // Pour pouvoir faire filtrer_L2:
  variables_internes().delta_u_interface.associer_eqn(*this);
  champs_compris.add(variables_internes().delta_u_interface.valeur());
  champs_compris_.ajoute_champ(variables_internes().delta_u_interface);
  dis.discretiser_champ("pression", mon_dom_dis,
                        "pression_laplacien_d", "",
                        1 /* composante */, temps,
                        variables_internes().laplacien_d);
  champs_compris.add(variables_internes().laplacien_d.valeur());
  champs_compris_.ajoute_champ(variables_internes().laplacien_d);
  dis.discretiser_champ("temperature", mon_dom_dis,
                        "temperature_mpoint", "",
                        1 /* composante */, temps,
                        variables_internes().mpoint);
  champs_compris.add(variables_internes().mpoint.valeur());
  champs_compris_.ajoute_champ(variables_internes().mpoint);
  dis.discretiser_champ("temperature", mon_dom_dis,
                        "temperature_mpointv", "",
                        1 /* composante */, temps,
                        variables_internes().mpoint_vap);
  champs_compris.add(variables_internes().mpoint_vap.valeur());
  champs_compris_.ajoute_champ(variables_internes().mpoint_vap);
  // Pour variation temporelle dI_dt
  dis.discretiser_champ("pression", mon_dom_dis,
                        "derivee_temporelle_indicatrice", "",
                        1 /* composante */, temps,
                        variables_internes().derivee_temporelle_indicatrice); // nombre de composantes invalide en VEF
  champs_compris.add(variables_internes().derivee_temporelle_indicatrice.valeur());
  champs_compris_.ajoute_champ(variables_internes().derivee_temporelle_indicatrice);
  //HMS
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "terme_source_collisions", "",
                        1 /* 1 nb composante */, temps,
                        variables_internes().terme_source_collisions);
  champs_compris.add(variables_internes().terme_source_collisions.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source_collisions);

  // Eulerian Interfacial area :
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        "interfacial_area", "m2",
                        1 /* composante */, temps,
                        variables_internes().ai);
  champs_compris.add(variables_internes().ai.valeur());
  champs_compris_.ajoute_champ(variables_internes().ai);

  dis.discretiser_champ("pression", mon_dom_dis,
                        "num_compo", "",
                        1 /* composante */ , temps,
                        variables_internes().num_compo);
  champs_compris.add(variables_internes().num_compo.valeur());
  champs_compris_.ajoute_champ(variables_internes().num_compo);
  //fin HMS


  // debut EB
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        "terme_correction_trainee", "",
                        1 , temps,
                        variables_internes().terme_correction_trainee);
  champs_compris.add(variables_internes().terme_correction_trainee.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_correction_trainee);

// fin EB
  // Velocity jump "u0" computed for phase 0 :
  Nom nom = Nom("vitesse_jump0_") + le_nom();
  dis.discretiser_champ("vitesse", mon_dom_dis, nom, "m/s", -1 /* nb composantes par defaut */,1,  temps,
                        variables_internes().vitesse_jump0_);
  variables_internes().vitesse_jump0_.associer_eqn(*this);
  champs_compris.add(variables_internes().vitesse_jump0_.valeur());
  champs_compris_.ajoute_champ(variables_internes().vitesse_jump0_);
}

/*! @brief Methode surchargee de Navier_Stokes_std, appelee par Navier_Stokes_std::discretiser().
 *
 *   L'assembleur pression est particulier pour le front-tracking
 *   en VEF (en attendant qu'on factorise tous ces assembleurs pression)
 *
 */
void Navier_Stokes_FT_Disc::discretiser_assembleur_pression()
{
  const Discretisation_base& dis = discretisation();
  if (dis.que_suis_je() == "VDF" || dis.que_suis_je()=="VDF+")
    {
      // Assembleur pression generique (prevu pour rho_variable)
      assembleur_pression_.typer("Assembleur_P_VDF");
    }
  else if (dis.que_suis_je() == "VEFPreP1B")
    {
      // Assembleur pression generique (prevu pour rho_variable)
      assembleur_pression_.typer("Assembleur_P_VEFPreP1B");
    }

  Assembleur_base& assembleur = assembleur_pression_.valeur();
  assembleur.associer_domaine_dis_base(domaine_dis().valeur());
  // B. Mathieu, 07_2004 : premier jet de la methode, on resout en pression.
  // Version actuelle : pas d'increment de pression
  assembleur.set_resoudre_increment_pression(0);
  assembleur.set_resoudre_en_u(1);
}

/*! @brief methode appelee par Navier_Stokes_std::preparer_calcul.
 *
 * ..
 *
 */
void Navier_Stokes_FT_Disc::projeter()
{
  if (Process::je_suis_maitre() && limpr())
    Cerr << "Navier_Stokes_FT_Disc::projeter does nothing" << finl;
}

/*! @brief methode appelee par Probleme_base::preparer_calcul()
 *
 */
int Navier_Stokes_FT_Disc::preparer_calcul()
{
  Cerr << "Navier_Stokes_FT_Disc::preparer_calcul()" << finl;
  Equation_base::preparer_calcul();

  le_modele_turbulence.preparer_calcul();


  {
    // Si le mot cle equation_interfaces_proprietes_fluide a ete specifie,
    // l'indicatrice de phase correspondant a l'equation sert a calculer
    // les proprietes physiques du milieu.
    // Sinon, le milieu doit etre de type Fluide_Incompressible.

    REF(Transport_Interfaces_FT_Disc) & ref_equation =
      variables_internes().ref_eq_interf_proprietes_fluide;

    // debut EB
    if (ref_equation.non_nul())
      {
        Transport_Interfaces_FT_Disc& eq_transport = ref_equation.valeur();
        // Pour utiliser l'indicatrice aux aretes, il faut renseigner les proprietes du fluide pour l'operateur (terme_diffusif->associer_proprietes_fluide)
        // Il faut absolument utiliser une discretisation VDF+ pour avoir le bon operateur de diffusion (Op_Dift_VDF_var_Face_FT)
        // On pourrait l'utiliser pour autre chose que des particules solides mais il faudrait le valider
        //
        if (eq_transport.is_solid_particle())
          {
            // TODO FIXME WILL NOT WORK FOR VEF .... should test subtype ... // EB : test not required because the baltik fluid_particle_interaction (is_solid_particle=1) should only be used in VDF
            //ref_cast(Op_Dift_VDF_Face_FT, terme_diffusif.valeur()).associer_indicatrices(eq_transport.inconnue().valeurs(),eq_transport.get_indicatrice_aretes());
            terme_diffusif->associer_indicatrices(eq_transport.inconnue().valeurs(),eq_transport.get_indicatrice_aretes());
            const int formule_mu=fluide_diphasique().formule_mu();
            const double mu_fluide  =  fluide_diphasique().fluide_phase(1).viscosite_dynamique().valeur().valeurs()(0, 0);
            const double mu_particule  = fluide_diphasique().fluide_phase(0).viscosite_dynamique().valeur().valeurs()(0, 0);
            //ref_cast(Op_Dift_VDF_Face_FT, terme_diffusif.valeur()).associer_proprietes_fluide(formule_mu,mu_particule,mu_fluide);
            terme_diffusif->associer_proprietes_fluide(formule_mu,mu_particule,mu_fluide);
            /*
            const Particule_Solide& particule_solide=ref_cast(Particule_Solide,fluide_diphasique().fluide_phase(0));  // a modifier pour des particules bidisperses ou de tailles differentes
            const double rayon_compo=eq_transport.get_rayons_compo()(0); // a modifier pour des particules bidisperses ou de tailles differentes
            // Les lignes suivantes seront utiles tant que l'on fera du monodisperse, cela evite de tout recalculer a chaque fois et evite des tableaux inutiles et lourds.
            particule_solide.set_rayon(rayon_compo);
            particule_solide.set_diametre(2*rayon_compo);
            particule_solide.set_volume_compo(4 * M_PI * pow(rayon_compo, 3) / 3);
            particule_solide.set_masse_compo((4 * M_PI * pow(rayon_compo, 3) / 3)*rho_solide);

            */ // EB : Voir Particule_Solide::set_param
          }
      }
    // fin EB

    if (ref_equation.non_nul())
      {

        // Couplage Navier-Stokes / Fluide : les interfaces determinent la position des phases.
        // On utilise les proprietes du fluide diphasique

        if (Process::je_suis_maitre())
          {
            Cerr << "Initialisation of the physical properties (rho, mu, ...)\n"
                 << " based on the indicatrice field of the equation " << ref_equation.valeur().le_nom() << finl;

          }
        const int calcul_precis_indic_face=ref_equation.valeur().calcul_precis_indic_faces();
        FT_disc_calculer_champs_rho_mu_nu_dipha(domaine_dis().valeur(),
                                                fluide_diphasique(),
                                                ref_equation.valeur().
                                                get_update_indicatrice().valeurs(), // indicatrice
                                                ref_equation.valeur(), // EB : eq_transport
                                                champ_rho_elem_.valeur().valeurs(),
                                                champ_nu_.valeur().valeurs(),
                                                champ_mu_.valeur().valeurs(),
                                                champ_rho_faces_.valeur().valeurs(), schema_temps().temps_courant(), calcul_precis_indic_face);
      }
    else
      {

        // Pas de couplage Navier-Stokes / Fluide : le fluide est monophasique.

        const Fluide_Incompressible& phase_0 = ref_cast(Fluide_Incompressible,milieu());
        const Domaine_dis_base& zdis = domaine_dis().valeur();
        FT_disc_calculer_champs_rho_mu_nu_mono(zdis,
                                               phase_0,
                                               champ_rho_elem_,
                                               champ_mu_,
                                               champ_nu_,
                                               champ_rho_faces_);
      }

    // Premiere version :les termes sources sont homogenes a rho*v
    //   sources().associer_champ_rho(champ_rho_elem_.valeur());
    // Nouvelle version: les termes "sources()" sont homogenes a v.
    //   on ne fait rien.
  }

  DoubleTab& tab_vitesse = inconnue().valeurs();

  // On assemble la matrice de pression pour la premiere fois.
  assembleur_pression().valeur().assembler_rho_variable(matrice_pression_,
                                                        champ_rho_faces_.valeur());
  // Informe le solveur que la matrice a change :
  solveur_pression().valeur().reinit();
  // Calcul du second membre :
  //  div(vpoint) a l'interieur du domaine,
  //  prise en compte des conditions aux limites de pression/vitesse
  DoubleTab& secmem = variables_internes().second_membre_projection.valeurs();
  divergence.calculer(tab_vitesse, secmem);
  secmem *= -1;
  // Il faut faire ceci car on ne resout pas en "increment de pression":
  assembleur_pression_.valeur().modifier_secmem(secmem);

  // Ajout pour la sauvegarde au premier pas de temps si reprise
  la_pression.changer_temps(schema_temps().temps_courant());

  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),
                                     secmem,
                                     la_pression.valeurs()
                                    );
  assembleur_pression_.modifier_solution(la_pression->valeurs());
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)
  DoubleTab& gradP = variables_internes().gradient_pression.valeurs();
  gradient.calculer(la_pression.valeur().valeurs(), gradP);
  solveur_masse.appliquer(gradP);
  // Correction de la vitesse :
  if (projection_a_faire()) // Temporaire pour permettre de ne pas resoudre NS avec mettant operateurs nuls et projection_initiale 0
    {
      int i, j;
      const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
      const int n = tab_vitesse.dimension(0);
      const int m = tab_vitesse.line_size();

      for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
          tab_vitesse(i,j) -= gradP(i,j) / rho_faces(i);

      tab_vitesse.echange_espace_virtuel();
    }

  if (le_traitement_particulier.non_nul())
    {
      if (le_traitement_particulier.valeur().support_ok())
        le_traitement_particulier.valeur().
        associer_champ_masse_volumique(champ_rho_faces_.valeur());
      le_traitement_particulier.preparer_calcul_particulier();
    }

// debut EB
  const Discretisation_base& dis = discretisation();
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  if (refeq_transport.non_nul())
    {
      const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
      const Postraitement_Forces_Interfaces_FT& les_post_interf=eq_transport.postraitement_forces_interf();
      const double temps = schema_temps().temps_courant();
      const Domaine_dis_base& mon_domaine_dis = domaine_dis().valeur();
      LIST(REF(Champ_base)) & champs_compris = variables_internes().liste_champs_compris;
      if (les_post_interf.calcul_forces_theoriques_stokes_)
        {
          dis.discretiser_champ("vitesse", mon_domaine_dis,
                                "vitesse_stokes_th", "m/s",
                                3 , temps,
                                variables_internes().vitesse_stokes_th_);
          champs_compris.add(variables_internes().vitesse_stokes_th_.valeur());
          champs_compris_.ajoute_champ(variables_internes().vitesse_stokes_th_);

          dis.discretiser_champ("pression", mon_domaine_dis,
                                "pression_stokes_th", "pa",
                                1 , temps,
                                variables_internes().pression_stokes_th_);
          champs_compris.add(variables_internes().pression_stokes_th_.valeur());
          champs_compris_.ajoute_champ(variables_internes().pression_stokes_th_);
        }
    }
  // fin EB
  return 1;
}

const int& Navier_Stokes_FT_Disc::get_is_penalized() const
{
  return variables_internes().is_penalized;
}

const int& Navier_Stokes_FT_Disc::get_new_mass_source() const
{
  return variables_internes().new_mass_source_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_interfacial_area() const
{
  return variables_internes().ai.valeur().valeurs();
}

DoubleTab& Navier_Stokes_FT_Disc::get_set_interfacial_area()
{
  return variables_internes().ai.valeur().valeurs();
}

const DoubleTab& Navier_Stokes_FT_Disc::get_mpoint() const
{
  return variables_internes().mpoint.valeur().valeurs();
}

DoubleTab& Navier_Stokes_FT_Disc::get_set_mpoint()
{
  return variables_internes().mpoint.valeur().valeurs();
}

void Navier_Stokes_FT_Disc::preparer_pas_de_temps()
{
}

void Navier_Stokes_FT_Disc::mettre_a_jour(double temps)
{
  Navier_Stokes_Turbulent::mettre_a_jour(temps);

  champ_rho_elem_.mettre_a_jour(temps);
  champ_rho_faces_.mettre_a_jour(temps);
  champ_mu_.mettre_a_jour(temps);
  champ_nu_.mettre_a_jour(temps);
  variables_internes().num_compo.mettre_a_jour(temps); // EB
}

const SolveurSys& Navier_Stokes_FT_Disc::get_solveur_pression() const
{
  return solveur_pression_;
}

/*! @brief Calcul des forces de tension superficielles (F_sigma) et de la partie a rotationnel non nul de la gravite (G_rot) (si GRAVITE_GRAD_I) :
 *
 *   F_sigma = INTEGRALE sur le volume de controle (
 *             sigma_aux_faces * courbure_aux_faces * gradient(indicatrice)
 *             + gradient_sigma )
 *   G_rot   = INTEGRALE sur le volume de controle (
 *                phi * gradient(rho) )   (avec phi = potentiel de pesanteur)
 *
 * @param (gradient_indicatrice) le gradient de l'indicatrice issu de l'operateur "gradient", donc homogene a l'integrale du gradient sur les volumes de controle de la vitesse.
 * @param (potentiel_faces) un champ aux faces a une composante, ou on stocke le "potentiel aux faces"
 * @param (champ) le champ aux faces (meme discretisation que la vitesse) ou on stocke le terme source des forces superficielles.
 */
//
//  The method is no inter const because it changes a member of variables_internes() to store the Eulerian interfacial Area.
void Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles(const Maillage_FT_Disc& maillage,
                                                                 const Champ_base& gradient_indicatrice,
                                                                 Champ_base& potentiel_elements,
                                                                 Champ_base& potentiel_faces,
                                                                 Champ_base& champ)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  // Nombre de faces
  const int nb_faces = domaine_vf.nb_faces();
  {
    // Le champ et le gradient de l'indicatrice doivent etre aux faces
    assert(champ.valeurs().dimension(0) == nb_faces);
    assert(potentiel_faces.valeurs().dimension(0) == nb_faces);
    assert(gradient_indicatrice.valeurs().dimension(0) == nb_faces);
  }
  const int dim = Objet_U::dimension;

  // Calcul du "potentiel aux sommets du maillage"
  ArrOfDouble potentiel_sommets;
  potentiel_sommets.resize_array(maillage.nb_sommets(), Array_base::NOCOPY_NOINIT);

  // (ce potentiel est constant sur chaque portion connexe d'interface si
  //  l'interface est a l'equilibre).
  //  courbure * sigma + potentiel_gravite_phase_1 - potentiel_gravite_phase_0
  {
    // Pour l'instant, la tension de surface est constante
    const Fluide_Diphasique& fluide_dipha = fluide_diphasique();
    const double sigma = fluide_dipha.sigma();

    const int n = maillage.nb_sommets();
    const ArrOfDouble& courbure_sommets = maillage.get_update_courbure_sommets();
    //Calcul du "potentiel aux sommets"
    int i;
    {
      const double clipping_courbure_max = variables_internes().clipping_courbure_interface;
      int clip_counter = 0;
      for (i = 0; i < n; i++)
        {
          double c = courbure_sommets[i];
          // Clipping de la courbure: si la courbure est superieure a la
          // valeur maxi autorisee, on limite (permet de ne pas plomber le
          // pas de temps s'il y a une singularite geometrique dans le maillage)
          if (std::fabs(c) > clipping_courbure_max)
            {
              clip_counter++;
              c = ((c > 0) ? 1. : -1.) * clipping_courbure_max;
            }
          potentiel_sommets[i] = c * sigma;
        }
      clip_counter = mp_sum(clip_counter);
      if (clip_counter > 0 && je_suis_maitre())
        {
          Cerr << "Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles : clip_count "
               << clip_counter << finl;
        }
    }

    //Ajout des termes de gravite
    if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_GRAD_I)
      {
        if (milieu().a_gravite())
          {
            const double rho_0 = fluide_dipha.fluide_phase(0).masse_volumique().valeurs()(0,0);
            const double rho_1 = fluide_dipha.fluide_phase(1).masse_volumique().valeurs()(0,0);
            const double delta_rho = rho_1 - rho_0;

            // Pour l'instant : gravite uniforme g => phi(s) = - x scalaire g
            const DoubleTab& gravite = milieu().gravite().valeurs();
            if (gravite.nb_dim() != 2 || gravite.line_size() != dim)
              {
                Cerr << "Error for the method calculer_champ_forces_superficielles\n";
                Cerr << "gravite.dimension(1) != Objet_U::dimension" << finl;
                Process::exit();
              }
            const DoubleTab& sommets = maillage.sommets();

            for (i = 0; i < n; i++)
              {
                double p = 0.;
                for (int j = 0; j < dim; j++)
                  p += sommets(i,j) * gravite(0,j);

                if(is_repulsion)
                  {
                    double dx=0.;
                    double px=sommets(i,0);
                    if(px>maxx)
                      dx=px-maxx;
                    else if (px<minx)
                      dx=minx-px;

                    double dy=0.;
                    double py=sommets(i,1);
                    if(py>maxx)
                      dy=py-maxx;
                    else if (py<minx)
                      dy=minx-py;

                    p+=sqrt(dx*dx+dy*dy)*pente;
                  }

                potentiel_sommets[i] -= p * delta_rho;
              }
          }
      }
  }
#if 0
  // Calcul du "potentiel aux faces" : interpolation sur les faces des elements
  // traverses par l'interface.
  {
    // Tableau des facettes du maillage interfaces
    const IntTab& facettes = maillage.facettes();
    // Tableau des faces des elements :
    const IntTab& elem_faces = domaine_vf.elem_faces();
    const int nb_faces_par_element = elem_faces.line_size();
    // Le tableau "potentiel aux faces" a remplir :
    DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    valeurs_potentiel_faces = 0.;
    // Somme des poids des contributions ajoutees aux faces
    ArrOfDouble poids(nb_faces);
    poids = 0.;

    // Boucle sur les faces de l'interface
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_facette_elem = intersections.index_facette();
    double pot[3] = {0., 0., 0.};

    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        // Potentiel aux sommets de la facette:
        pot[0] = potentiel_sommets[facettes(fa7,0)];
        pot[1] = potentiel_sommets[facettes(fa7,1)];
        if (dim == 3)
          pot[2] = potentiel_sommets[facettes(fa7,2)];
        // Boucle sur les elements traverses par la facette
        int index = index_facette_elem[fa7];
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            // Calcul du potentiel au centre de l'intersection
            double p = 0.;
            int i;
            for (i = 0; i < 3; i++)
              p += data.barycentre_[i] * pot[i];
            // La contribution aux faces est ponderee par la surface d'intersection
            p *= data.surface_intersection_;

            // Ajout de la contribution aux faces de l'element
            const int element = data.numero_element_;
            for (i = 0; i < nb_faces_par_element; i++)
              {
                const int face = elem_faces(element, i);
                poids[face] += data.surface_intersection_;
                valeurs_potentiel_faces(face) += p;
              }
            index = data.index_element_suivant_;
          }
      }

    // Il reste a diviser les valeurs aux face par la somme des poids
    for (int face = 0; face < nb_faces; face++)
      {
        double p = poids[face];
        if (p > 0.)
          valeurs_potentiel_faces(face) /= p;
      }
    valeurs_potentiel_faces.echange_espace_virtuel();
  }
#else
  // Calcul du "potentiel aux elements" :
  DoubleTab poids(potentiel_elements.valeurs());
  {
    // Tableau des facettes du maillage interfaces
    const IntTab& facettes = maillage.facettes();
    // Le tableau "potentiel aux faces" a remplir :
    DoubleTab& valeurs_potentiel_elements = potentiel_elements.valeurs();
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    valeurs_potentiel_elements = 0.;
    // Somme des poids des contributions ajoutees aux faces
    poids = 0.;

    // Boucle sur les faces de l'interface
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_facette_elem = intersections.index_facette();
    double pot[3] = {0., 0., 0.};

    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        // Potentiel aux sommets de la facette:
        pot[0] = potentiel_sommets[facettes(fa7,0)];
        pot[1] = potentiel_sommets[facettes(fa7,1)];
        if (dim == 3)
          pot[2] = potentiel_sommets[facettes(fa7,2)];
        // Boucle sur les elements traverses par la facette
        int index = index_facette_elem[fa7];
        const double surface_facette = surface_facettes[fa7];
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            // Calcul du potentiel au centre de l'intersection
            double p = 0.;
            int i;
            for (i = 0; i < 3; i++)
              p += data.barycentre_[i] * pot[i];
            // La contribution aux faces est ponderee par la surface d'intersection
#ifdef AVEC_BUG_SURFACES
            const double surf = data.surface_intersection_;
#else
            const double surf = data.fraction_surface_intersection_ * surface_facette;
#endif
            p *= surf;

            // Ajout de la contribution a l'element
            const int element = data.numero_element_;
            valeurs_potentiel_elements(element) += p;
            poids(element) += surf;

            index = data.index_element_suivant_;
          }
      }
    valeurs_potentiel_elements.echange_espace_virtuel();
    poids.echange_espace_virtuel();
    if (champ.valeurs().line_size() > 1)   // VEF ?
      {
        // Compute a potential on elements that are neighbours of
        // elements containing an interface :
        //  For VEF, the gradient of the indicator function can be
        //  non zero on the faces of these elements.
        int element;
        const int nb_elements = poids.dimension(0);
        const IntTab& face_voisins = le_dom_dis.valeur().valeur().face_voisins();
        const Domaine_VF& domainevf = ref_cast(Domaine_VF, le_dom_dis.valeur().valeur());
        const IntTab& elem_faces = domainevf.elem_faces();
        const int nb_faces_par_element = elem_faces.line_size();
        DoubleVect copie_valeurs_potentiel_elements(valeurs_potentiel_elements);
        DoubleVect copie_poids(poids);
        for (element = 0; element < nb_elements; element++)
          {
            double potential = 0.; // sum of potentials of neighbouring elements
            double p = 0.;   // sum of weights of neighbouring elements
            int i_face;
            for (i_face = 0; i_face < nb_faces_par_element; i_face++)
              {
                const int face = elem_faces(element, i_face);
                const int elem_voisin_0 = face_voisins(face, 0);
                const int elem_voisin_1 = face_voisins(face, 1);
                const int elem_voisin = elem_voisin_0 + elem_voisin_1 - element;
                if (elem_voisin >= 0)   // Not a boundary of the domain ?
                  {
                    potential += copie_valeurs_potentiel_elements(elem_voisin);
                    p += copie_poids(elem_voisin);
                  }
              }
            const double old_weight = copie_poids(element);
            // Do not change values for elements that contain an interface:
            if (p > 0. && old_weight == 0.)
              {
                // Decrease weight so that values on faces adjacent to elements
                // containing an interface are almost untouched.
                static double reduction_factor = 0.1;
                potential = potential * reduction_factor;
                p = p * reduction_factor;
                // Assign value
                valeurs_potentiel_elements(element) = potential;
                poids(element) = p;
              }
          }
        valeurs_potentiel_elements.echange_espace_virtuel();
        poids.echange_espace_virtuel();
        Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles poids:",poids);
      }
  }

  // I take the opportunity to store the Eulerian Interfacial Area...
  DoubleTab& interfacial_area = get_set_interfacial_area();
  interfacial_area=poids;
  {
    // Boucle sur les faces
    int face;
    DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    valeurs_potentiel_faces = 0.;
    const DoubleTab& valeurs_potentiel_elements = potentiel_elements.valeurs();
    const int nb_faces_pot = valeurs_potentiel_faces.dimension(0);
    const IntTab& face_voisins = le_dom_dis.valeur().valeur().face_voisins();
    for (face = 0; face < nb_faces_pot; face++)
      {
        double p = 0.; // Somme des poids des deux elements voisins
        double pot = 0.; // Somme des potentiels
        // Boucle sur les deux elements voisins de la face
        for (int i = 0; i < 2; i++)
          {
            int element = face_voisins(face, i);
            if (element >= 0)
              {
                // If neighbour exists (not a boundary face)
                p += poids(element);
                pot += valeurs_potentiel_elements(element);
              }
          }
        if (p > 0.)
          valeurs_potentiel_faces(face) = pot / p;
      }
    valeurs_potentiel_faces.echange_espace_virtuel();
    Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles valeurs_potentiel_faces:",valeurs_potentiel_faces);
  }
#endif
  // Derniere operation : calcul du champ
  //   champ = potentiel_faces * gradient_indicatrice
  {
    const DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    const DoubleTab& valeurs_gradient_i = gradient_indicatrice.valeurs();
    DoubleTab& valeurs_champ = champ.valeurs();
    const int nb_compo = valeurs_champ.line_size(); // 1 for VDF, 3 for VEF

    for (int face = 0; face < nb_faces; face++)
      {
        const double p = valeurs_potentiel_faces(face);
        for (int i = 0; i < nb_compo; i++)
          valeurs_champ(face, i) = valeurs_gradient_i(face, i) * p;
      }

    valeurs_champ.echange_espace_virtuel();
    Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles valeurs_champ:",valeurs_champ);
  }
}

double compute_indic_ghost(const int elem, const int num_face, const double indic,
                           const int normale_sortante_au_domaine,
                           const Domaine_VF& domVF, const Maillage_FT_Disc& maillage)
{
  double indic_ghost= 0.;
  ArrOfDouble normale(3), normale_face(3); // Normale sortante de I=0 vers I=1
  const double face_surface = domVF.face_surfaces(num_face); //  ==  0. ? 1. : domVF.face_surfaces(num_face);
  if (est_egal(face_surface, 0.))
    {
      // On est sur l'axe de rotation du bidim_axi :
      return indic;
    }
  for (int i=0; i< Objet_U::dimension; i++)
    normale_face[i] = domVF.face_normales(num_face, i)/face_surface; // pour normalisation
  const double norm = maillage.compute_normale_element(elem, true /* NORMALIZE */, normale);
  if (est_egal(norm, 0.))
    // L'element est pure, non traverse.
    return indic;
  // normale_sortante_au_domaine pour regler le signe du prodscal avec la normale_SORTANTE
  const double prodscal = normale_sortante_au_domaine * dotproduct_array(normale, normale_face);
  double val = 1./sqrt(2.); // Could be improved. Corresponds to cubic cell and Indic=0.5
  if (prodscal <-val)
    indic_ghost = 0.;
  else if (prodscal <0)
    // Linear interp to indic if 0
    indic_ghost=indic*(1+prodscal/val); // prodscal is negative
  else if (prodscal <val)
    indic_ghost=indic+(1.-indic)*(prodscal/val)  ;
  else
    indic_ghost = 1.;
  return indic_ghost;
}


// debut EB
// On envoie la liste des composantes reelles aux procs des zones_inf
IntTabFT echanger_recevoir_list_num_compo(ArrOfInt& liste_compo_reelles_to_send, const ArrOfInt& liste_zone_inf, const Schema_Comm_FT& comm)
{
  IntTabFT list_compo_recv;
  list_compo_recv.set_smart_resize(1);
  int nb_elem_recv=0;
  const int nb_compo_reelles_to_send=liste_compo_reelles_to_send.size_array();

  comm.begin_comm();

  for (int ind_pe_dest=0; ind_pe_dest<liste_zone_inf.size_array(); ind_pe_dest++)
    {
      int PE_dest=liste_zone_inf(ind_pe_dest);
      for(int ind_compo=0; ind_compo<nb_compo_reelles_to_send; ind_compo++)
        {
          const int num_compo=liste_compo_reelles_to_send(ind_compo);
          assert(PE_dest!=Process::me());
          comm.send_buffer(PE_dest)  << num_compo;
        }
    }
  comm.echange_taille_et_messages();

  list_compo_recv.resize_array(0);
  const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
  const int nb_recv_pe = recv_pe_list.size_array();
  for (int i=0; i<nb_recv_pe; i++)
    {
      const int pe_source = recv_pe_list[i];
      Entree& buffer = comm.recv_buffer(pe_source);
      while(1)
        {
          int num_compo_recv=-1;
          buffer >> num_compo_recv;
          if (buffer.eof())
            break;
          if (num_compo_recv<0) Process::exit();
          nb_elem_recv++;

          list_compo_recv.append_array(num_compo_recv);
        }
    }

  comm.end_comm();
  list_compo_recv.set_smart_resize(0);
  list_compo_recv.resize_array(nb_elem_recv);

  return list_compo_recv;

}
// fin EB

// HMS
// EB : modif de la fonction en cours pour integration dans Trio. Done ?
void Navier_Stokes_FT_Disc::calculer_champ_forces_collisions(const DoubleTab& indicatrice, DoubleTab& valeurs_champ, const Transport_Interfaces_FT_Disc& eq_transport, Transport_Interfaces_FT_Disc& eq_transport_non_const, REF(Transport_Interfaces_FT_Disc)& refeq_transport, const Maillage_FT_Disc& maillage)
{

  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_forces_collisions", 0);
  statistiques().begin_count(count);

//<editor-fold desc="Determiniations des vitesses et positions de chaques particule">
  Cerr << "Navier_Stokes_FT_Disc::calculer_champ_forces_collisions" << finl;
  Process::imprimer_ram_totale();
  /***********************************************/
  // ETAPE 0 : On recupere les grandeurs d'interet
  /***********************************************/
  Modele_Collision_FT& modele_collision_particule = eq_transport_non_const.collision_interface_particule();
  int& compteur_collisions=modele_collision_particule.compteur_collisions();
  compteur_collisions+=1;

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  // recuperation des vitesse compo et des centres de gravites

  DoubleTab& positions=refeq_transport.valeur().get_positions_compo();
  DoubleTab& vitesses=refeq_transport.valeur().get_vitesses_compo();
  const Fluide_Diphasique& mon_fluide = fluide_diphasique();
  const Particule_Solide& particule_solide=ref_cast(Particule_Solide,fluide_diphasique().fluide_phase(0));
  double ed = particule_solide.e_dry();
  const double& rayon_compo=particule_solide.rayon_collision(); // a modifier pour des particules bidisperses ou de tailles differentes
  const int nb_facettes = maillage.nb_facettes();
  IntVect compo_connexes_facettes(nb_facettes); // Init a zero

  int nb_compo_tot=positions.dimension(0);

  static const DoubleTab& positions_bords=modele_collision_particule.position_bords();
  int nb_bords = positions_bords.dimension(1);

  assert(nb_compo_tot != 0);
  ArrOfInt elem_cg;
  domaine_vf.domaine().chercher_elements_FT(positions, elem_cg);

  /***********************************************/
  // ETAPE 0 : Correspondance num eulerien au temps precedent et mauvais num lagrangien actuel
  /***********************************************/
  DoubleTab correct_vitesses(nb_compo_tot, dimension), correct_positions(nb_compo_tot, dimension);
  IntVect copy_elem_cg, correctNum(nb_compo_tot);
  DoubleTab& F_now=modele_collision_particule.get_F_now();
  F_now.resize(nb_compo_tot, nb_compo_tot + nb_bords);
  F_now=0;

  // EB : F_old et F_new sont necessaires pour savoir si la collision entre 2 particules est deja en cours ou pas
  // En effet, pour calculer la 2eme raideur du modele de collision,
  // on calcule un Stokes de collision en utilisant la vitesse d'impact. Il faut garder en memoire ce Stokes d'un pas de temps a l'autre de la collision.
  // Ces 3 tableaux de sont pas optimal car de taille Np^2, avec N_p le nombre de particules.
  DoubleTab& F_old=modele_collision_particule.get_F_old();
  DoubleTab& raideur=modele_collision_particule.get_raideur();
  DoubleTab& e_eff=modele_collision_particule.get_e_eff();

  /***********************************************/
  // ETAPE 1 : Initialisation des tableaux s'il n'y a pas encore eu de collisions.
  // Sinon, permutation des tableaux de positions et de vitesses pour assurer
  // la correspondance entre le numero eulerien au temps n-1 et le numero
  // lagrangien au temps n (suite a la renumerotation (voir search_connex_components_local_FT et compute_global_connex_components_FT)
  /***********************************************/
  // ATTENTION, il ne faut faire un resize qu'au premier pas de temps (t=0).
  // Lors d'une reprise, on recupere le "compteur collisions" pour ne pas ecraser les tableaux. Sinon, il y aura une erreur dans le calcul de la 2eme raideur a chaque
  // reprise de calcul.
  if (compteur_collisions==1)
    {
      F_old.resize(nb_compo_tot, nb_compo_tot + nb_bords);
      raideur.resize(nb_compo_tot, nb_compo_tot + nb_bords);
      e_eff.resize(nb_compo_tot, nb_compo_tot + nb_bords);

      F_old = 0;
      raideur=0;
      e_eff=0;
    }
  /***********************************************/
  // ETAPE 2 : Mise a jour des tables de Verlet
  // Methode inspiree du papier : X. Fang, J. Tang, H. Luo., Granular damping analysis using an improved discrete element approach, J. Sound Vib., 2007
  /***********************************************/
  static IntLists table_Verlet;
  static IntLists table_Verlet_bord;
  static ArrOfIntFT liste_composantes_reelles;
  ArrOfInt list_compo_to_check;
  static int nb_compo_reelles;
  static const int is_LC_on=modele_collision_particule.is_LC_activated();
  nb_compo_reelles= is_LC_on ? 0: nb_compo_tot;
  if (modele_collision_particule.is_detection_Verlet()==1)
    {
      double s_Verlet = modele_collision_particule.get_s_Verlet(); // initialise dans Transport_Interfaces_FT_Disc::preparer_calcul() par 0.3*2*rayon_compo
      int& nb_dt_Verlet = modele_collision_particule.get_nb_dt_Verlet();
      int& dt_compute_Verlet = modele_collision_particule.get_dt_compute_Verlet();
      int& nb_pas_dt_max_Verlet=modele_collision_particule.get_nb_pas_dt_max_Verlet();
      if (nb_dt_Verlet >= dt_compute_Verlet || schema_temps().nb_pas_dt()==0)
        {
          Cerr << "Mise a jour de la table de Verlet - " ;

          // ETAPE A : mise a jour des Linked Cells (methode LC).
          // Ici, chaque LC est une zone de calcul d'un processeur. On envoie aux LC des zones inf (pour recevoir les LC des zones sup) le numero des compos reelles
          // On identifie les elements (particules) reels pour chaque processeur
          if (is_LC_on)
            {
              nb_compo_reelles=0;
              liste_composantes_reelles.resize_array(0);
              liste_composantes_reelles.set_smart_resize(1);
              for (int compo=0; compo<elem_cg.size_array(); compo++)
                {
                  if (elem_cg(compo)>-1 && elem_cg(compo)<domaine_vf.nb_elem()) liste_composantes_reelles.append_array(compo), nb_compo_reelles++;
                }
              liste_composantes_reelles.set_smart_resize(0);
              liste_composantes_reelles.resize_array(nb_compo_reelles);

              Process::barrier(); // on est oblige d'attendre que tous les procs aient mis a jour leur liste avant de faire l'echange
              const ArrOfIntFT& liste_zone_inf=modele_collision_particule.get_liste_zone_inf();
              const Schema_Comm_FT& schema_com= maillage.get_schema_comm_FT();
              IntTabFT list_compo_virt=echanger_recevoir_list_num_compo(liste_composantes_reelles, liste_zone_inf, schema_com);

              if (nb_compo_reelles>0)
                {
                  list_compo_to_check.resize(liste_composantes_reelles.size_array()+list_compo_virt.size_array());
                  for (int k=0; k<liste_composantes_reelles.size_array()+list_compo_virt.size_array(); k++)
                    {
                      if (k<liste_composantes_reelles.size_array()) list_compo_to_check(k)=liste_composantes_reelles(k);
                      else list_compo_to_check(k)=list_compo_virt(k-liste_composantes_reelles.size_array());
                    }
                }
            }

          // ETAPE B : mise a jour des tables de Verlet
          // On calcule la distance entre chaque particule de list_compo_to_check. Si elle est inferieure a s (=0.3 * Dp), alors on garde en memoire la paire
          double max_vi_glob=0;
          double max_vi=0.;
          if (nb_compo_reelles>0)
            {
              table_Verlet=0;
              table_Verlet_bord=0;
              table_Verlet.dimensionner(nb_compo_reelles);
              table_Verlet_bord.dimensionner(nb_compo_reelles);

              if (is_LC_on)
                {
                  for (int ind_compo_i=0; ind_compo_i< nb_compo_reelles; ind_compo_i++)
                    {
                      int compo_i=liste_composantes_reelles(ind_compo_i);
                      double vi=fabs(sqrt( pow(vitesses(compo_i,0),2)+pow(vitesses(compo_i,1),2)+pow(vitesses(compo_i,2),2) ));
                      if (vi>max_vi) max_vi=vi;
                      // distance particule-particule
                      for (int ind_compo_j=ind_compo_i+1; ind_compo_j< list_compo_to_check.size_array(); ind_compo_j++)
                        {
                          int compo_j=list_compo_to_check(ind_compo_j);
                          double dij=sqrt(pow(positions(compo_j,0)-positions(compo_i,0),2)
                                          +pow(positions(compo_j,1)-positions(compo_i,1),2)
                                          +pow(positions(compo_j,2)-positions(compo_i,2),2));
                          if ((dij-2*rayon_compo)<=s_Verlet) table_Verlet[ind_compo_i].add(compo_j);
                        }

                      // distance particule-paroi
                      for (int ind_bord=0; ind_bord<nb_bords; ind_bord++)
                        {
                          int ori = ind_bord < dimension ? ind_bord : ind_bord - dimension;
                          double dij=fabs(positions_bords(compo_i,ind_bord)-positions(compo_i,ori));
                          if ((dij-2*rayon_compo)<=s_Verlet) table_Verlet_bord[ind_compo_i].add(ind_bord);
                        }
                    }
                }
              else // Methode Verlet "sequentiel" - tous les procs font les memes operations
                {
                  for (int compo_i=0; compo_i< nb_compo_reelles; compo_i++)
                    {
                      double vi=fabs(sqrt(pow(vitesses(compo_i,0),2)+pow(vitesses(compo_i,1),2)+pow(vitesses(compo_i,2),2)));
                      if (vi>max_vi) max_vi=vi;
                      // distance particule-particule
                      for (int compo_j=compo_i+1; compo_j< nb_compo_tot; compo_j++)
                        {
                          double dij=sqrt(pow(positions(compo_j,0)-positions(compo_i,0),2)
                                          +pow(positions(compo_j,1)-positions(compo_i,1),2)
                                          +pow(positions(compo_j,2)-positions(compo_i,2),2));
                          if ((dij-2*rayon_compo)<=s_Verlet) table_Verlet[compo_i].add(compo_j);
                        }
                      // distance particule-paroi
                      for (int ind_bord=0; ind_bord<nb_bords; ind_bord++)
                        {
                          int ori = ind_bord < dimension ? ind_bord : ind_bord - dimension;
                          double dij=fabs(positions_bords(compo_i,ind_bord)-positions(compo_i,ori));
                          if ((dij-2*rayon_compo)<=s_Verlet) table_Verlet_bord[compo_i].add(ind_bord);
                        }
                    }
                }
            }

          max_vi_glob=mp_max(max_vi); // si on fait du LC, alors on prend le max, sinon tous les procs ont la meme valeur donc mp_max(max_vi)=max_vi
          // calcul du prochain pas de temps ou il faudra remettre a jour les tables de Verlet
          int pas_dt_compute_min= (max_vi_glob>0) ? static_cast<int>(floor((s_Verlet/(2*max_vi_glob))/schema_temps().pas_de_temps())) : nb_pas_dt_max_Verlet;
          dt_compute_Verlet=std::min(pas_dt_compute_min,nb_pas_dt_max_Verlet);
          Cerr << "prochaine mise a jour dans  " <<dt_compute_Verlet << " pas de temps." << finl;
          nb_dt_Verlet=0;
        }
      nb_dt_Verlet++;
    }
  /***********************************************/
  // ETAPE 2 : calcul de la force de contact (discrete, par particule)
  /***********************************************/
  static const IntVect& orientation = domaine_vdf.orientation();
  double dt = schema_temps().pas_de_temps();
  static const double rho_solide = mon_fluide.fluide_phase(0).masse_volumique().valeurs()(0, 0);
  static const double rho_fluide = mon_fluide.fluide_phase(1).masse_volumique().valeurs()(0, 0);
  static const double mu_fluide  = rho_fluide * mon_fluide.fluide_phase(1).viscosite_cinematique().valeur().valeurs()(0, 0);

  static int indic_phase_fluide = 1; //, indic_phase_solide=0;
  static double d_act = modele_collision_particule.get_d_act_lub(); // a modifier dans un cas bidisperse ou particules de tailles differentes
  static double d_sat = modele_collision_particule.get_d_sat_lub(); // a modifier dans un cas bidisperse ou particules de tailles differentes

  //<editor-fold desc="Calcule des forces de collisions">
  DoubleTab& forces_solide=modele_collision_particule.get_forces_solide();
  forces_solide.resize(nb_compo_tot, dimension);
  forces_solide=0;
  DoubleVect& collision_detected=modele_collision_particule.get_collisions_detected();
  collision_detected.resize(nb_compo_tot) ;
  collision_detected = 0 ;
  static const double seuil=1e-10;
  static double volume_compo = particule_solide.volume_compo();
  static double masse_compo = particule_solide.masse_compo();
  //static double volume_compo_voisin = volume_compo;// a modifier pour un cas bidisperse
  //static double masse_compo_voisin = masse_compo;// a modifier pour un cas bidisperse
  static int isModeleLubrification =modele_collision_particule.modele_lubrification();
  DoubleTab dX(dimension), dU(dimension);

  IntTab Collision(nb_compo_tot,nb_compo_tot+nb_bords);
  if (modele_collision_particule.is_detection_Verlet()==1)
    {
      for (int ind_compo_i=0; ind_compo_i< nb_compo_reelles; ind_compo_i++)
        {
          int compo_i=is_LC_on ? liste_composantes_reelles(ind_compo_i) : ind_compo_i;
          // Premiere boucle : Collision Particule-Particule uniquement
          double rayon_eff = rayon_compo/2;  // a modifier dans un cas bidisperse ou particules de tailles differentes
          double masse_eff = masse_compo/2;
          for (int ind_compo_j=0; ind_compo_j< (table_Verlet[ind_compo_i].size()); ind_compo_j++)
            {
              int compo_j=table_Verlet[ind_compo_i][ind_compo_j];
              if (compo_i==compo_j) Process::exit("Navier_Stokes_FT_Disc::calculer_champ_forces_collisions compo_i=compo_j");
              dX = 0;
              dU = 0;

              for (int d = 0; d < dimension; d++)
                {
                  dX(d) = positions(compo_i, d) - positions(compo_j, d);
                  dU(d) = vitesses(compo_i, d) - vitesses(compo_j, d);
                }
              double dist_cg = sqrt(local_carre_norme_vect(dX));
              double dist_int = dist_cg - 2 * rayon_compo;
              F_now(compo_i, compo_j) = 0;
              //Cerr << "dist_int " << dist_int << finl;
              if (dist_int <= 0) // contact
                {
                  //<editor-fold desc="Calcule de la norme et de la vitesse relative normale">
                  DoubleTab norm(dimension);
                  for (int d = 0; d < dimension; d++) norm(d) = dX(d) / dist_cg;
                  double prod_sacl = local_prodscal(dX,dU);
                  DoubleTab dUn(dimension);
                  for (int d = 0; d < dimension; d++) dUn(d) = (prod_sacl / dist_cg) * norm(d);
                  double vitesseRelNorm =sqrt(local_carre_norme_vect(dUn));

                  //<editor-fold desc="Modele de lubrification">
                  // Calcul des forces de lubrifications

                  double d_int = fabs(dist_int) / rayon_compo; // a modifier dans un cas bidisperse ou particules de tailles differentes

                  if (isModeleLubrification && d_int <= d_act )
                    {
                      double lambda     =  0.5 / d_int - 9 * log(d_int) / 20 - 3 * d_int * log(d_int) / 56 ;
                      double lambda_act =  0.5 / d_act - 9 * log(d_act) / 20 - 3 * d_act * log(d_act) / 56 ;
                      double lambda_sat =  0.5 / d_sat - 9 * log(d_sat) / 20 - 3 * d_sat * log(d_sat) / 56 ;
                      double delta_lambda=0;
                      if (d_sat < d_int && d_int <= d_act)
                        delta_lambda= (lambda - lambda_act);
                      if (0 < d_int && d_int <= d_sat)
                        delta_lambda= (lambda_sat - lambda_act);
                      for (int d = 0; d < dimension; d++)
                        {
                          double force_lubrification = -6 * M_PI * mu_fluide * rayon_compo * dUn(d) * delta_lambda; // a modifier dans un cas bidisperse ou particules de tailles differentes
                          continue;
                          forces_solide(compo_i, d) += +force_lubrification / volume_compo;
                          forces_solide(compo_j, d) += -force_lubrification / volume_compo;
                        }
                    }
                  //remplisage de l'indicateur de collisions
                  collision_detected(compo_i) +=  1 ;
                  collision_detected(compo_j) +=  1 ;

                  F_now(compo_i, compo_j) = 1;
                  // EB : F_now et F_old : pour savoir dans quelle partie de la collision on est (pour le modele hybride)
                  int isFirstStepOfCollision = F_now(compo_i, compo_j) > F_old(compo_i, compo_j);
                  //<editor-fold desc="schema semi implicte">
                  // EB : A l'endroit de la collision : le mur apparait comme une sphere de rayon rayon_compo(compo) (methode de HMS)
                  DoubleTab next_dX(dimension);
                  for (int d = 0; d < dimension; d++) next_dX(d) = dX(d) + dt * dU(d);
                  double next_dist_cg = sqrt(local_carre_norme_vect(next_dX));
                  /*double next_dist_int = CollisionParticuleParticule ? next_dist_cg -
                  					 (rayon_compo + rayons_compo_(voisin)) :
                  					 next_dist_cg - 2 * rayons_compo_(compo); */
                  double next_dist_int = next_dist_cg - 2 * rayon_compo; // a modifier par la ligne precedente dans un cas bidisperse ou particules de tailles different
                  double Stb = rho_solide * 2 * rayon_eff * vitesseRelNorm / (9 * mu_fluide);
                  DoubleTab force_contact(dimension);
                  modele_collision_particule.calculer_force_contact(force_contact, isFirstStepOfCollision, dist_int, next_dist_int, norm, dUn, masse_eff, compo_i, compo_j, Stb, ed, vitesseRelNorm, dt, prod_sacl);

                  for (int d = 0; d < dimension; d++)
                    {
                      forces_solide(compo_i, d) += fabs(force_contact(d))<=seuil ? 0 : force_contact(d) / volume_compo;
                      forces_solide(compo_j, d) -= fabs(force_contact(d))<=seuil ? 0 :  force_contact(d) / volume_compo;
                    }
                  Collision(compo_i,compo_j)=1;
                }
              F_old(compo_i, compo_j) = F_now(compo_i, compo_j);
            }
          rayon_eff = rayon_compo;  // a modifier dans un cas bidisperse ou particules de tailles differentes
          masse_eff = masse_compo;

          // Deuxieme boucle : Collision Particule-Paroi
          for (int ind_bord=0; ind_bord < (table_Verlet_bord[ind_compo_i].size()); ind_bord++)
            {
              int bord = table_Verlet_bord[ind_compo_i][ind_bord];
              dU=0;
              dX=0;
              int ori = bord < dimension ? bord : bord - dimension;
              dX(ori) = positions(compo_i, ori) - positions_bords(compo_i,bord);

              for (int d = 0; d < dimension; d++) dU(d) = vitesses(compo_i, d);
              double dist_cg = sqrt(local_carre_norme_vect(dX));
              double dist_int = dist_cg - 2 * rayon_compo;

              F_now(compo_i,nb_compo_tot+bord) = 0;
              if (dist_int <= 0) //contact
                {
                  //<editor-fold desc="Calcule de la norme et de la vitesse relative normale">
                  DoubleTab norm(dimension);
                  for (int d = 0; d < dimension; d++) norm(d) = dX(d) / dist_cg;
                  double prod_sacl = local_prodscal(dX,dU);
                  DoubleTab dUn(dimension);
                  for (int d = 0; d < dimension; d++) dUn(d) = (prod_sacl / dist_cg) * norm(d);
                  double vitesseRelNorm =sqrt(local_carre_norme_vect(dUn));

                  //<editor-fold desc="Modele de lubrification">
                  // Calcul des forces de lubrifications
                  double d_int = fabs(dist_int) / rayon_compo; // a modifier dans un cas bidisperse ou particules de tailles differentes

                  if (isModeleLubrification && d_int <= d_act )
                    {
                      double lambda     = 1 / d_int - log(d_int) / 5 - d_int * log(d_int) / 21;
                      double lambda_act = 1 / d_act - log(d_act) / 5 - d_act * log(d_act) / 21;
                      double lambda_sat = 1 / d_sat - log(d_sat) / 5 - d_sat * log(d_sat) / 21;
                      double delta_lambda=0;
                      if (d_sat < d_int && d_int <= d_act)
                        delta_lambda= (lambda - lambda_act);

                      if (0 < d_int && d_int <= d_sat)
                        delta_lambda= (lambda_sat - lambda_act);

                      for (int d = 0; d < dimension; d++)
                        {
                          double force_lubrification = -6 * M_PI * mu_fluide * rayon_compo * dUn(d) * delta_lambda; // a modifier dans un cas bidisperse ou particules de tailles differentes
                          continue;
                          forces_solide(compo_i, d) += +force_lubrification / volume_compo;
                        }
                    }

                  // EB : Si dist_int <0 alors il faut appliquer la collision
                  collision_detected(compo_i) +=  0.1 ;
                  F_now(compo_i, nb_compo_tot+bord) = 1;
                  // EB : F_now et F_old : pour savoir dans quelle partie de la collision on est (pour le modele hybride)
                  int isFirstStepOfCollision = F_now(compo_i,nb_compo_tot+bord) > F_old(compo_i,nb_compo_tot+bord);
                  // EB : A l'endroit de la collision : le mur apparait comme une sphere de rayon rayon_compo(compo) (methode de HMS)
                  DoubleTab next_dX(dimension);
                  for (int d = 0; d < dimension; d++) next_dX(d) = dX(d) + dt * dU(d);
                  double next_dist_cg = sqrt(local_carre_norme_vect(next_dX));
                  /*double next_dist_int = CollisionParticuleParticule ? next_dist_cg -
                  					 (rayon_compo + rayons_compo_(voisin)) :
                  					 next_dist_cg - 2 * rayons_compo_(compo); */
                  double next_dist_int = next_dist_cg - 2 * rayon_compo; // a modifier par la ligne precedente dans un cas bidisperse ou particules de tailles differentes

                  double Stb = rho_solide * 2 * rayon_eff * vitesseRelNorm / (9 * mu_fluide);
                  DoubleTab force_contact(dimension);
                  int voisin=nb_compo_tot+bord;
                  modele_collision_particule.calculer_force_contact(force_contact, isFirstStepOfCollision, dist_int, next_dist_int, norm, dUn, masse_eff, compo_i, voisin, Stb, ed, vitesseRelNorm, dt, prod_sacl);
                  for (int d = 0; d < dimension; d++)
                    {
                      forces_solide(compo_i, d) +=  fabs(force_contact(d))<=seuil ? 0 : force_contact(d) / volume_compo;
                    }

                  Collision(compo_i,nb_compo_tot+bord)=1;
                }
              F_old(compo_i, nb_compo_tot+bord) = F_now(compo_i, nb_compo_tot+bord);
              if (F_old(compo_i, nb_compo_tot+bord)>1) F_old(compo_i, nb_compo_tot+bord) =1;
            }
        }
      if (is_LC_on)
        {
          mp_max_for_each_item(Collision);
          mp_sum_for_each_item(forces_solide);
          mp_sum_for_each_item(collision_detected);
          mp_max_for_each_item(F_old);
          mp_max_for_each_item(F_now);
          mp_max_for_each_item(raideur);
          mp_max_for_each_item(e_eff);
        }
    }
  else
    {
      //const DoubleVect& rayons_compo_=eq_transport.get_rayons_compo();
      double volume_compo_voisin=0;
      double delta_n = modele_collision_particule.delta_n();
      for (int compo = 0; compo < nb_compo_tot; compo++)
        {
          //volume_compo = 4 * M_PI * pow(rayons_compo_[compo], 3) / 3;
          //masse_compo = volume_compo* rho_solide;
          for (int voisin = compo + 1; voisin < nb_compo_tot + nb_bords; voisin++)
            {
              if (voisin<nb_compo_tot)
                {
                  volume_compo_voisin = volume_compo;//4 * M_PI * pow(rayons_compo_(voisin), 3) / 3;
                }
              //<editor-fold desc="Calcule de la distance entre deux interfaces">
              // EB : On calcule les distances entre les particules et entre particules et parois pour savoir si
              // on active la collision ou non
              dX = 0;
              dU = 0;
              int CollisionParticuleParticule = voisin < nb_compo_tot;
              if (CollisionParticuleParticule)
                {
                  for (int d = 0; d < dimension; d++)
                    {
                      dX(d) = positions(compo, d) - positions(voisin, d);
                      dU(d) = vitesses(compo, d) - vitesses(voisin, d);
                    }
                }
              else
                {
                  int bord = voisin - nb_compo_tot;
                  int ori = bord < dimension ? bord : bord - dimension;
                  dX(ori) = positions(compo, ori) - positions_bords(compo,bord);
                  for (int d = 0; d < dimension; d++) dU(d) = vitesses(compo, d);
                }
              double dist_cg = sqrt(local_carre_norme_vect(dX));
              if (dist_cg == 0)
                {
                  Cerr << "ERROR : dist_cg = 0 entre " << compo << " et " << voisin << finl;
                  exit();
                }

              //double dist_int = CollisionParticuleParticule ? dist_cg - (rayons_compo_(compo) + rayons_compo_(voisin)) :
              //                dist_cg - 2 * rayons_compo_(compo);
              double dist_int = dist_cg - 2 * rayon_compo;

              //Calcule de la norme et de la vitesse relative normale">
              DoubleTab norm(dimension);
              for (int d = 0; d < dimension; d++) norm(d) = dX(d) / dist_cg;
              double prod_sacl = local_prodscal(dX,dU);
              DoubleTab dUn(dimension);
              for (int d = 0; d < dimension; d++) dUn(d) = (prod_sacl / dist_cg) * norm(d);
              double vitesseRelNorm =sqrt(local_carre_norme_vect(dUn));

              //<editor-fold desc="Modele de lubrification">
              // Calcul des forces de lubrifications
              double d_int = fabs(dist_int) / rayon_compo;
              d_act = 2.0 * delta_n / rayon_compo;
              d_sat = 0.1 * delta_n /rayon_compo;

              if (isModeleLubrification && d_int <= d_act )
                {
                  double lambda     = CollisionParticuleParticule ? 0.5 / d_int - 9 * log(d_int) / 20 - 3 * d_int * log(d_int) / 56 : 1 / d_int - log(d_int) / 5 - d_int * log(d_int) / 21;
                  double lambda_act = CollisionParticuleParticule ? 0.5 / d_act - 9 * log(d_act) / 20 - 3 * d_act * log(d_act) / 56 : 1 / d_act - log(d_act) / 5 - d_act * log(d_act) / 21;
                  double lambda_sat = CollisionParticuleParticule ? 0.5 / d_sat - 9 * log(d_sat) / 20 - 3 * d_sat * log(d_sat) / 56 : 1 / d_sat - log(d_sat) / 5 - d_sat * log(d_sat) / 21;

                  double delta_lambda=0;

                  if (d_sat < d_int && d_int <= d_act)
                    delta_lambda= (lambda - lambda_act);

                  if (0 < d_int && d_int <= d_sat)
                    delta_lambda= (lambda_sat - lambda_act);

                  for (int d = 0; d < dimension; d++)
                    {
                      double force_lubrification = -6 * M_PI * mu_fluide * rayon_compo * dUn(d) * delta_lambda;;
                      continue;
                      forces_solide(compo, d) += +force_lubrification / volume_compo;
                      if (!CollisionParticuleParticule) continue; //collision avec un bord -> pas de force sur le bord
                      forces_solide(voisin, d) += -force_lubrification / volume_compo_voisin;
                    }
                }
              //</editor-fold>
              F_now(compo, voisin) = 0;
              // EB : Si dist_int <0 alors il faut appliquer la collision
              if (dist_int <= 0) //contact
                {
                  //remplisage de l'indicateur de collisions
                  if (CollisionParticuleParticule)
                    {
                      collision_detected(compo) +=  1 ;
                      collision_detected(voisin) +=  1 ;
                    }
                  else
                    {
                      collision_detected(compo) +=  0.1 ;
                    }

                  F_now(compo, voisin) = 1;
                  // EB : F_now et F_old : pour savoir dans quelle partie de la collision on est (pour le modele hybride)
                  int isFirstStepOfCollision = F_now(compo, voisin) > F_old(compo, voisin);
                  //<editor-fold desc="schema semi implicte">
                  // EB : A l'endroit de la collision : le mur apparait comme une sphere de rayon rayon_compo(compo) (methode de HMS)
                  DoubleTab next_dX(dimension);
                  for (int d = 0; d < dimension; d++) next_dX(d) = dX(d) + dt * dU(d);
                  double next_dist_cg = sqrt(local_carre_norme_vect(next_dX));
                  //double next_dist_int = CollisionParticuleParticule ? next_dist_cg -
                  //                     (rayons_compo_(compo) + rayons_compo_(voisin)) :
                  //                   next_dist_cg - 2 * rayons_compo_(compo);
                  double next_dist_int = next_dist_cg - 2 * rayon_compo; // a modifier par la ligne precedente dans un cas bidisperse ou particules de tailles differentes
                  //</editor-fold>
                  //new writing
                  if(1)
                    {
                      double rayon_eff = CollisionParticuleParticule ? rayon_compo/2
                                         : rayon_compo; //bord // a modifier dans un cas bidisperse ou particules de tailles differentes
                      double masse_eff = CollisionParticuleParticule ? masse_compo/2
                                         : masse_compo; //bord // a modifier dans un cas bidisperse ou particules de tailles differentes
                      double Stb = rho_solide * 2 * rayon_eff * vitesseRelNorm / (9 * mu_fluide);
                      DoubleTab force_contact(dimension);
                      modele_collision_particule.calculer_force_contact(force_contact, isFirstStepOfCollision, dist_int, next_dist_int, norm, dUn, masse_eff, compo, voisin, Stb, ed, vitesseRelNorm, dt, prod_sacl);

                      for (int d = 0; d < dimension; d++)
                        {
                          forces_solide(compo, d) += fabs(force_contact(d))<=seuil ? 0 : force_contact(d) / volume_compo;
                          if (!CollisionParticuleParticule) continue; // collision avec un bord -> pas de force sur le bord
                          forces_solide(voisin, d) -= fabs(force_contact(d))<=seuil ? 0 :  force_contact(d) / volume_compo_voisin;
                        }
                    }
                  Collision(compo,voisin)=1;
                }
              F_old(compo, voisin) = F_now(compo, voisin);
            } // fin boucle parti
        } // fin boucle compo
    }

  if (Process::je_suis_maitre()&& schema_temps().limpr_fpi() && nb_compo_tot>1)
    {
      SFichier Liste_Collision;
      const Navier_Stokes_FT_Disc& mon_eq = *this;
      ouvrir_fichier(Liste_Collision,"liste_collision",1,mon_eq);
      schema_temps().imprimer_temps_courant(Liste_Collision);
      for (int compo_i=0; compo_i<nb_compo_tot; compo_i++)
        {
          for (int compo_j=0; compo_j<nb_compo_tot+nb_bords; compo_j++)
            {
              if (Collision(compo_i,compo_j)>0) Liste_Collision << " " << compo_i << "-" << compo_j;
            }
        }
      Liste_Collision << finl;
    }

  /***********************************************/
  // ETAPE 3 : Correspondance entre le numero eulerien de la compo au temps n et le numero eulerien au temps n-1
  // Pour cela, on "marque" les elements fluide par -1 et les elements solide par 1 avant de rentrer dans search_connex_components_local
  // et compute_global_connex_components
  /***********************************************/

  //<editor-fold desc="Application des forces de collisions sur les faces du domaine">
  static const DoubleVect& volumes_entrelaces = domaine_vf.volumes_entrelaces();
  static const int nb_elem = domaine_vf.domaine().nb_elem();

  const DoubleTab& valeurs_v = inconnue().valeur().valeurs();

  DoubleTab& rms_vitesse=refeq_transport.valeur().get_rms_vitesses_compo();
  rms_vitesse.resize(nb_compo_tot,dimension);
  DoubleTab& moy_vitesse=refeq_transport.valeur().get_moy_vitesses_compo();
  moy_vitesse.resize(nb_compo_tot,dimension);
  DoubleTab& moy_vitesse_solide_carre=refeq_transport.valeur().get_moy_vitesses_carre_compo();
  moy_vitesse_solide_carre.resize(nb_compo_tot,dimension);

  double volume_solide=0;
  double vmoy_face=0;

  rms_vitesse=0.;
  moy_vitesse=0.;
  moy_vitesse_solide_carre=0.;

  DoubleTrav num_compo;
  domaine_vf.domaine().creer_tableau_elements(num_compo);

  IntVect num_lagrange(nb_compo_tot);

  // EB : on marque les elements eulerien avant renumerotation locale dans dans search_connex_components_local
  // puis globale dans compute_global_connex_components
  for (int elem = 0; elem < nb_elem; elem++)
    {
      if (modele_collision_particule.force_elem_diphasique())
        {
          //les elements diphasiques sont compris dans num compo
          num_compo[elem] = (indicatrice[elem] == indic_phase_fluide) ? -1 : 1;
        }
      else
        {
          //les elements diphasiques ne sont pas compris dans num compo
          num_compo[elem] = (indicatrice[elem] != 1 - indic_phase_fluide ) ? -1 : 1;
        }
    }
  num_compo.echange_espace_virtuel();

  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
  const int nb_local_connex_components = search_connex_components_local(elem_faces, faces_elem, num_compo);
  nb_compo_tot = compute_global_connex_components(num_compo, nb_local_connex_components);

  // on s'assure que le numero eulerien corresepond au bon numero lagrangien
  //remplissage du tableau de corespandance indice interface lagrange

  for (int compo = 0; compo < nb_compo_tot; compo++)
    {
      int elem = elem_cg[compo];
      if (elem != -1)
        {
          int num_euler = static_cast<int>(num_compo[elem]);
          num_lagrange[num_euler] = compo;
        }
    }
  mp_max_for_each_item(num_lagrange);

  //syncronisation entre les indice lagrangien et eulerien dans num_compo
  //on parcour toutes les cellules du maillage, on identifie son indice eulerien
  // on utilise l'indice eulerien  pour trouver l'indice lagrangien correspandant
  // on remplace l'indice eulerien par l'indice lagrengien

  /***********************************************/
  // ETAPE 4 : permutation de num compo suivant l'etape precedente et calcul de la vitesse moyenne et rms des compos
  /***********************************************/
  const DoubleVect& volumes_maille = domaine_vf.volumes();
  DoubleTab volumes_euler(nb_compo_tot);

  for (int elem = 0; elem < nb_elem; elem++)
    {
      if (num_compo[elem] != -1)
        {
          num_compo[elem] =  static_cast<double>(num_lagrange[static_cast<int>(num_compo[elem])]);
        }

      if (indicatrice(elem) != 1)
        {
          //volume au elements
          volumes_euler[0] += (1 - indicatrice[elem]) * volumes_maille[elem];
          //Cerr << finl << elem << " " << indicatrice[elem] << " " << volumes_maille[elem] << " ";
          //volume au face
        }
      if (indicatrice[elem]==0)
        {
          for (int dim=0; dim<dimension; dim++)
            {
              int compo= static_cast<int>(num_compo[elem]);
              vmoy_face=0.5*(valeurs_v(elem_faces(elem,dim))+valeurs_v(elem_faces(elem,dim+dimension))); // vmoy_face est plutot un vmoy_elem ici, revoir les notations
              moy_vitesse(compo,dim) +=vmoy_face*volumes_maille(elem);
              moy_vitesse_solide_carre(compo,dim) +=pow(vmoy_face,2)*volumes_maille(elem);
            }
          volume_solide+=volumes_maille(elem);
        }
    }

  double vol_solide=mp_sum(volume_solide);

  moy_vitesse/=vol_solide;
  moy_vitesse_solide_carre/=vol_solide;

  mp_sum_for_each_item(moy_vitesse);
  mp_sum_for_each_item(moy_vitesse_solide_carre);


  if (Process::je_suis_maitre())
    {
      for (int dim=0; dim<dimension; dim++)
        {
          for (int compo=0; compo<nb_compo_tot; compo++)
            {
              if (moy_vitesse_solide_carre(compo,dim)>0)
                {
                  rms_vitesse(compo,dim)=sqrt(fabs(pow(moy_vitesse(compo,dim),2)-moy_vitesse_solide_carre(compo,dim)));
                  //Cerr << "Vmoy au carre Vcarre moy " <<pow(moy_vitesse(dim),2) <<"\t"<< moy_vitesse_solide_carre(dim) << finl;
                }
              else
                {
                  rms_vitesse(compo,dim)=0;
                }
            }
        }
    }


  /***********************************************/
  // ETAPE 5 : Application de la force de contact sur les faces euleriennes. La force "discrete" est appliquee de maniere volumique
  // sur les faces euleriennes constituant les particules.
  /***********************************************/
  for (int elem = 0; elem < nb_elem; elem++)
    {
      int compo = static_cast<int>(num_compo[elem]);
      if (compo != -1)
        {
          for (int ifac = 0; ifac < 2 * dimension; ifac++)
            {
              int face = elem_faces(elem, ifac);
              int ori = orientation(face);
              valeurs_champ(face) =
                (1 - indicatrice_faces(face)) * volumes_entrelaces(face) * forces_solide(compo, ori);
            }
        }
    }
  valeurs_champ.echange_espace_virtuel();
  variables_internes().num_compo.valeur().valeurs() = num_compo;  // champ pret pour postraitement
// fin
  Cerr << "FIN Navier_Stokes_FT_Disc::calculer_champ_forces_collisions" << finl;

  statistiques().end_count(count);
}

// EB
/*! @brief calcule la force de correction hydrodynamique a appliquer pour compenser la sous resolution des gradients de vitesse et de pression a l'interface
 * fluide - particule. La correction est de la force Fc=F_stokes * alpha/(N^beta). F_stokes = 6pi*mu_f*r_p*u_inf : force theorique de stokes (egale a la force
 * hydrodynamique calcule sur maillage fin. Dans le cas de stokes, le calcul de la force est convergee a 40 mailles par diametre). alpha et beta, des coefficients
 * de correlation. N : nombre de mailles par diametre de particule. /!\ Pour le moment, la methode est uniquement validee pour une particule en sedimentation dans un milieu infini, pour Re<=0.1.
 * Pour plus de details, voir la publi :
 *  E.Butaye, A. Toutant, S. Mer and F. Bataille. Development of Particle Resolved - Subgrid Corrected Simulations: Hydrodynamics force calculation and flow sub-resolution corrections. Computers and Fluids, 2023.
 */
void Navier_Stokes_FT_Disc::calculer_correction_trainee( DoubleTab& valeurs_champ, const Transport_Interfaces_FT_Disc& eq_transport,Transport_Interfaces_FT_Disc& eq_transport_non_const, REF(Transport_Interfaces_FT_Disc)& refeq_transport, const Maillage_FT_Disc& maillage)
{
  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_correction_trainee", 0);
  statistiques().begin_count(count);

  Cerr << "Navier_Stokes_FT_Disc::calculer_correction_trainee" << finl;
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  //const Domaine& domaine = domaine_vf.domaine();
  const DoubleVect& volumes_entrelaces = domaine_vf.volumes_entrelaces();
  //const DoubleVect& volume = domaine_vf.volumes();
  const Fluide_Diphasique& mon_fluide = fluide_diphasique();
  const Particule_Solide& particule_solide=ref_cast(Particule_Solide,mon_fluide.fluide_phase(0));
  const DoubleVect& rayon_compo=eq_transport.get_rayons_compo();
  const double mu_f = mon_fluide.fluide_phase(1).viscosite_dynamique().valeur().valeurs()(0, 0);
  const double mu_p = particule_solide.viscosite_dynamique().valeur().valeurs()(0,0);
  const double phi_mu=mu_p/mu_f;
  const double rho_f = mon_fluide.fluide_phase(1).masse_volumique().valeurs()(0, 0);
  const double rho_p = particule_solide.masse_volumique().valeurs()(0, 0);
  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntVect& orientation = domaine_vdf.orientation();

  const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();

  const DoubleTab& vitesses =eq_transport.get_vitesses_compo();
  DoubleVect vitesse_compo(dimension);
  int nb_compo_tot = vitesses.dimension(0);

  const DoubleTab& force_pression= get_force_tot_pression_interf();
  const DoubleTab& force_frottements= get_force_tot_pression_interf();
  DoubleVect correction_trainee(nb_compo_tot);

  DoubleVect& longueurs = Modele_Collision_FT::get_longueurs();
  IntVect& nb_noeuds=Modele_Collision_FT::get_nb_noeuds();

  DoubleTab& gravite = milieu().gravite().valeurs();
  DoubleVect vect_g(dimension);
  for (int i=0; i<dimension; i++) vect_g(i)=gravite(0,i);
  const double norme_g = sqrt(local_carre_norme_vect(vect_g));

  // vitesse de sedimentation dans le cas d'un ecoulement de stokes

  // coeff definit par une fonction d'optimisation sur python
  const double alpha=variables_internes().alpha_correction_trainee_;
  const double beta=variables_internes().beta_correction_trainee_;


  double Re_p=0., Cd_p_Abraham=0.;

  DoubleVect& tab_num_compo=variables_internes().num_compo.valeur().valeurs();

  // On calcule la direction des particules pour application la correction de la trainee
  DoubleTab direction_compo(nb_compo_tot,dimension);
  const int correction_proportionnelle=variables_internes().proportionnel_;
  const int extension_reynolds=variables_internes().extension_reynolds_;

  // La condition qui suit separe deux methodes de calcul. On pourrait simplifier en regroupant les parties communes mais cela impliquerait
  // de faire deux boucles sur le nombre de composantes au lieu d'une.
  if (!correction_proportionnelle) // on calcule la direction de la particule avec sa vitesse
    {
      double norme_vitesse=0.;
      for (int compo=0; compo<nb_compo_tot; compo++)
        {
          double N=(nb_noeuds(0)-1)/(longueurs(0)/(2.*rayon_compo(compo)));
          double u_inf=2./3.*(pow(rayon_compo(compo),2)*norme_g/mu_f)*((1+phi_mu)/(2.+3.*phi_mu)*(rho_f-rho_p));
          for (int dim=0; dim<dimension; dim++) vitesse_compo(dim)=vitesses(compo,dim);
          norme_vitesse=sqrt(local_carre_norme_vect(vitesse_compo));
          for (int dim=0; dim<dimension; dim++)
            {
              if(norme_vitesse>0) direction_compo(compo,dim)=vitesses(compo,dim)/norme_vitesse;
              else direction_compo(compo,dim)=0.;
            }

          correction_trainee(compo)=fabs(6.*M_PI*mu_f*rayon_compo(compo)*u_inf*(alpha/pow(N,beta)));
          Cerr << "Correction trainee "<<correction_trainee(compo)<< finl;
          if (extension_reynolds)
            {
              Re_p=rho_f*norme_vitesse*2.*rayon_compo(compo)/mu_f;
              Cd_p_Abraham = (norme_vitesse>0) ? (24./(pow(9.06,2)))*pow(9.06/sqrt(Re_p)+1.,2) : 0.;
              correction_trainee(compo)*=Cd_p_Abraham;
            }

        }
    }
  else // on calcule direction d'application de la correction avec les composantes de la force de trainee
    {
      double norme_vitesse=0.;
      double norme_force_trainee=0.;
      for (int compo=0; compo<nb_compo_tot; compo++)
        {
          DoubleVect la_force_trainee(dimension);
          if (schema_temps().nb_pas_dt()>0)
            {
              for (int dim=0; dim<dimension; dim++) la_force_trainee(dim)=force_pression(compo,dim)+force_frottements(compo,dim);
            }
          else la_force_trainee=0;

          norme_force_trainee=sqrt(local_carre_norme_vect(la_force_trainee));
          for (int dim=0; dim<dimension; dim++)
            {
              if(norme_force_trainee>0) direction_compo(compo,dim)=-la_force_trainee(dim)/norme_force_trainee; // - car la trainee s'oppose a la vitesse
              else direction_compo(compo,dim)=0.;
            }
          double N=(nb_noeuds(0)-1)/(longueurs(0)/(2.*rayon_compo(compo)));
          correction_trainee(compo)=norme_force_trainee*(alpha/pow(N,beta));
          if (extension_reynolds)
            {
              for (int dim=0; dim<dimension; dim++) vitesse_compo(dim)=vitesses(compo,dim);
              norme_vitesse=sqrt(local_carre_norme_vect(vitesse_compo));
              Re_p=rho_f*norme_vitesse*2.*rayon_compo(compo)/mu_f;
              Cd_p_Abraham = (norme_vitesse>0) ? (24./(pow(9.06,2)))*pow(9.06/sqrt(Re_p)+1.,2) : 0.;
              correction_trainee(compo)*=Cd_p_Abraham;
            }
        }

    }
  const int faces_diphasiques= variables_internes().faces_diphasiques_;
  int nb_faces=domaine_vdf.nb_faces();

  const ArrOfInt& faces_doubles = domaine_vdf.faces_doubles();
  DoubleVect somme_volume_particule(dimension);
  somme_volume_particule= 0.;
  valeurs_champ=0;

  // On discretise la correction sur les faces solide et diphasiques
  for (int face=0; face<nb_faces; face++)
    {
      if ((indicatrice_faces(face)<1 && faces_diphasiques) || (indicatrice_faces(face)==1 && !faces_diphasiques) )
        {
          double contribution_face=(faces_doubles[face] == 1) ? 0.5 : 1.;

          int elem_gauche=face_voisins(face,0);
          int elem_droite=face_voisins(face,1);
          // si on a pas acces a l'element, ie elem=-1, alors on prend tab_num_compo(elem)=-1 car tab_num_compo(-1) non defini
          int num_compo_gauche=elem_gauche>=0 ? static_cast<int>(tab_num_compo(elem_gauche)) : -1;
          int num_compo_droite=elem_droite>=0 ? static_cast<int>(tab_num_compo(elem_droite)) : -1;
          int max_num_compo=std::max(num_compo_gauche,num_compo_droite);
          valeurs_champ(face)=-(1-indicatrice_faces(face))*volumes_entrelaces(face)*correction_trainee(max_num_compo)*direction_compo(max_num_compo,orientation(face))*contribution_face; // On met un signe "-" car la correction appliquee doit etre opposee a la direction de la particule

          somme_volume_particule(orientation(face))+=(1.-indicatrice_faces(face))*volumes_entrelaces(face)*contribution_face;

        }
    }
  mp_sum_for_each_item(somme_volume_particule);

  for (int face=0; face<nb_faces; face++)
    {
      if (indicatrice_faces(face)<1)
        {
          valeurs_champ(face)/=somme_volume_particule(orientation(face));
        }
    }
  valeurs_champ.echange_espace_virtuel();

  statistiques().end_count(count);

}


/*! @brief Calcul du gradient de l'indicatrice.
 *
 * Ce gradient est utilise pour calculer le second membre de l'equation de qdm,
 *   contenant les termes de tension de surface.
 *   En VEF, on commence par creer un champ P1B a partir du champ P0
 *   et on calcule le gradient.
 *   Design de classe a revoir pour separer VDF et VEF...
 *
 *  Parametre : indicatrice
 *  Signification : un champ aux elements (l'espace virtuel doit etre a jour)
 *  Parametre : gradient_i
 *  Signification : un champ discretise comme la vitesse dans lequel
 *   on met gradient(indicatrice).
 *
 */
#include<Neumann_sortie_libre.h>
#include<Dirichlet.h>
#include<Dirichlet_homogene.h>
#include<Periodique.h>
#include<Symetrie.h>
#include<Sortie_libre_rho_variable.h>
void Navier_Stokes_FT_Disc::calculer_gradient_indicatrice(
  const Champ_base& indicatrice,
  const DoubleTab& distance_interface_sommets,
  Champ_base& gradient_i)
{
  if (gradient_i.que_suis_je() == "Champ_Fonc_Face")
    {
      gradient.calculer(indicatrice.valeurs(), gradient_i.valeurs());
    }
  else
    {
      int i;
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
      const Domaine& domaine = domaine_vf.domaine();
      const int nb_elem_tot = domaine.nb_elem_tot();
      //const int nb_sommets = domaine.nb_som();
      //const IntTab & les_elems = domaine.les_elems();
      //const int nb_sommets_par_element = les_elems.dimension(1);

      // Calcul d'une indicatrice p1bulle
      DoubleTab& indic_p1b = variables_internes().indicatrice_p1b->valeurs();
      // Verification du support du champ indicatrice_p1b
      if (dimension==2 && indic_p1b.size_totale()!=nb_elem_tot+domaine.nb_som_tot())
        {
          Cerr << "The method Navier_Stokes_FT_Disc::calculer_gradient_indicatrice is developped" << finl;
          Cerr << "only for an indicatrice field discretized like P0+P1 for the 2D dimension case."<<finl;
          Cerr << "Please change the discretization." << finl;
          Process::exit();
        }
      if (dimension==3 && indic_p1b.size_totale()!=nb_elem_tot+domaine.nb_som_tot()+domaine.nb_aretes_tot() && indic_p1b.size_totale()!=nb_elem_tot+domaine.nb_som_tot())
        {
          Cerr << "The method Navier_Stokes_FT_Disc::calculer_gradient_indicatrice is developped" << finl;
          Cerr << "only for an indicatrice field discretized like P0+P1 or P0+P1+Pa for the 3D dimension case."<<finl;
          Cerr << "Please change the discretization." << finl;
          Process::exit();
        }

      // Le champ P1bulle contient
      //   nelements valeurs au centre des elements
      // suivies de
      //   nsommets valeurs aux sommets
      indic_p1b = 0.;
      // Valeurs aux centres des elements = indicatrice a l'element
      // On recopie des valeurs sur les elements virtuels car on en
      // a besoin pour le calcul des valeurs aux sommets.
      for (i = 0; i < nb_elem_tot; i++)
        {
          indic_p1b(i) = indicatrice(i);
        }
#if 0
      // Premiere version : pas terrible parce que le stencil est tres large
      //  et il faut calculer la "courbure de l'interface" sur une grande largeur
      //  autour de l'interface.

      // Valeurs aux sommets = moyenne des valeurs des elements adjacents
      //  ponderee par le volume des elements.
      // Il y a la un choix a faire qui peut etre important...
      // (autre idee : projection L2 du champ discontinu sur l'espace
      // P1bulle)
      const DoubleVect& volume = domaine_vf.volumes();

      // Somme des poids aux sommets:
      ArrOfDouble poids(nb_sommets);
      poids = 0.;
      // Boucle sur les elements reels et virtuels :
      for (i = 0; i < nb_elem_tot; i++)
        {
          // Ajout d'une contribution a chaque sommet reel de l'element
          double p = volume(i);
          double x = indicatrice(i);
          int j;
          for (j = 0; j < nb_sommets_par_element; j++)
            {
              const int sommet = les_elems(i, j);
              if (sommet < nb_sommets)   // sommet reel ?
                {
                  // Le tableau des ddl contient nb_elem_tot valeurs aux elements
                  // suivies de nb_sommets_tot valeurs aux sommets.
                  // L'inconnue associee au "sommet" se trouve a l'indice
                  // nb_elem_tot + sommet.
                  indic_p1b(nb_elem_tot + sommet) += p * x;
                  poids[sommet] += p;
                }
            }
        }
      // Division par le poids
      for (i = 0; i < nb_sommets; i++)
        {
          const double p = poids[i];
          indic_p1b(nb_elem_tot + i) /= p;
        }
#else
#if 0
      // Calcul de l'indicatrice aux sommets. Deuxieme version, pour
      // limiter le stencil de gradient_indicatrice :
      // l'indicatrice vaut 1 si le sommet est adjacent a un element
      // de phase 1, 0 si le sommet est adjacent a un element de phase 0

      // 1) on met une valeur invalide dans les inconnues:
      for (i = 0; i < nb_sommets; i++)
        indic_p1b(nb_elem_tot + i) = -1;
      // 2) On met 1. ou 0. pour les sommets adjacents a un element
      //    dont l'indicatrice vaut 1. ou 0.
      for (i = 0; i < nb_elem_tot; i++)
        {
          const double indic_element = indic_p1b(i);
          if (indic_element == 0. || indic_element == 1.)
            {
              int j;
              for (j = 0; j < nb_sommets_par_element; j++)
                {
                  const int sommet = les_elems(i,j);
                  indic_p1b(nb_elem_tot + sommet) = indic_element;
                }
            }
        }
      // 3) Pour les sommets indecis, on prend 0 ou 1 selon le signe
      //    de la fonction distance pour ce sommet.
      static const double valeur_distance_invalide = -1e30;
      int error_count = 0;
      for (i = 0; i < nb_sommets; i++)
        {
          const double valeur = indic_p1b(nb_elem_tot + i);
          if (valeur < 0.)
            {
              double newval;
              const double d = distance_interface_sommets(i);
              if (d < valeur_distance_invalide)
                {
                  error_count++;
                  newval = 0.;
                }
              else
                {
                  newval = (d >= 0.) ? 1. : 0.;
                }
              indic_p1b(nb_elem_tot + i) = newval;
            }
        }
      if (error_count > 0)
        {
          Process::Journal()
              << "Navier_Stokes_FT_Disc::calculer_gradient_indicatrice\n"
              << error_count << " sommets ont une valeur indeterminee" << finl;
        }
#else
      // On met l'indicatrice sur P0 et 0 sur P1
#endif
#endif
      indic_p1b.echange_espace_virtuel();

      gradient.calculer(indic_p1b, gradient_i.valeurs());
    }

  const bool ghost_correction = (variables_internes_->OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC);
  if (ghost_correction)
    {
      // Correction du gradient a la ligne de contact :
      const DoubleTab& inco=indicatrice.valeurs();
      DoubleTab& resu=gradient_i.valeurs();
      const int nbdimr = resu.line_size();
      //      assert_espace_virtuel_vect(inco);
      const Domaine_VF& zvf = ref_cast(Domaine_VF, domaine_dis().valeur());
      const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis().valeur();
      const DoubleVect& face_surfaces = zvf.face_surfaces();
      const IntTab& face_voisins = domaine_dis().valeur().face_voisins();

      double coef;
      int n0, n1;
      // Boucle sur les bords pour traiter les conditions aux limites
      int ndeb, nfin, num_face;
      for (int n_bord=0; n_bord<zvf.nb_front_Cl(); n_bord++)
        {
          // pour chaque Condition Limite on regarde son type
          // Si face de Dirichlet ou de Symetrie on ne fait rien
          // Si face de Neumann on calcule la contribution au terme source

          const Cond_lim& la_cl = zcldis.les_conditions_limites(n_bord);
          Cerr << que_suis_je() << "::calculer_gradient_indicatrice() correction du gradient a la CL : " <<  la_cl.valeur() << finl;
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          // Correction periodicite :
          if (sub_type(Periodique,la_cl.valeur()))
            {
              for (num_face=ndeb; num_face<nfin; num_face++)
                {
                  n0 = face_voisins(num_face,0);
                  n1 = face_voisins(num_face,1);
                  if (!est_egal(n0,n1))
                    {
                      Cerr << "Periodic boundary condition with FT is not supported yet." << finl;
                      Process::exit();
                    }
                  coef = face_surfaces(num_face);//*porosite_surf(num_face);
                  for (int k=0; k<nbdimr; k++)
                    {
                      const int normale_sortante_au_domaine = (n0 == -1) ? 1 : -1 ; // Si on a le ghost dans la case 0, la normale sortante pointe vers "x-", on met donc -1 dans normale_sortante_au_domaine
                      const double dSk = normale_sortante_au_domaine*zvf.face_normales(num_face, k); // its magnitude is the surface.
                      // dSk/coef should be either 0 or 1...
                      assert(dSk/coef>-10.*Objet_U::precision_geom*Objet_U::precision_geom);
                      Cerr << "Check if sign of nk is compatible with the expression" << finl;
                      const double nk = zvf.face_normales(num_face, k);
                      Cerr << "Check if sign of nk is compatible with the expression" << finl;
                      Cerr << "nk=" << nk << finl;
                      Cerr << "dSk=" << dSk << finl;
                      Cerr << "num_face=" << num_face << finl;
                      Cerr << "voisin 0/1 =" << zvf.face_voisins(num_face, 0)<< " " << zvf.face_voisins(num_face, 1) << finl;
                      Cerr << "code to validate" << finl;
                      Process::exit();
                      resu(num_face,k) = dSk*(inco(n1) - inco(n0));
                    }
                }
            }
          else if (sub_type(Symetrie,la_cl.valeur()))
            ;
          else if ( (sub_type(Dirichlet,la_cl.valeur()))
                    ||
                    (sub_type(Neumann_sortie_libre,la_cl.valeur()))
                    ||
                    (sub_type(Dirichlet_homogene,la_cl.valeur()))
                  )
            {
              REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
              const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
              const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
              const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
              const ArrOfInt& index_elem = intersections.index_elem();
              for (num_face=ndeb; num_face<nfin; num_face++)
                {
                  n0 = face_voisins(num_face,0);
                  n1 = face_voisins(num_face,1);
                  int elem = n0+n1+1;
                  //double surface_totale = 0.;
                  for (int k=0; k<nbdimr; k++)
                    resu(num_face,k) = 0.;
                  if (index_elem[elem]<0)
                    // element is pure
                    continue;
                  {
                    const double indic = inco(elem);
                    int normale_sortante_au_domaine = (n0 == -1) ? -1 : 1 ; // Si on a le ghost dans la case 0, la normale sortante pointe vers "x-", on met donc -1 dans normale_sortante_au_domaine
                    const double indic_ghost = compute_indic_ghost(elem, num_face,
                                                                   indic, normale_sortante_au_domaine,
                                                                   zvf, maillage);
                    //const double delta = zvdf.dim_elem(elem, orientation(num_face));
                    //const double volume = zvdf.volumes(elem);
                    coef = face_surfaces(num_face);//*porosite_surf(num_face);
                    if (nbdimr ==1)
                      {
                        resu(num_face)=coef*normale_sortante_au_domaine*(indic_ghost-indic);
                        assert(std::abs(coef-zvf.face_normales(num_face, zvf.orientation(num_face)))
                               < Objet_U::precision_geom*Objet_U::precision_geom);
                      }
                    else
                      {
                        for (int k=0; k<nbdimr; k++)
                          {
                            normale_sortante_au_domaine =1 ; // En VEF, la normale dSk est toujours sortante?
                            const double dSk = normale_sortante_au_domaine*zvf.face_normales(num_face, k); // its magnitude is the surface.
                            // dSk/coef should be either 0 or 1...
                            assert(dSk/coef>-10.*Objet_U::precision_geom*Objet_U::precision_geom);
                            Cerr << "Check if sign of nk is compatible with the expression" << finl;
                            Cerr << "nk=" << dSk/coef << finl;
                            Cerr << "num_face=" << num_face << finl;
                            Cerr << "voisin 0/1 =" << zvf.face_voisins(num_face, 0)<< " " << zvf.face_voisins(num_face, 1) << finl;
                            //Process::exit();
                            resu(num_face,k)=dSk*(indic_ghost-indic);
                          }
                      }
                  }
                }
            }
          // Fin de la boucle for
        }
    }
}

void Navier_Stokes_FT_Disc::correct_at_exit_bad_gradient(DoubleTab& u0) const
{
  // Correction du gradient a la sortie (car celui-ci ne doit pas se baser sur la CL de pression):
  // On prefere mettre la valeur d'en face qu'une valeur abherante. -> comme cela, la divergence (plus tard) tendra vers 0)
  const Domaine_VF& zvf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const IntTab& elem_faces = zvf.elem_faces();
  const IntTab& face_voisins = zvf.face_voisins();
  const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis().valeur();
  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord=0; n_bord<zvf.nb_front_Cl(); n_bord++)
    {
      // pour chaque Condition Limite on regarde son type
      // Si face de Dirichlet ou de Symetrie on ne fait rien
      // Si face de Neumann on calcule la contribution au terme source
      const Cond_lim& la_cl = zcldis.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if ( (sub_type(Dirichlet,la_cl.valeur()))
           ||
           (sub_type(Neumann_sortie_libre,la_cl.valeur()))
           ||
           (sub_type(Dirichlet_homogene,la_cl.valeur()))
         )
        {
          // const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, zvf);
          // Cerr << que_suis_je() << "::calculer_delta_u_interface() correction de u0 gradient(phi) a la CL : " <<  la_cl.valeur() << finl;
          const int nb_faces=elem_faces.dimension(1) ;
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              const int elem = face_voisins(num_face,0)+face_voisins(num_face,1)+1;
              int idx_face_de_lelem = 0;
              for (idx_face_de_lelem=0; idx_face_de_lelem<nb_faces; idx_face_de_lelem++)
                {
                  if (elem_faces(elem,idx_face_de_lelem) == num_face)
                    break; // Face found
                }
              if (nb_faces==idx_face_de_lelem)
                {
                  Cerr << "Face is not found!! " << finl;
                  Process::exit();
                }
              const int num_face_den_face=elem_faces(elem,(idx_face_de_lelem+Objet_U::dimension)%nb_faces);
              // const int elem_voisin = face_voisins(num_face_den_face,0)+face_voisins(num_face_den_face,1)-elem;
              // const double d_elem = dist[elem];
              // const double d_voisin = dist[elem_voisin];
              // const double delta = zvdf.dist_face(num_face, num_face_den_face , zvdf.orientation(num_face));
              // Il faudrait lineariser a l'ordre 2 car cela sert aussi au calcul de laplacien_d.
              // On extrapole le gradient mais on ne peut pas faire mieux car on a un petit stencil a 2 points sur dist.
              for (int c=0; c< u0.line_size(); c++)
                u0(num_face,c) = u0(num_face_den_face,c);//*(1.+xxxxx);
            }
        }
    }
}

int get_num_face_den_face(const int num_face, const int elem, const IntTab& elem_faces )
{
  const int nb_faces=elem_faces.dimension(1) ;
  int idx_face_de_lelem = 0;
  for (idx_face_de_lelem=0; idx_face_de_lelem<nb_faces; idx_face_de_lelem++)
    {
      if (elem_faces(elem,idx_face_de_lelem) == num_face)
        break; // Face found
    }
  if (nb_faces==idx_face_de_lelem)
    {
      Cerr << "Face is not found!! " << finl;
      Process::exit();
    }
  const int num_face_den_face=elem_faces(elem,(idx_face_de_lelem+Objet_U::dimension)%nb_faces);
  return num_face_den_face;
}

/*! @brief Calcul du saut de vitesse a l'interface du au changement de phase
 *
 *   phase_pilote = -1: u-u0 = champ de vitesse de deplacement de l'interface
 *   phase_pilote = 0 : u+u0 = champ de vitesse de la phase 0
 *   phase_pilote = 1 : u+u0 = champ de vitesse de la phase 1
 *   ordre = 0 : pas de prise en compte de la correction en courbure
 *   ordre = 1 : prise en compte de la correction en courbure a l'ordre 1
 *   ordre = 2 : prise en compte de la correction en courbure a l'ordre 2
 *
 */
void Navier_Stokes_FT_Disc::calculer_delta_u_interface(Champ_base& champ_u0,
                                                       int phase_pilote,
                                                       int ordre)
{
  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_delta_u_interface", 0);
  //  statistiques().begin_count(count);
  DoubleTab& u0 = champ_u0.valeurs();
  const Fluide_Diphasique& fluide_dipha = fluide_diphasique();
  const Fluide_Incompressible& phase_0 = fluide_dipha.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide_dipha.fluide_phase(1);
  const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
  const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
  const double rho_0 = tab_rho_phase_0(0,0);
  const double rho_1 = tab_rho_phase_1(0,0);
  //const double delta_un_sur_rho = 1. / rho_1 - 1. / rho_0;
  const double un_sur_rho_0 = 1. / rho_0;
  const double un_sur_rho_1 = 1. / rho_1;

  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

  // Pour permettre de calculer mpoint, mais sans l'utiliser pour le deplacement de l'interface,
  // il suffit de ne pas mettre le mot cle " equation_temperature_mpoint temp" dans le jdd.
  const int nn = variables_internes().second_membre_projection.valeurs().dimension(0); // nombre d'elems
  DoubleTab mpoint;
  mpoint.resize(nn);
  mpoint=0.; //pour initialiser
  if (variables_internes().ref_equation_mpoint_.non_nul())
    {
      const DoubleTab& mp = variables_internes().ref_equation_mpoint_.valeur().get_mpoint();
      // Si inactif, on ne prend pas en compte sa contribution dans le calcul de delta_u:
      if (!variables_internes().mpoint_inactif)
        mpoint = mp;
    }
  if (variables_internes().ref_equation_mpoint_vap_.non_nul())
    {
      const DoubleTab& mpv = variables_internes().ref_equation_mpoint_vap_.valeur().get_mpoint();
      // Si inactif, on ne prend pas en compte sa contribution dans le calcul de delta_u:
      if (!variables_internes().mpointv_inactif)
        mpoint +=mpv;
    }

  // GB2023.10.10 : I don't understand why I did a distinction only on phase_pilote == 1 (should be the same when it's phase_pilote == 0
  //                for instance in the convection of temperature in the vapour...
  //        BESIDES, the "switch (phase_pilote)" makes no sense if only one phase_pilote is used.
  if ((variables_internes().new_mass_source_) && (phase_pilote != 1))
    {
      const DoubleTab& normale_elements = eq_transport.get_update_normale_interface().valeurs();
      const DoubleTab& interfacial_area = variables_internes().ai.valeur().valeurs();

      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const int nb_faces = face_voisins.dimension(0);
      const int dim = Objet_U::dimension;
      const DoubleTab& xp = domaine_vf.xp(); // centres de gravite des elements
      const DoubleTab& xv = domaine_vf.xv(); // centres de gravite des faces.
      u0 = 0.;
      if (u0.line_size() == dim) // vef
        {
          Cerr << "Using option new_mass_source is not possible yet in VEF. Contact us. " << finl;
          Process::exit();
          for (int face = 0; face < nb_faces; face++)
            for (int j = 0; j < dim; j++)
              {
                double x = 0.; // a coder...
                u0(face,j) *= x;
              }
        }
      else
        {
          const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
          const IntVect& orientation = zvdf.orientation();
          for (int face = 0; face < nb_faces; face++)
            {
              const int dir = orientation[face];
              const double surface=domaine_vf.face_surfaces(face);
              const int e1 = face_voisins(face, 0);
              const int e2 = face_voisins(face, 1);
              const double xf = xv(face,dir);
              if (Objet_U::bidim_axi && fabs(xf)<=DMINFLOAT && (dir==0))
                {
                  // We are on the symmetry axis => surface = 0.;
                  u0(face) = 0.; // there is no normal velocity at a symmetry axis
                  continue;
                }
              //double x = 0.;
              double xx = 0.;
              // Si on n'est pas au bord...
              if (e1 >= 0)
                {
                  const double nx = normale_elements(e1, dir);
                  //x = c*secmem2[e1]*nx;
                  const double ai= interfacial_area(e1);
                  // nx pointe vers le liquide (sortant de phase 0)
                  if ((fabs(ai)>DMINFLOAT) && (fabs(nx)>DMINFLOAT))
                    {
                      // distance positive on the vapour side (chi_0 = 1, ie indicatrice = 0)
                      const double d = (xf-xp(e1,dir)) * nx;
                      switch (phase_pilote)
                        {
                        case -1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et egal a la
                            // vitesse de deplacement de l'interface
                            const double un_sur_rho = (d >0.) ? un_sur_rho_0 : un_sur_rho_1;
                            xx = ai/surface * un_sur_rho * mpoint[e1]*nx;
                            //Cerr << "diff " << x << " " << xx << finl;
                            break;
                          }
                        case 0:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 0 dans la phase 0
                            const double p = (d >0.) ? 0. : (un_sur_rho_1 - un_sur_rho_0);
                            xx = ai/surface * p * mpoint[e1]*nx;
                            //Cerr << "face " << face << " " << xx << finl;
                            break;
                          }
                        case 1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 1 dans la phase 1
                            const double p = (d >0.) ? (un_sur_rho_0 - un_sur_rho_1) : 0. ;
                            xx = ai/surface * p * mpoint[e1]*nx;
                            //Cerr << "face " << face << " " << xx << finl;
                            break;
                          }
                        default:
                          Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                          Process::exit();
                        }
                    }
                }
              // We ADD contribution of e2 if not a boundary to xx
              if (e2 >= 0)
                {
                  const double nx = normale_elements(e2, dir);
                  //x += c*secmem2[e2]*normale_elements(e2, dir);
                  const double ai= interfacial_area(e2);
                  if ((fabs(ai)>DMINFLOAT) && (fabs(nx)>DMINFLOAT))
                    {
                      // distance positive on the vapour side (chi_0 = 1, ie indicatrice = 0)
                      const double d = (xf-xp(e2,dir)) * nx;
                      switch (phase_pilote)
                        {
                        case -1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et egal a la
                            // vitesse de deplacement de l'interface
                            const double un_sur_rho = (d >0.) ? un_sur_rho_0 : un_sur_rho_1;
                            xx += ai/surface * un_sur_rho * mpoint[e2]*nx;
                            //Cerr << "diff2 " << x << " " << xx << finl;
                            break;
                          }
                        case 0:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 0 dans la phase 0
                            const double p = (d >0.) ? 0. : (un_sur_rho_1 - un_sur_rho_0);
                            xx += ai/surface * p * mpoint[e2]*nx;
                            //Cerr << "face2 " << face << " " << xx << finl;
                            break;
                          }
                        case 1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 1 dans la phase 1
                            const double p = (d >0.) ? (un_sur_rho_0 - un_sur_rho_1) : 0. ;
                            xx += ai/surface * p * mpoint[e2]*nx;
                            //Cerr << "face2 " << face << " " << xx << finl;
                            break;
                          }
                        default:
                          Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                          Process::exit();
                        }
                    }

                }
              u0(face) = xx;
            }

          // GB 2023.10.10. On inclined interfaces, we realised that some configurations (interface topologies)
          //                can lead to using faces where the velocity jump delta_u has not been extended
          //                for the velocity interpolation at the markers. This leads to a misprediction of delta_u_i
          u0.echange_espace_virtuel();

          //  Pour parcourir les elements qui sont coupes par une facette "facette":
          const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
          const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
          const ArrOfInt& index_facette = intersections.index_facette();
          const Domaine_VF& zvf = ref_cast(Domaine_VF, domaine_dis().valeur());
          const IntTab& elem_faces = zvf.elem_faces();
          const int nb_faces_per_elem = elem_faces.line_size();
          for (int facette=0; facette < mesh.nb_facettes(); facette++)
            {
              int index=index_facette[facette];
              if (index >= 0)
                {
                  const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                  const int elem = data.numero_element_;
                  // elem is mixed (crossed by an interface)
                  // loop on its faces to find a neighbour:
                  for (int idx_face_de_lelem=0; idx_face_de_lelem<nb_faces_per_elem; idx_face_de_lelem++)
                    {
                      const int face_of_elem = elem_faces(elem,idx_face_de_lelem);
                      // The neighbour :
                      const int e1 = face_voisins(face_of_elem, 0)+face_voisins(face_of_elem, 1)-elem;
                      if ((e1 >= 0) && (fabs(interfacial_area(e1)) <DMINFLOAT))
                        {
                          // e1 exists (ie face is not a boundary) and is pure
                          // face_of_elem is between elem and e1
                          const int num_face_den_face = get_num_face_den_face(face_of_elem, e1, elem_faces);
                          if (fabs(u0(num_face_den_face))<DMINFLOAT)
                            {
                              // There is no velocity on the face in front of "face_of_e1" (which is adjacent to a mixed element)
                              // we should extrapolate the other one there.
                              // (ie we assume that this other face was zero because the 2nd neighbour (rang 2) is pure too.
                              u0(num_face_den_face) = u0(face_of_elem);
                            }
                        }
                    }
                }
            }
        }
      u0.echange_espace_virtuel();
      return;
    }

  // Distance a l'interface discretisee aux elements:
  const DoubleTab& dist = eq_transport.get_update_distance_interface().valeurs();
  DoubleTab phi = calculer_div_normale_interface().valeurs();
  {
    const int n = phi.dimension(0);
    for (int i = 0; i < n; i++)
      {
        double d = dist(i);
        double p = 0.;
        if (d >= -1e20)
          {
            const double div_n = phi(i);
            // Distance calculee pour cet element ?
            // Calcul de la fonction phi pour cet element :
            const double mp = mpoint[i];
            switch (ordre)
              {
              case 0:
                // Pas de prise en compte de la correction en courbure
                p = d * mp;
                break;
              case 1:
                //  Prise en compte de la correction en courbure a l'ordre 1
                p = d * (1. - 0.5 * div_n * d) * mp;
                break;
              case 2:
                //  Prise en compte de la correction en courbure a l'ordre 2
                p = d * (1. - 0.5 * div_n * d + div_n * div_n * d * d / 6.) * mp;
                break;
              default:
                Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface ordre" << finl;
                Process::exit();
              }
            switch (phase_pilote)
              {
              case -1:
                // Champ de vitesse tel que u + u0 soit continu et egal a la
                // vitesse de deplacement de l'interface
                if (d < 0)
                  p *= un_sur_rho_0;
                else
                  p *= un_sur_rho_1;
                break;
              case 0:
                // Champ de vitesse tel que u + u0 soit continu et
                // u+u0 = la vitesse de la phase 0 dans la phase 0
                if (d < 0)
                  p = 0.; // dans la phase 0
                else
                  p *= (un_sur_rho_0 - un_sur_rho_1); // GB BugFix 2020/10/09
                break;
              case 1:
                // Champ de vitesse tel que u + u0 soit continu et
                // u+u0 = la vitesse de la phase 1 dans la phase 1
                if (d < 0)
                  p *= (un_sur_rho_1 - un_sur_rho_0); // GB BugFix 2020/10/09
                else
                  p = 0.; // dans la phase 1
                break;
              default:
                Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                Process::exit();
              }
          }
        phi(i) = p;
      }
    phi.echange_espace_virtuel();
  }

  // Gradient de phi:
  if (champ_u0.que_suis_je() == "Champ_Face")
    {
      gradient.calculer(phi, u0);
      correct_at_exit_bad_gradient(u0);
    }
  else
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface\n"
           << " Non code pour " << champ_u0.que_suis_je() << finl;
      Process::exit();
    }

  // On annule la vitesse calculee pour les faces adjacentes a un element
  // invalide.
  {
    const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
    const IntTab&   face_voisins = domaine_vf.face_voisins();
    const int nb_faces = domaine_vf.nb_faces();;
    for (int i = 0; i < nb_faces; i++)
      {
        for (int j = 0; j < 2; j++)
          {
            const int elem = face_voisins(i, j);
            if (elem >= 0 && dist(elem) < -1e20)
              {
                u0(i) = 0.;
                break;
              }
          }
      }
  }
  u0.echange_espace_virtuel();
  solveur_masse.appliquer(u0);
}

double calculer_indicatrice_face_privilegie_pure(const DoubleTab& indicatrice, const IntTab& face_voisins, const int num_face)
{
  const int elem0 = face_voisins(num_face, 0);
  const int elem1 = face_voisins(num_face, 1);
  double indic_0 = (elem0 >= 0) ? indicatrice[elem0] : indicatrice[elem1];
  double indic_1 = (elem1 >= 0) ? indicatrice[elem1] : indicatrice[elem0];
  // Decentrage : si la face est adjacente a une face monophasique,
  //  prendre la phase pure pour toute la face:
  if (indic_0 == 0. || indic_0 == 1.)
    indic_1 = indic_0;
  if (indic_1 == 0. || indic_1 == 1.)
    indic_0 = indic_1;
  const double indic_face = (indic_0 + indic_1) * 0.5;
  return indic_face;
}

double calculer_indicatrice_face_based_on_ai(const DoubleTab& indicatrice,
                                             const DoubleTab& indicatrice_faces,
                                             const IntTab& face_voisins,
                                             const Domaine_VF& domaine_vf,
                                             const DoubleTab& interfacial_area,
                                             const DoubleTab& normale_elements,
                                             const int face,
                                             const int dim)
{
  double indic_face = 0.;
  int v;
  for (v = 0; v < 2; v++)
    {
      const int elem = face_voisins(face, v);
      if (elem >=0)
        {
          // If a neighbour is pure, we use that value at the face and stop further calculation.
          const double indic = indicatrice[elem]; // This is the value of chi_1 (ie =1 in phase 1!)
          //if (indic == 0. || indic == 1.)
          if (indic <=5e-3 || indic >= 1.-5e-3)
            {
              indic_face = 1-indic; // We want chi of phase_0
              break;
            }
          else
            {
              const double surface=domaine_vf.face_surfaces(face);
              //const DoubleVect& volumes = domaine_vf.volumes();
              const double ai= interfacial_area(elem); // nx pointe vers le liquide (sortant de phase 0)
              if (fabs(ai)>DMINFLOAT)
                {
                  double x = 0.;
                  if (dim>1) // dim==1 en VDF; sinon VEF
                    {
                      for (int j = 0; j < dim; j++)
                        {
                          const double dSj = domaine_vf.face_normales(face , j);
                          const double nx = normale_elements(elem, j);
                          // produit scalaire :
                          x +=  dSj*nx;
                          x *= ai/surface;
                          // Que/comment Choisir?
                          indic_face += 1-x;
                          // indic_face += x;
                          Cerr << "Never tested. To be verified. It should depend on a scalar product with the vect (xp-xv)" << finl;
                          Process::exit();
                        }
                    }
                  else
                    {
                      // En VDF, l'acces a orientation permet d'eviter le calcul du produit scalaire.
                      const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_vf);
                      const IntVect& orientation = zvdf.orientation();
                      const int dir = orientation[face];
                      const double nx = normale_elements(elem, dir);
                      // Assumes a cube, nx larger than diag means we can use the method rather safely
                      if (nx>0.707)
                        {
                          x = ai/surface*nx;
                          // On suppose que v0 est a gauche et v1 a droite!!!
                          if (v==0)
                            indic_face += 1-x; // This way, we build chi_0 because normale points towards chi_1
                          else
                            indic_face += x;
                        }
                      else
                        {
                          // L'interface croise probablement la face d'en face et la methode ne marche plus.
                          // on revient a la methode classique :
                          indic_face += (1. - indicatrice_faces[face]);
                        }
                    }
                }
              else
                {
                  Cerr <<" WTF, c'est impossible" << finl;
                  Process::exit();
                }
            }
        }
      else
        {
          // The only neighbour to the face :
          const int elem_voisin = face_voisins(face, 1-v); // The other one is accessed by 1-v
          const double indic = indicatrice[elem_voisin]; // This is the value of chi_1 (ie =1 in phase 1!)
          indic_face = 1-indic; // We want chi of phase_0
          break; // c'est important pour le if d'apres.
        }
    }
  if (v==2)
    // On n'a pas touche le break, on est donc passe 2 fois. donc :
    indic_face*=0.5;

  // assert((indic_face >=0) && (indic_face<=1.));
  // ca arrive des petits derapages..
  if (indic_face <0)
    indic_face=0.;
  if (indic_face >1.)
    indic_face=1.;
  return indic_face;
}

void correct_indicatrice_face_bord(const int num_face,
                                   const Maillage_FT_Disc& maillage,
                                   const Domaine_VF& zvf,
                                   const IntTab& face_voisins,
                                   const DoubleTab& indicatrice,
                                   const bool privilegie_pure,
                                   double& indic_face)
{
  // Correction de l'indicatrice face :
  const int n0 = face_voisins(num_face, 0);
  const int n1 = face_voisins(num_face, 1);
  if ((n0==-1) or (n1==-1))
    {
      // On a boundary face
      const int outward_normal = (n0 == -1) ? -1 : 1 ;
      const int elem = n0+n1+1;
      const double indic_ghost = compute_indic_ghost(elem, num_face, indicatrice(elem), outward_normal, zvf,
                                                     maillage);
      indic_face = indic_ghost;
      if ((privilegie_pure) && (est_egal( indicatrice(elem), 1.) || est_egal( indicatrice(elem), 0.)))
        {
          indic_face =  indicatrice(elem);
        }
    }
}

// Calcul de l'integrale de dI_dt sur chaque element du maillage.
// Le tableau dI_dt doit avoir la bonne structure. L'espace virtuel est
// mis a jour. La method n'est plus const a cause des options
// INTERP_MODIFIEE et AI_BASED qui recalculent indicatrice_faces.
void Navier_Stokes_FT_Disc::calculer_dI_dt(DoubleVect& dI_dt) //const
{
  const double rho_0 = fluide_diphasique().fluide_phase(0).masse_volumique().valeurs()(0,0);
  const double rho_1 = fluide_diphasique().fluide_phase(1).masse_volumique().valeurs()(0,0);
  const double delta_rho = rho_0-rho_1;

  double rho_0_sur_delta_rho = 0.;
  if (delta_rho != 0)
    rho_0_sur_delta_rho = rho_0 / delta_rho;

  const DoubleTab& tab_vitesse = inconnue().valeurs();
  const IntTab& face_voisins = domaine_dis().valeur().face_voisins();

  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const DoubleTab& indicatrice = eq_transport.inconnue().valeur().valeurs();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const Domaine_VF& domVF = ref_cast(Domaine_VF, domaine_dis().valeur());
  //const IntVect& orientation = ref_cast(Domaine_VF, domaine_dis().valeur()).orientation();
  //  const DoubleTab& indicatrice = variables_internes().ref_eq_interf_proprietes_fluide.valeur().inconnue().valeur().valeurs();
  DoubleTab tmp(tab_vitesse); // copie du tableau des vitesses de ns
  const int dim = tab_vitesse.line_size();
  const int n = tab_vitesse.dimension(0);

  // On cree un tableau avec la meme structure que la pression
  DoubleTab resu;
  resu.copy(variables_internes().second_membre_projection.valeurs(), Array_base::NOCOPY_NOINIT);
  resu=0.;

  //On utilise un operateur de divergence temporaire et pas celui porte par l equation
  //pour ne pas modifier les flux_bords_ rempli au cours de Navier_Stokes_std::mettre_a_jour
  Operateur_Div div_tmp;
  div_tmp.associer_eqn(*this);
  div_tmp.typer();
  div_tmp.l_op_base().associer_eqn(*this);
  div_tmp->completer();

  if ((variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT)
      || (variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT)
      || (variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT))
    {
      // Avec changement de phase, on veut reconstruire u_vap (ie phase 0)
      // Prise en compte du terme source div(u) du changement de phase
      if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
        {
          calculer_delta_u_interface(variables_internes().vitesse_jump0_, 0,  variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);
          // u+u0 = champ de vitesse de la phase 0
          tmp += variables_internes().vitesse_jump0_.valeurs();
        }
    }
  if ((variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT)
      || (variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT)
      || (variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT))
    {
      // Reconstruction d'un champ de vitesse interfaciale (ie phase 1)
      // Prise en compte du terme source div(u) du changement de phase
      if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
        {
          calculer_delta_u_interface(variables_internes().delta_u_interface, -1,  variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);
          // u-dui = champ de vitesse d'interface
          tmp -= variables_internes().delta_u_interface.valeurs();

          // Question: il y a un assert_espace_virtuel_vect dans divergence.calculer,
          //  mais l'operateur n'a normalement pas besoin de l'espace virtuel !
          //  La ligne suivante devrait pouvoir etre retiree:
          tmp.echange_espace_virtuel();

          div_tmp->calculer(tmp,resu);
          for (int i = 0; i < resu.size_array(); i++)
            resu[i] *= indicatrice[i];
        }
    }

  const bool ghost_correction = ( variables_internes_->OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC);
  switch(variables_internes_->type_interpol_indic_pour_dI_dt_)
    {
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD:
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT:
      {
        for (int i = 0; i < n; i++)
          {
            double indic_face = calculer_indicatrice_face_privilegie_pure(indicatrice, face_voisins, i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice,
                                            true /* privilegie_pure */, indic_face);
            const double x = rho_0_sur_delta_rho - indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i,j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT:
      {
        for (int i = 0; i < n; i++)
          {
            double indic_face = calculer_indicatrice_face_privilegie_pure(indicatrice, face_voisins, i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice,
                                            true /* privilegie_pure */, indic_face);
            const double x = -indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i,j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE:
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT:
      {
        const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
        for (int i = 0; i < n; i++)
          {
            double indic_face = indicatrice_faces(i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice, false, indic_face);
            const double x = rho_0_sur_delta_rho - indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i,j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT:
      {
        const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
        for (int i = 0; i < n; i++)
          {
            double indic_face = indicatrice_faces(i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice, false, indic_face);
            const double x = -indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i,j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT:
      {
        const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
        const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
        if (Process::je_suis_maitre())
          Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area"
               << " and on the normal to the interface." << finl;

        const DoubleTab& normale_elements = eq_transport.get_update_normale_interface().valeurs();
        const DoubleTab& interfacial_area = variables_internes().ai.valeur().valeurs();

#if NS_VERBOSE
        const DoubleTab& xp = domaine_vf.xp(); // centres de gravite des elements
        const DoubleTab& xv = domaine_vf.xv(); // centres de gravite des faces.
#endif
        // On fait la moyenne des 2 valeurs calculees sur les voisins
        // ATTENTION, ici on veut la valeur de chiv (cad chi_0) a la face.
        for (int face = 0; face < n; face++)
          {
            double indic_face = calculer_indicatrice_face_based_on_ai(indicatrice,
                                                                      indicatrice_faces,
                                                                      face_voisins,
                                                                      domaine_vf,
                                                                      interfacial_area,
                                                                      normale_elements,
                                                                      face,
                                                                      dim);
            if (ghost_correction)
              correct_indicatrice_face_bord(face, maillage, domVF, face_voisins, indicatrice, false, indic_face);

#if NS_VERBOSE
            const double val = 1.-indicatrice_faces[face]; // indicatrice_faces=chi_1  whereas indic_face=chi_0 !!! Hence, the "1.-"
            if (fabs(indic_face-val)>DMINFLOAT)
              {
                const int elem0 = face_voisins(face, 0);
                const int elem1 = face_voisins(face, 1);
                double indic_0 = (elem0 >= 0) ? indicatrice[elem0] : indicatrice[elem1];
                double indic_1 = (elem1 >= 0) ? indicatrice[elem1] : indicatrice[elem0];
                if (elem0>=0)
                  Cerr << "xp0["<< elem0<< "]: " << xp(elem0, 0) << " "  << xp(elem0, 1) << " "  << finl;
                else
                  Cerr << "xp0: bord!" <<finl;
                if (elem1>=0)
                  Cerr << "xp1["<< elem1<< "]: " << xp(elem1, 0) << " "  << xp(elem1, 1) << " "  << finl;
                else
                  Cerr << "xp1: bord!" <<finl;

                Cerr << "xv: " << xv(face, 0) << " "  << xv(face, 1) << " "  << finl;
                Cerr << "voisins (ou ghost): " << indic_0 << " " << indic_1 << finl;
                Cerr << "GB whats up?face="<< face <<" indic / val / diff " << indic_face << " " << val << " " << indic_face-val << finl;
              }
#endif
            // chi_v * u_v_ext
            for (int j = 0; j < dim; j++)
              tmp(face,j) *= indic_face;
          }
        break;
      }
    default:
      Cerr << " Navier_Stokes_FT_Disc::calculer_dI_dt \n"
           << " unknown case?" << finl;
      Process::exit();
    }

  if (variables_internes_->type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT)
    tmp *=-1;

  // Question: il y a un assert_espace_virtuel_vect dans divergence.calculer,
  //  mais l'operateur n'a normalement pas besoin de l'espace virtuel !
  //  La ligne suivante devrait pouvoir etre retiree:
  tmp.echange_espace_virtuel();

  div_tmp->ajouter(tmp,resu);

  // Correction des flux bords :
  // resu = int_V div(tmp) dv    avec:   tmp = (rho_0_sur_delta_rho - indic_face)*vitesse_ns
  // 		(a voir qui est indic_face -> celle de phase1?)
  // Les flux aux bords dans les mailles diphasiques avec sortie libre sont mal calcules si l'indic_phase l'est mal.
  // (car il est difficile de savoir si l'indic_face est en faite encore pure et qu'il ne sort que d'une phase
  //  dans certaines mailles mixtes).
  // Or, les flux_bord ont deja affecte resu par l'operation :
  //   resu += flux // ...
  // Plusieurs options sont testees :
  //   - On fait rien (NO_CORRECTION)
  //   - On fait la correction de l'indic_face avant (avant le calcul de tmp qui est en entree de la div) (CORRECTION_GHOST_INDIC)
  //   - On suppose que globalement, les cellules touchant le bord sortie sont
  // a flux de masse phasique nulle (ce qui rentre sort). (ZERO_NET_FLUX_ON_MIXED_CELLS)
  // Cela signifie que le flux sur la face de bord est l'oppose des flux internes.
  // Concretement, il suffit de faire "resu(elem) =0." dans les mailles concernees.
  //   - On recalcule le flux sur la face de bord des mailles mixtes et on retranche celui qui avait ete ajoute (ZERO_OUT_FLUX_ON_MIXED_CELLS)
  // Cela ne laisse pas sortir la vapeur, ce n'est donc pas bon.
  switch(variables_internes_->OutletCorrection_pour_dI_dt_)
    {
    case Navier_Stokes_FT_Disc_interne::NO_CORRECTION:
    case Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC:
      {
        break;
      }
    case Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS:
    case Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS:
      {
        const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
        const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis().valeur();
        const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
        const DoubleVect& face_surfaces = zvdf.face_surfaces();
        // Boucle sur les bords pour traiter les conditions aux limites
        int ndeb, nfin;
        for (int n_bord=0; n_bord<zvdf.nb_front_Cl(); n_bord++)
          {
            // pour chaque Condition Limite on regarde son type
            // Si face de Dirichlet ou de Symetrie on ne fait rien
            // Si face de Neumann on calcule la contribution au terme source
            const Cond_lim& la_cl = zclvdf.les_conditions_limites(n_bord);
            //Cerr << que_suis_je() << "::calculer_dI_dt() correction du dIdt a la CL : " <<  la_cl.valeur() << finl;
            const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
            ndeb = le_bord.num_premiere_face();
            nfin = ndeb + le_bord.nb_faces();

            // Correction sortie libre :
            // On ne sait pas bien calculer la correction de volume liee a dIdt a la sortie libre.
            // On prefere donc l'annuler dans cet element:
            if ( (sub_type(Dirichlet,la_cl.valeur()))
                 ||
                 (sub_type(Neumann_sortie_libre,la_cl.valeur()))
                 ||
                 (sub_type(Dirichlet_homogene,la_cl.valeur()))
               )
              {
                for (int num_face=ndeb; num_face<nfin; num_face++)
                  {
                    const int n0 = face_voisins(num_face,0);
                    const int n1 = face_voisins(num_face,1);
                    const int elem = n0+n1+1;
                    const double indic = indicatrice(elem);
                    if (indic*(1-indic)> 1e-6) // In a mixed cell!
                      {
                        const double coef = face_surfaces(num_face);//*porosite_surf(num_face);
                        if ( variables_internes_->OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS)
                          resu(elem) =0.;
                        else if ( variables_internes_->OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS)
                          {
                            for (int j = 0; j < dim; j++)
                              {
                                const int outward_normal = (n0 == -1) ? -1 : 1 ;
                                double flux_bord_calcule_par_operateur_div = 0.;
                                flux_bord_calcule_par_operateur_div =outward_normal*coef*tmp(num_face,j);
                                resu(elem) -= flux_bord_calcule_par_operateur_div; // On RETRANCHE le flux qui avait ete pris a l'etape de calcul du div
                              }
                          }
                        else
                          {
                            Process::exit();
                          }

                      }
                  }
              }
          }
        // Fin de la boucle for
        break;
      }
    default:
      Cerr << "unexpected" <<finl;
      Process::exit();
    }

  // Extraction des valeurs
  const int nb_elem = domaine_dis().valeur().nb_elem();
  assert(nb_elem == dI_dt.size());

  // Simple copie
  dI_dt.inject_array(resu, nb_elem);

  // L'integrale de div sur l'element est dimension * la valeur aux elements
  //  renvoyee par l'operateur.
  // Les valeurs aux elements sont au debut du tableau resu.
  if(tab_vitesse.line_size() > 1) // i.e. VEF
    dI_dt *= dimension;

#if NS_VERBOSE
  {
    Cerr << "[BEFORE-PCH] Locally, the maximum of dI_dt is : " << dI_dt.mp_max_abs_vect() << finl;
    const double temps = schema_temps().temps_courant();
    double sum = 0.;
    for (int i=0; i< nb_elem; i++)
      sum += dI_dt[i];
    Cerr << "[BEFORE-PCH] " << temps << " The sum is : " << sum << " [not valid in //]" << finl;
  }
#endif

  switch(variables_internes_->type_interpol_indic_pour_dI_dt_)
    {
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED:
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT:
      {
        // dI_dt contains div(chi_v * u_v^ext) * vol_cell (sign to be checked!!)
        // Integrated, this is equal to 0 or to the integral of chi_f entering or leaving the domain through boundaries...
        // It is lacking the phase-change contribution
#define CHECK 0
#if CHECK
        {
          Cerr << "Checking if secmem2 is correctly filled in that case. Check it if you need it..." <<finl;
          if (0) Process::exit();
          DoubleTab& secmem2 = variables_internes().second_membre_projection_jump_.valeurs();

          //if (variables_internes().ref_equation_mpoint_.non_nul())
          //  variables_internes().ref_equation_mpoint_.valeur().calculer_mpoint(variables_internes().mpoint.valeur());
          //if (variables_internes().ref_equation_mpoint_vap_.non_nul())
          //  variables_internes().ref_equation_mpoint_vap_.valeur().calculer_mpoint(variables_internes().mpoint_vap.valeur());

          const Fluide_Diphasique& fluide = fluide_diphasique();
          const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
          const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
          const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
          const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
          const double rho_phase_0 = tab_rho_phase_0(0,0);
          const double rho_phase_1 = tab_rho_phase_1(0,0);
          const double jump_inv_rho = 1./rho_phase_1 - 1./rho_phase_0;

          const DoubleTab& interfacial_area = variables_internes().ai.valeur().valeurs();
          // Pas une ref, mais un tableau de travail local dans lequel on peut ajouter mointv
          // DoubleTab mpoint = variables_internes().mpoint.valeur().valeurs();
          DoubleTab mpoint = variables_internes().ref_equation_mpoint_.valeur().get_mpoint();
          if (variables_internes().ref_equation_mpoint_vap_.non_nul())
            {
              //const DoubleTab& mpointv = variables_internes().mpoint_vap.valeur().valeurs();
              const DoubleTab& mpointv = variables_internes().ref_equation_mpoint_vap_.valeur().get_mpoint();
              for (int elem = 0; elem < nb_elem; elem++)
                mpoint[elem] += mpointv[elem];
            }
          for (int elem = 0; elem < nb_elem; elem++)
            {
              double diff = (secmem2[elem] - jump_inv_rho*interfacial_area[elem]*mpoint[elem]);
              if (fabs(diff)>DMINFLOAT)
                {
                  Cerr << "Problem, secmem2 is not filled properly in that case? "
                  << "elem= " << elem << " diff=" << diff <<finl;
                  Process::exit();
                }
            }
        }
#endif
        if ((variables_internes().ref_equation_mpoint_.non_nul())&& !variables_internes().mpoint_inactif)
          {
            // Is it necessary to recompute them? no
            //variables_internes().ref_equation_mpoint_.valeur().calculer_mpoint(variables_internes().mpoint.valeur());
            const double un_sur_rho_0 =  1./rho_0;

            // We can't use secmem2 because we can't undo jump_inv_rho as it may be 0. !!!
            // DoubleTab& secmem2 = variables_internes().second_membre_projection_jump_.valeurs();
            const DoubleTab& interfacial_area = variables_internes().ai.valeur().valeurs();
            // const DoubleVect& volumes = domaine_vf.volumes();
            // Pas une ref, mais un tableau de travail local dans lequel on peut ajouter mointv
            DoubleTab mpoint = variables_internes().ref_equation_mpoint_.valeur().get_mpoint();
            if (variables_internes().ref_equation_mpoint_vap_.non_nul())
              {
                // Is it necessary to recompute them? no
                //variables_internes().ref_equation_mpoint_vap_.valeur().calculer_mpoint(variables_internes().mpoint_vap.valeur());
                const DoubleTab& mpointv = variables_internes().ref_equation_mpoint_vap_.valeur().get_mpoint();
                for (int elem = 0; elem < nb_elem; elem++)
                  mpoint[elem] += mpointv[elem];
              }
#if TCL_MODEL
            // At this point in the algorithm, mpoint has not been augmented by the CL contribution yet, even though
            // the TCL contribution (micro and meso) are already computed and associated to elements.
            // It is done so because just earlier in this method, we call calculer_delta_u_interface
            // to compute a continuous extension of velocity (to displace the interface markers).

            // Adding the TCL contribution to the **local** DoubleTab mpoint
            // (it is thus temporary for the moment and will be done on the real table later at the second call to
            // Triple_Line_Model_FT_Disc::corriger_mpoint)
            // The difference if we correct it here directly is subtile. It is probably fine for the extension of
            // the interface velocity, but the extension for the convective velocity in the temperature has not been done yet.
            // It may cause trouble in the convection then, but the GFM method kind of erase those cells.
            if (probleme_ft().tcl().is_activated())
              {
                Cerr << "[TCL] Contact line model activated in volume correction" << finl;
                probleme_ft().tcl().corriger_mpoint(mpoint);
              }
#endif
            for (int elem=0; elem< nb_elem; elem++)
              {
                // By convention, mpoint is positive in condensation. Hence, mpoint >0 is responsible for dIv_dt < 0  => a minus sign!
                //                But, \nabla \chi_v = -n_v \delta_i. => another minus sign!
                //                ==> Consequently, it's a "+"
                // Besides, ai = \int_cell delta^i dv => It's homogeneous to the integral, there's no need for an additional "*volumes[elem]"
                //                                       It can be directly summed to the divergence computed before.
                const double x = mpoint[elem]*interfacial_area[elem]*un_sur_rho_0;
                dI_dt[elem] += x;
              }
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT:
      {
        if (probleme_ft().tcl().is_activated())
          {
            const double un_sur_rho_0 =  1./rho_0;
            const DoubleTab& interfacial_area = variables_internes().ai.valeur().valeurs();
            Cerr << "[TCL] Contact line model activated in volume correction" << finl;
            ArrOfInt& tcl_elems =  probleme_ft().tcl().elems();
            ArrOfDouble& tcl_mp =  probleme_ft().tcl().mp();
            // We accounted for the contribution that is in Ui but not the contribution from TCL yet:
            for (int idx=0; idx< tcl_elems.size_array(); idx++)
              {
                const int elem = tcl_elems[idx];
                // By convention, mpoint is positive in condensation. Hence, mpoint >0 is responsible for dIv_dt < 0  => a minus sign!
                //                But, \nabla \chi_v = -n_v \delta_i. => another minus sign!
                //                ==> Consequently, it's a "+"
                // Besides, ai = \int_cell delta^i dv => It's homogeneous to the integral, there's no need for an additional "*volumes[elem]"
                //                                       It can be directly summed to the divergence computed before.
                const double x = tcl_mp[idx]*interfacial_area[elem]*un_sur_rho_0;
                dI_dt[elem] += x;
              }
          }
        break;
      }
    default:
      break;
    }

#if NS_VERBOSE
  {
    Cerr << "[AFTER-PCH] Locally, the maximum of dI_dt is : " << dI_dt.mp_max_abs_vect() << finl;
    const double temps = schema_temps().temps_courant();
    double sum = 0.;
    for (int i=0; i< nb_elem; i++)
      sum += dI_dt[i];
    Cerr << "[AFTER-PCH] " << temps << " The sum is : " << sum << " [not valid in //]" << finl;
  }
#endif
  dI_dt.echange_espace_virtuel();
}

// Description:
//  Compute in one Eulerian cell, the average of Front properties in it :
//  normale, bary_facettes_dans_elem and surface_tot are area-weighted averages in the given cell.
//  Warning : normale is not a UNIT vector!
void compute_normale_barycenter_area_in_cell(const int elem,
                                             const  Maillage_FT_Disc& mesh,
                                             Vecteur3& normale,
                                             Vecteur3& bary_facettes_dans_elem,
                                             double& surface_tot)
{
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normale_facettes = mesh.get_update_normale_facettes();

  surface_tot = 0.;
  normale = 0.;
  bary_facettes_dans_elem = 0.;
  // Get the begining index defining the position in the list index_elem that contains information
  // relative to the current element :
  int index=intersections.index_elem()[elem];
  if (index < 0)
    return; // No facette in this element.

  // Loop over the facettes crossing the element
  while (index >= 0)
    {
      // Accessing the structure containing all the relevant information for facette number fa7
      // Beware, fraction_surface_intersection_ gives the fraction of the facette that is within the current element.
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
      const int fa7 = data.numero_facette_;
      const double surface_facette = surface_facettes[fa7];
      const double surf = data.fraction_surface_intersection_ * surface_facette;
      // We compute the (real) coordinates of barycenter of the fraction of facette within the elem
      //    from the Barycentric Coordinates stored in data.barycentre_
      //    (see the web for more information on Barycentric or Areal coordinates)
      //    coord_barycentre_fraction will contain real coordinates (x,y,z) of the barycenter
      Vecteur3 coord_barycentre_fraction(0., 0., 0.) ;
      for (int dir = 0; dir< 3; dir++)
        {
          const double nx = normale_facettes(fa7,dir);
          normale[dir] += nx * surf;
        }
      for (int isom = 0; isom< 3; isom++)
        {
          const int num_som = facettes(fa7, isom); // numero du sommet dans le tableau sommets
          const double bary_som = data.barycentre_[isom];
          for (int dir = 0; dir< 3; dir++)
            coord_barycentre_fraction[dir] += bary_som * sommets(num_som,dir);

        }
      coord_barycentre_fraction *= surf;
      surface_tot +=surf;
      bary_facettes_dans_elem += coord_barycentre_fraction; // This is done for all the 3 components.

      index = data.index_facette_suivante_;
    }

  if (surface_tot > 0.)
    {
      normale *= 1./surface_tot;
      bary_facettes_dans_elem *= 1./surface_tot;
    }
  else
    {
      normale = 0.;
      Cerr << " Error in compute_normale_barycenter_area_in_cell (Navier_Stokes_FT_Disc.cpp)." << finl;
      Cerr << "The element " << elem << " only contains facettes of surface=0, so that surface_totale is zero!" << finl;
      Cerr << "What a mess for the barycentre? ..." << finl;
      //    assert(0);
      Process::exit();
      bary_facettes_dans_elem = 0.;
    }
#if NS_VERBOSE
  const double norm =  normale[0]*normale[0] +  normale[1]*normale[1] +  normale[2]*normale[2];
  if (norm<0.9)
    {
      Cerr << " In Navier_Stokes_FT_Disc.cpp compute_normale_barycenter_area_in_cell." << finl;
      Cerr << "Small normal : " << count << " facettes dans l'element " << elem
           << ". surface_tot = "<< surface_tot << "Norm**2 = " << norm << finl;
    }
#endif

}

void Navier_Stokes_FT_Disc::compute_boussinesq_additional_gravity(
  const Convection_Diffusion_Temperature_FT_Disc& eq,
  const Fluide_Diphasique& fluide_dipha,
  const IntTab& face_voisins,
  const DoubleVect& volumes_entrelaces,
  const IntVect& orientation,
  const DoubleTab& indicatrice,
  const ArrOfDouble& g, // Vect3
  DoubleTab& gravite_face) const
{
  const int phase_eq = eq.get_phase();
  const DoubleTab& temperature_eq = eq.inconnue().valeur().valeurs();
  const Fluide_Incompressible& fluide_phase_eq = fluide_dipha.fluide_phase(phase_eq);
  const DoubleTab& tab_beta_th_phase_eq = fluide_phase_eq.beta_t().valeurs();
  const double beta_th_phase_eq = tab_beta_th_phase_eq(0,0);

  for (int face=0; face<gravite_face.dimension(0); face++)
    {
      const int elem0 = face_voisins(face, 0);
      const int elem1 = face_voisins(face, 1);
      double coef = 0.;
      // On suppose la ref T0 egale a Tsat
      //
      // Pour les mailles monophasiques, on peut faire simplement l'hypothese que T = chi_k T_k
      // Dans les mailles diphasiques, T = chi_k T_k est une hypothese discutable.
      // Pour les mailles diphasiques, on pourrait envisager une reconstruction plus precise de la
      // temperature monofluide a partir du gradient (cad de mpoint).
      // Neglected in first approximation. we simply compute T = chi_k T_k
      if (elem0 >= 0)
        {
          double chi = (2*phase_eq-1)*indicatrice[elem0]+1-phase_eq;
          double T_eq = temperature_eq[elem0];
          coef = chi*T_eq;
        }
      if (elem1 >= 0)
        {
          double chi = (2*phase_eq-1)*indicatrice[elem1]+1-phase_eq;
          double T_eq = temperature_eq[elem1];
          coef += chi*T_eq;
        }
      if (elem0 >= 0 && elem1 >= 0) // Not a boundary of the domain ?
        coef *= 0.5;
      gravite_face(face)-=volumes_entrelaces(face)*g(orientation[face])*coef*beta_th_phase_eq;
    }
}


// EB
/*! @brief Cette fonction calcule les composantes du tenseur gradU aux points P1 pour chaque facette du maillage lagrangien.
 *  En entree: indicatrice_face, valeurs_champs (tableau de vitesse), coord : coordonnes xyz des points P1 pour chaque facette lagrangienne, resu : composantes
 *  du tenseur interpole pour chaque facette lagrangienne.
 *  Boucle sur les facetttes lagrangiennes :
 *  	Dans un premier temps, identification des 8 faces voisines de normale x (4 en 2D) les plus proches.
 *  	Ensuite, pour chaque face, calcul du tenseur des contraintes
 *  	Pour chaque composante du tenseur, interpolation trilineaire au point de coordonnees (coord(fa7,0), coord(fa7,1), cooord(fa7,2)).
 *
 */
int Navier_Stokes_FT_Disc::trilinear_interpolation_gradU_face(const DoubleTab& indicatrice_face, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu)
{
  // On identifie l'element dans lequel appartient le point de coordonnees coord
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  const Domaine& domaine = domaine_vdf.domaine();
  //const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();


  IntVect faces_elem_interp(2*dimension);
  int nb_voisins=8;

  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!maillage.facette_virtuelle(fa7))
        {
          const int elem=domaine.chercher_elements(coord(fa7,0), coord(fa7,1), coord(fa7,2));
          for (int dim=0; dim<dimension; dim++)
            {
              faces_elem_interp(dim)=domaine_vdf.elem_faces_pour_interp(elem,dim);
              faces_elem_interp(dimension+dim)=domaine_vdf.elem_faces_pour_interp(elem,dimension+dim);
              if (faces_elem_interp(dim)<0 || faces_elem_interp(dimension+dim)<0) return 0;
            }

          DoubleTab coord_face(2*dimension,dimension);
          for (int i=0; i<2*dimension; i++)
            {
              for (int dim=0; dim<dimension; dim++)
                {
                  coord_face(i,dim)=domaine_vdf.xv(faces_elem_interp(i),dim);
                }
            }

          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect faces_voisines(nb_voisins); // 8 elements voisins au point de coordonnees coord. L'element elem est inclu dedans.
          chercher_faces_voisines(coord_elem_interp,faces_voisines,0);
          for (int i=0; i<nb_voisins; i++)
            {
              if (faces_voisines(i)<0) return 0;
            }
          DoubleTab gradUx(nb_voisins, dimension, dimension); // le x signifie que l'on calcule le tenseur en une facette dont la composante de la vitesse est en x
          for (int face=0; face<8; face++) //on calcule le tenseur gradient de la vitesse pour chaque face voisine
            {
              // 														SCHEMA EN 2D
              //  										 ---- --- --- --- --- --- --- --- --- --- -
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //					 										 3
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  										 ---- --- --- --- -7- -8- --- --- --- --- -
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //											  	  	  	 1 	 x 	 2
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  									     ---- --- --- --- -5- -6- --- --- --- --- -	  y
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |	  ^
              //					 										 4						  |
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |	  |
              //  										 ---- --- --- --- --- --- --- --- --- --- -   o----->x
              //																					  z
              // Evaluation du gradient sur les faces de normale x (choix arbitraire mais peut tuer une certaine symetrie)
              // Les termes + designent les faces juxtaposees selon +z
              // Les termes - designent les faces juxtaposees selon -z
              // Les termes 57 et 86 pour w designent les faces de normale z secantes aux faces 5,7 et 8,6
              // 		   --					   													   	   										   													   --
              //		   |  u_2-u_1	 		   				  	  u_3-u_4				 						 		  u_x+ - u_x- 	 								   							|
              //		   | ----------								------------	        								---------------																|
              //		   | (x_2-x1)			 					(y_3-y_4)												  (z_x+ - z_x-)	  															|
              //		   |												    																														|
              //  grad U = | 1    v_7-v_8	 v_5 - v_6    			1	v_5-v_7	  v_6-v_8								  1	  v_7+ - v_7- 	     v_8+ - v_8-	   v_5+ - v_5-	      v_6+ - v_6-		|
              // 	  	   | - ( --------- + --------- ) 			- ( ------ + -------)								  - ( --------------- + --------------- + --------------- + --------------- )	|
              //		   | 2    x_7-x_8	 x_5 - x_6   			2	y_5-y_7	  y_6-y_8								  4	  (z_7+ - z_7-)	    (z_8+ - z_8-)	   (z_5+ - z_8-)	  (z_6+ - z_6-)	|
              //		   |																																											|
              // 	  	   | 1	 w57+ - w86+	 w57- - w86-		1	w27-w29	 	w28-w30	    w23-w25      w24-w26	  1   w_57+ - w_57-   w_86+ - w_86-												|
              //		   | - ( ----------- +	 ----------- )		- (	-------- + 	-------- + --------- +  --------- )   - ( ------------ + --------------)	    									|
              //  		   | 2	 x57+ - x86+	 x57- - x86-		4	y27-y29	 	y28-y30	    y23-y25		 y24-y26 	  2   z_57+ - z_57-   z_86+ - z_86-												|
              //		   --					   		   									              									   	   													   --
              //
              // Pour chaque face, il faut donc 30 faces voisines : 1-8, x+, x-, 5+, 5-, 6+, 6-, 7+, 7-, 8+, 8-, 57+, 57-, 86+, 86-
              // Pour plus de facilites, on numerote :  x+:9, x-:10, 5+:11, 5-:12, 6+:13, 6-:14, 7+:15, 7-:16, 8+:17, 8-:18, 57+:19, 57-:20, 86+:21, 86-:22
              // 23-30 : les faces de normales z autour de x. Numerotation de bas en haut, de devant a derriere, de gauche a droite
              //Cerr << "Debut de la definition des elements voisins " << finl;

              int nb_voisinsx=30;
              IntVect voisinsx(nb_voisinsx);
              int num_facex=faces_voisines(face);
              int elem_gauche=domaine_vdf.face_voisins_pour_interp(num_facex, 0);
              int elem_droite=domaine_vdf.face_voisins_pour_interp(num_facex, 1);
              int elem_haut_gauche=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_gauche,2+dimension), 1);
              int elem_haut_droite=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_droite,2+dimension), 1);
              int elem_bas_gauche=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_gauche,2), 0);
              int elem_bas_droite=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_droite,2), 0);
              int elem_avant_gauche=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_gauche,1+dimension),1);
              int elem_avant_droite=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_droite,1+dimension),1);
              int elem_arriere_gauche=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_gauche,1),0);
              int elem_arriere_droite=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_droite,1),0);

              voisinsx(0)  = domaine_vdf.elem_faces_pour_interp(elem_gauche,0);
              voisinsx(1)  = domaine_vdf.elem_faces_pour_interp(elem_droite,0+dimension);
              voisinsx(2)  = domaine_vdf.elem_faces_pour_interp(elem_avant_gauche,0+dimension);
              voisinsx(3)  = domaine_vdf.elem_faces_pour_interp(elem_arriere_gauche,0+dimension);
              voisinsx(4)  = domaine_vdf.elem_faces_pour_interp(elem_gauche,1);
              voisinsx(5)  = domaine_vdf.elem_faces_pour_interp(elem_droite,1);
              voisinsx(6)  = domaine_vdf.elem_faces_pour_interp(elem_gauche,1+dimension);
              voisinsx(7)  = domaine_vdf.elem_faces_pour_interp(elem_droite,1+dimension);
              voisinsx(8)  = domaine_vdf.elem_faces_pour_interp(elem_haut_gauche,0+dimension);
              voisinsx(9)  = domaine_vdf.elem_faces_pour_interp(elem_bas_gauche,0+dimension);
              voisinsx(10) = domaine_vdf.elem_faces_pour_interp(elem_haut_gauche,1);
              voisinsx(11) = domaine_vdf.elem_faces_pour_interp(elem_bas_gauche,1);
              voisinsx(12) = domaine_vdf.elem_faces_pour_interp(elem_haut_droite,1);
              voisinsx(13) = domaine_vdf.elem_faces_pour_interp(elem_bas_droite,1);
              voisinsx(14) = domaine_vdf.elem_faces_pour_interp(elem_haut_gauche,1+dimension);
              voisinsx(15) = domaine_vdf.elem_faces_pour_interp(elem_bas_gauche,1+dimension);
              voisinsx(16) = domaine_vdf.elem_faces_pour_interp(elem_haut_droite,1+dimension);
              voisinsx(17) = domaine_vdf.elem_faces_pour_interp(elem_bas_droite,1+dimension);
              voisinsx(18) = domaine_vdf.elem_faces_pour_interp(elem_gauche,2+dimension);
              voisinsx(19) = domaine_vdf.elem_faces_pour_interp(elem_gauche,2);
              voisinsx(20) = domaine_vdf.elem_faces_pour_interp(elem_droite,2+dimension);
              voisinsx(21) = domaine_vdf.elem_faces_pour_interp(elem_droite,2);
              voisinsx(22) = domaine_vdf.elem_faces_pour_interp(elem_arriere_gauche,2);
              voisinsx(23) = domaine_vdf.elem_faces_pour_interp(elem_arriere_droite,2);
              voisinsx(24) = domaine_vdf.elem_faces_pour_interp(elem_avant_gauche,2);
              voisinsx(25) = domaine_vdf.elem_faces_pour_interp(elem_avant_droite,2);
              voisinsx(26) = domaine_vdf.elem_faces_pour_interp(elem_arriere_gauche,2+dimension);
              voisinsx(27) = domaine_vdf.elem_faces_pour_interp(elem_arriere_droite,2+dimension);
              voisinsx(28) = domaine_vdf.elem_faces_pour_interp(elem_avant_gauche,2+dimension);
              voisinsx(29) = domaine_vdf.elem_faces_pour_interp(elem_avant_droite,2+dimension);
              for (int i=0; i<nb_voisinsx; i++)
                {
                  if (voisinsx(i)<0) return 0;
                }
              gradUx(face,0,0) = 	  (valeurs_champ(voisinsx(1))  - valeurs_champ(voisinsx(0)))  / (domaine_vdf.xv(voisinsx(1),0)  - domaine_vdf.xv(voisinsx(0),0));
              gradUx(face,0,1) =      (valeurs_champ(voisinsx(2))  - valeurs_champ(voisinsx(3)))  / (domaine_vdf.xv(voisinsx(2),1)  - domaine_vdf.xv(voisinsx(3),1));
              gradUx(face,0,2) =      (valeurs_champ(voisinsx(8))  - valeurs_champ(voisinsx(9)))  / (domaine_vdf.xv(voisinsx(8),2)  - domaine_vdf.xv(voisinsx(9),2));
              gradUx(face,1,0) = 1./2.*((valeurs_champ(voisinsx(6))  - valeurs_champ(voisinsx(7)))  / (domaine_vdf.xv(voisinsx(6),0)  - domaine_vdf.xv(voisinsx(7),0))  + (valeurs_champ(voisinsx(4))  - valeurs_champ(voisinsx(5)))  / (domaine_vdf.xv(voisinsx(4),0)  - domaine_vdf.xv(voisinsx(5),0)));
              gradUx(face,1,1) = 1./2.*((valeurs_champ(voisinsx(4))  - valeurs_champ(voisinsx(6)))  / (domaine_vdf.xv(voisinsx(4),1)  - domaine_vdf.xv(voisinsx(6),1))  + (valeurs_champ(voisinsx(5))  - valeurs_champ(voisinsx(7)))  / (domaine_vdf.xv(voisinsx(5),1)  - domaine_vdf.xv(voisinsx(7),1)));
              gradUx(face,1,2) = 1./4.*((valeurs_champ(voisinsx(14)) - valeurs_champ(voisinsx(15))) / (domaine_vdf.xv(voisinsx(14),2) - domaine_vdf.xv(voisinsx(15),2)) + (valeurs_champ(voisinsx(16)) - valeurs_champ(voisinsx(17))) / (domaine_vdf.xv(voisinsx(16),2) - domaine_vdf.xv(voisinsx(17),2))  + (valeurs_champ(voisinsx(10)) - valeurs_champ(voisinsx(11))) / (domaine_vdf.xv(voisinsx(10),2) - domaine_vdf.xv(voisinsx(11),2)) + (valeurs_champ(voisinsx(12)) - valeurs_champ(voisinsx(13))) / (domaine_vdf.xv(voisinsx(12),2) - domaine_vdf.xv(voisinsx(13),2)));
              gradUx(face,2,0) = 1./2.*((valeurs_champ(voisinsx(18)) - valeurs_champ(voisinsx(20))) / (domaine_vdf.xv(voisinsx(18),0) - domaine_vdf.xv(voisinsx(20),0)) + (valeurs_champ(voisinsx(19)) - valeurs_champ(voisinsx(21))) / (domaine_vdf.xv(voisinsx(19),0) - domaine_vdf.xv(voisinsx(21),0)));
              gradUx(face,2,1) = 1./4.*((valeurs_champ(voisinsx(26)) - valeurs_champ(voisinsx(28))) / (domaine_vdf.xv(voisinsx(26),1) - domaine_vdf.xv(voisinsx(28),1)) + (valeurs_champ(voisinsx(27)) - valeurs_champ(voisinsx(29))) / (domaine_vdf.xv(voisinsx(27),1) - domaine_vdf.xv(voisinsx(29),1))  + (valeurs_champ(voisinsx(22)) - valeurs_champ(voisinsx(24))) / (domaine_vdf.xv(voisinsx(22),1) - domaine_vdf.xv(voisinsx(24),1)) + (valeurs_champ(voisinsx(23)) - valeurs_champ(voisinsx(25))) / (domaine_vdf.xv(voisinsx(23),1) - domaine_vdf.xv(voisinsx(25),1)));
              gradUx(face,2,2) = 1./2.*((valeurs_champ(voisinsx(18)) - valeurs_champ(voisinsx(19))) / (domaine_vdf.xv(voisinsx(18),2) - domaine_vdf.xv(voisinsx(19),2)) + (valeurs_champ(voisinsx(20)) - valeurs_champ(voisinsx(21))) / (domaine_vdf.xv(voisinsx(20),2) - domaine_vdf.xv(voisinsx(21),2)));
            }

          double xfact;
          double yfact;
          double zfact;

          DoubleTab Delta_x(dimension);

          Delta_x(0)=fabs(domaine_vdf.dist_face(faces_voisines(0),faces_voisines(1),0));
          Delta_x(1)=fabs(domaine_vdf.dist_face(faces_voisines(0),faces_voisines(2),1));
          Delta_x(2)=fabs(domaine_vdf.dist_face(faces_voisines(0),faces_voisines(4),2));

          xfact=fabs((coord(fa7,0)-coord_face(0,0))/Delta_x(0));
          yfact=fabs((coord(fa7,1)-coord_face(0,1))/Delta_x(1));
          zfact=fabs((coord(fa7,2)-coord_face(0,2))/Delta_x(2));

          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  resu(fa7,i,j)=(1-zfact)*((1-yfact)*((1-xfact)*(gradUx(0,i,j)) + xfact*(gradUx(1,i,j))) +
                                           yfact*((1-xfact)*(gradUx(2,i,j)) + xfact*(gradUx(3,i,j)))) +
                                zfact*((1-yfact)*((1-xfact)*(gradUx(4,i,j)) + xfact*(gradUx(5,i,j))) +
                                       yfact*((1-xfact)*(gradUx(6,i,j)) + xfact*(gradUx(7,i,j))));
                }
            }
        }
    }
  return 1;
}
// EB
/*! @brief Cette fonction calcule les composantes du tenseur gradU aux points P1 pour chaque facette du maillage lagrangien.
 *  En entree: indicatrice_face, indicatrice, valeurs_champs (tableau de vitesse), coord : coordonnes xyz des points P1 pour chaque facette lagrangienne, resu : composantes
 *  du tenseur interpole pour chaque facette lagrangienne.
 *  Boucle sur les facetttes lagrangiennes :
 *  	Dans un premier temps, identification des 8 elements voisins (4 en 2D) les plus proches.
 *  	Ensuite, pour chaque element, calcul du tenseur des contraintes
 *  	Pour chaque composante du tenseur, interpolation trilineaire au point de coordonnees (coord(fa7,0), coord(fa7,1), cooord(fa7,2)).
 *
 */
int Navier_Stokes_FT_Disc::trilinear_interpolation_gradU_elem(const DoubleTab& indicatrice_face, const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu)
{
  // On identifie l'element dans lequel appartient le point de coordonnees coord
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  int nb_voisins=8;
  const double rho_fluide = fluide_diphasique().fluide_phase(1).masse_volumique().valeurs()(0, 0);
  const double mu_fluide =rho_fluide * fluide_diphasique().fluide_phase(1).viscosite_cinematique().valeur().valeurs()(0, 0);

  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!maillage.facette_virtuelle(fa7))
        {
          // On recupere les elements 8 voisins ou le tenseur gradient sera calcule
          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect elem_voisins(nb_voisins);
          chercher_elem_voisins(indicatrice,coord_elem_interp,elem_voisins);

          for (int i=0; i<nb_voisins; i++)
            {
              if (elem_voisins(i)<0) return 0;
              //if (indicatrice(elem_voisins(i))!=1) Cerr << "L'interpolation du gradient de la vitesse utilise un element diphasique"<< finl;
            }
          // On calcule les distances entre mailles voisines dans chaque direction pour calculer les coefficients d'interpolation
          DoubleVect delta_i(dimension);
          delta_i(0) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(1), 0));
          delta_i(1) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(3), 1));
          delta_i(2) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(5), 2));
          DoubleVect coord_elem_0(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_0(dim)=domaine_vdf.xp(elem_voisins(0),dim);

          double xfact=fabs((coord_elem_interp(0)-coord_elem_0(0))/delta_i(0));
          double yfact=fabs((coord_elem_interp(1)-coord_elem_0(1))/delta_i(1));
          double zfact=fabs((coord_elem_interp(2)-coord_elem_0(2))/delta_i(2));

          //on calcule le tenseur gradient de la vitesse pour chaque element voisin
          DoubleTab gradU(nb_voisins, dimension, dimension);
          for (int elem=0; elem<nb_voisins; elem++)
            {

              // 														SCHEMA EN 2D
              //  					AVANT				 ---- ---- --- --- --- --- --- --- --- ---   AVANT
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //														 9	 10
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  										 ---- --- --- -5- -3- -6- --- --- ---
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //	  	  	  	  	GAUCHE  	  	  	  	  	  	     1 x 2              		DROITE
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  									     ---- --- --- -7- -4- -8- --- --- --- 	      			y
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			^
              //												         11  12									|
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			|
              //  					ARRIERE					 ---- --- --- --- --- --- --- --- ---    ARRIERE    o----->x
              //																					 			 z
              // Evaluation du gradient sur les faces de normale x (choix arbitraire mais peut tuer une certaine symetrie)
              // Les termes + designent les faces juxtaposees selon +z
              // Les termes - designent les faces juxtaposees selon -z
              // Les termes 57 et 86 pour w designent les faces de normale z secantes aux faces 5,7 et 8,6
              // 		   --					   													   	   									--
              //		   |   u_2-u_1	 		   				  1	  u_9-u_11	   u_10-u_12		          1   u_1+ - u_1- 	 u_2+ - u_2-	|
              //		   | ----------							  -	(---------- + ----------- )	        	  -	( -----------  + ----------- )	|
              //		   |   x_2-x1			 				  2	  y_9-y_11	   y_10-y_12				  2   z_1+ - z_1-	 z_2+ - z_2- 	|
              //		   |												    																|
              //  grad U = | 1    v_5-v_6	 v_7 - v_8    			 	v_3-v_4								  1	  v_3+ -v_3-	v_4+ - v_4-     |
              // 	  	   | - ( --------- + --------- ) 			 	------- 						      - ( ----------- + ----------- )	|
              //		   | 2    x_5-x_6	 x_7 - x_8   			 	y_3-y_4						     	  2	  z_3+ - z_3-	z_4+ - z_4-		|
              //		   |																													|
              // 	  	   | 1	 w57+ - w86+	 w57- - w86-	  1  	w1112+ - w910+    w1112- -w910-  	      w_34+ - w_34-					|
              //		   | - ( ----------- +	 ----------- )	  - (	-------------- + -------------- )		  -------------	    			|
              //  		   | 2	 x57+ - x86+	 x57- - x86-  	  2	    y1112+ - y910+	  y1112--y910-	  		  z_34+ - z_34-					|
              //		   --					   		   									              									   --
              //
              // Pour chaque face, il faut donc 30 faces voisines : 1-12, 1+, 1-, 2+, 2-, 3+, 3-, 4+, 4-, 34+, 34-, 57+, 57-, 86+, 86-, 1112+, 1112-, 910+, 910-
              // Pour plus de facilites, on numerote :  1+:13, 1-:14, 2+:15, 2-:16, 3+:17, 3-:18, 4+:19, 4-:20, 34+:21, 34-:22, 57+:23, 57-:24, 86+:25, 86-:26, 1112+:27, 1112-:28, 910+:29, 910-:30

              // On identifie les faces voisines pour le calcul des composantes du tenseur
              int nb_faces_voisines=30;
              IntVect les_faces_voisines(nb_faces_voisines);
              int elem_=elem_voisins(elem);
              int elem_gauche = domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,0),0);
              int elem_droite = domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,0+dimension),1);
              int elem_haut=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,2+dimension), 1);
              int elem_bas=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,2), 0);
              int elem_avant=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,1+dimension),1);
              int elem_arriere=domaine_vf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,1),0);

              les_faces_voisines(0)=domaine_vdf.elem_faces_pour_interp(elem_,0);
              les_faces_voisines(1)=domaine_vdf.elem_faces_pour_interp(elem_,0+dimension);
              les_faces_voisines(2)=domaine_vdf.elem_faces_pour_interp(elem_,1);
              les_faces_voisines(3)=domaine_vdf.elem_faces_pour_interp(elem_,1+dimension);
              les_faces_voisines(4)=domaine_vdf.elem_faces_pour_interp(elem_gauche,1+dimension);
              les_faces_voisines(5)=domaine_vdf.elem_faces_pour_interp(elem_droite,1+dimension);
              les_faces_voisines(6)=domaine_vdf.elem_faces_pour_interp(elem_gauche,1);
              les_faces_voisines(7)=domaine_vdf.elem_faces_pour_interp(elem_droite,1);
              les_faces_voisines(8)=domaine_vdf.elem_faces_pour_interp(elem_avant,0);
              les_faces_voisines(9)=domaine_vdf.elem_faces_pour_interp(elem_avant,0+dimension);
              les_faces_voisines(10)=domaine_vdf.elem_faces_pour_interp(elem_arriere,0);
              les_faces_voisines(11)=domaine_vdf.elem_faces_pour_interp(elem_arriere,0+dimension);

              //Cerr << "Elem de coord " << domaine_vdf.xp(elem_,0) << " " << domaine_vdf.xp(elem_,1) << " " << domaine_vdf.xp(elem_,2) << finl;
              //for (int ii=4; ii<12; ii++) Cerr << "Face voisine " << ii << " de coord " << domaine_vdf.xv(les_faces_voisines(ii),0) << " " << domaine_vdf.xv(les_faces_voisines(ii),1) << " " << domaine_vdf.xv(les_faces_voisines(ii),2) << finl;
              les_faces_voisines(12)=domaine_vdf.elem_faces_pour_interp(elem_haut,0);
              les_faces_voisines(13)=domaine_vdf.elem_faces_pour_interp(elem_bas,0);
              les_faces_voisines(14)=domaine_vdf.elem_faces_pour_interp(elem_haut,0+dimension);
              les_faces_voisines(15)=domaine_vdf.elem_faces_pour_interp(elem_bas,0+dimension);
              les_faces_voisines(16)=domaine_vdf.elem_faces_pour_interp(elem_haut,1+dimension);
              les_faces_voisines(17)=domaine_vdf.elem_faces_pour_interp(elem_bas,1+dimension);
              les_faces_voisines(18)=domaine_vdf.elem_faces_pour_interp(elem_haut,1);
              les_faces_voisines(19)=domaine_vdf.elem_faces_pour_interp(elem_bas,1);
              les_faces_voisines(20)=domaine_vdf.elem_faces_pour_interp(elem_,2);
              les_faces_voisines(21)=domaine_vdf.elem_faces_pour_interp(elem_,2+dimension);
              les_faces_voisines(22)=domaine_vdf.elem_faces_pour_interp(elem_gauche,2);
              les_faces_voisines(23)=domaine_vdf.elem_faces_pour_interp(elem_gauche,2+dimension);
              les_faces_voisines(24)=domaine_vdf.elem_faces_pour_interp(elem_droite,2);
              les_faces_voisines(25)=domaine_vdf.elem_faces_pour_interp(elem_droite,2+dimension);
              les_faces_voisines(26)=domaine_vdf.elem_faces_pour_interp(elem_arriere,2);
              les_faces_voisines(27)=domaine_vdf.elem_faces_pour_interp(elem_arriere,2+dimension);
              les_faces_voisines(28)=domaine_vdf.elem_faces_pour_interp(elem_avant,2);
              les_faces_voisines(29)=domaine_vdf.elem_faces_pour_interp(elem_avant,2+dimension);

              for (int i=0; i<nb_faces_voisines; i++)
                {
                  if (les_faces_voisines(i)<0) return 0;
                  /*if (indicatrice_face(les_faces_voisines(i))<1)
                    {
                      if (indicatrice(domaine_vdf.face_voisins_pour_interp(les_faces_voisines(i),0))<1 && indicatrice(domaine_vdf.face_voisins_pour_interp(les_faces_voisines(i),1))) Cerr << "Une des face voisines pour le calcul du gradient est solide, coord : " <<domaine_vdf.xv(les_faces_voisines(i),0) << " " << domaine_vdf.xv(les_faces_voisines(i),1) << " " <<domaine_vdf.xv(les_faces_voisines(i),2)  << finl;
                    } */
                }
              gradU(elem,0,0) = 	    ( (valeurs_champ(les_faces_voisines(1))  - valeurs_champ(les_faces_voisines(0)))  / (domaine_vdf.xv(les_faces_voisines(1),0)  - domaine_vdf.xv(les_faces_voisines(0),0)));
              gradU(elem,0,1) = 1./2.*( (valeurs_champ(les_faces_voisines(8))  - valeurs_champ(les_faces_voisines(10))) / (domaine_vdf.xv(les_faces_voisines(8),1)  - domaine_vdf.xv(les_faces_voisines(10),1))   + (valeurs_champ(les_faces_voisines(9))  - valeurs_champ(les_faces_voisines(11)))  / (domaine_vdf.xv(les_faces_voisines(9),1)   - domaine_vdf.xv(les_faces_voisines(11),1)));
              gradU(elem,0,2) = 1./2.*( (valeurs_champ(les_faces_voisines(12)) - valeurs_champ(les_faces_voisines(13))) / (domaine_vdf.xv(les_faces_voisines(12),2) - domaine_vdf.xv(les_faces_voisines(13),2))   + (valeurs_champ(les_faces_voisines(14)) - valeurs_champ(les_faces_voisines(15)))  / (domaine_vdf.xv(les_faces_voisines(14),2)  - domaine_vdf.xv(les_faces_voisines(15),2)));
              gradU(elem,1,0) = 1./2.*( (valeurs_champ(les_faces_voisines(4))  - valeurs_champ(les_faces_voisines(5)))  / (domaine_vdf.xv(les_faces_voisines(4),0)  - domaine_vdf.xv(les_faces_voisines(5),0))    + (valeurs_champ(les_faces_voisines(6))  - valeurs_champ(les_faces_voisines(7)))   / (domaine_vdf.xv(les_faces_voisines(6),0)   - domaine_vdf.xv(les_faces_voisines(7),0)));
              gradU(elem,1,1) = 	    ( (valeurs_champ(les_faces_voisines(2))  - valeurs_champ(les_faces_voisines(3)))  / (domaine_vdf.xv(les_faces_voisines(2),1)  - domaine_vdf.xv(les_faces_voisines(3),1)));
              gradU(elem,1,2) = 1./2.*( (valeurs_champ(les_faces_voisines(16)) - valeurs_champ(les_faces_voisines(17))) / (domaine_vdf.xv(les_faces_voisines(16),2) - domaine_vdf.xv(les_faces_voisines(17),2))   + (valeurs_champ(les_faces_voisines(18)) - valeurs_champ(les_faces_voisines(19)))  / (domaine_vdf.xv(les_faces_voisines(18),2)  - domaine_vdf.xv(les_faces_voisines(19),2)));
              gradU(elem,2,0) = 1./2.*( (valeurs_champ(les_faces_voisines(22)) - valeurs_champ(les_faces_voisines(24))) / (domaine_vdf.xv(les_faces_voisines(22),0) - domaine_vdf.xv(les_faces_voisines(24),0))   + (valeurs_champ(les_faces_voisines(23)) - valeurs_champ(les_faces_voisines(25)))  / (domaine_vdf.xv(les_faces_voisines(23),0)  - domaine_vdf.xv(les_faces_voisines(25),0)));
              gradU(elem,2,1) = 1./2.*( (valeurs_champ(les_faces_voisines(26)) - valeurs_champ(les_faces_voisines(28))) / (domaine_vdf.xv(les_faces_voisines(26),1) - domaine_vdf.xv(les_faces_voisines(28),1))   + (valeurs_champ(les_faces_voisines(27)) - valeurs_champ(les_faces_voisines(29)))  / (domaine_vdf.xv(les_faces_voisines(27),1)  - domaine_vdf.xv(les_faces_voisines(29),1)));
              gradU(elem,2,2) = 	  ( (valeurs_champ(les_faces_voisines(20)) - valeurs_champ(les_faces_voisines(21))) / (domaine_vdf.xv(les_faces_voisines(20),2) - domaine_vdf.xv(les_faces_voisines(21),2)));


            }

          // On fait une interpolation trilineaire de chaque composante du tenseur
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  resu(fa7,i,j)=mu_fluide*((1-zfact)*((1-yfact)*((1-xfact)*(gradU(0,i,j)) + xfact*(gradU(1,i,j))) +
                                                      yfact*((1-xfact)*(gradU(2,i,j)) + xfact*(gradU(3,i,j)))) +
                                           zfact*((1-yfact)*((1-xfact)*(gradU(4,i,j)) + xfact*(gradU(5,i,j))) +
                                                  yfact*((1-xfact)*(gradU(6,i,j)) + xfact*(gradU(7,i,j)))));
                }
            }
        }
    }



  return 1;
}
// EB
/*! @brief Voir Navier_Stokes_FT_Disc::trilinear_interpolation_gradU_elem. La seule difference avec la methode citee est que la viscosite aux aretes est utilisee pour calculer
 *  les composantes du tenseur gradU. On ne calcule donc pas gradU mais mu*gradU.
 *
 */
int Navier_Stokes_FT_Disc::trilinear_interpolation_gradU_elem_P1(const DoubleTab& indicatrice_face, const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu)
{
  // On identifie l'element dans lequel appartient le point de coordonnees coord
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  //const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  int nb_voisins=8;
  const double mu_f =fluide_diphasique().fluide_phase(1).viscosite_dynamique().valeur().valeurs()(0, 0);
  const double mu_p = fluide_diphasique().fluide_phase(0).viscosite_dynamique().valeur().valeurs()(0, 0);

  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!maillage.facette_virtuelle(fa7))
        {
          // On recupere les elements 8 voisins ou le tenseur gradient sera calcule
          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect elem_voisins(nb_voisins);
          chercher_elem_voisins(indicatrice,coord_elem_interp,elem_voisins);

          for (int i=0; i<nb_voisins; i++)
            {
              if (elem_voisins(i)<0) return 0;
              //if (indicatrice(elem_voisins(i))!=1) Cerr << "L'interpolation du gradient de la vitesse utilise un element diphasique"<< finl;
            }
          // On calcule les distances entre mailles voisines dans chaque direction pour calculer les coefficients d'interpolation
          DoubleVect delta_i(dimension);
          delta_i(0) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(1), 0));
          delta_i(1) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(3), 1));
          delta_i(2) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(5), 2));
          DoubleVect coord_elem_0(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_0(dim)=domaine_vdf.xp(elem_voisins(0),dim);

          double xfact=fabs((coord_elem_interp(0)-coord_elem_0(0))/delta_i(0));
          double yfact=fabs((coord_elem_interp(1)-coord_elem_0(1))/delta_i(1));
          double zfact=fabs((coord_elem_interp(2)-coord_elem_0(2))/delta_i(2));

          //on calcule le tenseur gradient de la vitesse pour chaque element voisin
          DoubleTab gradU(nb_voisins, dimension, dimension);
          for (int elem=0; elem<nb_voisins; elem++)
            {

              // 														SCHEMA EN 2D
              //  					AVANT				 ---- ---- --- --- --- --- --- --- ---     AVANT
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //														 9	 10
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  										 ---- --- --- -5- -3- -6- --- --- ---
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //	  	  	  	  	GAUCHE  	  	  	  	  	  	     1 x 2              		DROITE
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  									     ---- --- --- -7- -4- -8- --- --- --- 	      			y
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			^
              //												         11  12									|
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			|
              //  					ARRIERE				 ---- --- --- --- --- --- --- --- --- ---    ARRIERE    o----->x
              //																					 			 z
              // Evaluation du gradient sur les faces de normale x (choix arbitraire mais peut tuer une certaine symetrie)
              // Les termes + designent les faces juxtaposees selon +z
              // Les termes - designent les faces juxtaposees selon -z
              // Les termes 57 et 86 pour w designent les faces de normale z secantes aux faces 5,7 et 8,6
              //
              //   ON CALCULE mu*gradU et pas gradU
              //   Pour plus de lisibilite, mu n'est pas reecrit pour chaque composante (mais sans ce mu_h local, la decomposition n'as plus d'interet)
              //   ON UTILISE LE MU HARMONIQUE LOCAL, IE : LE MU_HARMONIQUE AUX ARETES
              //
              // 		   --					   													   	   																															   								   --
              //		   |   u_2-u_1	 		   				                                   1  1  u_9-u_1	 u_1-u_11   1  u_10-u_2	   u_2-u_12	          								1     1   u_1+ - u_1    u_1 - u_1-	   1   u_2+ - u_2   u_2 - u_2-      |
              //		   | ----------							                                   -( - (--------- + --------)+ - (--------- + ---------) )	      								-	( - ( ----------- + ---------- ) + - ( ----------- + ---------- ) ) |
              //		   |   x_2-x1			 				                                   2  2  y_9-y_1	 y_1-y_11   2  y_10-y_2	   y_2-y_12			  								2     2   z_1+ - z_1	z_1 - z_1- 	   2   z_2+ - z_2   z_2 - z_2-      |
              //		   |												    																																														|
              //  grad U = | 1   1   v_5-v_3	v_3 - v_6    1  v_7-v_4	   v_4-v_8 	                       v_3-v_4								   				  								1     1   v_3+ - v_3    v_3 - v_3-	   1   v_4+ - v_4    v_4 - v_4-		|
              // 	  	   | - ( - ( -------- + ---------)+ - (-------- + --------) ) 			 	       ------- 						    					  								-	( - ( ----------- + ---------- ) + - ( ----------- + ---------- ) )	|
              //		   | 2   2   x_5-x_3    x_3 - x_6    2  x_7-x_4	   x_4-x_8	                       y_3-y_4						     					  								2     2   z_3+ - z_3	z_3 - z_3- 	   2   z_4+ - z_4    z_4 - z_4- 	|
              //		   |																																																											|
              // 	  	   | 1	 1  w57+ - w34+	 w34+ - w86+	 1 w57- - w34-	 w34- - w86-       1   1   w_1112+ - w_34+   w_34+ - w_910+	  1     w_1112- - w_34-   w_34- - w_910-	      									w_34+ - w_34-							|
              //		   | - ( - (---------- + ----------- ) + -(---------- + ------------) )    - ( - ( --------------- + -------------- ) + - ( --------------- + -------------- ) )					    				-------------	    					|
              //  		   | 2	 2  x57+ - x34+	 x34+ - x86+     2 x57- - x34-	 x34- - x86-       2   2   y_1112+ - y_34+	 y_34+ - y_910+   2     y_1112- - y_34-   y_34- - y_910- 	  		  	     					    z_34+ - z_34-							|
              //		   --					   		   									              									 																						   								   --
              //
              // Pour chaque face, il faut donc 30 faces voisines : 1-12, 1+, 1-, 2+, 2-, 3+, 3-, 4+, 4-, 34+, 34-, 57+, 57-, 86+, 86-, 1112+, 1112-, 910+, 910-
              // Pour plus de facilites, on numerote :  1+:13, 1-:14, 2+:15, 2-:16, 3+:17, 3-:18, 4+:19, 4-:20, 34+:21, 34-:22, 57+:23, 57-:24, 86+:25, 86-:26, 1112+:27, 1112-:28, 910+:29, 910-:30

              // On identifie les faces voisines pour le calcul des composantes du tenseur
              int nb_faces_voisines=30;
              IntVect les_faces_voisines(nb_faces_voisines);
              int elem_=elem_voisins(elem);
              int elem_gauche = domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,0),0);
              int elem_droite = domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,0+dimension),1);
              int elem_haut=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,2+dimension), 1);
              int elem_bas=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,2), 0);
              int elem_avant=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,1+dimension),1);
              int elem_arriere=domaine_vdf.face_voisins_pour_interp(domaine_vdf.elem_faces_pour_interp(elem_,1),0);

              les_faces_voisines(0)=domaine_vdf.elem_faces_pour_interp(elem_,0);
              les_faces_voisines(1)=domaine_vdf.elem_faces_pour_interp(elem_,0+dimension);
              les_faces_voisines(2)=domaine_vdf.elem_faces_pour_interp(elem_,1);
              les_faces_voisines(3)=domaine_vdf.elem_faces_pour_interp(elem_,1+dimension);
              les_faces_voisines(4)=domaine_vdf.elem_faces_pour_interp(elem_gauche,1+dimension);
              les_faces_voisines(5)=domaine_vdf.elem_faces_pour_interp(elem_droite,1+dimension);
              les_faces_voisines(6)=domaine_vdf.elem_faces_pour_interp(elem_gauche,1);
              les_faces_voisines(7)=domaine_vdf.elem_faces_pour_interp(elem_droite,1);
              les_faces_voisines(8)=domaine_vdf.elem_faces_pour_interp(elem_avant,0);
              les_faces_voisines(9)=domaine_vdf.elem_faces_pour_interp(elem_avant,0+dimension);
              les_faces_voisines(10)=domaine_vdf.elem_faces_pour_interp(elem_arriere,0);
              les_faces_voisines(11)=domaine_vdf.elem_faces_pour_interp(elem_arriere,0+dimension);

              //Cerr << "Elem de coord " << zone_vdf.xp(elem_,0) << " " << zone_vdf.xp(elem_,1) << " " << zone_vdf.xp(elem_,2) << finl;
              //for (int ii=4; ii<12; ii++) Cerr << "Face voisine " << ii << " de coord " << domaine_vdf.xv(les_faces_voisines(ii),0) << " " << domaine_vdf.xv(les_faces_voisines(ii),1) << " " << domaine_vdf.xv(les_faces_voisines(ii),2) << finl;
              les_faces_voisines(12)=domaine_vdf.elem_faces_pour_interp(elem_haut,0);
              les_faces_voisines(13)=domaine_vdf.elem_faces_pour_interp(elem_bas,0);
              les_faces_voisines(14)=domaine_vdf.elem_faces_pour_interp(elem_haut,0+dimension);
              les_faces_voisines(15)=domaine_vdf.elem_faces_pour_interp(elem_bas,0+dimension);
              les_faces_voisines(16)=domaine_vdf.elem_faces_pour_interp(elem_haut,1+dimension);
              les_faces_voisines(17)=domaine_vdf.elem_faces_pour_interp(elem_bas,1+dimension);
              les_faces_voisines(18)=domaine_vdf.elem_faces_pour_interp(elem_haut,1);
              les_faces_voisines(19)=domaine_vdf.elem_faces_pour_interp(elem_bas,1);
              les_faces_voisines(20)=domaine_vdf.elem_faces_pour_interp(elem_,2);
              les_faces_voisines(21)=domaine_vdf.elem_faces_pour_interp(elem_,2+dimension);
              les_faces_voisines(22)=domaine_vdf.elem_faces_pour_interp(elem_gauche,2);
              les_faces_voisines(23)=domaine_vdf.elem_faces_pour_interp(elem_gauche,2+dimension);
              les_faces_voisines(24)=domaine_vdf.elem_faces_pour_interp(elem_droite,2);
              les_faces_voisines(25)=domaine_vdf.elem_faces_pour_interp(elem_droite,2+dimension);
              les_faces_voisines(26)=domaine_vdf.elem_faces_pour_interp(elem_arriere,2);
              les_faces_voisines(27)=domaine_vdf.elem_faces_pour_interp(elem_arriere,2+dimension);
              les_faces_voisines(28)=domaine_vdf.elem_faces_pour_interp(elem_avant,2);
              les_faces_voisines(29)=domaine_vdf.elem_faces_pour_interp(elem_avant,2+dimension);


              for (int i=0; i<nb_faces_voisines; i++)
                {
                  if (les_faces_voisines(i)<0) return 0;
                  /*if (indicatrice_face(les_faces_voisines(i))<1)
                    {
                      if (indicatrice(domaine_vdf.face_voisins_pour_interp(les_faces_voisines(i),0))<1 && indicatrice(domaine_vdf.face_voisins_pour_interp(les_faces_voisines(i),1))) Cerr << "Une des face voisines pour le calcul du gradient est solide, coord : " <<zone_vdf.xv(les_faces_voisines(i),0) << " " << zone_vdf.xv(les_faces_voisines(i),1) << " " <<zone_vdf.xv(les_faces_voisines(i),2)  << finl;
                    } */
                }
              int compo= static_cast<int>(get_num_compo().valeur().valeurs()(elem));
              gradU(elem,0,0) = 	    ( (mu_p*mu_f/(mu_f-indicatrice(elem_voisins(elem))*(mu_f-mu_p)))*(valeurs_champ(les_faces_voisines(1))  - valeurs_champ(les_faces_voisines(0)))  / (domaine_vdf.xv(les_faces_voisines(1),0)  - domaine_vdf.xv(les_faces_voisines(0),0)));
              gradU(elem,0,1) = 1./2.*  ( 1./2.*( calculer_viscosite_arete(les_faces_voisines(8),les_faces_voisines(0),compo)*(valeurs_champ(les_faces_voisines(8))  - valeurs_champ(les_faces_voisines(0))) / (domaine_vdf.xv(les_faces_voisines(8),1)  - domaine_vdf.xv(les_faces_voisines(0),1))   + calculer_viscosite_arete(les_faces_voisines(0),les_faces_voisines(10),compo)*(valeurs_champ(les_faces_voisines(0))  - valeurs_champ(les_faces_voisines(10)))  / (domaine_vdf.xv(les_faces_voisines(0),1) - domaine_vdf.xv(les_faces_voisines(10),1)))
                                          + 1./2.*( calculer_viscosite_arete(les_faces_voisines(9),les_faces_voisines(1),compo)*(valeurs_champ(les_faces_voisines(9))  - valeurs_champ(les_faces_voisines(1))) / (domaine_vdf.xv(les_faces_voisines(9),1)  - domaine_vdf.xv(les_faces_voisines(1),1))   + calculer_viscosite_arete(les_faces_voisines(1),les_faces_voisines(11),compo)*(valeurs_champ(les_faces_voisines(1))  - valeurs_champ(les_faces_voisines(11)))  / (domaine_vdf.xv(les_faces_voisines(1),1) - domaine_vdf.xv(les_faces_voisines(11),1))) );
              gradU(elem,0,2) = 1./2.*  ( 1./2.*( calculer_viscosite_arete(les_faces_voisines(12),les_faces_voisines(0),compo)*(valeurs_champ(les_faces_voisines(12)) - valeurs_champ(les_faces_voisines(0))) / (domaine_vdf.xv(les_faces_voisines(12),2) - domaine_vdf.xv(les_faces_voisines(0),2))   + calculer_viscosite_arete(les_faces_voisines(0),les_faces_voisines(13),compo)*(valeurs_champ(les_faces_voisines(0)) - valeurs_champ(les_faces_voisines(13)))  / (domaine_vdf.xv(les_faces_voisines(0),2)  - domaine_vdf.xv(les_faces_voisines(13),2)))
                                          + 1./2.*( calculer_viscosite_arete(les_faces_voisines(14),les_faces_voisines(1),compo)*(valeurs_champ(les_faces_voisines(14)) - valeurs_champ(les_faces_voisines(1))) / (domaine_vdf.xv(les_faces_voisines(14),2) - domaine_vdf.xv(les_faces_voisines(1),2))   + calculer_viscosite_arete(les_faces_voisines(1),les_faces_voisines(15),compo)*(valeurs_champ(les_faces_voisines(1)) - valeurs_champ(les_faces_voisines(15)))  / (domaine_vdf.xv(les_faces_voisines(1),2)  - domaine_vdf.xv(les_faces_voisines(15),2))) );
              gradU(elem,1,0) = 1./2.*  ( 1./2.*( calculer_viscosite_arete(les_faces_voisines(4),les_faces_voisines(2),compo)*(valeurs_champ(les_faces_voisines(4))  - valeurs_champ(les_faces_voisines(2)))  / (domaine_vdf.xv(les_faces_voisines(4),0)  - domaine_vdf.xv(les_faces_voisines(2),0))    + calculer_viscosite_arete(les_faces_voisines(3),les_faces_voisines(5),compo)*(valeurs_champ(les_faces_voisines(2))  - valeurs_champ(les_faces_voisines(5)))   / (domaine_vdf.xv(les_faces_voisines(2),0)   - domaine_vdf.xv(les_faces_voisines(5),0)))
                                          + 1./2.*( calculer_viscosite_arete(les_faces_voisines(6),les_faces_voisines(3),compo)*(valeurs_champ(les_faces_voisines(6))  - valeurs_champ(les_faces_voisines(3)))  / (domaine_vdf.xv(les_faces_voisines(6),0)  - domaine_vdf.xv(les_faces_voisines(3),0))    + calculer_viscosite_arete(les_faces_voisines(3),les_faces_voisines(7),compo)*(valeurs_champ(les_faces_voisines(3))  - valeurs_champ(les_faces_voisines(7)))   / (domaine_vdf.xv(les_faces_voisines(3),0)   - domaine_vdf.xv(les_faces_voisines(7),0))) );
              gradU(elem,1,1) = 	    ( (mu_p*mu_f/(mu_f-indicatrice(elem_voisins(elem))*(mu_f-mu_p)))*(valeurs_champ(les_faces_voisines(2))  - valeurs_champ(les_faces_voisines(3)))  / (domaine_vdf.xv(les_faces_voisines(2),1)  - domaine_vdf.xv(les_faces_voisines(3),1)));
              gradU(elem,1,2) = 1./2.*  ( 1./2.*( calculer_viscosite_arete(les_faces_voisines(16),les_faces_voisines(2),compo)*(valeurs_champ(les_faces_voisines(16)) - valeurs_champ(les_faces_voisines(2))) / (domaine_vdf.xv(les_faces_voisines(16),2) - domaine_vdf.xv(les_faces_voisines(2),2))   + calculer_viscosite_arete(les_faces_voisines(2),les_faces_voisines(17),compo)*(valeurs_champ(les_faces_voisines(2)) - valeurs_champ(les_faces_voisines(17)))  / (domaine_vdf.xv(les_faces_voisines(2),2)  - domaine_vdf.xv(les_faces_voisines(17),2)))
                                          + 1./2.*( calculer_viscosite_arete(les_faces_voisines(18),les_faces_voisines(3),compo)*(valeurs_champ(les_faces_voisines(18)) - valeurs_champ(les_faces_voisines(3))) / (domaine_vdf.xv(les_faces_voisines(18),2) - domaine_vdf.xv(les_faces_voisines(3),2))   + calculer_viscosite_arete(les_faces_voisines(3),les_faces_voisines(19),compo)*(valeurs_champ(les_faces_voisines(3)) - valeurs_champ(les_faces_voisines(19)))  / (domaine_vdf.xv(les_faces_voisines(3),2)  - domaine_vdf.xv(les_faces_voisines(19),2))));
              gradU(elem,2,0) = 1./2.*  ( 1./2.*( calculer_viscosite_arete(les_faces_voisines(22),les_faces_voisines(20),compo)*(valeurs_champ(les_faces_voisines(22)) - valeurs_champ(les_faces_voisines(20))) / (domaine_vdf.xv(les_faces_voisines(22),0) - domaine_vdf.xv(les_faces_voisines(20),0))   + calculer_viscosite_arete(les_faces_voisines(20),les_faces_voisines(24),compo)*(valeurs_champ(les_faces_voisines(20)) - valeurs_champ(les_faces_voisines(24)))  / (domaine_vdf.xv(les_faces_voisines(20),0)  - domaine_vdf.xv(les_faces_voisines(24),0)))
                                          + 1./2.*( calculer_viscosite_arete(les_faces_voisines(23),les_faces_voisines(21),compo)*(valeurs_champ(les_faces_voisines(23)) - valeurs_champ(les_faces_voisines(21))) / (domaine_vdf.xv(les_faces_voisines(23),0) - domaine_vdf.xv(les_faces_voisines(21),0))   + calculer_viscosite_arete(les_faces_voisines(21),les_faces_voisines(25),compo)*(valeurs_champ(les_faces_voisines(21)) - valeurs_champ(les_faces_voisines(25)))  / (domaine_vdf.xv(les_faces_voisines(21),0)  - domaine_vdf.xv(les_faces_voisines(25),0))) );
              gradU(elem,2,1) = 1./2.*  ( 1./2.*( calculer_viscosite_arete(les_faces_voisines(26),les_faces_voisines(20),compo)*(valeurs_champ(les_faces_voisines(26)) - valeurs_champ(les_faces_voisines(20))) / (domaine_vdf.xv(les_faces_voisines(26),1) - domaine_vdf.xv(les_faces_voisines(20),1))   + calculer_viscosite_arete(les_faces_voisines(20),les_faces_voisines(28),compo)*(valeurs_champ(les_faces_voisines(20)) - valeurs_champ(les_faces_voisines(28)))  / (domaine_vdf.xv(les_faces_voisines(20),1)  - domaine_vdf.xv(les_faces_voisines(28),1)))
                                          + 1./2.*( calculer_viscosite_arete(les_faces_voisines(27),les_faces_voisines(21),compo)*(valeurs_champ(les_faces_voisines(27)) - valeurs_champ(les_faces_voisines(21))) / (domaine_vdf.xv(les_faces_voisines(27),1) - domaine_vdf.xv(les_faces_voisines(21),1))   + calculer_viscosite_arete(les_faces_voisines(21),les_faces_voisines(29),compo)*(valeurs_champ(les_faces_voisines(21)) - valeurs_champ(les_faces_voisines(29)))  / (domaine_vdf.xv(les_faces_voisines(21),1)  - domaine_vdf.xv(les_faces_voisines(29),1))));
              gradU(elem,2,2) = 	  ( (mu_p*mu_f/(mu_f-indicatrice(elem_voisins(elem))*(mu_f-mu_p)))*(valeurs_champ(les_faces_voisines(20)) - valeurs_champ(les_faces_voisines(21))) / (domaine_vdf.xv(les_faces_voisines(20),2) - domaine_vdf.xv(les_faces_voisines(21),2)));

            }

          // On fait une interpolation trilineaire de chaque composante du tenseur
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  resu(fa7,i,j)=((1-zfact)*((1-yfact)*((1-xfact)*(gradU(0,i,j)) + xfact*(gradU(1,i,j))) +
                                            yfact*((1-xfact)*(gradU(2,i,j)) + xfact*(gradU(3,i,j)))) +
                                 zfact*((1-yfact)*((1-xfact)*(gradU(4,i,j)) + xfact*(gradU(5,i,j))) +
                                        yfact*((1-xfact)*(gradU(6,i,j)) + xfact*(gradU(7,i,j)))));
                }
            }
        }
    }
  return 1;
}

// EB
/*! @brief Calcul de la viscosite aux aretes par la methode de Vincent et al (2014). Pour l'arete commune aux faces 1 et 2, on repartit uniformement dans le volume 10 points par direction.
 *  Pour chaque point, on calcule la distance au centre de gravite de la particule poour savoir si le point est solide ou fluide. On compte ensuite le nombre de points fluide pour connaitre
 *  l'indicatrice a l'arete.
 *  Cette methode est utilisee dans Navier_Stokes_FT_Disc::trilinear_interpolation_gradU_elem_P1.
 */
double Navier_Stokes_FT_Disc::calculer_viscosite_arete(int face1, int face2, int compo)
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, Domaine_dis().valeur());
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const DoubleTab& positions_compo=eq_transport.get_positions_compo();

  int elem1 = domaine_vdf.face_voisins_pour_interp(face1,0);
  int elem2 = domaine_vdf.face_voisins_pour_interp(face1,1);
  int elem3 = domaine_vdf.face_voisins_pour_interp(face2,0);
  int elem4 = domaine_vdf.face_voisins_pour_interp(face2,1);

  const double mu_f =fluide_diphasique().fluide_phase(1).viscosite_dynamique().valeur().valeurs()(0, 0);
  const double mu_p =fluide_diphasique().fluide_phase(0).viscosite_dynamique().valeur().valeurs()(0, 0);
  const double rayon_particule = eq_transport.get_rayons_compo()(compo);
  double indicatrice_arete=0;
  int pvx=10;
  int pvy=10;
  int pvz=10;
  double x_min,y_min,z_min,x_max;
  double x_cg=positions_compo(compo,0);
  double y_cg=positions_compo(compo,1);
  double z_cg=positions_compo(compo,2);
  x_min=std::min(domaine_vdf.xp(elem1,0),domaine_vdf.xp(elem2,0));
  x_min=std::min(x_min,domaine_vdf.xp(elem3,0));
  x_min=std::min(x_min,domaine_vdf.xp(elem4,0));
  x_max=std::max(domaine_vdf.xp(elem1,0),domaine_vdf.xp(elem2,0));
  x_max=std::max(x_min,domaine_vdf.xp(elem3,0));
  x_max=std::max(x_min,domaine_vdf.xp(elem4,0));
  y_min=std::min(domaine_vdf.xp(elem1,1),domaine_vdf.xp(elem2,1));
  y_min=std::min(y_min,domaine_vdf.xp(elem3,1));
  y_min=std::min(y_min,domaine_vdf.xp(elem4,1));
  z_min=std::min(domaine_vdf.xp(elem1,2),domaine_vdf.xp(elem2,2));
  z_min=std::min(z_min,domaine_vdf.xp(elem3,2));
  z_min=std::min(z_min,domaine_vdf.xp(elem4,2));

  double delta_x=x_max-x_min;

  for (int i=0; i<pvx; i++)
    {
      for (int j=0; j<pvy; j++)
        {
          for (int k=0; k<pvz; k++)
            {
              double x=x_min+i*delta_x/pvx;
              double y=y_min+i*delta_x/pvy;
              double z=z_min+i*delta_x/pvz;
              if (sqrt(pow(x-x_cg,2)+pow(y-y_cg,2)+pow(z-z_cg,2))>rayon_particule) indicatrice_arete+=1;
            }
        }
    }
  indicatrice_arete/=(pvx*pvy*pvz);

  return(mu_p*mu_f/(mu_f-indicatrice_arete*(mu_f-mu_p)));
}

/*! @brief Interpolation trilineaire d'un champs aux coordonnees du tableau coord
 * Il y a  3 interpolations : une utilisant les faces de normale x --> resu(fa7,0), la faces de nomale y --> resu(fa7,1) et les faces de normale z --> resu(fa7,2)
 */
int Navier_Stokes_FT_Disc::trilinear_interpolation_face(const DoubleTab& indicatrice_faces, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu)
{
  // On identifie l'element dans lequel appartient le point de coordonnees coord
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  int nb_voisins=8;
  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!maillage.facette_virtuelle(fa7))
        {
          // On recupere les faces voisines : les faces euleriennes qui vont servir a l'interpolation
          DoubleVect coord_elem_interp(dimension); // coordonnee du cg de la fa7
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntTab faces_voisines(dimension,nb_voisins); // 8 elements voisins au point de coordonnees coord. L'element elem est inclu dedans.
          chercher_faces_voisines_xyz(coord_elem_interp,faces_voisines);
          int acces_faces=1;
          for (int dim=0; dim<dimension; dim++)
            {
              for (int i=0; i<nb_voisins; i++)
                {
                  if (faces_voisines(dim,i)<0) acces_faces= 0; // s'il y a eu un bug dans la recuperation des faces (pbm d'acces : zone de joint, bord), on renvoie 0 et la fa7 n'est pas prise en compte dans le calcul
                }
            }
          if (acces_faces)
            {
              DoubleVect xfact(dimension);
              DoubleVect yfact(dimension);
              DoubleVect zfact(dimension);
              DoubleTab Delta_(dimension);
              Delta_(0)=fabs(domaine_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,1),0));
              Delta_(1)=fabs(domaine_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,2),1));
              Delta_(2)=fabs(domaine_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,4),2));
              for (int dim=0; dim<dimension; dim++)
                {
                  xfact(dim)=fabs((coord(fa7,0)-domaine_vdf.xv(faces_voisines(dim,0),0))/Delta_(0));
                  yfact(dim)=fabs((coord(fa7,1)-domaine_vdf.xv(faces_voisines(dim,0),1))/Delta_(1));
                  zfact(dim)=fabs((coord(fa7,2)-domaine_vdf.xv(faces_voisines(dim,0),2))/Delta_(2));
                }
              for (int dim=0; dim<dimension; dim++)
                {
                  if (xfact(dim)>1) Cerr << "xfact > 1 pour le point " <<coord(fa7,0)<< " " << coord(fa7,1)<< " " << coord(fa7,2) <<finl;
                  if (yfact(dim)>1) Cerr << "yfact > 1 pour le point " <<coord(fa7,0)<< " " << coord(fa7,1)<< " " << coord(fa7,2) <<finl;
                  if (zfact(dim)>1) Cerr << "zfact > 1 pour le point " <<coord(fa7,0)<< " " << coord(fa7,1)<< " " << coord(fa7,2) <<finl;
                }

              for (int dim=0; dim<dimension; dim++)
                {
                  resu(fa7,dim)=(1.-zfact(dim))*((1.-yfact(dim))*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,0))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,1)))) +
                                                 yfact(dim)*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,2))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,3))))) +
                                zfact(dim)*((1.-yfact(dim))*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,4))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,5)))) +
                                            yfact(dim)*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,6))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,7)))));
                }
            }
          else
            {
              for (int dim=0; dim<dimension; dim++)
                {
                  resu(fa7,dim)=-1e15; // de cette maniere, on ne calcule pas la force pour la fa7 pour laquelle on n'a pas acces a P2
                }
            }
        }
    }
  return 1;
}

/*! @brief Interpolation trilineaire d'un champs aux coordonnees du tableau coord
 * idem trilinear_interpolation_face mais le tableau coord contient les coordonnees des points P1 a partir des sommets lagrangien et non a partir des facettes lagrangiennes
 */
int Navier_Stokes_FT_Disc::trilinear_interpolation_face_sommets(const DoubleTab& indicatrice_faces, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu)
{
  // On identifie l'element dans lequel appartient le point de coordonnees coord
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const ArrOfInt& elem = maillage.sommet_elem();
  const DoubleTab& pos = maillage.sommets();
  const int nb_pos_tot = pos.dimension(0);
  int nb_voisins=8;
  for (int som=0; som<nb_pos_tot; som++)
    {
      const int element = elem[som];
      if (element >= 0)   // sommet reel ?
        {

          // On recupere les faces voisines : les faces euleriennes qui vont servir a l'interpolation
          DoubleVect coord_elem_interp(dimension); // coordonnee du cg de la fa7
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(som,dim);
          IntTab faces_voisines(dimension,nb_voisins); // 8 elements voisins au point de coordonnees coord. L'element elem est inclu dedans.
          chercher_faces_voisines_xyz(coord_elem_interp,faces_voisines);
          for (int dim=0; dim<dimension; dim++)
            {
              for (int i=0; i<nb_voisins; i++)
                {
                  if (faces_voisines(dim,i)<0) return 0; // s'il y a eu un bug dans la recuperation des faces (pbm d'acces : zone de joint, bord), on renvoie 0 et la fa7 n'est pas prise en compte dans le calcul
                }
            }


          DoubleVect xfact(dimension);
          DoubleVect yfact(dimension);
          DoubleVect zfact(dimension);

          DoubleTab Delta_(dimension);

          Delta_(0)=fabs(domaine_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,1),0));
          Delta_(1)=fabs(domaine_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,2),1));
          Delta_(2)=fabs(domaine_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,4),2));

          for (int dim=0; dim<dimension; dim++)
            {
              xfact(dim)=fabs((coord(som,0)-domaine_vdf.xv(faces_voisines(dim,0),0))/Delta_(0));
              yfact(dim)=fabs((coord(som,1)-domaine_vdf.xv(faces_voisines(dim,0),1))/Delta_(1));
              zfact(dim)=fabs((coord(som,2)-domaine_vdf.xv(faces_voisines(dim,0),2))/Delta_(2));
            }

          for (int dim=0; dim<dimension; dim++)
            {
              if (xfact(dim)>1) Cerr << "xfact > 1 pour le point " <<coord(som,0)<< " " << coord(som,1)<< " " << coord(som,2) <<finl;
              if (yfact(dim)>1) Cerr << "yfact > 1 pour le point " <<coord(som,0)<< " " << coord(som,1)<< " " << coord(som,2) <<finl;
              if (zfact(dim)>1) Cerr << "zfact > 1 pour le point " <<coord(som,0)<< " " << coord(som,1)<< " " << coord(som,2) <<finl;
            }
          //Cerr << "Delta_ " << Delta_(0) << " " << Delta_(1) << " " << Delta_(2) << finl;

          for (int dim=0; dim<dimension; dim++)
            {
              resu(som,dim)=(1.-zfact(dim))*((1.-yfact(dim))*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,0))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,1)))) +
                                             yfact(dim)*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,2))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,3))))) +
                            zfact(dim)*((1.-yfact(dim))*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,4))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,5)))) +
                                        yfact(dim)*((1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,6))) + xfact(dim)*(valeurs_champ(faces_voisines(dim,7)))));
            }
        }
    }

  return 1;
}

int Navier_Stokes_FT_Disc::trilinear_interpolation_elem(const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu)
{
  return trilinear_interpolation_elem(indicatrice, valeurs_champ, coord, resu, 0, 0);
}
// EB
/*! @brief Interpolation trilineaire d'un champs "valeurs_champs" aux coordonnees coord a partir des 8 elements (4 en 2D) voisins les plus proches. Valeurs_chaamps contient
 * typiquement le champ de presssion ou le champ de temperature.
 * On utilise egalement cette fonction pour sauvegarder la liste des elments P1 ainsi que l'indiatrice des elements P2
 */
int Navier_Stokes_FT_Disc::trilinear_interpolation_elem(const DoubleTab& indicatrice, const DoubleTab& valeurs_champ, DoubleTab& coord, DoubleTab& resu, const int is_P2, const int discr)
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  int nb_voisins=8; // 8 elements voisins au point de coordonnees coord.
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const int nb_fa7=coord.dimension(0);

  DoubleVect& tab_num_compo=variables_internes().num_compo.valeur().valeurs();
  IntTab& list_elem_P1=variables_internes().list_elem_P1_;
  IntTab& list_elem_P1_all=variables_internes().list_elem_P1_all_;
  const int sauv_list_P1 = (discr==0) ? 1 : 0; // si discr ==0 alors on discretize en P1
  if (sauv_list_P1)
    {
      list_elem_P1.resize(nb_fa7,2); // En (j,0) on stocke le num de l'element, en (j,1) on stocke le num de la particule
      list_elem_P1=-1;
    }
  const int sauv_list_P1_all = (discr==2) ? 1 : 0; // si discr ==2 alors on discretize sur tous les elements utilises pour l'interpolation en P1
  int nb_elem_P1_all;


  if (sauv_list_P1_all)
    {
      list_elem_P1_all.resize(0,2); // En (j,0) on stocke le num de l'element, en (j,1) on stocke le num de la particule
      list_elem_P1_all.set_smart_resize(1);
      nb_elem_P1_all=0;
    }
  const DoubleTab& cg_fa7=maillage.cg_fa7();

  for (int fa7=0; fa7<nb_fa7; fa7++)
    {
      if (!maillage.facette_virtuelle(fa7))
        {
          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect elem_voisins(nb_voisins);
          DoubleTab& indic_elem_P2=variables_internes().Indic_elem_P2_;
          if (is_P2)
            {
              indic_elem_P2(fa7)=chercher_elem_voisins(indicatrice,coord_elem_interp,elem_voisins);
            }
          else
            {
              chercher_elem_voisins(indicatrice,coord_elem_interp,elem_voisins,sauv_list_P1,fa7);
              if (sauv_list_P1)
                {
                  int elem_diph=domaine_vdf.domaine().chercher_elements(cg_fa7(fa7,0),cg_fa7(fa7,1),cg_fa7(fa7,2));
                  int num_compo = static_cast<int>(tab_num_compo(elem_diph));
                  list_elem_P1(fa7,1)= num_compo;
                }
            }
          int acces_elems=1;
          for (int i=0; i<nb_voisins; i++)
            {
              if (elem_voisins(i)<0) acces_elems= 0;
              if (sauv_list_P1_all)
                {
                  int elem_voisin=elem_voisins(i);
                  int elem_diph=domaine_vdf.domaine().chercher_elements(cg_fa7(fa7,0),cg_fa7(fa7,1),cg_fa7(fa7,2));
                  int num_compo = static_cast<int>(tab_num_compo(elem_diph));
                  if (elem_voisin>=0 && indicatrice(elem_voisin)==1)
                    {
                      list_elem_P1_all.append_line(elem_voisin,num_compo);
                      nb_elem_P1_all++;
                    }
                }
            }
          if (acces_elems)
            {
              DoubleVect delta_i(dimension);
              delta_i(0) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(1), 0));
              delta_i(1) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(2), 1));
              delta_i(2) = fabs(domaine_vdf.dist_elem(elem_voisins(0), elem_voisins(4), 2));
              DoubleVect coord_elem_0(dimension);
              for (int dim=0; dim<dimension; dim++) coord_elem_0(dim)=domaine_vdf.xp(elem_voisins(0),dim);

              double xfact=fabs((coord_elem_interp(0)-coord_elem_0(0))/delta_i(0));
              double yfact=fabs((coord_elem_interp(1)-coord_elem_0(1))/delta_i(1));
              double zfact=fabs((coord_elem_interp(2)-coord_elem_0(2))/delta_i(2));

              resu(fa7)=(1-zfact)*((1-yfact)*((1-xfact)*(valeurs_champ(elem_voisins(0))) + xfact*(valeurs_champ(elem_voisins(1)))) +
                                   yfact*((1-xfact)*(valeurs_champ(elem_voisins(2))) + xfact*(valeurs_champ(elem_voisins(3))))) +
                        zfact*((1-yfact)*((1-xfact)*(valeurs_champ(elem_voisins(4))) + xfact*(valeurs_champ(elem_voisins(5)))) +
                               yfact*((1-xfact)*(valeurs_champ(elem_voisins(6))) + xfact*(valeurs_champ(elem_voisins(7)))));
            }
          else
            {
              for (int dim=0; dim<dimension; dim++)
                {
                  resu(fa7)=-1e15; // de cette maniere, on ne calcule pas la force pour la fa7 pour laquelle on n'a pas acces a P2
                }
            }
        }
    }
  if (sauv_list_P1_all)
    {
      list_elem_P1_all.set_smart_resize(0);
      list_elem_P1_all.resize(nb_elem_P1_all,2);
    }
  return 1;
}

// EB
/*! @brief calcul des efforts hydrodynamiques a l'interface fluide - particule. Les composantes de frottements et de pression sont
 * calculees pour chaque facette du maillage lagrangien. Ces composantes peuvent etre ecrites dans les fichiers latas (en Newton, on a
 * multiplie l'effort locale exerce sur la facette par la surface de la facette). Les valeurs integrales, par particule, sont egalement imprimees dans les fichiers.out
 * Valable uniquement en 3D.
 */
void Navier_Stokes_FT_Disc::calcul_forces_interface()
{
  /***************************************************************
  *	  	  	  	  	  	RECUPERATION DES GRANDEURS		         *
  ***************************************************************/

  Cerr << "Navier_Stokes_FT_Disc::calcul_forces_interface"  <<  finl;
  if (dimension!=3) Process::exit("Navier_Stokes_FT_Disc::calcul_forces_interface code uniquement pour un cas 3D.");
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  const Domaine& domaine = domaine_vdf.domaine();

  DoubleVect& surface_tot_interf=variables_internes().surface_tot_interf_;
  DoubleTab& force_pression_tot_interf=variables_internes().force_pression_tot_interf_;
  DoubleTab& force_frottements_tot_interf=variables_internes().force_frottements_tot_interf_;

  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface_pour_post();

  const int nb_fa7 = maillage.nb_facettes();
  int nb_fa7_reelle=0;
  for (int i=0; i<nb_fa7; i++)
    if (!maillage.facette_virtuelle(i)) nb_fa7_reelle++;
  IntVect compo_connexes_fa7(nb_fa7); // Init a zero
  int n = search_connex_components_local_FT(maillage, compo_connexes_fa7);
  int nb_compo_tot=compute_global_connex_components_FT(maillage, compo_connexes_fa7, n);
  const Fluide_Diphasique& mon_fluide = fluide_diphasique();
  double mu_f=mon_fluide.fluide_phase(1).viscosite_dynamique().valeurs()(0, 0);

  surface_tot_interf.resize(nb_compo_tot);
  force_pression_tot_interf.resize(nb_compo_tot,dimension);
  force_frottements_tot_interf.resize(nb_compo_tot,dimension);
  surface_tot_interf=0;
  force_pression_tot_interf=0;
  force_frottements_tot_interf=0;
  variables_internes().Prop_P2_fluide_compo_.resize(nb_compo_tot);

  const Postraitement_Forces_Interfaces_FT& les_post_interf=eq_transport.postraitement_forces_interf();
  const double& distance_interpolation_pression_P1=les_post_interf.get_distance_interpolation_pression_P1();
  const double& distance_interpolation_pression_P2=les_post_interf.get_distance_interpolation_pression_P2();
  const double& distance_interpolation_gradU_P1=les_post_interf.get_distance_interpolation_gradU_P1();
  const double& distance_interpolation_gradU_P2=les_post_interf.get_distance_interpolation_gradU_P2();
  const DoubleTab& indicatrice = refeq_transport.valeur().get_update_indicatrice().valeurs();
  const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();

  DoubleTab& Nb_fa7_ok_prop = variables_internes().Proportion_fa7_ok_UP2_;
  Nb_fa7_ok_prop.resize(nb_compo_tot,dimension);
  Nb_fa7_ok_prop=0;
  IntTab Nb_fa7_tot_par_compo;
  Nb_fa7_tot_par_compo.resize(nb_compo_tot);
  Nb_fa7_tot_par_compo=0;
  DoubleTab& U_P2_moy = variables_internes().U_P2_moy_;
  U_P2_moy.resize(nb_compo_tot, dimension);
  U_P2_moy=0;
  DoubleTab& Prop_P2_fluide_compo=variables_internes().Prop_P2_fluide_compo_;
  if (nb_fa7>0)
    {
      variables_internes().Indic_elem_P2_.resize(nb_fa7);
      const ArrOfDouble& les_surfaces_fa7 = maillage.get_update_surface_facettes();
      const DoubleTab& les_normales_fa7 = maillage.get_update_normale_facettes();
      DoubleTab& pression_interf = variables_internes().pression_interf_;
      if (les_post_interf.flag_pression_facettes_)
        {
          pression_interf.resize(nb_fa7);
          pression_interf=3e15;
        }
      DoubleTab& force_frottements_interf = variables_internes().force_frottements_interf_;
      if (les_post_interf.flag_force_frottements_facettes_)
        {
          force_frottements_interf.resize(nb_fa7,dimension);
          force_frottements_interf=3e15;

        }
      DoubleTab& force_pression_interf = variables_internes().force_pression_interf_;
      if (les_post_interf.flag_force_pression_facettes_)
        {
          force_pression_interf.resize(nb_fa7,dimension);

          force_pression_interf=3e15;
        }

      DoubleTab& sigma_xx_interf = variables_internes().sigma_xx_interf_;
      DoubleTab& sigma_xy_interf = variables_internes().sigma_xy_interf_;
      DoubleTab& sigma_xz_interf = variables_internes().sigma_xz_interf_;
      DoubleTab& sigma_yx_interf = variables_internes().sigma_yx_interf_;
      DoubleTab& sigma_yy_interf = variables_internes().sigma_yy_interf_;
      DoubleTab& sigma_yz_interf = variables_internes().sigma_yz_interf_;
      DoubleTab& sigma_zx_interf = variables_internes().sigma_zx_interf_;
      DoubleTab& sigma_zy_interf = variables_internes().sigma_zy_interf_;
      DoubleTab& sigma_zz_interf = variables_internes().sigma_zz_interf_;
      DoubleTab& dUdx_P1 = variables_internes().dUdx_P1_;
      DoubleTab& dUdy_P1 = variables_internes().dUdy_P1_;
      DoubleTab& dUdz_P1 = variables_internes().dUdz_P1_;
      DoubleTab& dVdx_P1 = variables_internes().dVdx_P1_;
      DoubleTab& dVdy_P1 = variables_internes().dVdy_P1_;
      DoubleTab& dVdz_P1 = variables_internes().dVdz_P1_;
      DoubleTab& dWdx_P1 = variables_internes().dWdx_P1_;
      DoubleTab& dWdy_P1 = variables_internes().dWdy_P1_;
      DoubleTab& dWdz_P1 = variables_internes().dWdz_P1_;
      DoubleTab& dUdx_P2 = variables_internes().dUdx_P2_;
      DoubleTab& dUdy_P2 = variables_internes().dUdy_P2_;
      DoubleTab& dUdz_P2 = variables_internes().dUdz_P2_;
      DoubleTab& dVdx_P2 = variables_internes().dVdx_P2_;
      DoubleTab& dVdy_P2 = variables_internes().dVdy_P2_;
      DoubleTab& dVdz_P2 = variables_internes().dVdz_P2_;
      DoubleTab& dWdx_P2 = variables_internes().dWdx_P2_;
      DoubleTab& dWdy_P2 = variables_internes().dWdy_P2_;
      DoubleTab& dWdz_P2 = variables_internes().dWdz_P2_;

      if (les_post_interf.flag_tenseur_contraintes_facettes_)
        {
          sigma_xx_interf.resize(nb_fa7);
          sigma_xy_interf.resize(nb_fa7);
          sigma_xz_interf.resize(nb_fa7);
          sigma_yx_interf.resize(nb_fa7);
          sigma_yy_interf.resize(nb_fa7);
          sigma_yz_interf.resize(nb_fa7);
          sigma_zx_interf.resize(nb_fa7);
          sigma_zy_interf.resize(nb_fa7);
          sigma_zz_interf.resize(nb_fa7);

          dUdx_P1.resize(nb_fa7);
          dUdy_P1.resize(nb_fa7);
          dUdz_P1.resize(nb_fa7);
          dVdx_P1.resize(nb_fa7);
          dVdy_P1.resize(nb_fa7);
          dVdz_P1.resize(nb_fa7);
          dWdx_P1.resize(nb_fa7);
          dWdy_P1.resize(nb_fa7);
          dWdz_P1.resize(nb_fa7);

          dUdx_P2.resize(nb_fa7);
          dUdy_P2.resize(nb_fa7);
          dUdz_P2.resize(nb_fa7);
          dVdx_P2.resize(nb_fa7);
          dVdy_P2.resize(nb_fa7);
          dVdz_P2.resize(nb_fa7);
          dWdx_P2.resize(nb_fa7);
          dWdy_P2.resize(nb_fa7);
          dWdz_P2.resize(nb_fa7);
        }

      const DoubleTab& les_cg_fa7=maillage.cg_fa7();
      DoubleTab coord_voisin_fluide_fa7_pression_1(nb_fa7,dimension);
      DoubleTab coord_voisin_fluide_fa7_pression_2(nb_fa7,dimension);
      DoubleTab coord_voisin_fluide_fa7_gradU_1(nb_fa7,dimension);
      DoubleTab coord_voisin_fluide_fa7_gradU_2(nb_fa7,dimension);

      const int discr= variables_internes().ref_equation_mpoint_.non_nul() ? variables_internes().ref_equation_mpoint_.valeur().get_discretization_correction() : 0; // 0 si P1, 1 si elem_diph

      IntTab& list_elem_diph=variables_internes().list_elem_diph_; // EB
      if (discr==1) list_elem_diph.resize(nb_fa7,2); // EB (j,0) : num_elem_diph, (j,1) : num compo associee

      const DoubleVect& tab_num_compo=variables_internes().num_compo.valeur().valeurs();
      for (int fa7 =0 ; fa7<nb_fa7 ; fa7++)
        {
          if (!maillage.facette_virtuelle(fa7))
            {
              DoubleVect normale_fa7(dimension);
              int elem_diph=domaine.chercher_elements(les_cg_fa7(fa7,0), les_cg_fa7(fa7,1), les_cg_fa7(fa7,2));
              if (discr==1)
                {
                  list_elem_diph(fa7,0) = elem_diph;
                  list_elem_diph(fa7,1) = static_cast<int>(tab_num_compo(elem_diph));
                }
              if (elem_diph>=0)
                {
                  DoubleVect delta_i(dimension);
                  // On calcule les epaisseurs des mailles euleriennes  dans lesquelles se trouvent les facettes
                  // Si on y a acces, on prend l'epaisseur a l'exterieur de la particule
                  // Sinon, on prend l'epaisseur dans la particule
                  // Cela revient simplement a choisir la maille juxtaposee a la maille diphasique
                  int acces=1;
                  for (int dim=0; dim<dimension; dim++)
                    {
                      int elem_haut=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, dim+dimension),1);
                      int elem_bas=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, dim),0);
                      if (elem_bas<0 && elem_haut<0) acces=0;
                      if (les_normales_fa7(fa7,dim)>0) delta_i(dim) =  (elem_haut>=0) ? fabs(domaine_vdf.dist_elem(elem_diph,elem_haut, dim)) : fabs(domaine_vdf.dist_elem(elem_diph,elem_bas, dim));
                      else delta_i(dim) =  (elem_bas>=0) ? fabs(domaine_vdf.dist_elem(elem_diph,elem_bas, dim)) : fabs(domaine_vdf.dist_elem(elem_diph,elem_haut, dim)) ;
                    }
                  if (acces==0)
                    {
                      for (int dim=0; dim<dimension; dim++)
                        {
                          coord_voisin_fluide_fa7_pression_1(fa7,dim)=-1e15;
                          coord_voisin_fluide_fa7_pression_2(fa7,dim)=-1e15;
                          coord_voisin_fluide_fa7_gradU_1(fa7,dim)=-1e15;
                          coord_voisin_fluide_fa7_gradU_2(fa7,dim)=-1e15;
                        }
                    }
                  else
                    {
                      double epsilon=0;
                      for (int dim=0; dim<dimension; dim++)
                        {
                          epsilon+= fabs(delta_i(dim)*fabs(les_normales_fa7(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage
                        }
                      for (int dim=0; dim<dimension; dim++)
                        {
                          normale_fa7(dim)=les_normales_fa7(fa7,dim);
                          coord_voisin_fluide_fa7_pression_1(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_pression_P1*epsilon*normale_fa7(dim);
                          coord_voisin_fluide_fa7_pression_2(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_pression_P2*epsilon*normale_fa7(dim);
                          coord_voisin_fluide_fa7_gradU_1(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_gradU_P1*epsilon*normale_fa7(dim);
                          coord_voisin_fluide_fa7_gradU_2(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_gradU_P2*epsilon*normale_fa7(dim);
                        }
                    }
                }
              else
                {
                  for (int dim=0; dim<dimension; dim++)
                    {
                      coord_voisin_fluide_fa7_pression_1(fa7,dim)=-1e15;
                      coord_voisin_fluide_fa7_pression_2(fa7,dim)=-1e15;
                      coord_voisin_fluide_fa7_gradU_1(fa7,dim)=-1e15;
                      coord_voisin_fluide_fa7_gradU_2(fa7,dim)=-1e15;
                    }
                }
            }
        }

      // ----------------------------- Force de pression -----------------------------
      DoubleTab pression_P1(nb_fa7);
      DoubleTab pression_P2(nb_fa7);

      if (trilinear_interpolation_elem(indicatrice, la_pression.valeurs(), coord_voisin_fluide_fa7_pression_1, pression_P1, 0, discr) && trilinear_interpolation_elem(indicatrice, la_pression.valeurs(), coord_voisin_fluide_fa7_pression_2, pression_P2,1,0)) // soit on est capable d'interpoler en P1 et P1 et on calcule la force de p
        {
          DoubleTab pression_extrapolee(nb_fa7);
          Prop_P2_fluide_compo=0;
          for (int fa7=0; fa7<nb_fa7; fa7++)
            {
              if (!maillage.facette_virtuelle(fa7))
                {
                  int compo = compo_connexes_fa7(fa7);
                  if (variables_internes().Indic_elem_P2_(fa7)==1) Prop_P2_fluide_compo(compo)+=1;
                  if (pression_P2(fa7)<-1e10)
                    {
                      if (les_post_interf.flag_pression_facettes_) pression_interf(fa7)=1e15;
                      pression_extrapolee(fa7)=1e15;
                      continue;
                    }

                  pression_extrapolee(fa7) = pression_P2(fa7)-distance_interpolation_pression_P2*(pression_P2(fa7)-pression_P1(fa7))/(distance_interpolation_pression_P2-distance_interpolation_pression_P1); //3*pression_P2(compo,fa7)-2*pression_P1(compo,fa7); // Si on n'est pas capable d'interpoler en P1 ET P2, alors on ne calcule pas la force de pression
                  if (les_post_interf.flag_pression_facettes_) pression_interf(fa7)=pression_extrapolee(fa7);
                }
              else
                {
                  if (les_post_interf.flag_pression_facettes_) pression_interf(fa7)=1e15;
                  pression_extrapolee(fa7)=-1e15;
                }
            }
          for (int fa7=0; fa7<nb_fa7; fa7++)
            {
              int compo = compo_connexes_fa7(fa7);
              if (pression_extrapolee(fa7)>1e10) continue;
              double coeff=-pression_extrapolee(fa7)*les_surfaces_fa7(fa7);
              DoubleVect pressure_force_fa7(dimension);
              for (int dim=0; dim<dimension; dim++) pressure_force_fa7(dim)=coeff*les_normales_fa7(fa7,dim);
              if (!maillage.facette_virtuelle(fa7))
                {
                  surface_tot_interf(compo)+=les_surfaces_fa7(fa7);
                  if (les_post_interf.flag_force_pression_facettes_)
                    {
                      for (int dim=0; dim<dimension; dim++) force_pression_interf(fa7,dim)=pressure_force_fa7(dim);
                    }
                  for (int dim=0; dim<dimension; dim++) force_pression_tot_interf(compo,dim)+=pressure_force_fa7(dim);
                }
            }
        }
      else
        {
          for (int compo=0; compo<nb_compo_tot; compo++)
            {
              surface_tot_interf(compo)=1e18;
              for (int dim=0; dim<dimension; dim++) force_pression_tot_interf(compo,dim)+=1e18;
            }
        }
      // ----------------------------- Force de frottements -----------------------------
      if (les_post_interf.methode_calcul_force_frottements_==Postraitement_Forces_Interfaces_FT::TRILINEAIRE_LINEAIRE_TENSEUR_COMPLET)
        {
          DoubleTab grad_U_P1(nb_fa7, dimension, dimension);
          DoubleTab grad_U_P2(nb_fa7, dimension, dimension);
          grad_U_P1=-1e15;
          grad_U_P2=-1e30;
          int interp_gradU_P1_ok=0;
          int interp_gradU_P2_ok=0;
          if (les_post_interf.localisation_tenseur_contraintes_== Postraitement_Forces_Interfaces_FT::FACES_NORMALE_X)
            {
              interp_gradU_P1_ok=trilinear_interpolation_gradU_face(indicatrice_faces, la_vitesse.valeurs(), coord_voisin_fluide_fa7_gradU_1, grad_U_P1);
              interp_gradU_P2_ok=trilinear_interpolation_gradU_face(indicatrice_faces, la_vitesse.valeurs(), coord_voisin_fluide_fa7_gradU_2, grad_U_P2);
            }
          else if (les_post_interf.localisation_tenseur_contraintes_== Postraitement_Forces_Interfaces_FT::ELEMENTS)
            {
              interp_gradU_P1_ok=trilinear_interpolation_gradU_elem(indicatrice_faces, indicatrice, la_vitesse.valeurs(), coord_voisin_fluide_fa7_gradU_1, grad_U_P1);
              interp_gradU_P2_ok=trilinear_interpolation_gradU_elem(indicatrice_faces, indicatrice, la_vitesse.valeurs(), coord_voisin_fluide_fa7_gradU_2, grad_U_P2);
            }
          if ( interp_gradU_P1_ok==1 && interp_gradU_P2_ok==1 )
            {
              DoubleTab grad_U_extrapole(nb_fa7, dimension, dimension);
              grad_U_extrapole=1e20;

              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  dUdx_P1(fa7)=grad_U_P1(fa7,0,0);
                  dUdy_P1(fa7)=grad_U_P1(fa7,0,1);
                  dUdz_P1(fa7)=grad_U_P1(fa7,0,2);
                  dVdx_P1(fa7)=grad_U_P1(fa7,1,0);
                  dVdy_P1(fa7)=grad_U_P1(fa7,1,1);
                  dVdz_P1(fa7)=grad_U_P1(fa7,1,2);
                  dWdx_P1(fa7)=grad_U_P1(fa7,2,0);
                  dWdy_P1(fa7)=grad_U_P1(fa7,2,1);
                  dWdz_P1(fa7)=grad_U_P1(fa7,2,2);

                  dUdx_P2(fa7)=grad_U_P2(fa7,0,0);
                  dUdy_P2(fa7)=grad_U_P2(fa7,0,1);
                  dUdz_P2(fa7)=grad_U_P2(fa7,0,2);
                  dVdx_P2(fa7)=grad_U_P2(fa7,1,0);
                  dVdy_P2(fa7)=grad_U_P2(fa7,1,1);
                  dVdz_P2(fa7)=grad_U_P2(fa7,1,2);
                  dWdx_P2(fa7)=grad_U_P2(fa7,2,0);
                  dWdy_P2(fa7)=grad_U_P2(fa7,2,1);
                  dWdz_P2(fa7)=grad_U_P2(fa7,2,2);
                  for (int i=0; i<dimension; i++)
                    {
                      for (int j=0; j<dimension; j++)
                        {
                          grad_U_extrapole(fa7,i,j) = grad_U_P2(fa7,i,j)-distance_interpolation_gradU_P2*(grad_U_P2(fa7,i,j)-grad_U_P1(fa7,i,j))/(distance_interpolation_gradU_P2-distance_interpolation_gradU_P1);
                        }
                    }
                }
              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  int compo=compo_connexes_fa7(fa7);
                  Matrice_Dense tenseur_contrainte(dimension,dimension);
                  for (int i=0; i<dimension; i++)
                    {
                      for (int j=0; j<dimension; j++)
                        {
                          tenseur_contrainte(i,j) = (grad_U_extrapole(fa7,i,j) + grad_U_extrapole(fa7,j,i));
                        }
                    }

                  if (les_post_interf.flag_tenseur_contraintes_facettes_)
                    {
                      sigma_xx_interf(fa7)=tenseur_contrainte(0,0);
                      sigma_xy_interf(fa7)=tenseur_contrainte(0,1);
                      sigma_xz_interf(fa7)=tenseur_contrainte(0,2);
                      sigma_yx_interf(fa7)=tenseur_contrainte(1,0);
                      sigma_yy_interf(fa7)=tenseur_contrainte(1,1);
                      sigma_yz_interf(fa7)=tenseur_contrainte(1,2);
                      sigma_zx_interf(fa7)=tenseur_contrainte(2,0);
                      sigma_zy_interf(fa7)=tenseur_contrainte(2,1);
                      sigma_zz_interf(fa7)=tenseur_contrainte(2,2);
                    }
                  DoubleTab la_normale_fa7_x_surface(dimension);
                  for (int dim=0; dim<dimension; dim++) la_normale_fa7_x_surface(dim) =les_surfaces_fa7(fa7)*les_normales_fa7(fa7,dim);
                  DoubleVect friction_force_fa7=tenseur_contrainte*la_normale_fa7_x_surface;

                  if (!maillage.facette_virtuelle(fa7))
                    {
                      if (les_post_interf.flag_force_frottements_facettes_)
                        {
                          for (int dim=0; dim<dimension; dim++) force_frottements_interf(fa7,dim)=friction_force_fa7(dim);
                        }
                      for (int dim=0; dim<dimension; dim++) force_frottements_tot_interf(compo,dim)+=friction_force_fa7(dim);
                    }
                }

            }
          else
            {
              for (int compo=0; compo<nb_compo_tot; compo++)
                {
                  for (int dim=0; dim<dimension; dim++) force_frottements_tot_interf(compo,dim)+=0;
                }
            }
        }
      else if (les_post_interf.methode_calcul_force_frottements_==Postraitement_Forces_Interfaces_FT::TRILINENAIRE_TENSEUR_PROJETE) // Voir Butaye et al., Computers and Fluids, 2023.
        {
          int interp_U_P1_ok=0;
          int interp_U_P2_ok=0;
          DoubleTab& U_P1 = variables_internes().U_P1_;
          U_P1.resize(nb_fa7, dimension);
          DoubleTab& U_P2 = variables_internes().U_P2_;
          U_P2.resize(nb_fa7, dimension);

          DoubleTab U_P1_spherique(nb_fa7, dimension);
          DoubleTab U_P2_spherique(nb_fa7, dimension);
          DoubleTab U_cg_spherique(nb_fa7, dimension);
          DoubleTab Urr(nb_fa7);
          DoubleTab Uthetar(nb_fa7);
          DoubleTab Uphir(nb_fa7);

          U_P1=-1e15;
          U_P2=-1e30;
          U_P1_spherique=-1e15;
          U_P2_spherique=-1e30;
          U_cg_spherique=-1e20;
          Urr=1e8;
          Uthetar=1e12;
          Uphir=1e15;

          double theta=0;
          double phi=0;
          double distance_au_cg=0;
          const DoubleTab& positions_compo=eq_transport.get_positions_compo();
          const DoubleTab& vitesses_compo = eq_transport.get_vitesses_compo();

          //Cerr << "debut calcul composante frottements" << finl;
          if (les_post_interf.localisation_tenseur_contraintes_== Postraitement_Forces_Interfaces_FT::ELEMENTS)
            {
              // 1. On calcule (interpolation trilineaire) la vitesse en P1 et P2 en coord cartesiennes
              interp_U_P1_ok=trilinear_interpolation_face(indicatrice_faces, la_vitesse.valeurs(), coord_voisin_fluide_fa7_gradU_1, U_P1);
              interp_U_P2_ok=trilinear_interpolation_face(indicatrice_faces, la_vitesse.valeurs(), coord_voisin_fluide_fa7_gradU_2, U_P2);
            }
          if ( interp_U_P1_ok && interp_U_P2_ok )
            {
              // 2. On passe en coordonnees spheriques

              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  int compo=compo_connexes_fa7(fa7);
                  if (!maillage.facette_virtuelle(fa7))
                    {
                      Nb_fa7_tot_par_compo(compo)++;

                      int cont=0;
                      for (int dim=0; dim<dimension; dim++)
                        {
                          if (U_P2(fa7,dim)>-1e10)
                            {
                              //Cerr << "avant ajout U_P2_moy" << finl;
                              U_P2_moy(compo,dim)+= U_P2(fa7,dim); // calcul de la vitesse moyenne en P2
                              //Cerr << "apres" << finl;
                              Nb_fa7_ok_prop(compo,dim)+=1;
                            }
                          else cont=1;
                        }
                      if (cont) continue;
                      DoubleVect distance_cg_vect(dimension);
                      for (int i=0; i<dimension; i++) distance_cg_vect(i)=coord_voisin_fluide_fa7_gradU_1(fa7,i)-positions_compo(compo,i);

                      distance_au_cg=sqrt(local_carre_norme_vect(distance_cg_vect));
                      if (fabs((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg)<=1) theta=acos((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg);
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg>1) theta=0;
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg<-1) theta=M_PI;
                      if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))>0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))>=0)
                        {
                          phi=atan((coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))/(coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0)));
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))>0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))<0)
                        {
                          phi=atan((coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))/(coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0)))+2*M_PI;
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))<0)
                        {
                          phi=atan((coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))/(coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0)))+M_PI;
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))==0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))>0)
                        {
                          phi=M_PI/2.;
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))==0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))<0)
                        {
                          phi=3.*M_PI/2.;
                        }

                      U_P1_spherique(fa7,0)=sin(theta)*cos(phi)*U_P1(fa7,0)+sin(theta)*sin(phi)*U_P1(fa7,1)+cos(theta)*U_P1(fa7,2);
                      U_P1_spherique(fa7,1)=cos(theta)*cos(phi)*U_P1(fa7,0)+cos(theta)*sin(phi)*U_P1(fa7,1)-sin(theta)*U_P1(fa7,2);
                      U_P1_spherique(fa7,2)=sin(phi)*U_P1(fa7,0)+cos(phi)*U_P1(fa7,1);

                      U_P2_spherique(fa7,0)=sin(theta)*cos(phi)*U_P2(fa7,0)+sin(theta)*sin(phi)*U_P2(fa7,1)+cos(theta)*U_P2(fa7,2);
                      U_P2_spherique(fa7,1)=cos(theta)*cos(phi)*U_P2(fa7,0)+cos(theta)*sin(phi)*U_P2(fa7,1)-sin(theta)*U_P2(fa7,2);
                      U_P2_spherique(fa7,2)=sin(phi)*U_P2(fa7,0)+cos(phi)*U_P2(fa7,1);

                      U_cg_spherique(fa7,0)=sin(theta)*cos(phi)*vitesses_compo(compo,0)+sin(theta)*sin(phi)*vitesses_compo(compo,1)+cos(theta)*vitesses_compo(compo,2);
                      U_cg_spherique(fa7,1)=cos(theta)*cos(phi)*vitesses_compo(compo,0)+cos(theta)*sin(phi)*vitesses_compo(compo,1)-sin(theta)*vitesses_compo(compo,2);
                      U_cg_spherique(fa7,2)=sin(phi)*vitesses_compo(compo,0)+cos(phi)*vitesses_compo(compo,1);

                      // On recalcule delta --> epsilon
                      int elem_diph=domaine.chercher_elements(les_cg_fa7(fa7,0), les_cg_fa7(fa7,1),les_cg_fa7(fa7,2));
                      DoubleVect delta_i(dimension);
                      int elem00=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, 0+dimension),1);
                      int elem11=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, 1+dimension),1);
                      int elem22=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, 2+dimension),1);
                      int elem33=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, 2),0);

                      delta_i(0) = elem00>=0 ? fabs(domaine_vdf.dist_elem(elem_diph, elem00, 0)) : -1e15;
                      delta_i(1) = elem11>=0 ? fabs(domaine_vdf.dist_elem(elem_diph, elem11, 1)) : -1e15;

                      if (les_normales_fa7(fa7,2)>0) delta_i(2) = elem22>=0 ? fabs(domaine_vdf.dist_elem(elem_diph, elem22, 2)) : -1e15;
                      else delta_i(2) = elem33>=0 ? fabs(domaine_vdf.dist_elem(elem_diph, elem33 , 2)) : -1e15;
                      double epsilon=0;
                      for (int dim=0; dim<dimension; dim++) epsilon+= fabs(delta_i(dim)*fabs(les_normales_fa7(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage
                      // On calcule les composantes de la force de frottements en coord spherique (apres simplifications) : ff=mu*(2*Urr, Uthetar, Uphir)
                      Urr(fa7)=(-U_P2_spherique(fa7,0)+4.*U_P1_spherique(fa7,0)-3.*U_cg_spherique(fa7,0))/(2.*epsilon);
                      Uthetar(fa7)=(-U_P2_spherique(fa7,1)+4.*U_P1_spherique(fa7,1)-3.*U_cg_spherique(fa7,1))/(2.*epsilon);
                      Uphir(fa7)=(-U_P2_spherique(fa7,2)+4.*U_P1_spherique(fa7,2)-3.*U_cg_spherique(fa7,2))/(2.*epsilon);

                      DoubleVect ff(dimension);
                      ff(0)=mu_f*les_surfaces_fa7(fa7)*(2.*sin(theta)*cos(phi)*Urr(fa7)+cos(theta)*cos(phi)*Uthetar(fa7)-sin(phi)*Uphir(fa7));
                      ff(1)=mu_f*les_surfaces_fa7(fa7)*(2.*sin(theta)*sin(phi)*Urr(fa7)+cos(theta)*sin(phi)*Uthetar(fa7)+cos(phi)*Uphir(fa7));
                      ff(2)=mu_f*les_surfaces_fa7(fa7)*(2.*cos(theta)*Urr(fa7)-sin(theta)*Uthetar(fa7));

                      for (int dim=0; dim<dimension; dim++) force_frottements_tot_interf(compo,dim)+=ff(dim);

                      if (les_post_interf.flag_force_frottements_facettes_)
                        for (int dim=0; dim<dimension; dim++) force_frottements_interf(fa7,dim)=ff(dim);
                    }
                }

            }
          else
            {
              for (int compo=0; compo<nb_compo_tot; compo++)
                {
                  for (int dim=0; dim<dimension; dim++) force_frottements_tot_interf(compo,dim)+=0;
                }
            }
        }
    }
  mp_sum_for_each_item(force_pression_tot_interf);
  mp_sum_for_each_item(surface_tot_interf);
  mp_sum_for_each_item(force_frottements_tot_interf);
  mp_sum_for_each_item(Nb_fa7_ok_prop);
  mp_sum_for_each_item(U_P2_moy);
  mp_sum_for_each_item(Prop_P2_fluide_compo);

  for (int compo=0; compo<nb_compo_tot; compo++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          if (Nb_fa7_ok_prop(compo,dim)>0) U_P2_moy(compo,dim)/=Nb_fa7_ok_prop(compo,dim);
        }
    }
  mp_sum_for_each_item(Nb_fa7_tot_par_compo);
  for (int compo=0; compo<nb_compo_tot; compo++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          if (Nb_fa7_tot_par_compo(compo)>0) Nb_fa7_ok_prop(compo,dim)/=Nb_fa7_tot_par_compo(compo);
        }
      if (Nb_fa7_tot_par_compo(compo)>0) Prop_P2_fluide_compo(compo)/=Nb_fa7_tot_par_compo(compo);

    }
  if (les_post_interf.calcul_forces_theoriques_stokes_ && schema_temps().nb_pas_dt()==0) calcul_forces_interface_stokes_th();
}
// EB
/*! @brief idem que Navier_Stokes_FT_Disc::calcul_forces_interface_stokes mais utilise les champs theoriques de vitesse et de pression
 * discretisee, respectivement, aux faces et aux elements du maillage eulerien. Valable pour une simulation avec une seule particule,
 * centre en 0 dans le domaine. Le calcul n'est effectue qu'au premier pas de temps. Les resultats sont disponibles dans les latas (pour chaque facette),
 * et dans les .out (grandeurs integrales). Les efforts sont calcules a partir d'une discretisation des champs theoriques et l'utilisation de la methode de calcul --> th_dis.
 * On calcule egalement les efforts theoriques exerces sur chaque facette --> th (Le suffixe "dis" est mal choisi car "th" est aussi discretise aux facettes mais...trop tard).
 * En effet, on connait la solution theorique des contraintes exerces sur l'interface en fonction de theta. Ces contraintes sont donnees par la solution de Stokes.
 * Comparer th_dis et th_dis permet de calculer l'erreur de la methode, en prenant en compte la meme erreur de discretisation. En effet,
 * dans les deux cas, on calcule les grandeurs locales sur la surface des facettes.
 * Valable uniquement en 3D.
 */
void Navier_Stokes_FT_Disc::calcul_forces_interface_stokes_th()
{

  /***************************************************************
  *	  	  	  	  	  	RECUPERATION DES GRANDEURS		         *
  ***************************************************************/
  Cerr << "Navier_Stokes_FT_Disc::calcul_forces_interface_stokes_th" << finl;
  REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
  const Domaine& domaine = domaine_vdf.domaine();

  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface_pour_post();

  const int nb_fa7 = maillage.nb_facettes();
  int nb_fa7_reelle=0;
  for (int i=0; i<nb_fa7; i++)
    if (!maillage.facette_virtuelle(i)) nb_fa7_reelle++;
  IntVect compo_connexes_fa7(nb_fa7); // Init a zero
  int n = search_connex_components_local_FT(maillage, compo_connexes_fa7);
  int nb_compo_tot=compute_global_connex_components_FT(maillage, compo_connexes_fa7, n);
  const Fluide_Diphasique& mon_fluide = fluide_diphasique();
  const Particule_Solide& particule_solide=ref_cast(Particule_Solide,mon_fluide.fluide_phase(0));
  const double rayon_particule =eq_transport.get_rayons_compo()(0); // On ne se sert de cette fonction que lorsqu'il y a une seule particule dans le domaine
  double mu_p=particule_solide.viscosite_dynamique().valeurs()(0, 0);
  double mu_f=mon_fluide.fluide_phase(1).viscosite_dynamique().valeurs()(0, 0);
  double phi_mu=mu_p/mu_f;

  DoubleTab& force_pression_tot_interf_stokes_th=variables_internes().force_pression_tot_interf_stokes_th_;
  DoubleTab& force_frottements_tot_interf_stokes_th=variables_internes().force_frottements_tot_interf_stokes_th_;
  DoubleTab& force_pression_tot_interf_stokes_th_dis=variables_internes().force_pression_tot_interf_stokes_th_dis_;
  DoubleTab& force_frottements_tot_interf_stokes_th_dis=variables_internes().force_frottements_tot_interf_stokes_th_dis_;

  force_pression_tot_interf_stokes_th.resize(nb_compo_tot,dimension);
  force_frottements_tot_interf_stokes_th.resize(nb_compo_tot,dimension);
  force_pression_tot_interf_stokes_th_dis.resize(nb_compo_tot,dimension);
  force_frottements_tot_interf_stokes_th_dis.resize(nb_compo_tot,dimension);
  force_pression_tot_interf_stokes_th=0;
  force_frottements_tot_interf_stokes_th=0;
  force_pression_tot_interf_stokes_th_dis=0;
  force_frottements_tot_interf_stokes_th_dis=0;

  const Postraitement_Forces_Interfaces_FT& les_post_interf=eq_transport.postraitement_forces_interf();
  const double& distance_interpolation_pression_P1=les_post_interf.get_distance_interpolation_pression_P1();
  const double& distance_interpolation_pression_P2=les_post_interf.get_distance_interpolation_pression_P2();
  const double& distance_interpolation_gradU_P1=les_post_interf.get_distance_interpolation_gradU_P1();
  const double& distance_interpolation_gradU_P2=les_post_interf.get_distance_interpolation_gradU_P2();
  DoubleTab& gravite = milieu().gravite().valeurs();
  DoubleVect vect_g(dimension);
  for (int i=0; i<dimension; i++) vect_g(i)=gravite(0,i);
  const double norme_g = sqrt(local_carre_norme_vect(vect_g));
  const double rho_f = mon_fluide.fluide_phase(1).masse_volumique().valeurs()(0, 0);
  const double rho_p = particule_solide.masse_volumique().valeurs()(0, 0);
  const double& v_inf_stokes=2./3.*(pow(rayon_particule,2)*norme_g/mu_f)*((1+phi_mu)/(2.+3.*phi_mu)*(rho_f-rho_p));
  const DoubleTab& indicatrice = refeq_transport.valeur().get_update_indicatrice().valeurs();
  const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();

  // Idealement il faudrait plutot mettre le calcul des champs theoriques de stokes sur le champs eulerien dans une fonction inclue a la fin de derivee_en_temps_inco mais comme on va quasi jamais s'en servir...
  DoubleTab& vitesse_stokes_th=variables_internes().vitesse_stokes_th_.valeur().valeurs();
  vitesse_stokes_th=0;
  //vitesse_stokes_th.resize(la_vitesse.valeurs().dimension_tot(0),dimension);

  for (int num_face=0; num_face<vitesse_stokes_th.dimension_tot(0); num_face++)
    {
      double x=domaine_vdf.xv(num_face,0);
      double y=domaine_vdf.xv(num_face,1);
      double z=domaine_vdf.xv(num_face,2);
      if (sqrt(x*x+y*y+z*z)>=rayon_particule)
        {
          if (domaine_vf.orientation(num_face)==0) vitesse_stokes_th(num_face)=(v_inf_stokes*z*x*((2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(4*(pow(x*x+y*y+z*z,1.5)))-3*phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,2.5))));
          if (domaine_vf.orientation(num_face)==1) vitesse_stokes_th(num_face)=(v_inf_stokes*z*y*((2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(4*(pow(x*x+y*y+z*z,1.5)))-3*phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,2.5))));
          if (domaine_vf.orientation(num_face)==2) vitesse_stokes_th(num_face)=-v_inf_stokes*(+(z*z)*(1/(x*x+y*y+z*z)-(2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(2*pow(x*x+y*y+z*z,1.5))+phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(2*pow(x*x+y*y+z*z,2.5)))+(1-z*z/(x*x+y*y+z*z))*(1-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*sqrt(x*x+y*y+z*z))-phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,1.5))));
        }
      else
        {
          if (domaine_vf.orientation(num_face)==0) vitesse_stokes_th(num_face)=v_inf_stokes/(2*pow(rayon_particule,2))*(1/(1+phi_mu))*x*z;
          if (domaine_vf.orientation(num_face)==1) vitesse_stokes_th(num_face)=v_inf_stokes/(2*pow(rayon_particule,2))*(1/(1+phi_mu))*y*z;
          if (domaine_vf.orientation(num_face)==2) vitesse_stokes_th(num_face)=v_inf_stokes/(2*(1+phi_mu))*(1+z*z/(pow(rayon_particule,2))-2*(x*x+y*y+z*z)/(pow(rayon_particule,2)));
        }
    }
  vitesse_stokes_th.echange_espace_virtuel();

  DoubleTab& pression_stokes_th=variables_internes().pression_stokes_th_.valeur().valeurs();

  for (int num_elem=0; num_elem<pression_stokes_th.dimension_tot(0); num_elem++)
    {
      double x=domaine_vdf.xp(num_elem,0);
      double y=domaine_vdf.xp(num_elem,1);
      double z=domaine_vdf.xp(num_elem,2);
      if (sqrt(x*x+y*y+z*z)>=rayon_particule)
        {
          pression_stokes_th(num_elem)=mu_f*v_inf_stokes*((2+3*phi_mu)/(1+phi_mu))*1/2*rayon_particule*z/(pow(x*x+y*y+z*z,1.5));
        }
      else
        {
          pression_stokes_th(num_elem)=mu_p*v_inf_stokes*5/((1+phi_mu)*(pow(rayon_particule,2)))*z;
        }
    }

  pression_stokes_th.echange_espace_virtuel();

  if (nb_fa7>0)
    {
      const ArrOfDouble& les_surfaces_fa7 = maillage.get_update_surface_facettes();
      const DoubleTab& les_normales_fa7 = maillage.get_update_normale_facettes();
      DoubleTab& pression_interf_stokes_th_dis = variables_internes().pression_interf_stokes_th_dis_;
      if (les_post_interf.flag_pression_facettes_)
        {
          pression_interf_stokes_th_dis.resize(nb_fa7);
          pression_interf_stokes_th_dis=3e15;
        }
      DoubleTab& force_pression_th = variables_internes().force_pression_stokes_th_;
      force_pression_th.resize(nb_fa7,dimension);
      force_pression_th=3e15;

      DoubleTab& force_frottements_th = variables_internes().force_frottements_stokes_th_;
      force_frottements_th.resize(nb_fa7,dimension);
      force_frottements_th=3e15;

      DoubleTab& force_pression_th_dis = variables_internes().force_pression_stokes_th_dis_;
      force_pression_th_dis.resize(nb_fa7,dimension);
      force_pression_th_dis=3e15;

      DoubleTab& force_frottements_th_dis = variables_internes().force_frottements_stokes_th_dis_;
      force_frottements_th_dis.resize(nb_fa7,dimension);
      force_frottements_th_dis=3e15;

      const DoubleTab& la_vitesse_stokes_th=variables_internes().vitesse_stokes_th_.valeurs();
      const DoubleTab& la_pression_stokes_th=variables_internes().pression_stokes_th_.valeurs();

      DoubleTab& sigma_xx_interf_stokes_th_dis = variables_internes().sigma_xx_interf_stokes_th_dis_;
      DoubleTab& sigma_xy_interf_stokes_th_dis = variables_internes().sigma_xy_interf_stokes_th_dis_;
      DoubleTab& sigma_xz_interf_stokes_th_dis = variables_internes().sigma_xz_interf_stokes_th_dis_;
      DoubleTab& sigma_yx_interf_stokes_th_dis = variables_internes().sigma_yx_interf_stokes_th_dis_;
      DoubleTab& sigma_yy_interf_stokes_th_dis = variables_internes().sigma_yy_interf_stokes_th_dis_;
      DoubleTab& sigma_yz_interf_stokes_th_dis = variables_internes().sigma_yz_interf_stokes_th_dis_;
      DoubleTab& sigma_zx_interf_stokes_th_dis = variables_internes().sigma_zx_interf_stokes_th_dis_;
      DoubleTab& sigma_zy_interf_stokes_th_dis = variables_internes().sigma_zy_interf_stokes_th_dis_;
      DoubleTab& sigma_zz_interf_stokes_th_dis = variables_internes().sigma_zz_interf_stokes_th_dis_;
      DoubleTab& sigma_xx_interf_stokes_th = variables_internes().sigma_xx_interf_stokes_th_;
      DoubleTab& sigma_xy_interf_stokes_th = variables_internes().sigma_xy_interf_stokes_th_;
      DoubleTab& sigma_xz_interf_stokes_th = variables_internes().sigma_xz_interf_stokes_th_;
      DoubleTab& sigma_yy_interf_stokes_th = variables_internes().sigma_yy_interf_stokes_th_;
      DoubleTab& sigma_yz_interf_stokes_th = variables_internes().sigma_yz_interf_stokes_th_;
      DoubleTab& sigma_zz_interf_stokes_th = variables_internes().sigma_zz_interf_stokes_th_;
      DoubleTab& dUdx_P1_th_dis = variables_internes().dUdx_P1_th_dis_;
      DoubleTab& dUdy_P1_th_dis = variables_internes().dUdy_P1_th_dis_;
      DoubleTab& dUdz_P1_th_dis = variables_internes().dUdz_P1_th_dis_;
      DoubleTab& dVdx_P1_th_dis = variables_internes().dVdx_P1_th_dis_;
      DoubleTab& dVdy_P1_th_dis = variables_internes().dVdy_P1_th_dis_;
      DoubleTab& dVdz_P1_th_dis = variables_internes().dVdz_P1_th_dis_;
      DoubleTab& dWdx_P1_th_dis = variables_internes().dWdx_P1_th_dis_;
      DoubleTab& dWdy_P1_th_dis = variables_internes().dWdy_P1_th_dis_;
      DoubleTab& dWdz_P1_th_dis = variables_internes().dWdz_P1_th_dis_;
      DoubleTab& dUdx_P2_th_dis = variables_internes().dUdx_P2_th_dis_;
      DoubleTab& dUdy_P2_th_dis = variables_internes().dUdy_P2_th_dis_;
      DoubleTab& dUdz_P2_th_dis = variables_internes().dUdz_P2_th_dis_;
      DoubleTab& dVdx_P2_th_dis = variables_internes().dVdx_P2_th_dis_;
      DoubleTab& dVdy_P2_th_dis = variables_internes().dVdy_P2_th_dis_;
      DoubleTab& dVdz_P2_th_dis = variables_internes().dVdz_P2_th_dis_;
      DoubleTab& dWdx_P2_th_dis = variables_internes().dWdx_P2_th_dis_;
      DoubleTab& dWdy_P2_th_dis = variables_internes().dWdy_P2_th_dis_;
      DoubleTab& dWdz_P2_th_dis = variables_internes().dWdz_P2_th_dis_;
      DoubleTab& dUdx_P1_th = variables_internes().dUdx_P1_th_;
      DoubleTab& dUdy_P1_th = variables_internes().dUdy_P1_th_;
      DoubleTab& dUdz_P1_th = variables_internes().dUdz_P1_th_;
      DoubleTab& dVdx_P1_th = variables_internes().dVdx_P1_th_;
      DoubleTab& dVdy_P1_th = variables_internes().dVdy_P1_th_;
      DoubleTab& dVdz_P1_th = variables_internes().dVdz_P1_th_;
      DoubleTab& dWdx_P1_th = variables_internes().dWdx_P1_th_;
      DoubleTab& dWdy_P1_th = variables_internes().dWdy_P1_th_;
      DoubleTab& dWdz_P1_th = variables_internes().dWdz_P1_th_;
      DoubleTab& dUdx_P2_th = variables_internes().dUdx_P2_th_;
      DoubleTab& dUdy_P2_th = variables_internes().dUdy_P2_th_;
      DoubleTab& dUdz_P2_th = variables_internes().dUdz_P2_th_;
      DoubleTab& dVdx_P2_th = variables_internes().dVdx_P2_th_;
      DoubleTab& dVdy_P2_th = variables_internes().dVdy_P2_th_;
      DoubleTab& dVdz_P2_th = variables_internes().dVdz_P2_th_;
      DoubleTab& dWdx_P2_th = variables_internes().dWdx_P2_th_;
      DoubleTab& dWdy_P2_th = variables_internes().dWdy_P2_th_;
      DoubleTab& dWdz_P2_th = variables_internes().dWdz_P2_th_;

      if (les_post_interf.flag_tenseur_contraintes_facettes_)
        {
          sigma_xx_interf_stokes_th_dis.resize(nb_fa7);
          sigma_xy_interf_stokes_th_dis.resize(nb_fa7);
          sigma_xz_interf_stokes_th_dis.resize(nb_fa7);
          sigma_yx_interf_stokes_th_dis.resize(nb_fa7);
          sigma_yy_interf_stokes_th_dis.resize(nb_fa7);
          sigma_yz_interf_stokes_th_dis.resize(nb_fa7);
          sigma_zx_interf_stokes_th_dis.resize(nb_fa7);
          sigma_zy_interf_stokes_th_dis.resize(nb_fa7);
          sigma_zz_interf_stokes_th_dis.resize(nb_fa7);

          sigma_xx_interf_stokes_th.resize(nb_fa7);
          sigma_xy_interf_stokes_th.resize(nb_fa7);
          sigma_xz_interf_stokes_th.resize(nb_fa7);
          sigma_yy_interf_stokes_th.resize(nb_fa7);
          sigma_yz_interf_stokes_th.resize(nb_fa7);
          sigma_zz_interf_stokes_th.resize(nb_fa7);

          dUdx_P1_th_dis.resize(nb_fa7);
          dUdy_P1_th_dis.resize(nb_fa7);
          dUdz_P1_th_dis.resize(nb_fa7);
          dVdx_P1_th_dis.resize(nb_fa7);
          dVdy_P1_th_dis.resize(nb_fa7);
          dVdz_P1_th_dis.resize(nb_fa7);
          dWdx_P1_th_dis.resize(nb_fa7);
          dWdy_P1_th_dis.resize(nb_fa7);
          dWdz_P1_th_dis.resize(nb_fa7);

          dUdx_P2_th_dis.resize(nb_fa7);
          dUdy_P2_th_dis.resize(nb_fa7);
          dUdz_P2_th_dis.resize(nb_fa7);
          dVdx_P2_th_dis.resize(nb_fa7);
          dVdy_P2_th_dis.resize(nb_fa7);
          dVdz_P2_th_dis.resize(nb_fa7);
          dWdx_P2_th_dis.resize(nb_fa7);
          dWdy_P2_th_dis.resize(nb_fa7);
          dWdz_P2_th_dis.resize(nb_fa7);

          dUdx_P1_th.resize(nb_fa7);
          dUdy_P1_th.resize(nb_fa7);
          dUdz_P1_th.resize(nb_fa7);
          dVdx_P1_th.resize(nb_fa7);
          dVdy_P1_th.resize(nb_fa7);
          dVdz_P1_th.resize(nb_fa7);
          dWdx_P1_th.resize(nb_fa7);
          dWdy_P1_th.resize(nb_fa7);
          dWdz_P1_th.resize(nb_fa7);

          dUdx_P2_th.resize(nb_fa7);
          dUdy_P2_th.resize(nb_fa7);
          dUdz_P2_th.resize(nb_fa7);
          dVdx_P2_th.resize(nb_fa7);
          dVdy_P2_th.resize(nb_fa7);
          dVdz_P2_th.resize(nb_fa7);
          dWdx_P2_th.resize(nb_fa7);
          dWdy_P2_th.resize(nb_fa7);
          dWdz_P2_th.resize(nb_fa7);
        }

      DoubleTab& U_P1_th=variables_internes().U_P1_th_;
      U_P1_th.resize(nb_fa7,dimension);
      DoubleTab& U_P2_th=variables_internes().U_P2_th_;
      U_P2_th.resize(nb_fa7,dimension);

      const DoubleTab& les_cg_fa7=maillage.cg_fa7();
      DoubleTab coord_voisin_fluide_fa7_pression_1(nb_fa7,dimension);
      DoubleTab coord_voisin_fluide_fa7_pression_2(nb_fa7,dimension);
      DoubleTab coord_voisin_fluide_fa7_gradU_1(nb_fa7,dimension);
      DoubleTab coord_voisin_fluide_fa7_gradU_2(nb_fa7,dimension);

      for (int fa7 =0 ; fa7<nb_fa7 ; fa7++)
        {
          if (!maillage.facette_virtuelle(fa7))
            {
              DoubleVect normale_fa7(dimension);
              int elem_diph=domaine.chercher_elements(les_cg_fa7(fa7,0), les_cg_fa7(fa7,1),les_cg_fa7(fa7,2));
              DoubleVect delta_i(dimension);

              // On calcule les epaisseurs des mailles euleriennes  dans lesquelles se trouvent les facettes
              // Si on y a acces, on prend l'epaisseur a l'exterieur de la particule
              // Sinon, on prend l'epaisseur dans la particule
              // Cela revient simplement a choisir la maille juxtaposee a la maille diphasique
              for (int dim=0; dim<dimension; dim++)
                {
                  int elem_haut=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, dim+dimension),1);
                  int elem_bas=domaine_vdf.face_voisins_pour_interp(domaine_vf.elem_faces_pour_interp(elem_diph, dim),0);
                  if (les_normales_fa7(fa7,dim)>0) delta_i(dim) =  (elem_haut>=0) ? fabs(domaine_vdf.dist_elem(elem_diph,elem_haut, dim)) : fabs(domaine_vdf.dist_elem(elem_diph,elem_bas, dim));
                  else delta_i(dim) =  (elem_bas>=0) ? fabs(domaine_vdf.dist_elem(elem_diph,elem_bas, dim)) : fabs(domaine_vdf.dist_elem(elem_diph,elem_haut, dim));
                }

              double epsilon=0;
              for (int dim=0; dim<dimension; dim++)
                {
                  epsilon+= fabs(delta_i(dim)*fabs(les_normales_fa7(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage
                }

              for (int dim=0; dim<dimension; dim++)
                {
                  normale_fa7(dim)=les_normales_fa7(fa7,dim);
                  coord_voisin_fluide_fa7_pression_1(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_pression_P1*epsilon*normale_fa7(dim);
                  coord_voisin_fluide_fa7_pression_2(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_pression_P2*epsilon*normale_fa7(dim);
                  coord_voisin_fluide_fa7_gradU_1(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_gradU_P1*epsilon*normale_fa7(dim);
                  coord_voisin_fluide_fa7_gradU_2(fa7,dim)=les_cg_fa7(fa7,dim)+distance_interpolation_gradU_P2*epsilon*normale_fa7(dim);
                }

              double x=coord_voisin_fluide_fa7_pression_1(fa7,0);
              double y=coord_voisin_fluide_fa7_pression_1(fa7,1);
              double z=coord_voisin_fluide_fa7_pression_1(fa7,2);
              U_P1_th(fa7,0) = (v_inf_stokes*z*x*((2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(4*(pow(x*x+y*y+z*z,1.5)))-3*phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,2.5))));
              U_P1_th(fa7,1) = (v_inf_stokes*z*y*((2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(4*(pow(x*x+y*y+z*z,1.5)))-3*phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,2.5))));
              U_P1_th(fa7,2) = -v_inf_stokes*(+(z*z)*(1/(x*x+y*y+z*z)-(2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(2*pow(x*x+y*y+z*z,1.5))+phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(2*pow(x*x+y*y+z*z,2.5)))+(1-z*z/(x*x+y*y+z*z))*(1-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*sqrt(x*x+y*y+z*z))-phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,1.5))));

              x=coord_voisin_fluide_fa7_pression_2(fa7,0);
              y=coord_voisin_fluide_fa7_pression_2(fa7,1);
              z=coord_voisin_fluide_fa7_pression_2(fa7,2);
              U_P2_th(fa7,0) = (v_inf_stokes*z*x*((2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(4*(pow(x*x+y*y+z*z,1.5)))-3*phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,2.5))));
              U_P2_th(fa7,1) = (v_inf_stokes*z*y*((2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(4*(pow(x*x+y*y+z*z,1.5)))-3*phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,2.5))));
              U_P2_th(fa7,2) = -v_inf_stokes*(+(z*z)*(1/(x*x+y*y+z*z)-(2+3*phi_mu)/(1+phi_mu)*(rayon_particule)/(2*pow(x*x+y*y+z*z,1.5))+phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(2*pow(x*x+y*y+z*z,2.5)))+(1-z*z/(x*x+y*y+z*z))*(1-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*sqrt(x*x+y*y+z*z))-phi_mu/(1+phi_mu)*pow(rayon_particule,3)/(4*pow(x*x+y*y+z*z,1.5))));

            }
        }

      // ----------------------------- Force de pression -----------------------------
      DoubleTab pression_P1(nb_fa7);
      DoubleTab pression_P2(nb_fa7);

      if (trilinear_interpolation_elem(indicatrice, la_pression_stokes_th, coord_voisin_fluide_fa7_pression_1, pression_P1) && trilinear_interpolation_elem(indicatrice, la_pression_stokes_th, coord_voisin_fluide_fa7_pression_2, pression_P2)) // soit on est capable d'interpoler en P1 et P1 et on calcule la force de p
        {
          DoubleTab pression_extrapolee(nb_fa7);

          for (int fa7=0; fa7<nb_fa7; fa7++)
            {
              if (!maillage.facette_virtuelle(fa7))
                {
                  pression_extrapolee(fa7) = pression_P2(fa7)-distance_interpolation_pression_P2*(pression_P2(fa7)-pression_P1(fa7))/(distance_interpolation_pression_P2-distance_interpolation_pression_P1); //3*pression_P2(compo,fa7)-2*pression_P1(compo,fa7); // Si on n'est pas capable d'interpoler en P1 ET P2, alors on ne calcule pas la force de pression
                  pression_interf_stokes_th_dis(fa7)=pression_extrapolee(fa7);
                }
              else
                {
                  pression_interf_stokes_th_dis(fa7)=1e15;
                  pression_extrapolee(fa7)=-1e15;
                }
            }


          for (int fa7=0; fa7<nb_fa7; fa7++)
            {
              int compo=compo_connexes_fa7(fa7);
              double coeff=-pression_extrapolee(fa7)*les_surfaces_fa7(fa7);
              DoubleVect pressure_force_fa7(dimension);
              for (int dim=0; dim<dimension; dim++) pressure_force_fa7(dim)=coeff*les_normales_fa7(fa7,dim);
              if (!maillage.facette_virtuelle(fa7))
                {

                  for (int dim=0; dim<dimension; dim++) force_pression_th_dis(fa7,dim)=pressure_force_fa7(dim);
                  force_pression_th(fa7,0)=-(mu_f*v_inf_stokes*(2.+3.*phi_mu)/(1.+phi_mu)*1./(2.*rayon_particule)*les_cg_fa7(fa7,2)/rayon_particule)*les_normales_fa7(fa7,0)*les_surfaces_fa7(fa7);
                  force_pression_th(fa7,1)=-(mu_f*v_inf_stokes*(2.+3.*phi_mu)/(1.+phi_mu)*1./(2.*rayon_particule)*les_cg_fa7(fa7,2)/rayon_particule)*les_normales_fa7(fa7,1)*les_surfaces_fa7(fa7);
                  force_pression_th(fa7,2)=-(mu_f*v_inf_stokes*(2.+3.*phi_mu)/(1.+phi_mu)*1./(2.*rayon_particule)*les_cg_fa7(fa7,2)/rayon_particule)*les_normales_fa7(fa7,2)*les_surfaces_fa7(fa7);

                  for (int dim=0; dim<dimension; dim++)
                    {
                      force_pression_tot_interf_stokes_th_dis(compo,dim)+=pressure_force_fa7(dim);
                      force_pression_tot_interf_stokes_th(compo,dim)+=force_pression_th(fa7,dim);
                    }
                }
            }

        }
      else
        {
          for (int compo=0; compo<nb_compo_tot; compo++)
            {
              for (int dim=0; dim<dimension; dim++) force_pression_tot_interf_stokes_th(compo,dim)+=1e18;
            }
        }

      // ----------------------------- Force de frottements -----------------------------
      DoubleTab grad_U_P1(nb_fa7, dimension, dimension);
      DoubleTab grad_U_P2(nb_fa7, dimension, dimension);
      grad_U_P1=-1e15;
      grad_U_P2=-1e30;

      if (les_post_interf.methode_calcul_force_frottements_==Postraitement_Forces_Interfaces_FT::TRILINEAIRE_LINEAIRE_TENSEUR_COMPLET)
        {
          int interp_gradU_P1_ok=0;
          int interp_gradU_P2_ok=0;
          if (les_post_interf.localisation_tenseur_contraintes_==Postraitement_Forces_Interfaces_FT::FACES_NORMALE_X)
            {
              interp_gradU_P1_ok=trilinear_interpolation_gradU_face(indicatrice_faces, la_vitesse_stokes_th, coord_voisin_fluide_fa7_gradU_1, grad_U_P1);
              interp_gradU_P2_ok=trilinear_interpolation_gradU_face(indicatrice_faces, la_vitesse_stokes_th, coord_voisin_fluide_fa7_gradU_2, grad_U_P2);
            }
          else if (les_post_interf.localisation_tenseur_contraintes_==Postraitement_Forces_Interfaces_FT::ELEMENTS)
            {
              interp_gradU_P1_ok=trilinear_interpolation_gradU_elem(indicatrice_faces, indicatrice, la_vitesse_stokes_th, coord_voisin_fluide_fa7_gradU_1, grad_U_P1);
              interp_gradU_P2_ok=trilinear_interpolation_gradU_elem(indicatrice_faces, indicatrice, la_vitesse_stokes_th, coord_voisin_fluide_fa7_gradU_2, grad_U_P2);
            }
          if ( interp_gradU_P1_ok==1 && interp_gradU_P2_ok==1 )
            {
              DoubleTab grad_U_extrapole(nb_fa7, dimension, dimension);
              grad_U_extrapole=1e20;

              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  dUdx_P1_th_dis(fa7)=grad_U_P1(fa7,0,0);
                  dUdy_P1_th_dis(fa7)=grad_U_P1(fa7,0,1);
                  dUdz_P1_th_dis(fa7)=grad_U_P1(fa7,0,2);
                  dVdx_P1_th_dis(fa7)=grad_U_P1(fa7,1,0);
                  dVdy_P1_th_dis(fa7)=grad_U_P1(fa7,1,1);
                  dVdz_P1_th_dis(fa7)=grad_U_P1(fa7,1,2);
                  dWdx_P1_th_dis(fa7)=grad_U_P1(fa7,2,0);
                  dWdy_P1_th_dis(fa7)=grad_U_P1(fa7,2,1);
                  dWdz_P1_th_dis(fa7)=grad_U_P1(fa7,2,2);

                  dUdx_P2_th_dis(fa7)=grad_U_P2(fa7,0,0);
                  dUdy_P2_th_dis(fa7)=grad_U_P2(fa7,0,1);
                  dUdz_P2_th_dis(fa7)=grad_U_P2(fa7,0,2);
                  dVdx_P2_th_dis(fa7)=grad_U_P2(fa7,1,0);
                  dVdy_P2_th_dis(fa7)=grad_U_P2(fa7,1,1);
                  dVdz_P2_th_dis(fa7)=grad_U_P2(fa7,1,2);
                  dWdx_P2_th_dis(fa7)=grad_U_P2(fa7,2,0);
                  dWdy_P2_th_dis(fa7)=grad_U_P2(fa7,2,1);
                  dWdz_P2_th_dis(fa7)=grad_U_P2(fa7,2,2);
                  if (!maillage.facette_virtuelle(fa7))
                    {

                      double x_fa7_P1=coord_voisin_fluide_fa7_gradU_1(fa7,0);
                      double y_fa7_P1=coord_voisin_fluide_fa7_gradU_1(fa7,1);
                      double z_fa7_P1=coord_voisin_fluide_fa7_gradU_1(fa7,2);
                      double r_fa7_P1=sqrt(x_fa7_P1*x_fa7_P1+y_fa7_P1*y_fa7_P1+z_fa7_P1*z_fa7_P1);

                      dUdx_P1_th(fa7)=v_inf_stokes*(pow(x_fa7_P1,2)*z_fa7_P1/pow(r_fa7_P1,3)*( -(2+3*phi_mu)/(1+phi_mu)*3*rayon_particule/(pow(2*r_fa7_P1,2)) + phi_mu/(1+phi_mu)*9*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)) )  +
                                                    pow(x_fa7_P1,2)/(pow(r_fa7_P1,2)-pow(z_fa7_P1,2))*z_fa7_P1/r_fa7_P1*( (1-pow(z_fa7_P1/r_fa7_P1,2))*( (2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P1,2)) + phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)) )
                                                                                                                          + pow(z_fa7_P1/r_fa7_P1,2)*( (2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P1,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)) ) )
                                                    + pow(y_fa7_P1,2)*z_fa7_P1/(r_fa7_P1*(pow(r_fa7_P1,2)-pow(z_fa7_P1,2)))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(pow(2*r_fa7_P1,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4))));

                      dUdy_P1_th(fa7)=1e8; // valeur abritraire car pas besoin de cette composante

                      dUdz_P1_th(fa7)=v_inf_stokes*(x_fa7_P1/r_fa7_P1*(pow(z_fa7_P1/r_fa7_P1,2)*(-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(2*pow(r_fa7_P1,2))+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P1,4)))
                                                                       +(1-pow(z_fa7_P1/r_fa7_P1,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P1,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4))))+
                                                    pow(z_fa7_P1,2)*x_fa7_P1/pow(r_fa7_P1,3)*phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P1,4)));

                      dVdx_P1_th(fa7)=1e8;

                      dVdy_P1_th(fa7)=-1e8;

                      dVdz_P1_th(fa7)=v_inf_stokes*(y_fa7_P1/r_fa7_P1*(pow(z_fa7_P1/r_fa7_P1,2)*(-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(2*pow(r_fa7_P1,2))+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P1,4)))
                                                                       +(1-pow(z_fa7_P1/r_fa7_P1,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P1,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4))))+
                                                    pow(z_fa7_P1,2)*y_fa7_P1/pow(r_fa7_P1,3)*phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P1,4)));

                      dWdx_P1_th(fa7)=v_inf_stokes*(pow(z_fa7_P1,2)*x_fa7_P1/pow(r_fa7_P1,3)*(-(2+3*phi_mu)/(1+phi_mu)*3*rayon_particule/(pow(2*r_fa7_P1,2))+phi_mu/(1+phi_mu)*9*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)))-
                                                    x_fa7_P1/r_fa7_P1*((1-pow(z_fa7_P1/r_fa7_P1,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P1,2)+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)))
                                                                       +pow(z_fa7_P1/r_fa7_P1,2)*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P1,2)-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)))));

                      dWdy_P1_th(fa7)=v_inf_stokes*(pow(z_fa7_P1,2)*y_fa7_P1/pow(r_fa7_P1,3)*(-(2+3*phi_mu)/(1+phi_mu)*3*rayon_particule/(pow(2*r_fa7_P1,2))+phi_mu/(1+phi_mu)*9*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)))-
                                                    y_fa7_P1/r_fa7_P1*((1-pow(z_fa7_P1/r_fa7_P1,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P1,2)+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)))
                                                                       +pow(z_fa7_P1/r_fa7_P1,2)*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P1,2)-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4)))));

                      dWdz_P1_th(fa7)=v_inf_stokes*(z_fa7_P1/r_fa7_P1*(pow(z_fa7_P1/r_fa7_P1,2)*(-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(2*pow(r_fa7_P1,2))+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P1,4)))+
                                                                       (1-pow(z_fa7_P1/r_fa7_P1,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(pow(2*r_fa7_P1,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P1,4))))-
                                                    z_fa7_P1/r_fa7_P1*(1-pow(z_fa7_P1/r_fa7_P1,2))*phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P1,4)));



                      double x_fa7_P2=coord_voisin_fluide_fa7_gradU_2(fa7,0);
                      double y_fa7_P2=coord_voisin_fluide_fa7_gradU_2(fa7,1);
                      double z_fa7_P2=coord_voisin_fluide_fa7_gradU_2(fa7,2);
                      double r_fa7_P2=sqrt(x_fa7_P2*x_fa7_P2+y_fa7_P2*y_fa7_P2+z_fa7_P2*z_fa7_P2);

                      dUdx_P2_th(fa7)=v_inf_stokes*(pow(x_fa7_P2,2)*z_fa7_P2/pow(r_fa7_P2,3)*( -(2+3*phi_mu)/(1+phi_mu)*3*rayon_particule/(pow(2*r_fa7_P2,2)) + phi_mu/(1+phi_mu)*9*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)) )  +
                                                    pow(x_fa7_P2,2)/(pow(r_fa7_P2,2)-pow(z_fa7_P2,2))*z_fa7_P2/r_fa7_P2*( (1-pow(z_fa7_P2/r_fa7_P2,2))*( (2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P2,2)) + phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)) )
                                                                                                                          + pow(z_fa7_P2/r_fa7_P2,2)*( (2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P2,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)) ) )
                                                    + pow(y_fa7_P2,2)*z_fa7_P2/(r_fa7_P2*(pow(r_fa7_P2,2)-pow(z_fa7_P2,2)))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(pow(2*r_fa7_P2,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4))));

                      dUdy_P2_th(fa7)=1e8;

                      dUdz_P2_th(fa7)=v_inf_stokes*(x_fa7_P2/r_fa7_P2*(pow(z_fa7_P2/r_fa7_P2,2)*(-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(2*pow(r_fa7_P2,2))+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P2,4)))
                                                                       +(1-pow(z_fa7_P2/r_fa7_P2,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P2,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4))))+
                                                    pow(z_fa7_P2,2)*x_fa7_P2/pow(r_fa7_P2,3)*phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P2,4)));

                      dVdx_P2_th(fa7)=1e8;

                      dVdy_P2_th(fa7)=-1e8;

                      dVdz_P2_th(fa7)=v_inf_stokes*(y_fa7_P2/r_fa7_P2*(pow(z_fa7_P2/r_fa7_P2,2)*(-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(2*pow(r_fa7_P2,2))+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P2,4)))
                                                                       +(1-pow(z_fa7_P2/r_fa7_P2,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(4*pow(r_fa7_P2,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4))))+
                                                    pow(z_fa7_P2,2)*y_fa7_P2/pow(r_fa7_P2,3)*phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P2,4)));

                      dWdx_P2_th(fa7)=v_inf_stokes*(pow(z_fa7_P2,2)*x_fa7_P2/pow(r_fa7_P2,3)*(-(2+3*phi_mu)/(1+phi_mu)*3*rayon_particule/(pow(2*r_fa7_P2,2))+phi_mu/(1+phi_mu)*9*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)))-
                                                    x_fa7_P2/r_fa7_P2*((1-pow(z_fa7_P2/r_fa7_P2,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P2,2)+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)))
                                                                       +pow(z_fa7_P2/r_fa7_P2,2)*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P2,2)-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)))));

                      dWdy_P2_th(fa7)=v_inf_stokes*(pow(z_fa7_P2,2)*y_fa7_P2/pow(r_fa7_P2,3)*(-(2+3*phi_mu)/(1+phi_mu)*3*rayon_particule/(pow(2*r_fa7_P2,2))+phi_mu/(1+phi_mu)*9*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)))-
                                                    y_fa7_P2/r_fa7_P2*((1-pow(z_fa7_P2/r_fa7_P2,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P2,2)+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)))
                                                                       +pow(z_fa7_P2/r_fa7_P2,2)*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/pow(2*r_fa7_P2,2)-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4)))));

                      dWdz_P2_th(fa7)=v_inf_stokes*(z_fa7_P2/r_fa7_P2*(pow(z_fa7_P2/r_fa7_P2,2)*(-(2+3*phi_mu)/(1+phi_mu)*rayon_particule/(2*pow(r_fa7_P2,2))+phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P2,4)))+
                                                                       (1-pow(z_fa7_P2/r_fa7_P2,2))*((2+3*phi_mu)/(1+phi_mu)*rayon_particule/(pow(2*r_fa7_P2,2))-phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(4*pow(r_fa7_P2,4))))-
                                                    z_fa7_P2/r_fa7_P2*(1-pow(z_fa7_P2/r_fa7_P2,2))*phi_mu/(1+phi_mu)*3*pow(rayon_particule,3)/(2*pow(r_fa7_P2,4)));

                    }
                  for (int i=0; i<dimension; i++)
                    {
                      for (int j=0; j<dimension; j++)
                        {
                          grad_U_extrapole(fa7,i,j) = grad_U_P2(fa7,i,j)-distance_interpolation_gradU_P2*(grad_U_P2(fa7,i,j)-grad_U_P1(fa7,i,j))/(distance_interpolation_gradU_P2-distance_interpolation_gradU_P1);
                        }
                    }
                }

              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  int compo=compo_connexes_fa7(fa7);
                  Matrice_Dense tenseur_contrainte(dimension,dimension);
                  for (int i=0; i<dimension; i++)
                    {
                      for (int j=0; j<dimension; j++)
                        {
                          tenseur_contrainte(i,j) = (grad_U_extrapole(fa7,i,j) + grad_U_extrapole(fa7,j,i));
                        }
                    }

                  if (les_post_interf.flag_tenseur_contraintes_facettes_)
                    {
                      sigma_xx_interf_stokes_th_dis(fa7)=tenseur_contrainte(0,0);
                      sigma_xy_interf_stokes_th_dis(fa7)=tenseur_contrainte(0,1);
                      sigma_xz_interf_stokes_th_dis(fa7)=tenseur_contrainte(0,2);
                      sigma_yx_interf_stokes_th_dis(fa7)=tenseur_contrainte(1,0);
                      sigma_yy_interf_stokes_th_dis(fa7)=tenseur_contrainte(1,1);
                      sigma_yz_interf_stokes_th_dis(fa7)=tenseur_contrainte(1,2);
                      sigma_zx_interf_stokes_th_dis(fa7)=tenseur_contrainte(2,0);
                      sigma_zy_interf_stokes_th_dis(fa7)=tenseur_contrainte(2,1);
                      sigma_zz_interf_stokes_th_dis(fa7)=tenseur_contrainte(2,2);

                      double X_p=les_cg_fa7(fa7,0);
                      double Y_p=les_cg_fa7(fa7,1);
                      double Z_p=les_cg_fa7(fa7,2);

                      sigma_xx_interf_stokes_th(fa7)=(mu_f*v_inf_stokes/(rayon_particule))*(Z_p/(rayon_particule*(1.+phi_mu)))*(pow(X_p/rayon_particule,2)*(Z_p*Z_p/(pow(rayon_particule,2)*(1.-pow(Z_p/rayon_particule,2)))+3.*phi_mu-2.)+(pow(Y_p/rayon_particule,2)*(1.-pow(Z_p/rayon_particule,2))));
                      sigma_xy_interf_stokes_th(fa7)=(mu_f*v_inf_stokes/(rayon_particule))*(Z_p/(rayon_particule*(1.+phi_mu)))*X_p*Y_p/(rayon_particule*rayon_particule)*(3.*phi_mu-3.);
                      sigma_xz_interf_stokes_th(fa7)=(mu_f*v_inf_stokes/(rayon_particule))*X_p/(rayon_particule*(1.+phi_mu))*(-3./2.*phi_mu*(1.-pow(Z_p/rayon_particule,2))+3.*pow(Z_p/rayon_particule,2)*(phi_mu/2.-1.));
                      sigma_yy_interf_stokes_th(fa7)=(mu_f*v_inf_stokes/(rayon_particule))*(Z_p/(rayon_particule*(1.+phi_mu)))*(Y_p*Y_p/pow(rayon_particule,2)*(3.*phi_mu-2.+(Z_p*Z_p/pow(rayon_particule,2))/(1.-Z_p*Z_p/pow(rayon_particule,2)))+X_p*X_p/(pow(rayon_particule,2)*(1.-Z_p*Z_p/pow(rayon_particule,2))));
                      sigma_yz_interf_stokes_th(fa7)=(mu_f*v_inf_stokes/(rayon_particule))*Y_p/(rayon_particule*(1.+phi_mu))*(-3./2.*phi_mu*(1.-Z_p*Z_p/pow(rayon_particule,2))+3.*Z_p*Z_p/pow(rayon_particule,2)*(phi_mu/2.-1.));
                      sigma_zz_interf_stokes_th(fa7)=(mu_f*v_inf_stokes/(rayon_particule))*(Z_p/(rayon_particule*(1.+phi_mu)))*(-3.*phi_mu/2.*(1.-Z_p*Z_p/pow(rayon_particule,2))+1.-3.*Z_p*Z_p/pow(rayon_particule,2));
                    }
                  DoubleTab la_normale_fa7_x_surface(dimension);
                  for (int dim=0; dim<dimension; dim++) la_normale_fa7_x_surface(dim) =les_surfaces_fa7(fa7)*les_normales_fa7(fa7,dim);
                  DoubleVect friction_force_fa7=tenseur_contrainte*la_normale_fa7_x_surface;

                  if (!maillage.facette_virtuelle(fa7))
                    {

                      for (int dim=0; dim<dimension; dim++) force_frottements_th_dis(fa7,dim)=friction_force_fa7(dim);

                      force_frottements_th(fa7,0)=-mu_f*v_inf_stokes*les_cg_fa7(fa7,0)*les_cg_fa7(fa7,2)/(pow(rayon_particule,3))*(9.*phi_mu+8.)/(1.+phi_mu)*les_surfaces_fa7(fa7);
                      force_frottements_th(fa7,1)=-mu_f*v_inf_stokes*les_cg_fa7(fa7,1)*les_cg_fa7(fa7,2)/(pow(rayon_particule,3))*(9.*phi_mu+8.)/(1.+phi_mu)*les_surfaces_fa7(fa7);
                      force_frottements_th(fa7,2)=mu_f*v_inf_stokes*1./(2.*rayon_particule*(1.+phi_mu))*(-3.*phi_mu+pow(les_cg_fa7(fa7,2)/rayon_particule,2)*(3.*phi_mu-4.))*les_surfaces_fa7(fa7);

                      for (int dim=0; dim<dimension; dim++)
                        {
                          force_frottements_tot_interf_stokes_th_dis(compo,dim)+=friction_force_fa7(dim);
                          force_frottements_tot_interf_stokes_th(compo,dim)+=force_frottements_th(fa7,dim);
                        }
                    }
                }
            }
          else
            {
              for (int compo=0; compo<nb_compo_tot; compo++)
                {
                  for (int dim=0; dim<dimension; dim++) force_frottements_tot_interf_stokes_th(compo,dim)+=0;
                }
            }
        }
      else if (les_post_interf.methode_calcul_force_frottements_==Postraitement_Forces_Interfaces_FT::TRILINENAIRE_TENSEUR_PROJETE)
        {
          int interp_U_P1_ok=0;
          int interp_U_P2_ok=0;
          DoubleTab& U_P1_th_dis = variables_internes().U_P1_th_dis_;
          U_P1_th_dis.resize(nb_fa7, dimension);
          DoubleTab& U_P2_th_dis = variables_internes().U_P2_th_dis_;
          U_P2_th_dis.resize(nb_fa7, dimension);

          DoubleTab U_P1_spherique(nb_fa7, dimension);
          DoubleTab U_P2_spherique(nb_fa7, dimension);
          DoubleTab U_cg_spherique(nb_fa7, dimension);
          DoubleTab Urr(nb_fa7);
          DoubleTab Uthetar(nb_fa7);
          DoubleTab Uphir(nb_fa7);

          U_P1_th_dis=-1e15;
          U_P2_th_dis=-1e30;
          U_P1_spherique=-1e15;
          U_P2_spherique=-1e30;
          U_cg_spherique=-1e20;
          Urr=1e8;
          Uthetar=1e12;
          Uphir=1e15;
          const DoubleTab& positions_compo=eq_transport.get_positions_compo();
          double theta=0;
          double phi=0;
          double distance_au_cg=0;

          if (les_post_interf.localisation_tenseur_contraintes_== Postraitement_Forces_Interfaces_FT::ELEMENTS)
            {
              // 1. On calcule (interpolation trilineaire) la vitesse en P1 et P2 en coord cartesiennes
              interp_U_P1_ok=trilinear_interpolation_face(indicatrice_faces, la_vitesse_stokes_th, coord_voisin_fluide_fa7_gradU_1, U_P1_th_dis);
              interp_U_P2_ok=trilinear_interpolation_face(indicatrice_faces, la_vitesse_stokes_th, coord_voisin_fluide_fa7_gradU_2, U_P2_th_dis);
            }
          if ( interp_U_P1_ok && interp_U_P2_ok )
            {
              // 2. On passe en coordonnees spheriques

              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  int compo=compo_connexes_fa7(fa7);
                  if (!maillage.facette_virtuelle(fa7))
                    {
                      DoubleVect distance_cg_vect(dimension);
                      for (int i=0; i<dimension; i++) distance_cg_vect(i)=coord_voisin_fluide_fa7_gradU_1(fa7,i)-positions_compo(compo,i);

                      distance_au_cg=sqrt(local_carre_norme_vect(distance_cg_vect));
                      if (fabs((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg)<=1) theta=acos((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg);
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg>1) theta=0;
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,2)-positions_compo(compo,2))/distance_au_cg<-1) theta=M_PI;
                      if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))>0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))>=0)
                        {
                          phi=atan((coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))/(coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0)));
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))>0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))<0)
                        {
                          phi=atan((coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))/(coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0)))+2*M_PI;
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))<0)
                        {
                          phi=atan((coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))/(coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0)))+M_PI;
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))==0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))>0)
                        {
                          phi=M_PI/2.;
                        }
                      else if ((coord_voisin_fluide_fa7_gradU_1(fa7,0)-positions_compo(compo,0))==0 && (coord_voisin_fluide_fa7_gradU_1(fa7,1)-positions_compo(compo,1))<0)
                        {
                          phi=3.*M_PI/2.;
                        }

                      U_P1_spherique(fa7,0)=sin(theta)*cos(phi)*U_P1_th_dis(fa7,0)+sin(theta)*sin(phi)*U_P1_th_dis(fa7,1)+cos(theta)*U_P1_th_dis(fa7,2);
                      U_P1_spherique(fa7,1)=cos(theta)*cos(phi)*U_P1_th_dis(fa7,0)+cos(theta)*sin(phi)*U_P1_th_dis(fa7,1)-sin(theta)*U_P1_th_dis(fa7,2);
                      U_P1_spherique(fa7,2)=-sin(phi)*U_P1_th_dis(fa7,0)+cos(phi)*U_P1_th_dis(fa7,1);

                      U_P2_spherique(fa7,0)=sin(theta)*cos(phi)*U_P2_th_dis(fa7,0)+sin(theta)*sin(phi)*U_P2_th_dis(fa7,1)+cos(theta)*U_P2_th_dis(fa7,2);
                      U_P2_spherique(fa7,1)=cos(theta)*cos(phi)*U_P2_th_dis(fa7,0)+cos(theta)*sin(phi)*U_P2_th_dis(fa7,1)-sin(theta)*U_P2_th_dis(fa7,2);
                      U_P2_spherique(fa7,2)=-sin(phi)*U_P2_th_dis(fa7,0)+cos(phi)*U_P2_th_dis(fa7,1);

                      U_cg_spherique(fa7,0)=0;//sin(theta)*cos(phi)*vitesses_compo(compo,0)+sin(theta)*sin(phi)*vitesses_compo(compo,1)+cos(theta)*vitesses_compo(compo,2);
                      U_cg_spherique(fa7,1)=0;//cos(theta)*cos(phi)*vitesses_compo(compo,0)+cos(theta)*sin(phi)*vitesses_compo(compo,1)-sin(theta)*vitesses_compo(compo,2);
                      U_cg_spherique(fa7,2)=0;//-sin(phi)*vitesses_compo(compo,0)+cos(phi)*vitesses_compo(compo,1);

                      // On recalcule delta --> epsilon
                      int elem_diph=domaine.chercher_elements(les_cg_fa7(fa7,0), les_cg_fa7(fa7,1),les_cg_fa7(fa7,2));
                      DoubleVect delta_i(dimension);
                      delta_i(0) = fabs(domaine_vdf.dist_elem(elem_diph, domaine_vdf.face_voisins(domaine_vf.elem_faces(elem_diph, 0+dimension),1), 0));
                      delta_i(1) = fabs(domaine_vdf.dist_elem(elem_diph, domaine_vdf.face_voisins(domaine_vf.elem_faces(elem_diph, 1+dimension),1), 1));
                      if (les_normales_fa7(fa7,2)>0) delta_i(2) = fabs(domaine_vdf.dist_elem(elem_diph, domaine_vdf.face_voisins(domaine_vf.elem_faces(elem_diph, 2+dimension),1), 2));
                      else delta_i(2) = fabs(domaine_vdf.dist_elem(elem_diph, domaine_vdf.face_voisins(domaine_vf.elem_faces(elem_diph, 2),0), 2));
                      double epsilon=0;
                      for (int dim=0; dim<dimension; dim++) epsilon+= fabs(delta_i(dim)*fabs(les_normales_fa7(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage

                      // On calcule les composantes de la force de frottements en coord spherique (apres simplifications) : ff=mu*(2*Urr, Uthetar, Uphir)
                      Urr(fa7)=(-U_P2_spherique(fa7,0)+4.*U_P1_spherique(fa7,0)-3.*U_cg_spherique(fa7,0))/(2.*epsilon);
                      Uthetar(fa7)=(-U_P2_spherique(fa7,1)+4.*U_P1_spherique(fa7,1)-3.*U_cg_spherique(fa7,1))/(2.*epsilon);
                      Uphir(fa7)=(-U_P2_spherique(fa7,2)+4.*U_P1_spherique(fa7,2)-3.*U_cg_spherique(fa7,2))/(2.*epsilon);

                      force_frottements_th_dis(fa7,0)=mu_f*les_surfaces_fa7(fa7)*(2.*sin(theta)*cos(phi)*Urr(fa7)+cos(theta)*cos(phi)*Uthetar(fa7)-sin(phi)*Uphir(fa7));
                      force_frottements_th_dis(fa7,1)=mu_f*les_surfaces_fa7(fa7)*(2.*sin(theta)*sin(phi)*Urr(fa7)+cos(theta)*sin(phi)*Uthetar(fa7)+cos(phi)*Uphir(fa7));
                      force_frottements_th_dis(fa7,2)=mu_f*les_surfaces_fa7(fa7)*(2.*cos(theta)*Urr(fa7)-sin(theta)*Uthetar(fa7));


                      force_frottements_th(fa7,0)=-mu_f*v_inf_stokes*les_cg_fa7(fa7,0)*les_cg_fa7(fa7,2)/(pow(rayon_particule,3))*(9.*phi_mu+8.)/(1.+phi_mu)*les_surfaces_fa7(fa7);
                      force_frottements_th(fa7,1)=-mu_f*v_inf_stokes*les_cg_fa7(fa7,1)*les_cg_fa7(fa7,2)/(pow(rayon_particule,3))*(9.*phi_mu+8.)/(1.+phi_mu)*les_surfaces_fa7(fa7);
                      force_frottements_th(fa7,2)=mu_f*v_inf_stokes*1./(2.*rayon_particule*(1.+phi_mu))*(-3.*phi_mu+pow(les_cg_fa7(fa7,2)/rayon_particule,2)*(3.*phi_mu-4.))*les_surfaces_fa7(fa7);

                      for (int dim=0; dim<dimension; dim++)
                        {
                          force_frottements_tot_interf_stokes_th_dis(compo,dim)+=force_frottements_th_dis(fa7,dim);
                          force_frottements_tot_interf_stokes_th(compo,dim)+=force_frottements_th(fa7,dim);
                        }
                    }
                }

            }
          else
            {
              for (int compo=0; compo<nb_compo_tot; compo++)
                {
                  for (int dim=0; dim<dimension; dim++) force_frottements_tot_interf_stokes_th(compo,dim)+=0;
                }
            }
        }
    }

  mp_sum_for_each_item(force_pression_tot_interf_stokes_th);
  mp_sum_for_each_item(force_frottements_tot_interf_stokes_th);
  mp_sum_for_each_item(force_pression_tot_interf_stokes_th_dis);
  mp_sum_for_each_item(force_frottements_tot_interf_stokes_th_dis);
}
// fin EB

/*! @brief Calcul de la derivee en temps de la vitesse.
 *
 */
DoubleTab& Navier_Stokes_FT_Disc::derivee_en_temps_inco(DoubleTab& vpoint)
{
  // Preparation des champs utilises pour le calcul des derivees en temps
  // S'il n'y a pas d'equation de transport des interfaces entre phases fluides,
  // on ne recalcule pas les proprietes.
  {
    REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;

    if (refeq_transport.non_nul())
      {
        const int calcul_precis_indic_face=refeq_transport.valeur().calcul_precis_indic_faces();
        FT_disc_calculer_champs_rho_mu_nu_dipha(domaine_dis().valeur(),
                                                fluide_diphasique(),
                                                refeq_transport.valeur().inconnue().valeur().valeurs(),
                                                // (indicatrice)
                                                refeq_transport.valeur(), // indicatrice_face
                                                champ_rho_elem_.valeur().valeurs(),
                                                champ_nu_.valeur().valeurs(),
                                                champ_mu_.valeur().valeurs(),
                                                champ_rho_faces_.valeur().valeurs(), schema_temps().temps_courant(),calcul_precis_indic_face);
      }
    else
      {

        if (sub_type(Fluide_Incompressible,milieu()))
          {
            const Domaine_dis_base& zdis = domaine_dis().valeur();
            const Fluide_Incompressible& phase_0 = ref_cast(Fluide_Incompressible,milieu());
            FT_disc_calculer_champs_rho_mu_nu_mono(zdis,
                                                   phase_0,
                                                   champ_rho_elem_,
                                                   champ_mu_,
                                                   champ_nu_,
                                                   champ_rho_faces_);
          }
        else if (sub_type(Particule_Solide,milieu()))
          {
            const Particule_Solide& phase_0 = ref_cast(Particule_Solide,milieu());//EB
            const Domaine_dis_base& zdis = domaine_dis().valeur();
            FT_disc_calculer_champs_rho_mu_nu_mono(zdis,
                                                   phase_0,
                                                   champ_rho_elem_,
                                                   champ_mu_,
                                                   champ_nu_,
                                                   champ_rho_faces_);
          }

      }
  }

  vpoint = 0.;

  // =====================================================================
  // Methode de projection :
  // Premiere etape : calcul de u_etoile (tous les termes de N.S. sauf la pression)

  // Contribution des operateurs diffusion et convection :
  // Operateur de diffusion : valeurs discretes homogenes a
  //                                             / d             \    //
  //                INTEGRALE                    | -- (rho * v)  |    //
  //                (sur le volume de controle)  \ dt            /    //
  // B.M. 08/2004 : on envoie la vitesse "v" a l'operateur qui calcule
  //                div (mu * (grad(v)+tr(grad(v))))
  //                (on a associe "mu" a la "diffusivite" de l'operateur,
  //                 voir Navier_Stokes_FT_Disc::lire)

  // EB
  {
    REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide; // EB
    if (refeq_transport.non_nul())
      {
        const int calcul_precis_indic_arete=refeq_transport.valeur().calcul_precis_indic_aretes(); // EB
        if (calcul_precis_indic_arete) refeq_transport.valeur().get_compute_indicatrice_aretes_internes();
      }
  }
  terme_diffusif.calculer(la_vitesse.valeurs(),
                          variables_internes().terme_diffusion.valeur().valeurs());
  solveur_masse.appliquer(variables_internes().terme_diffusion.valeur().valeurs());

  // Termes sources : gravite et tension de surface,
  // valeurs discretes homogenes a
  //                                             / d             \    //
  //                INTEGRALE                    | -- (rho * v)  |    //
  //                (sur le volume de controle)  \ dt            /    //

  // HMS // EB
  int flag_correction_trainee,is_solid_particle=0;
  DoubleTab& terme_source_collisions=variables_internes().terme_source_collisions.valeur().valeurs() ;
  terme_source_collisions=0;

  DoubleTab& terme_correction_trainee=variables_internes().terme_correction_trainee.valeur().valeurs() ;
  terme_correction_trainee=0;
  flag_correction_trainee=variables_internes().flag_correction_trainee_; // EB
  {
    REF(Transport_Interfaces_FT_Disc) & refeq_transport =
      variables_internes().ref_eq_interf_proprietes_fluide;
    if (refeq_transport.non_nul())
      {
        const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
        is_solid_particle=eq_transport.is_solid_particle(); // EB
      }
  }

  {
    // Si une equation de transport est associee aux proprietes du fluide,
    // on ajoute le terme de tension de surface.
    REF(Transport_Interfaces_FT_Disc) & refeq_transport =
      variables_internes().ref_eq_interf_proprietes_fluide;

    if (refeq_transport.non_nul())
      {
        const Champ_base& indicatrice = refeq_transport.valeur().get_update_indicatrice();
        //const Champ_base& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces();
        Champ_base& gradient_i = variables_internes().gradient_indicatrice.valeur();
        // Note:
        // On appelle la version const de maillage_interface() (qui est publique) car
        // on passe par const Transport_Interfaces_FT_Disc :
        const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
        Transport_Interfaces_FT_Disc& eq_transport_non_const = refeq_transport.valeur();
        const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
        const DoubleTab& distance_interface_sommets =
          eq_transport.get_update_distance_interface_sommets();

        calculer_gradient_indicatrice(indicatrice,
                                      distance_interface_sommets,
                                      gradient_i);
        calculer_champ_forces_superficielles(maillage,
                                             gradient_i,
                                             variables_internes().potentiel_elements,
                                             variables_internes().potentiel_faces,
                                             variables_internes().terme_source_interfaces);

        if (is_solid_particle) calculer_champ_forces_collisions(indicatrice.valeurs(), terme_source_collisions, eq_transport, eq_transport_non_const, refeq_transport, maillage); //HMS
        if (flag_correction_trainee) calculer_correction_trainee(terme_correction_trainee, eq_transport, eq_transport_non_const, refeq_transport, maillage);
      }
    else
      {
        variables_internes().terme_source_interfaces.valeurs() = 0;
      }
  }
  solveur_masse.appliquer(variables_internes().terme_source_interfaces.valeur().valeurs());
  if (is_solid_particle) solveur_masse.appliquer(terme_source_collisions); //HMS
  if (flag_correction_trainee)  solveur_masse.appliquer(terme_correction_trainee); // EB
  // Autres termes sources (acceleration / repere mobile)
  //  Valeurs homogenes a
  //                                             / d       \          //
  //                INTEGRALE                    | -- (v)  |          //
  //                (sur le volume de controle)  \ dt      /          //
  // (voir "preparer_calcul", commentaire sources().associer_rho...)
  variables_internes().terme_source.valeur().valeurs() = 0.;
  les_sources.ajouter(variables_internes().terme_source.valeur().valeurs());
  solveur_masse.appliquer(variables_internes().terme_source.valeur().valeurs());

  // Operateur de convection : valeurs discretes homogenes a
  //                                             / d       \          //
  //                INTEGRALE                    | -- (v)  |          //
  //                (sur le volume de controle)  \ dt      /          //
  // B.M. 08/2004 : on transporte "v" et non "rho_v"...

  DoubleTab& terme_convection_valeurs = variables_internes().terme_convection.valeur().valeurs();
  bool calcul_explicite = false;
  if (parametre_equation_.non_nul() && sub_type(Parametre_implicite, parametre_equation_.valeur()))
    {
      Parametre_implicite& param2 = ref_cast(Parametre_implicite, parametre_equation_.valeur());
      calcul_explicite = param2.calcul_explicite();
    }
  if (schema_temps().diffusion_implicite() && !calcul_explicite)
    {
      terme_convection_valeurs=0;
      derivee_en_temps_conv(terme_convection_valeurs,la_vitesse.valeurs());
    }
  else
    {
      terme_convectif.calculer(la_vitesse.valeurs(),terme_convection_valeurs);
    }
  solveur_masse.appliquer(variables_internes().terme_convection.valeur().valeurs());

  // Ajout des differentes contributions a vpoint :
  const DoubleTab& tab_rho_faces = champ_rho_faces_.valeur().valeurs();
  const DoubleVect& volumes_entrelaces = ref_cast(Domaine_VF, domaine_dis().valeur()).volumes_entrelaces();
  const DoubleTab& tab_diffusion = variables_internes().terme_diffusion.valeur().valeurs();
  const DoubleTab& termes_sources_interf = variables_internes().terme_source_interfaces.valeur().valeurs();
  const DoubleTab& termes_sources = variables_internes().terme_source.valeur().valeurs();
  const DoubleTab& tab_convection = variables_internes().terme_convection.valeur().valeurs();
  const int n = vpoint.dimension(0);
  const int nbdim1 = (vpoint.line_size() == 1);
  const int m =  vpoint.line_size();

  DoubleTab gravite_face(inconnue().valeurs());
  if (milieu().a_gravite())
    {
      ArrOfDouble g(dimension);
      // Pour l'instant : gravite uniforme g => phi(s) = - x scalaire g
      const DoubleTab& gravite = milieu().gravite().valeurs();
      if (gravite.nb_dim() != 2 || gravite.line_size() != dimension)
        {
          Cerr << "Error for calculer_champ_forces_superficielles\n";
          Cerr << " gravite.line_size() != Objet_U::dimension" << finl;
          Process::exit();
        }
      for (int j = 0; j < dimension; j++)
        g[j] =  gravite(0,j);

      // On multiplie par les volumes entrelaces et on applique ensuite le solveur masse
      //  (traitement special des CL de Dirichlet)
      if (nbdim1)
        {
          const IntTab& face_voisins = le_dom_dis.valeur().valeur().face_voisins();
          const IntVect& orientation = ref_cast(Domaine_VDF, domaine_dis().valeur()).orientation();
          if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
            {
              for (int face=0; face<n; face++)
                gravite_face(face,0)=volumes_entrelaces(face)*g[orientation[face]];
            }
          else
            {
              gravite_face = 0.; // En gradI, on ne prend pas directement la gravite
            }
          // Boussinesq Approximation in use :
          if (variables_internes().is_boussinesq_)
            {
              REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
              if (!refeq_transport.non_nul())
                {
                  Cerr << "Trying to use Boussinesq approximation on a 2phase flow when the transport equation is not specified" << finl;
                  Process::exit();
                }
              const DoubleTab& indicatrice = refeq_transport.valeur().get_update_indicatrice().valeurs();
              // First phase with temperature :
              if (variables_internes().ref_equation_mpoint_.non_nul())
                {
                  compute_boussinesq_additional_gravity(variables_internes().ref_equation_mpoint_.valeur(),
                                                        fluide_diphasique(),
                                                        face_voisins,volumes_entrelaces,orientation,
                                                        indicatrice,
                                                        g,
                                                        gravite_face);
                }

              // Second phase with temperature :
              if (variables_internes().ref_equation_mpoint_vap_.non_nul())
                {
                  compute_boussinesq_additional_gravity(variables_internes().ref_equation_mpoint_vap_.valeur(),
                                                        fluide_diphasique(),
                                                        face_voisins,volumes_entrelaces,orientation,
                                                        indicatrice,
                                                        g,
                                                        gravite_face);
                }
              // The end of boussinesq force source term.
            }
          // End of if nbdim1 that is VDF.
        }
      else
        {
          // VEF case:
          if (variables_internes().is_boussinesq_)
            {
              Cerr << "Trying to use Boussinesq approximation on a 2phase flow in VEF? Not yet available. Ask TRUST support." << finl;
              Process::exit();
            }
          if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
            {
              for (int face=0; face<n; face++)
                for (int dim=0; dim<m; dim++)
                  gravite_face(face,dim)=volumes_entrelaces(face)*g[dim];
            }
          else
            {
              gravite_face = 0.; // En gradI, on ne prend pas directement la gravite
            }
        }
      solveur_masse.appliquer(gravite_face);
    }
  else
    {
      // Pas de gravite :
      gravite_face = 0.;
    }

  IntTab flag_gradP(n);
  IntTab coef_TSF(n);
  coef_TSF = 1;
  int flag_diff;
  DoubleTab& gradP = variables_internes().gradient_pression.valeurs();
  if (schema_temps().diffusion_implicite())
    {
      //on calcule gradP (pour le qdm)
      gradient.calculer(la_pression.valeur().valeurs(), gradP);
      solveur_masse.appliquer(gradP);
      // si variables_internes().is_penalized=1 alors variables_internes().is_explicite=0
      if( variables_internes().is_penalized )
        flag_gradP = 0 ; //(grad P raide)
      else
        flag_gradP = 1 ;
      // on ajoute pas la diffusion cela sera fait par Gradient_conjugue_diff_impl
      // sauf si !is_explicite (terme forcage vitesse implicite; resolution conjointe forcage diffusion)
      // car dans ce cas la vitesse a imposer a besoin d'une bonne approximation de vpoint
      if( !variables_internes().is_explicite && !calcul_explicite)
        flag_diff = 1;
      else
        flag_diff = 0;
    }
  else
    {
      flag_gradP = 0;
      flag_diff = 1;
    }

  bool interf_vitesse_imposee_ok = false;
  int nb_eqs = variables_internes().ref_eq_interf_vitesse_imposee.size();
  int nb_eq_non_nul = 0;
  for (int k=0; k<nb_eqs; k++)
    {
      REF(Transport_Interfaces_FT_Disc) & refeq_transport =
        variables_internes().ref_eq_interf_vitesse_imposee[k];

      if (refeq_transport.non_nul()) nb_eq_non_nul +=1;
    }
  DoubleTab terme_mul;
  if ( nb_eq_non_nul == nb_eqs && nb_eqs != 0 )
    {
      interf_vitesse_imposee_ok = true;
      terme_mul.copy(champ_rho_faces_.valeur().valeurs(), Array_base::COPY_INIT);
      terme_mul = 0.;
    }

  REF(Transport_Interfaces_FT_Disc) & refeq_transport_2pha =
    variables_internes().ref_eq_interf_proprietes_fluide;
  if (refeq_transport_2pha.non_nul() && interf_vitesse_imposee_ok && variables_internes().is_penalized)
    {
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
      const Domaine_VF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const IntTab& elem_faces = domaine_vf.elem_faces();
      const IntVect& orientation = domaine_vdf.orientation();
      const int   nb_faces_elem = elem_faces.line_size();
      for (int k=0; k<nb_eqs; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];
          const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
          for (int i = 0; i < face_voisins.dimension(0) ; i++)
            {
              if (indicatrice_faces(i) > 0.)
                {
                  flag_gradP(i) = 0;
                  coef_TSF(i) = 0; //annulation local terme source interf
                  const int ori=orientation(i);
                  const int voisin0 = face_voisins(i,0);
                  if (voisin0 >= 0)
                    {
                      int face_visavi = elem_faces(voisin0, ori) + elem_faces(voisin0, ori+Objet_U::dimension) - i;
                      for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                        {
                          const int face = elem_faces(voisin0, i_face);
                          if (indicatrice_faces(face) == 0. && face == face_visavi)
                            {
                              flag_gradP(face) = 0;
                              coef_TSF(face) = 0; //annulation local terme source interf
                            }
                        }
                    }
                  const int voisin1 = face_voisins(i,1);
                  if (voisin1 >= 0)
                    {
                      int face_visavi = elem_faces(voisin1, ori) + elem_faces(voisin1, ori+Objet_U::dimension) - i;
                      for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                        {
                          const int face = elem_faces(voisin1, i_face);
                          if (indicatrice_faces(face) == 0. && face == face_visavi)
                            {
                              flag_gradP(face) = 0;
                              coef_TSF(face) = 0; //annulation local terme source interf
                            }
                        }
                    }
                }
            }
        }
    }

  // Ajout des differentes contributions a vpoint :
  for (int i = 0; i < n; i++)
    {
      const double rho_face = tab_rho_faces(i);

      for (int j = 0; j < m; j++)
        vpoint(i, j) = ( - flag_gradP(i) * gradP(i,j) + flag_diff * tab_diffusion(i,j) + coef_TSF(i) *termes_sources_interf(i,j) + is_solid_particle*terme_source_collisions(i) + flag_correction_trainee*terme_correction_trainee(i)) / rho_face
                       + tab_convection(i,j) + termes_sources(i,j) + gravite_face(i,j);
    }
  vpoint.echange_espace_virtuel();

  //  si terme forcage vitesse explicite => J'ai tout, je peux resoudre (Gradient_conjugue_diff_impl)
  if (schema_temps().diffusion_implicite() && !calcul_explicite && variables_internes().is_explicite )
    {
      DoubleTab derivee(la_vitesse.valeurs());
      // on indique au solveur masse de diviser par rho en plus du volume car l'operateur de diffusion renvoit qqqc en rho d u/dt
      solveur_masse->set_name_of_coefficient_temporel(champ_rho_faces_.valeur().le_nom());

      DoubleTrav tt(vpoint);
      tt=vpoint;
      derivee=inconnue().valeurs();
      Equation_base::Gradient_conjugue_diff_impl( tt, derivee ) ;

      solveur_masse->set_name_of_coefficient_temporel("no_coeff");

      vpoint=derivee;
      // on retire le gradient si on ne penalise pas:
      for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
          vpoint(i, j) += gradP(i,j) / tab_rho_faces(i) ;
    }

  const int nfaces = vpoint.dimension_tot(0);

  if( interf_vitesse_imposee_ok )
    {
      int compteur_vimp_regul =0;
      DoubleTrav vpoint0(vpoint) ;
      vpoint0 = vpoint ;
      terme_mul = 1.0;
      // S'il y a une equation de transport avec vitesse imposee, on impose:
      DoubleTrav forces_tot(vpoint);
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k=0; k<nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];

          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

          const DoubleTab& inco_val = inconnue().valeur().valeurs();
          const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
          DoubleTab& source_val = variables_internes().terme_source.valeur().valeurs();
          const double temps = schema_temps().temps_courant();
          const double dt = schema_temps().pas_de_temps();

          //On ajoute un terme source a vpoint pour imposer au fluide la vitesse de l interface
          //source_val est rempli (peut etre postraite)
          eq_transport.modifier_vpoint_pour_imposer_vit(inco_val,vpoint0,vpoint,rho_faces,source_val,temps,dt,variables_internes().is_explicite,variables_internes().eta);
          source_val.echange_espace_virtuel();
          forces_tot += source_val;

          // Afin de savoir s il existe des interfaces IBC/fluide regularisees.
          // Si oui, on ne penalise que pour indic=1.
          if (eq_transport.get_vimp_regul()) compteur_vimp_regul++;
          //calcul du terme 1 + somme Xs / eta
          if ( !variables_internes().is_explicite )
            {
              const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_indicatrice_faces().valeurs();
              for (int i = 0; i < nfaces; i++)
                {
                  if ((eq_transport.get_vimp_regul()==0 && indicatrice_faces(i) > 0.) ||
                      (eq_transport.get_vimp_regul()==1 && indicatrice_faces(i)== 1.)) terme_mul(i) += 1. / variables_internes().eta;
                }
            }
        }
      terme_mul.echange_espace_virtuel();
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco terme_mul:",terme_mul);

      if (schema_temps().diffusion_implicite())
        {
          // si !variables_internes().is_explicite (terme forcage implicite) on retire la diffusion explicite
          // de vpoint (implicitee dans Gradient_conjugue_diff_impl_IBC)
          if( !variables_internes().is_explicite && !calcul_explicite)
            {
              const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
              const DoubleTab& diffusion = variables_internes().terme_diffusion.valeur().valeurs();

              for (int i = 0; i < vpoint.dimension(0); i++)
                for (int j = 0; j < vpoint.line_size(); j++)
                  vpoint(i,j) -= diffusion(i,j) / rho_faces(i) ;

              vpoint.echange_espace_virtuel() ;
              DoubleTab derivee(la_vitesse.valeurs());
              // on indique au solveur masse de diviser par rho en plus du volume car l'operateur de diffusion renvoit qqqc en rho d u/dt
              {
                solveur_masse->set_name_of_coefficient_temporel(champ_rho_faces_.valeur().le_nom());

                DoubleTrav tt(vpoint);
                tt=vpoint;
                Equation_base::Gradient_conjugue_diff_impl( tt, derivee, terme_mul ) ;

                solveur_masse->set_name_of_coefficient_temporel("no_coeff");
              }
              vpoint=derivee;

              // on retire le gradient
              const int nbis = vpoint.dimension(0);
              const int mbis = vpoint.line_size();
              for (int i = 0; i < nbis; i++)
                for (int j = 0; j < mbis; j++)
                  vpoint(i, j) += ( flag_gradP(i)*gradP(i,j)  ) / rho_faces(i) ;
            }
          vpoint.echange_espace_virtuel() ;
        }
      else
        {
          // si on implicite le calcul de vpoint
          // On divise vpoint par le terme multiplicatif calcule avant
          if ( !variables_internes().is_explicite )
            {
              const int mbis = vpoint.line_size();
              // calcul de vpoint : vpoint / (1 + somme Xs/eta )
              for (int i = 0; i < nfaces; i++)
                for (int j = 0; j < mbis; j++)
                  vpoint(i,j) /= terme_mul(i);

              vpoint.echange_espace_virtuel();
            }

        }
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco vpoint:",vpoint);

      // Dans le cas penalize + vitesse imposee regularisee,
      // seuls les ddl tels que indic_faces=1 sont penalisees, les
      // autres sont forces avec un DF :
      if( !variables_internes().is_explicite && compteur_vimp_regul)
        {
          vpoint0=vpoint;
          // S'il y a une equation de transport avec vitesse imposee, on impose:
          DoubleTrav forces_totbis(vpoint);
          for (int k=0; k<nb_eqs_bis; k++)
            {
              REF(Transport_Interfaces_FT_Disc) & refeq_transport =
                variables_internes().ref_eq_interf_vitesse_imposee[k];

              Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

              const DoubleTab& inco_val = inconnue().valeur().valeurs();
              const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
              DoubleTab& source_val = variables_internes().terme_source.valeur().valeurs();
              const double temps = schema_temps().temps_courant();
              const double dt = schema_temps().pas_de_temps();

              //On ajoute un terme source a vpoint pour imposer au fluide la vitesse de l interface
              //source_val est rempli (peut etre postraite)
              eq_transport.modifier_vpoint_pour_imposer_vit(inco_val,vpoint0,vpoint,rho_faces,source_val,temps,dt,/* is_explicite */ 1, /* eta */ 1.);
              source_val.echange_espace_virtuel();
              forces_totbis += source_val;
            }
        }
      vpoint.echange_espace_virtuel() ;

      //Calcul des efforts exerces par le fluide sur chaque interface
      //Attention valable si les ibc ne se chevauchent pas
      if (limpr())
        {
          const DoubleTab& tab_rho_facesbis = champ_rho_faces_.valeur().valeurs();
          for (int k=0; k<nb_eqs_bis; k++)
            {
              REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
              Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
              eq_transport.calcul_effort_fluide_interface(vpoint,tab_rho_facesbis,forces_tot,variables_internes().is_explicite,variables_internes().eta);
              forces_tot.echange_espace_virtuel();
            }
        }
    }

  // Assemblage de la matrice INTEGRALE            ( div(1/rho * grad(P)) ) :
  //                          (volume de controle pression)
  // Si l'option "matrice_pression_invariante" est activee, on ne recalcule
  // pas la matrice :
  if ( !interf_vitesse_imposee_ok )
    {
      if ( !variables_internes().matrice_pression_invariante )
        {
          assembleur_pression().valeur().assembler_rho_variable(matrice_pression_,
                                                                champ_rho_faces_.valeur());
          // On a modifie la matrice, il faut reinitialiser le solveur
          //  (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
          solveur_pression().valeur().reinit();
        }
    }

  // ====================================================================
  // Methode de projection :
  // Deuxieme etape : projection du champ de vitesse sur le sous-espace
  //  a divergence nulle.
  // Trouver la pression P telle que
  //   d(u)/dt = vpoint - 1/rho * grad(P)
  //   div( d(u)/dt ) = 0
  // Soit :
  //   div(1/rho * grad(P)) = div(vpoint)

  // Calcul du second membre :
  //  div(vpoint) a l'interieur du domaine,
  //  prise en compte des conditions aux limites de pression/vitesse

  DoubleTab& secmem = variables_internes().second_membre_projection.valeurs();
  DoubleTab& secmem2 = variables_internes().second_membre_projection_jump_.valeurs();
  const int nb_elem = secmem2.dimension(0);
  const double dt = schema_temps().pas_de_temps();
  const DoubleTab& inco = inconnue().valeur().valeurs();
  // secmem = div(U/dt+vpoint) = div(U(n+1)/dt)
  DoubleTab du(inco);
  du /= dt;
  du += vpoint;
  divergence.calculer(du, secmem);
  secmem *= -1;
#if NS_VERBOSE
  double int_sec_mem = 0;
  for (int elem = 0; elem < secmem.dimension(0); elem++)
    int_sec_mem +=secmem(elem);
  Cerr << "Secmem before tcl= " << int_sec_mem << finl;
#endif
  // Prise en compte des sources de constituant (terme source dans une equation concentration)
  // Pour ne pas faire figurer la source deux fois, je la mets uniquement dans l'equation
  // de concentration... elle produit alors un terme source de div_u dans Navier_Stokes:
  const Noms& noms_eq = variables_internes().equations_concentration_source_fluide_;
  for (int i_eq = 0; i_eq < noms_eq.size(); i_eq++)
    {
      const Equation_base& eq = probleme().get_equation_by_name(noms_eq[i_eq]);
      for (int i_source = 0; i_source < eq.sources().size(); i_source++)
        {
          if (!sub_type(Terme_Source_Constituant_Vortex_VEF_Face, eq.sources()[i_source].valeur()))
            continue;

          const Terme_Source_Constituant_Vortex_VEF_Face& src = ref_cast(Terme_Source_Constituant_Vortex_VEF_Face, eq.sources()[i_source].valeur());
          src.ajouter_terme_div_u(secmem, schema_temps().pas_de_temps());
        }

    }
  secmem.echange_espace_virtuel();

  // Prise en compte du terme source div(u) du changement de phase
  if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
    {
      // GB2016 : Le calcul de mpoint ci-dessous me semble inutile car il est fait au debut de calculer_delta_u_interface:
      // GB2016 : Mais si je ne le fais pas, j'ai parfois : 'vx.get_md_vector() == md' dans calculer_delta_u_interface

      // GB2022 : Je ne comprend toujours pas bien pourquoi, mais il faut calculer les mpoints
      //          pour avoir un cas FTD_Boiling_bubble avec une extension correcte
      //          (sinon le champ postraite de T est moche dans l'extension)
      if (variables_internes().ref_equation_mpoint_.non_nul())
        variables_internes().ref_equation_mpoint_.valeur().calculer_mpoint(variables_internes().mpoint.valeur());
      if (variables_internes().ref_equation_mpoint_vap_.non_nul())
        variables_internes().ref_equation_mpoint_vap_.valeur().calculer_mpoint(variables_internes().mpoint_vap.valeur());

      // Pas une ref, mais un tableau de travail local dans lequel on peut ajouter mointv
      DoubleTab mpoint = variables_internes().ref_equation_mpoint_.valeur().get_mpoint();
      if (variables_internes().ref_equation_mpoint_vap_.non_nul())
        {
          const DoubleTab& mpointv = variables_internes().ref_equation_mpoint_vap_.valeur().get_mpoint();
          for (int elem = 0; elem < nb_elem; elem++)
            mpoint[elem] += mpointv[elem];
        }
      // We can compute delta_u_interface:
      // depending on the option, either historical or new, the calculation may be based on the values filled in secmem2
      calculer_delta_u_interface(variables_internes().delta_u_interface, -1 /* vitesse de l'interface */, variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);

      const Fluide_Diphasique& fluide_diph = fluide_diphasique();
      const Fluide_Incompressible& phase_0 = fluide_diph.fluide_phase(0);
      const Fluide_Incompressible& phase_1 = fluide_diph.fluide_phase(1);
      const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
      const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
      const double rho_phase_0 = tab_rho_phase_0(0,0);
      const double rho_phase_1 = tab_rho_phase_1(0,0);
      const double jump_inv_rho = 1./rho_phase_1 - 1./rho_phase_0;
      if (variables_internes().new_mass_source_)
        {

          const DoubleTab& interfacial_area = variables_internes().ai.valeur().valeurs();
          for (int elem = 0; elem < nb_elem; elem++)
            secmem2[elem] = jump_inv_rho*interfacial_area[elem]*mpoint[elem];
        }
      else
        {
          const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
          const IntTab& face_voisins = domaine_vf.face_voisins();
          const IntTab& elem_faces = domaine_vf.elem_faces();
          const int nb_faces_elem = elem_faces.line_size();
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_proprietes_fluide;
          const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
#if NS_VERBOSE
          const DoubleTab& indicatrice = eq_transport.inconnue().valeurs();
#endif
          // Distance a l'interface discretisee aux elements:
          const DoubleTab& distance = eq_transport.get_update_distance_interface().valeurs();
          divergence.calculer(variables_internes().delta_u_interface->valeurs(), secmem2);
          // On ne conserve que la divergence des elements traverses par l'interface
          for (int elem = 0; elem < nb_elem; elem++)
            {
              const double dist = distance(elem);
              int i_face = -1;
              if (dist < -1e20)
                {
                  // Distance invalide: on est loin de l'interface
                  i_face = nb_faces_elem;
                }
              else
                {
                  // Y a-t-il un voisin pour lequel la distance est de signe oppose
                  for (i_face = 0; i_face < nb_faces_elem; i_face++)
                    {
                      const int face = elem_faces(elem, i_face);
                      const int voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                      if (voisin >= 0)
                        {
                          const double d = distance(voisin);
                          if (d > -1e20 && d * dist < 0.)
                            {
#if NS_VERBOSE
                              Cerr << "Compa "<< secmem2[elem] << " " << jump_inv_rho*interfacial_area[elem]*mpoint[elem] << finl;
#endif
                              break; // Changement de signe
                            }
                        }
                    }
                }
              if (i_face == nb_faces_elem)
                {
                  // Tous les voisins sont du meme cote de l'interface
                  secmem2(elem) = 0.;
#if NS_VERBOSE
                  if(interfacial_area[elem]>DMINFLOAT)
                    {
                      Cerr << "[WARNING] secmem2 is set to zero in element whereas phase is not pure (indic= "
                           << indicatrice[elem] << "). This is because the choice is based on the signs of distance." << finl;
                      // Indeed, in a diagonal case, for the cell with "x", indic can be close to (but not) pure and all neighbouring
                      // cells can still have the same distance.
                      // __________
                      // |  | /|  |
                      // |__|/_|__|
                      // | /| x|  |
                      // |/_|__|__|
                    }
#endif
                }
            }
        }
#if TCL_MODEL

      /* double int_sec_mem2 = 0.;
      for (int elem = 0; elem < nb_elem; elem++)
        {
          int_sec_mem2 +=secmem2(elem);
          int_sec_mem +=secmem(elem);
        }
      Cerr << "Integral of secmem2 before TCL and /DT : " << int_sec_mem2 << finl; */

      // Now that the correction "corriger_mpoint" is performed directly into Convection_Diffusion_Temperature_FT_Disc,
      // it is no inter required to compute it here. The correction should propagate to calculer_delta_vitesse (in the discreete sense),
      // and subsequently to secmem2. So that in the end, the whole TCL contribution should be accounted for naturally.
      // However, in the discretized version, It is still required to correct secmem2 even though mpoint was corrected itself.
      // It is because of the way secmem is computed (as the div( delta u)) that is using part of cells not crossed by the interface
      // (in the wall-normal direction). It results in an underestimation of the real TCL contribution...
      if (probleme_ft().tcl().is_activated())
        {
          Cerr << "[TCL] Contact line model is activated" << finl;


          const Triple_Line_Model_FT_Disc& tcl = ref_cast(Triple_Line_Model_FT_Disc, probleme_ft().tcl());

          const ArrOfInt& elems_with_CL_contrib = tcl.elems();
          // const ArrOfInt& faces_with_CL_contrib = probleme_ft().tctl().boundary_faces();
          const ArrOfDouble& Q_from_CL = tcl.Q();
          // const ArrOfDouble& mpoint_from_CL = probleme_ft().tcl().mp();

          // ArrOfInt elems_with_CL_contrib;
          // ArrOfInt faces_with_CL_contrib;
          // ArrOfDouble mpoint_from_CL;
          // ArrOfDouble Q_from_CL;
          // GB. 18/12/19. This call is actually the one filling the TCL tables (elems_, mp_ and Q_);
          // probleme_ft().tcl().compute_TCL_fluxes_in_all_boundary_cells(elems_with_CL_contrib,
          //                                                             faces_with_CL_contrib,
          //                                                             mpoint_from_CL,
          //                                                             Q_from_CL);

          // Correct the field mpoint in wall-adjacent cells to account for TCL model:
          // ---> It's not added to mpoint now. Its contribution is added in
          //      Convection_Diffusion_Temperature_FT_Disc::derivee_en_temps_inco
          //      to be after the evaluation of the extended velocities (interfacial and liquid)
          //      and their interpolation. That way, interpolation can still operate on a smooth field.
          // DoubleTab& mpoint = variables_internes().mpoint.valeur().valeurs();
          // probleme_ft().tcl().corriger_mpoint(elems_with_CL_contrib,mpoint_from_CL,mpoint);

          const double Lvap = fluide_diph.chaleur_latente();
          const double coef = jump_inv_rho/Lvap;
          // Correct the secmem2 contribution due to TCL :
          if (!variables_internes().mpoint_inactif)
            probleme_ft().tcl().corriger_secmem(coef, secmem2);

          const int check_consistency = 1 ; // local option to check that secmem2 in near-wall cell is actually well calculated
          if (check_consistency)
            {
              Cerr << "Verifying Contact line model consistency" << finl;
              double error = 0.;
              const int nb_contact_line_contribution = elems_with_CL_contrib.size_array();
              for (int idx = 0; idx < nb_contact_line_contribution; idx++)
                {
                  const int elem = elems_with_CL_contrib[idx];
                  const double sec = secmem2(elem);
                  double Q = 0.;
                  // Go through the list to find all occurences of elem;
                  for (int idx2 = 0; idx2 < nb_contact_line_contribution; idx2++)
                    {
                      if (elem == elems_with_CL_contrib[idx2])
                        {
                          Q +=Q_from_CL[idx2];
                        }
                    }
                  const double value = coef*Q;

                  // sec and value should be the same:
                  error +=fabs(sec - value);
                  if (fabs(sec - value) > 1.e-12) // changed from 1.e-12 to 1.e-7 ---- for test
                    {
                      Cerr << "local difference sec-value=" << sec <<" - " << value << " = " << (sec - value) << finl;
                    }
                }

              if (error > 1.e-8)
                {
                  Cerr << "Final error : " << error << " is fatal!" << finl;
                  Process::exit();
                }

            }
        }
#endif
      secmem2 /= schema_temps().pas_de_temps();
      secmem += secmem2;
      secmem.echange_espace_virtuel();
#if NS_VERBOSE
      double int_sec_mem2 = 0;
      double int_sec_mem = 0;
      for (int elem = 0; elem < nb_elem; elem++)
        {
          int_sec_mem2 +=secmem2(elem);
          int_sec_mem +=secmem(elem);
        }
      Cerr << "Integral of secmem2 after TCL and /DT : " << int_sec_mem2 << finl;
      Cerr << "Integral of secmem after TCL and /DT : " << int_sec_mem << finl;
#endif
    }

  Champ_Fonc champ_rho_faces_modifie(champ_rho_faces_);
  DoubleTab& rho_faces_modifie = champ_rho_faces_modifie.valeur().valeurs();

  if ( interf_vitesse_imposee_ok )
    {

      // On verifie s il existe des interfaces IBC/fluide regularisees.
      int compteur_vimp_regul =0;
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k=0; k<nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];
          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

          if (eq_transport.get_vimp_regul()==1)
            compteur_vimp_regul++;
        }

      // Creation d'un nouveau rho qui prend en compte les vitesses imposees dans
      // le terme de forcage
      // Si des interfaces IBC/fluide sont regularisees, la projection devient classique
      // (a ameliorer potentiellement).
      // Mais dans le cas ou on a une interface diphasique, on bloque toutes les vitesses impossees
      // en PDF (regularise ou non)
      int modif_rho_true = 0;
      if (variables_internes().is_penalized && (compteur_vimp_regul==0 || refeq_transport_2pha.non_nul()))
        {
          for (int i = 0; i < nfaces; i++)
            {
              if (compteur_vimp_regul!=0)
                {
                  for (int k = 0; k < nb_eqs ; k++)
                    {
                      REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
                      Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
                      const DoubleTab& indicatrice_faces = eq_transport.get_indicatrice_faces().valeurs();
                      if (eq_transport.get_vimp_regul()==1 &&
                          indicatrice_faces(i) > 0.0 && indicatrice_faces(i)!= 1.) terme_mul(i) += 1.0/variables_internes().eta;
                    }
                }
              rho_faces_modifie(i) *= terme_mul(i);
            }
          rho_faces_modifie.echange_espace_virtuel();
          modif_rho_true = 1;
        }
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco rho_faces_modifie:",rho_faces_modifie);

      // Assemblage de la matrice INTEGRALE            ( div(1/rho * grad(P)) ) :
      //                          (volume de controle pression)

      assembleur_pression().valeur().assembler_rho_variable(matrice_pression_,champ_rho_faces_modifie.valeur());

      //Penalization L2 de la pression si necessaire
      if (modif_rho_true == 1 && variables_internes().p_ref_pena != -1.e40)
        {
          // On se base sur le nombre de composantes par faces pour la discretisation
          Matrice_Morse_Sym& matrice_valeurs = (vpoint.line_size() == 1
                                                ? ref_cast(Matrice_Morse_Sym, (ref_cast(Matrice_Bloc, matrice_pression_.valeur())).get_bloc(0,0).valeur()) // VDF (1)
                                                : ref_cast(Matrice_Morse_Sym, matrice_pression_.valeur())) ;                                              // VEF (>1)
          DoubleTab& pressu = la_pression.valeurs();
          assert(nb_elem == champ_rho_elem_.valeur().valeurs().dimension(0));
          const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
          const Domaine& mon_dom = domaine_dis().domaine ();
          const IntTab& elem_faces = domaine_vf.elem_faces();
          const int   nb_faces_elem = elem_faces.line_size();
          const int   nb_sommet = mon_dom.nb_som_elem();
          int numero_global_som, ligne_mat;
          int point_fluide_dirichlet=-1;
          if (variables_internes().is_pfl_flottant)
            {
              if( Objet_U::dimension == 3 )
                {
                  point_fluide_dirichlet = mon_dom.chercher_elements( variables_internes().x_pfl_imp , variables_internes().y_pfl_imp , variables_internes().z_pfl_imp );
                }
              else
                {
                  point_fluide_dirichlet = mon_dom.chercher_elements( variables_internes().x_pfl_imp , variables_internes().y_pfl_imp );
                }
              if (mp_max(point_fluide_dirichlet)==-1)
                {
                  Cerr << "Point de reference pression fluide situe en dehors du domaine !" << finl;
                  Process::exit();
                }
            }
          for (int e = 0; e < nb_elem; e++)
            {
              int nbfglob = 0;
              int nbfpena = 0;
              for ( int f = 0; f < nb_faces_elem; f++ )
                {
                  int fglob = elem_faces(e,f);
                  if (fglob >= 0)
                    {
                      nbfglob += 1;
                      int dejafait =0 ;
                      for (int k = 0; k < nb_eqs_bis && dejafait == 0 ; k++)
                        {
                          REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
                          const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_indicatrice_faces().valeurs();
                          if (indicatrice_faces(fglob) > 0.0)
                            {
                              dejafait = 1 ;
                              nbfpena += 1 ;
                            }
                        }
                    }
                }
              if (nbfpena == nbfglob && nbfpena != 0)
                {
                  matrice_valeurs(e,e) += 1.0 / variables_internes().eta;
                  secmem(e) +=  variables_internes().p_ref_pena / variables_internes().eta;
                  pressu(e) = variables_internes().p_ref_pena ;
                  if (vpoint.line_size() > 1) // VEF
                    {
                      for (int somloc = 0; somloc < nb_sommet; somloc ++)
                        {
                          numero_global_som = mon_dom.sommet_elem(e,somloc);
                          ligne_mat = nb_elem+numero_global_som;
                          //                         matrice_valeurs(ligne_mat,ligne_mat) += 1.0 / variables_internes().eta;
                          pressu(ligne_mat) = 0.0 ;
                        }
                    }
                }
              if (point_fluide_dirichlet == e)
                {
                  if (nbfpena != nbfglob && nbfpena != 0)
                    {
                      // On impose une reference de pression (p_ref)sur un bord
                      secmem(e) +=  matrice_valeurs(e,e) * variables_internes().p_ref_pena / (float(nbfglob-nbfpena));
                      matrice_valeurs(e,e) *= float(nbfglob-nbfpena+1)/(float(nbfglob-nbfpena));
                    }
                  else
                    {
                      Cerr<<"Point de reference pression fluide non situe dans une cellule fluide voisin d'une IBC !"<<finl;
                      Cerr<<"Nb faces IBC = "<<nbfpena<<finl;
                      Process::exit();
                    }
                }
            }
          secmem.echange_espace_virtuel();
          pressu.echange_espace_virtuel();
        }

      // On a modifie la matrice, il faut reinitialiser le solveur
      //  (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
      solveur_pression().valeur().reinit();
    }

  assembleur_pression_.valeur().modifier_secmem(secmem);

  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),
                                     secmem,
                                     la_pression.valeurs()
                                    );
  assembleur_pression_.modifier_solution(la_pression->valeurs());
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)
  gradient.calculer(la_pression.valeur().valeurs(), gradP);
  solveur_masse.appliquer(gradP);

  // Correction de vpoint :
  if (projection_a_faire()) // Temporaire pour permettre de ne pas resoudre NS avec mettant operateurs nuls et projection_initiale 0
    {
      const int nbis = vpoint.dimension(0);
      const int mbis = vpoint.line_size();
      for (int i = 0; i < nbis; i++)
        for (int j = 0; j < mbis; j++)
          vpoint(i,j) -= gradP(i,j) / rho_faces_modifie(i);

      vpoint.echange_espace_virtuel();
    }

  // Calcul des efforts sur les IBCs et impression dans un fichier
  if (interf_vitesse_imposee_ok && limpr())
    {
      const DoubleTab& rho = champ_rho_faces_.valeur().valeurs();

      DoubleTrav forces_tot_2(vpoint) ;
      DoubleTrav pressure_part(vpoint) ;
      DoubleTrav diffusion_part(vpoint) ;

      //Calcul de la vitesse au temps n+1
      DoubleTab vv(vpoint) ;
      vv *= schema_temps().pas_de_temps() ;
      vv += inconnue().valeur().valeurs() ;

      // Terme de diffusion
      terme_diffusif.calculer(vv, variables_internes().terme_diffusion.valeur().valeurs());
      variables_internes().terme_diffusion.valeur().valeurs().echange_espace_virtuel() ;
      solveur_masse.appliquer(variables_internes().terme_diffusion.valeur().valeurs());
      const DoubleTab& diffusion = variables_internes().terme_diffusion.valeur().valeurs();
      // Terme de convection
      DoubleTrav trav(variables_internes().terme_convection.valeur().valeurs());
      derivee_en_temps_conv( trav, la_vitesse.valeurs());
      variables_internes().terme_convection.valeur().valeurs()=trav;
      solveur_masse.appliquer(variables_internes().terme_convection.valeur().valeurs());
      const DoubleTab& convection = variables_internes().terme_convection.valeur().valeurs();
      const int nbis = vpoint.dimension(0);
      const int mbis = vpoint.line_size();

      for (int i = 0; i<nbis; i++)
        for (int j = 0; j < mbis; j++)
          {
            pressure_part(i,j) = gradP(i,j) ;
            diffusion_part(i,j) = -rho(i) * convection(i,j) - diffusion(i,j) ;
            forces_tot_2(i,j) = pressure_part(i,j) + diffusion_part(i,j) ;
          }

      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k=0; k<nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
          // Impression des efforts
          eq_transport.impr_effort_fluide_interface( forces_tot_2, pressure_part, diffusion_part ) ;
        }
    }

  return vpoint;
}


const Probleme_FT_Disc_gen& Navier_Stokes_FT_Disc::probleme_ft() const
{
  return probleme_ft_.valeur();
}

Probleme_FT_Disc_gen& Navier_Stokes_FT_Disc::probleme_ft()
{
  return probleme_ft_.valeur();
}

/*! @brief In Front Tracking, pression is in Pa and so pression_pa field <=> pression field
 *
 */
void Navier_Stokes_FT_Disc::calculer_la_pression_en_pa()
{
  la_pression_en_pa.valeurs()=la_pression.valeurs();
}

const Navier_Stokes_FT_Disc_interne& Navier_Stokes_FT_Disc::variables_internes() const
{
  return *variables_internes_;
}
Navier_Stokes_FT_Disc_interne& Navier_Stokes_FT_Disc::variables_internes()
{
  return *variables_internes_;
}

/*! @brief Si le champ de vitesse est discontinu (calcul avec changement de phase), renvoie un pointeur vers le champ delta_v de "discontinuite", tel que
 *
 *   inconnue - delta_v = vitesse de deplacement des interfaces
 *   (voir Transport_Interfaces_FT_Disc::deplacer_maillage_v_fluide())
 *   Si pas de changement de phase, renvoie un pointeur nul.
 *
 */
const Champ_base *  Navier_Stokes_FT_Disc::get_delta_vitesse_interface() const
{
  if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
    return & (variables_internes().delta_u_interface.valeur());
  else
    return 0;
}
// Description : calcul de div(n) (la courbure discretisee sur le maillage volumique)
//  Faudrait deplacer cette methode dans transport interfaces...
const Champ_base& Navier_Stokes_FT_Disc::calculer_div_normale_interface()
{
  REF(Transport_Interfaces_FT_Disc) & refeq_transport =  variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  // Distance a l'interface discretisee aux elements:
  const DoubleTab& dist = eq_transport.get_update_distance_interface().valeurs();

  DoubleTab& phi = variables_internes().laplacien_d.valeur().valeurs();
  DoubleTab u0(inconnue().valeurs());

  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_div_normale", 0);
  //  statistiques().begin_count(count);

  phi = dist;

  const int n = phi.dimension(0);
  for (int i = 0; i < n; i++)
    {
      if (phi(i) < -1e20)
        phi(i) = 0;
    }
  phi.echange_espace_virtuel();
  //  static const Stat_Counter_Id count2 = statistiques().new_counter(1, "calculer_gradient", 0);
  //  statistiques().begin_count(count2);
  gradient.calculer(phi, u0);
  correct_at_exit_bad_gradient(u0);
  //  statistiques().end_count(count2);
  //  static const Stat_Counter_Id count4 = statistiques().new_counter(1, "calculer_solveur_masse", 0);
  //  statistiques().begin_count(count4);
  solveur_masse.appliquer(u0);
  //  statistiques().end_count(count4);

  // Calcul de integrale(div(u0)) sur les mailles:
  //  static const Stat_Counter_Id count3 = statistiques().new_counter(1, "calculer_div", 0);
  //  statistiques().begin_count(count3);
  divergence.calculer(u0, phi);
  //  statistiques().end_count(count3);
  // Division par le volume des mailles:
  const DoubleVect& volumes = ref_cast(Domaine_VF, domaine_dis().valeur()).volumes();
  for (int i = 0; i < n; i++)
    {
      const double p = phi(i);
      if (p != 0.)
        {
          const double v = volumes[i];
          phi(i) = p / v;
        }
    }

  //  statistiques().end_count(count);

  return variables_internes().laplacien_d.valeur();
}

const Champ_Fonc& Navier_Stokes_FT_Disc::champ_rho_faces() const
{
  return champ_rho_faces_;
}

//Renvoie 1 si l option GRAVITE_RHO_G est activee 0 sinon
int Navier_Stokes_FT_Disc::is_terme_gravite_rhog() const
{
  if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
    return 1;
  else
    return 0;
}

const Champ_Fonc& Navier_Stokes_FT_Disc::get_num_compo() const
{
  return variables_internes().num_compo;
}


// Debut EB
const DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_interf() const
{
  return variables_internes().force_pression_interf_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_interf() const
{
  return variables_internes().force_frottements_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_interf()
{
  return variables_internes().force_pression_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_interf()
{
  return variables_internes().force_frottements_interf_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_pression_interf() const
{
  return variables_internes().pression_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_pression_interf()
{
  return variables_internes().pression_interf_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_force_tot_pression_interf() const
{
  return variables_internes().force_pression_tot_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_tot_frottements_interf() const
{
  return variables_internes().force_frottements_tot_interf_;
}
const DoubleVect& Navier_Stokes_FT_Disc::get_surface_tot_interf() const
{
  return variables_internes().surface_tot_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_tot_interf_stokes_th() const
{
  return variables_internes().force_pression_tot_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_tot_interf_stokes_th() const
{
  return variables_internes().force_frottements_tot_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_tot_interf_stokes_th_dis() const
{
  return variables_internes().force_pression_tot_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_tot_interf_stokes_th_dis() const
{
  return variables_internes().force_frottements_tot_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_force_tot_pression_interf()
{
  return variables_internes().force_pression_tot_interf_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_tot_frottements_interf()
{
  return variables_internes().force_frottements_tot_interf_;
}
DoubleVect& Navier_Stokes_FT_Disc::get_surface_tot_interf()
{
  return variables_internes().surface_tot_interf_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_tot_interf_stokes_th()
{
  return variables_internes().force_pression_tot_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_tot_interf_stokes_th()
{
  return variables_internes().force_frottements_tot_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_tot_interf_stokes_th_dis()
{
  return variables_internes().force_pression_tot_interf_stokes_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_tot_interf_stokes_th_dis()
{
  return variables_internes().force_frottements_tot_interf_stokes_th_dis_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xx_interf() const
{
  return variables_internes().sigma_xx_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xx_interf()
{
  return variables_internes().sigma_xx_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xy_interf() const
{
  return variables_internes().sigma_xy_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xy_interf()
{
  return variables_internes().sigma_xy_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xz_interf() const
{
  return variables_internes().sigma_xz_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xz_interf()
{
  return variables_internes().sigma_xz_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yx_interf() const
{
  return variables_internes().sigma_yx_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yx_interf()
{
  return variables_internes().sigma_yx_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yy_interf() const
{
  return variables_internes().sigma_yy_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yy_interf()
{
  return variables_internes().sigma_yy_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yz_interf() const
{
  return variables_internes().sigma_yz_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yz_interf()
{
  return variables_internes().sigma_yz_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zx_interf() const
{
  return variables_internes().sigma_zx_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zx_interf()
{
  return variables_internes().sigma_zx_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zy_interf() const
{
  return variables_internes().sigma_zy_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zy_interf()
{
  return variables_internes().sigma_zy_interf_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zz_interf() const
{
  return variables_internes().sigma_zz_interf_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zz_interf()
{
  return variables_internes().sigma_zz_interf_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_stokes_th() const
{
  return variables_internes().force_pression_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_stokes_th()
{
  return variables_internes().force_pression_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc:: get_force_frottements_stokes_th() const
{
  return variables_internes().force_frottements_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc:: get_force_frottements_stokes_th()
{
  return variables_internes().force_frottements_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_stokes_th_dis() const
{
  return variables_internes().force_pression_stokes_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_force_pression_stokes_th_dis()
{
  return variables_internes().force_pression_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_force_frottements_stokes_th_dis() const
{
  return variables_internes().force_frottements_stokes_th_dis_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_pression_interf_stokes_th_dis() const
{
  return variables_internes().pression_interf_stokes_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_pression_interf_stokes_th_dis()
{
  return variables_internes().pression_interf_stokes_th_dis_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xx_interf_stokes_th_dis() const
{
  return variables_internes().sigma_xx_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xx_interf_stokes_th_dis()
{
  return variables_internes().sigma_xx_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xy_interf_stokes_th_dis() const
{
  return variables_internes().sigma_xy_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xy_interf_stokes_th_dis()
{
  return variables_internes().sigma_xy_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xz_interf_stokes_th_dis() const
{
  return variables_internes().sigma_xz_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xz_interf_stokes_th_dis()
{
  return variables_internes().sigma_xz_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yx_interf_stokes_th_dis() const
{
  return variables_internes().sigma_yx_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yx_interf_stokes_th_dis()
{
  return variables_internes().sigma_yx_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yy_interf_stokes_th_dis() const
{
  return variables_internes().sigma_yy_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yy_interf_stokes_th_dis()
{
  return variables_internes().sigma_yy_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yz_interf_stokes_th_dis() const
{
  return variables_internes().sigma_yz_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yz_interf_stokes_th_dis()
{
  return variables_internes().sigma_yz_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zx_interf_stokes_th_dis() const
{
  return variables_internes().sigma_zx_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zx_interf_stokes_th_dis()
{
  return variables_internes().sigma_zx_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zy_interf_stokes_th_dis() const
{
  return variables_internes().sigma_zy_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zy_interf_stokes_th_dis()
{
  return variables_internes().sigma_zy_interf_stokes_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zz_interf_stokes_th_dis() const
{
  return variables_internes().sigma_zz_interf_stokes_th_dis_;
}

DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zz_interf_stokes_th_dis()
{
  return variables_internes().sigma_zz_interf_stokes_th_dis_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xx_interf_stokes_th() const
{
  return variables_internes().sigma_xx_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xx_interf_stokes_th()
{
  return variables_internes().sigma_xx_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xy_interf_stokes_th() const
{
  return variables_internes().sigma_xy_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xy_interf_stokes_th()
{
  return variables_internes().sigma_xy_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xz_interf_stokes_th() const
{
  return variables_internes().sigma_xz_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_sigma_xz_interf_stokes_th()
{
  return variables_internes().sigma_xz_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yy_interf_stokes_th() const
{
  return variables_internes().sigma_yy_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yy_interf_stokes_th()
{
  return variables_internes().sigma_yy_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yz_interf_stokes_th() const
{
  return variables_internes().sigma_yz_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_sigma_yz_interf_stokes_th()
{
  return variables_internes().sigma_yz_interf_stokes_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zz_interf_stokes_th() const
{
  return variables_internes().sigma_zz_interf_stokes_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_sigma_zz_interf_stokes_th()
{
  return variables_internes().sigma_zz_interf_stokes_th_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P1() const
{
  return variables_internes().dUdx_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P1()
{
  return variables_internes().dUdx_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P1() const
{
  return variables_internes().dUdy_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P1()
{
  return variables_internes().dUdy_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P1() const
{
  return variables_internes().dUdz_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P1()
{
  return variables_internes().dUdz_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P1() const
{
  return variables_internes().dVdx_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P1()
{
  return variables_internes().dVdx_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P1() const
{
  return variables_internes().dVdy_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P1()
{
  return variables_internes().dVdy_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P1() const
{
  return variables_internes().dVdz_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P1()
{
  return variables_internes().dVdz_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P1() const
{
  return variables_internes().dWdx_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P1()
{
  return variables_internes().dWdx_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P1() const
{
  return variables_internes().dWdy_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P1()
{
  return variables_internes().dWdy_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P1() const
{
  return variables_internes().dWdz_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P1()
{
  return variables_internes().dWdz_P1_;
}


const DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P2() const
{
  return variables_internes().dUdx_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P2()
{
  return variables_internes().dUdx_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P2() const
{
  return variables_internes().dUdy_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P2()
{
  return variables_internes().dUdy_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P2() const
{
  return variables_internes().dUdz_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P2()
{
  return variables_internes().dUdz_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P2() const
{
  return variables_internes().dVdx_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P2()
{
  return variables_internes().dVdx_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P2() const
{
  return variables_internes().dVdy_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P2()
{
  return variables_internes().dVdy_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P2() const
{
  return variables_internes().dVdz_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P2()
{
  return variables_internes().dVdz_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P2() const
{
  return variables_internes().dWdx_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P2()
{
  return variables_internes().dWdx_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P2() const
{
  return variables_internes().dWdy_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P2()
{
  return variables_internes().dWdy_P2_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P2() const
{
  return variables_internes().dWdz_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P2()
{
  return variables_internes().dWdz_P2_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P1_th_dis() const
{
  return variables_internes().dUdx_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P1_th_dis()
{
  return variables_internes().dUdx_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P1_th_dis() const
{
  return variables_internes().dUdy_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P1_th_dis()
{
  return variables_internes().dUdy_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P1_th_dis() const
{
  return variables_internes().dUdz_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P1_th_dis()
{
  return variables_internes().dUdz_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P1_th_dis() const
{
  return variables_internes().dVdx_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P1_th_dis()
{
  return variables_internes().dVdx_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P1_th_dis() const
{
  return variables_internes().dVdy_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P1_th_dis()
{
  return variables_internes().dVdy_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P1_th_dis() const
{
  return variables_internes().dVdz_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P1_th_dis()
{
  return variables_internes().dVdz_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P1_th_dis() const
{
  return variables_internes().dWdx_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P1_th_dis()
{
  return variables_internes().dWdx_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P1_th_dis() const
{
  return variables_internes().dWdy_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P1_th_dis()
{
  return variables_internes().dWdy_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P1_th_dis() const
{
  return variables_internes().dWdz_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P1_th_dis()
{
  return variables_internes().dWdz_P1_th_dis_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P2_th_dis() const
{
  return variables_internes().dUdx_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P2_th_dis()
{
  return variables_internes().dUdx_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P2_th_dis() const
{
  return variables_internes().dUdy_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P2_th_dis()
{
  return variables_internes().dUdy_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P2_th_dis() const
{
  return variables_internes().dUdz_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P2_th_dis()
{
  return variables_internes().dUdz_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P2_th_dis() const
{
  return variables_internes().dVdx_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P2_th_dis()
{
  return variables_internes().dVdx_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P2_th_dis() const
{
  return variables_internes().dVdy_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P2_th_dis()
{
  return variables_internes().dVdy_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P2_th_dis() const
{
  return variables_internes().dVdz_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P2_th_dis()
{
  return variables_internes().dVdz_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P2_th_dis() const
{
  return variables_internes().dWdx_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P2_th_dis()
{
  return variables_internes().dWdx_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P2_th_dis() const
{
  return variables_internes().dWdy_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P2_th_dis()
{
  return variables_internes().dWdy_P2_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P2_th_dis() const
{
  return variables_internes().dWdz_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P2_th_dis()
{
  return variables_internes().dWdz_P2_th_dis_;
}


const DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P1_th() const
{
  return variables_internes().dUdx_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P1_th()
{
  return variables_internes().dUdx_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P1_th() const
{
  return variables_internes().dUdy_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P1_th()
{
  return variables_internes().dUdy_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P1_th() const
{
  return variables_internes().dUdz_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P1_th()
{
  return variables_internes().dUdz_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P1_th() const
{
  return variables_internes().dVdx_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P1_th()
{
  return variables_internes().dVdx_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P1_th() const
{
  return variables_internes().dVdy_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P1_th()
{
  return variables_internes().dVdy_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P1_th() const
{
  return variables_internes().dVdz_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P1_th()
{
  return variables_internes().dVdz_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P1_th() const
{
  return variables_internes().dWdx_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P1_th()
{
  return variables_internes().dWdx_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P1_th() const
{
  return variables_internes().dWdy_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P1_th()
{
  return variables_internes().dWdy_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P1_th() const
{
  return variables_internes().dWdz_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P1_th()
{
  return variables_internes().dWdz_P1_th_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P2_th() const
{
  return variables_internes().dUdx_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdx_P2_th()
{
  return variables_internes().dUdx_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P2_th() const
{
  return variables_internes().dUdy_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdy_P2_th()
{
  return variables_internes().dUdy_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P2_th() const
{
  return variables_internes().dUdz_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dUdz_P2_th()
{
  return variables_internes().dUdz_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P2_th() const
{
  return variables_internes().dVdx_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdx_P2_th()
{
  return variables_internes().dVdx_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P2_th() const
{
  return variables_internes().dVdy_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdy_P2_th()
{
  return variables_internes().dVdy_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P2_th() const
{
  return variables_internes().dVdz_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dVdz_P2_th()
{
  return variables_internes().dVdz_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P2_th() const
{
  return variables_internes().dWdx_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdx_P2_th()
{
  return variables_internes().dWdx_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P2_th() const
{
  return variables_internes().dWdy_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdy_P2_th()
{
  return variables_internes().dWdy_P2_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P2_th() const
{
  return variables_internes().dWdz_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_dWdz_P2_th()
{
  return variables_internes().dWdz_P2_th_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_U_P1() const
{
  return variables_internes().U_P1_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_U_P1()
{
  return variables_internes().U_P1_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_U_P2() const
{
  return variables_internes().U_P2_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_U_P2()
{
  return variables_internes().U_P2_;
}


const DoubleTab& Navier_Stokes_FT_Disc::get_U_P1_th() const
{
  return variables_internes().U_P1_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_U_P1_th()
{
  return variables_internes().U_P1_th_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_U_P2_th() const
{
  return variables_internes().U_P2_th_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_U_P2_th()
{
  return variables_internes().U_P2_th_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_U_P1_th_dis() const
{
  return variables_internes().U_P1_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_U_P1_th_dis()
{
  return variables_internes().U_P1_th_dis_;
}
const DoubleTab& Navier_Stokes_FT_Disc::get_U_P2_th_dis() const
{
  return variables_internes().U_P2_th_dis_;
}
DoubleTab& Navier_Stokes_FT_Disc::get_U_P2_th_dis()
{
  return variables_internes().U_P2_th_dis_;
}

const IntTab& Navier_Stokes_FT_Disc::get_list_elem_P1() const
{
  return variables_internes().list_elem_P1_;
}
IntTab& Navier_Stokes_FT_Disc::get_list_elem_P1()
{
  return variables_internes().list_elem_P1_;
}
const IntTab& Navier_Stokes_FT_Disc::get_list_elem_diph() const
{
  return variables_internes().list_elem_diph_;
}
IntTab& Navier_Stokes_FT_Disc::get_list_elem_diph()
{
  return variables_internes().list_elem_diph_;
}
const IntTab& Navier_Stokes_FT_Disc::get_list_elem_P1_all() const
{
  return variables_internes().list_elem_P1_all_;
}
IntTab& Navier_Stokes_FT_Disc::get_list_elem_P1_all()
{
  return variables_internes().list_elem_P1_all_;
}

// debut EB

// On initialise les champs a 0 pour le premier postraitement
void Navier_Stokes_FT_Disc::init_champs_forces_interf()
{
  REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const int nb_fa7 = maillage.nb_facettes();
  if (eq_transport.postraitement_forces_interf().flag_force_pression_facettes_)
    {
      variables_internes().force_pression_interf_.resize(nb_fa7,dimension);
      variables_internes().force_pression_interf_=3.14e31;
    }

  if (eq_transport.postraitement_forces_interf().flag_force_frottements_facettes_)
    {
      variables_internes().force_frottements_interf_.resize(nb_fa7,dimension);
      variables_internes().force_frottements_interf_=3.14e31;
    }

  if (eq_transport.postraitement_forces_interf().flag_pression_facettes_)
    {
      variables_internes().pression_interf_.resize(nb_fa7);
      variables_internes().pression_interf_=3.14e31;
    }

  if (eq_transport.postraitement_forces_interf().flag_tenseur_contraintes_facettes_)
    {
      variables_internes().sigma_xx_interf_.resize(nb_fa7);
      variables_internes().sigma_xx_interf_=0;
      variables_internes().sigma_xy_interf_.resize(nb_fa7);
      variables_internes().sigma_xy_interf_=0;
      variables_internes().sigma_xz_interf_.resize(nb_fa7);
      variables_internes().sigma_xz_interf_=0;
      variables_internes().sigma_yx_interf_.resize(nb_fa7);
      variables_internes().sigma_yz_interf_=0;
      variables_internes().sigma_yy_interf_.resize(nb_fa7);
      variables_internes().sigma_yy_interf_=0;
      variables_internes().sigma_yz_interf_.resize(nb_fa7);
      variables_internes().sigma_yz_interf_=0;
      variables_internes().sigma_zx_interf_.resize(nb_fa7);
      variables_internes().sigma_zx_interf_=0;
      variables_internes().sigma_zy_interf_.resize(nb_fa7);
      variables_internes().sigma_zy_interf_=0;
      variables_internes().sigma_zz_interf_.resize(nb_fa7);
      variables_internes().sigma_zz_interf_=0;
    }
  if (eq_transport.postraitement_forces_interf().calcul_forces_theoriques_stokes_)
    {
      variables_internes().force_pression_stokes_th_.resize(nb_fa7,dimension);
      variables_internes().force_pression_stokes_th_=3.14e31;
      variables_internes().force_frottements_stokes_th_.resize(nb_fa7,dimension);
      variables_internes().force_frottements_stokes_th_=3.14e31;
      variables_internes().force_pression_stokes_th_dis_.resize(nb_fa7,dimension);
      variables_internes().force_pression_stokes_th_dis_=3.14e31;
      variables_internes().force_frottements_stokes_th_dis_.resize(nb_fa7,dimension);
      variables_internes().force_frottements_stokes_th_dis_=3.14e31;
      variables_internes().pression_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().pression_interf_stokes_th_dis_=3.14e31;
    }
  if (eq_transport.postraitement_forces_interf().flag_tenseur_contraintes_facettes_ &&
      eq_transport.postraitement_forces_interf().calcul_forces_theoriques_stokes_)
    {
      variables_internes().sigma_xx_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_xx_interf_stokes_th_dis_=0;
      variables_internes().sigma_xy_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_xy_interf_stokes_th_dis_=0;
      variables_internes().sigma_xz_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_xz_interf_stokes_th_dis_=0;
      variables_internes().sigma_yx_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_yx_interf_stokes_th_dis_=0;
      variables_internes().sigma_yy_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_yy_interf_stokes_th_dis_=0;
      variables_internes().sigma_yz_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_yz_interf_stokes_th_dis_=0;
      variables_internes().sigma_zx_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_zx_interf_stokes_th_dis_=0;
      variables_internes().sigma_zy_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_zy_interf_stokes_th_dis_=0;
      variables_internes().sigma_zz_interf_stokes_th_dis_.resize(nb_fa7);
      variables_internes().sigma_zz_interf_stokes_th_dis_=0;

      variables_internes().sigma_xx_interf_stokes_th_.resize(nb_fa7);
      variables_internes().sigma_xx_interf_stokes_th_=0;
      variables_internes().sigma_xy_interf_stokes_th_.resize(nb_fa7);
      variables_internes().sigma_xy_interf_stokes_th_=0;
      variables_internes().sigma_xz_interf_stokes_th_.resize(nb_fa7);
      variables_internes().sigma_xz_interf_stokes_th_=0;
      variables_internes().sigma_yy_interf_stokes_th_.resize(nb_fa7);
      variables_internes().sigma_yy_interf_stokes_th_=0;
      variables_internes().sigma_yz_interf_stokes_th_.resize(nb_fa7);
      variables_internes().sigma_yz_interf_stokes_th_=0;
      variables_internes().sigma_zz_interf_stokes_th_.resize(nb_fa7);
      variables_internes().sigma_zz_interf_stokes_th_=0;

      variables_internes().dUdx_P1_.resize(nb_fa7);
      variables_internes().dUdx_P1_=0;
      variables_internes().dUdy_P1_.resize(nb_fa7);
      variables_internes().dUdy_P1_=0;
      variables_internes().dUdz_P1_.resize(nb_fa7);
      variables_internes().dUdz_P1_=0;
      variables_internes().dVdx_P1_.resize(nb_fa7);
      variables_internes().dVdx_P1_=0;
      variables_internes().dVdy_P1_.resize(nb_fa7);
      variables_internes().dVdy_P1_=0;
      variables_internes().dVdz_P1_.resize(nb_fa7);
      variables_internes().dVdz_P1_=0;
      variables_internes().dWdx_P1_.resize(nb_fa7);
      variables_internes().dWdx_P1_=0;
      variables_internes().dWdy_P1_.resize(nb_fa7);
      variables_internes().dWdy_P1_=0;
      variables_internes().dWdz_P1_.resize(nb_fa7);
      variables_internes().dWdz_P1_=0;

      variables_internes().dUdx_P2_.resize(nb_fa7);
      variables_internes().dUdx_P2_=0;
      variables_internes().dUdy_P2_.resize(nb_fa7);
      variables_internes().dUdy_P2_=0;
      variables_internes().dUdz_P2_.resize(nb_fa7);
      variables_internes().dUdz_P2_=0;
      variables_internes().dVdx_P2_.resize(nb_fa7);
      variables_internes().dVdx_P2_=0;
      variables_internes().dVdy_P2_.resize(nb_fa7);
      variables_internes().dVdy_P2_=0;
      variables_internes().dVdz_P2_.resize(nb_fa7);
      variables_internes().dVdz_P2_=0;
      variables_internes().dWdx_P2_.resize(nb_fa7);
      variables_internes().dWdx_P2_=0;
      variables_internes().dWdy_P2_.resize(nb_fa7);
      variables_internes().dWdy_P2_=0;
      variables_internes().dWdz_P2_.resize(nb_fa7);
      variables_internes().dWdz_P2_=0;

      variables_internes().dUdx_P1_th_dis_.resize(nb_fa7);
      variables_internes().dUdx_P1_th_dis_=0;
      variables_internes().dUdy_P1_th_dis_.resize(nb_fa7);
      variables_internes().dUdy_P1_th_dis_=0;
      variables_internes().dUdz_P1_th_dis_.resize(nb_fa7);
      variables_internes().dUdz_P1_th_dis_=0;
      variables_internes().dVdx_P1_th_dis_.resize(nb_fa7);
      variables_internes().dVdx_P1_th_dis_=0;
      variables_internes().dVdy_P1_th_dis_.resize(nb_fa7);
      variables_internes().dVdy_P1_th_dis_=0;
      variables_internes().dVdz_P1_th_dis_.resize(nb_fa7);
      variables_internes().dVdz_P1_th_dis_=0;
      variables_internes().dWdx_P1_th_dis_.resize(nb_fa7);
      variables_internes().dWdx_P1_th_dis_=0;
      variables_internes().dWdy_P1_th_dis_.resize(nb_fa7);
      variables_internes().dWdy_P1_th_dis_=0;
      variables_internes().dWdz_P1_th_dis_.resize(nb_fa7);
      variables_internes().dWdz_P1_th_dis_=0;

      variables_internes().dUdx_P2_th_dis_.resize(nb_fa7);
      variables_internes().dUdx_P2_th_dis_=0;
      variables_internes().dUdy_P2_th_dis_.resize(nb_fa7);
      variables_internes().dUdy_P2_th_dis_=0;
      variables_internes().dUdz_P2_th_dis_.resize(nb_fa7);
      variables_internes().dUdz_P2_th_dis_=0;
      variables_internes().dVdx_P2_th_dis_.resize(nb_fa7);
      variables_internes().dVdx_P2_th_dis_=0;
      variables_internes().dVdy_P2_th_dis_.resize(nb_fa7);
      variables_internes().dVdy_P2_th_dis_=0;
      variables_internes().dVdz_P2_th_dis_.resize(nb_fa7);
      variables_internes().dVdz_P2_th_dis_=0;
      variables_internes().dWdx_P2_th_dis_.resize(nb_fa7);
      variables_internes().dWdx_P2_th_dis_=0;
      variables_internes().dWdy_P2_th_dis_.resize(nb_fa7);
      variables_internes().dWdy_P2_th_dis_=0;
      variables_internes().dWdz_P2_th_dis_.resize(nb_fa7);
      variables_internes().dWdz_P2_th_dis_=0;

      variables_internes().dUdx_P1_th_.resize(nb_fa7);
      variables_internes().dUdx_P1_th_=0;
      variables_internes().dUdy_P1_th_.resize(nb_fa7);
      variables_internes().dUdy_P1_th_=0;
      variables_internes().dUdz_P1_th_.resize(nb_fa7);
      variables_internes().dUdz_P1_th_=0;
      variables_internes().dVdx_P1_th_.resize(nb_fa7);
      variables_internes().dVdx_P1_th_=0;
      variables_internes().dVdy_P1_th_.resize(nb_fa7);
      variables_internes().dVdy_P1_th_=0;
      variables_internes().dVdz_P1_th_.resize(nb_fa7);
      variables_internes().dVdz_P1_th_=0;
      variables_internes().dWdx_P1_th_.resize(nb_fa7);
      variables_internes().dWdx_P1_th_=0;
      variables_internes().dWdy_P1_th_.resize(nb_fa7);
      variables_internes().dWdy_P1_th_=0;
      variables_internes().dWdz_P1_th_.resize(nb_fa7);
      variables_internes().dWdz_P1_th_=0;

      variables_internes().dUdx_P2_th_.resize(nb_fa7);
      variables_internes().dUdx_P2_th_=0;
      variables_internes().dUdy_P2_th_.resize(nb_fa7);
      variables_internes().dUdy_P2_th_=0;
      variables_internes().dUdz_P2_th_.resize(nb_fa7);
      variables_internes().dUdz_P2_th_=0;
      variables_internes().dVdx_P2_th_.resize(nb_fa7);
      variables_internes().dVdx_P2_th_=0;
      variables_internes().dVdy_P2_th_.resize(nb_fa7);
      variables_internes().dVdy_P2_th_=0;
      variables_internes().dVdz_P2_th_.resize(nb_fa7);
      variables_internes().dVdz_P2_th_=0;
      variables_internes().dWdx_P2_th_.resize(nb_fa7);
      variables_internes().dWdx_P2_th_=0;
      variables_internes().dWdy_P2_th_.resize(nb_fa7);
      variables_internes().dWdy_P2_th_=0;
      variables_internes().dWdz_P2_th_.resize(nb_fa7);
      variables_internes().dWdz_P2_th_=0;
    }
  if (eq_transport.postraitement_forces_interf().calcul_forces_)
    {
      variables_internes().U_P1_.resize(nb_fa7,dimension);
      variables_internes().U_P1_=0;
      variables_internes().U_P2_.resize(nb_fa7,dimension);
      variables_internes().U_P2_=0;
    }
  if (eq_transport.postraitement_forces_interf().calcul_forces_theoriques_stokes_)
    {
      variables_internes().U_P1_th_.resize(nb_fa7,dimension);
      variables_internes().U_P1_th_=0;
      variables_internes().U_P2_th_.resize(nb_fa7,dimension);
      variables_internes().U_P2_th_=0;
      variables_internes().U_P1_th_dis_.resize(nb_fa7,dimension);
      variables_internes().U_P1_th_dis_=0;
      variables_internes().U_P2_th_dis_.resize(nb_fa7,dimension);
      variables_internes().U_P2_th_dis_=0;
    }
}

/*
void Navier_Stokes_FT_Disc::corriger_mpoint()
{
  if (probleme_ft().tcl().is_activated())
    {
      DoubleTab& mpoint = variables_internes().mpoint.valeur().valeurs();
      probleme_ft().tcl().corriger_mpoint(mpoint);
    }
}
*/
