/*
 * PostraitementForcesInterfaces.h
 *
 *  Created on: Jun 29, 2022
 *      Author: edouard
 */
/////////////////////////////////////////////////////////////////////////////
//
// File      : Postraitement_Forces_Interfaces_FT.h
// Directory : $FPI184_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef POSTRAITEMENT_FORCES_INTERFACES_FT_included
#define POSTRAITEMENT_FORCES_INTERFACES_FT_included

#include <Objet_U.h>
#include <FTd_tools.h>
//#include <Ref_Transport_Interfaces_FT_Disc.h>
#include <MD_Vector.h>
#include <TRUSTList.h>
#include <TRUSTTabFT.h>

class Param;


class Postraitement_Forces_Interfaces_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Postraitement_Forces_Interfaces_FT);
public:
  Postraitement_Forces_Interfaces_FT();
  void set_param(Param& p);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int postraiter_forces();
  int postraiter_flux();
  int postraiter_forces() const;
  int postraiter_flux() const;
  double get_distance_interpolation_pression_P1() const;
  double get_distance_interpolation_pression_P2() const;
  double get_distance_interpolation_temperature_P1() const;
  double get_distance_interpolation_temperature_P2() const;
  double get_distance_interpolation_gradU_P1() const;
  double get_distance_interpolation_gradU_P2() const;

  int calcul_forces_;
  int calcul_forces_theoriques_stokes_;
  int calcul_flux_;
  int flag_pression_facettes_;
  int flag_force_pression_facettes_;
  int flag_force_frottements_facettes_;
  int flag_tenseur_contraintes_facettes_;


  enum Methode_calcul_force_pression { TRILINEAIRE_LINEAIRE, TAYLOR_LAGRANGE_ORDRE3_PRESSION }; // Le developpement de la methode TAYLOR_LAGRANGE_ORDRE3 n'est pas termine (et doit surement etre repris depuis le debut)
  Methode_calcul_force_pression methode_calcul_force_pression_; // pas utilise (pour le moment ?)
  enum Methode_calcul_force_frottements { TRILINEAIRE_LINEAIRE_TENSEUR_COMPLET, TRILINENAIRE_TENSEUR_PROJETE, TAYLOR_LAGRANGE_ORDRE3_FROTTEMENTS }; // Le developpement de la methode TAYLOR_LAGRANGE_ORDRE3 n'est pas termine (et doit surement etre repris depuis le debut)
  Methode_calcul_force_frottements methode_calcul_force_frottements_;
  enum Localisation_tenseur_contraintes { FACES_NORMALE_X, ELEMENTS };
  Localisation_tenseur_contraintes localisation_tenseur_contraintes_;

  //IntList liste_flag_posts_interf;
  //void set_liste_flag_posts_interf();

private:
  double distance_interpolation_pression_P1_;
  double distance_interpolation_pression_P2_;
  double distance_interpolation_temp_P1_;
  double distance_interpolation_temp_P2_;
  double distance_interpolation_gradU_P1_;
  double distance_interpolation_gradU_P2_;
};

#endif
