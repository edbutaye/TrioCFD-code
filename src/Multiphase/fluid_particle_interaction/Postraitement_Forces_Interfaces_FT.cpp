/*
 * PostraitementForcesInterfaces.cpp
 *
 *  Created on: Jun 29, 2022
 *      Author: edouard
 */

/////////////////////////////////////////////////////////////////////////////
//
// File      : Postraitement_Forces_Interfaces_FT.cpp
// Directory : $FPI184_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Postraitement_Forces_Interfaces_FT.h>
#include <Param.h>


Implemente_instanciable_sans_constructeur(Postraitement_Forces_Interfaces_FT,"Postraitement_Forces_Interfaces_FT",Objet_U);

Postraitement_Forces_Interfaces_FT::Postraitement_Forces_Interfaces_FT() :
  calcul_forces_(0),
  calcul_forces_theoriques_stokes_(0),
  calcul_flux_(0),
  flag_pression_facettes_(0),
  flag_force_pression_facettes_(0),
  flag_force_frottements_facettes_(0),
  flag_tenseur_contraintes_facettes_(0),
  distance_interpolation_pression_P1_(0.),
  distance_interpolation_pression_P2_(0.),
  distance_interpolation_temp_P1_(0.),
  distance_interpolation_temp_P2_(0.),
  distance_interpolation_gradU_P1_(0.),
  distance_interpolation_gradU_P2_(0.)
{

}

Entree& Postraitement_Forces_Interfaces_FT::readOn (Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  return is;
}

Sortie& Postraitement_Forces_Interfaces_FT::printOn(Sortie& os) const
{
  Cerr << "Erreur Postraitement_Forces_Interfaces_FT::printOn n'est pas code." << finl;
  assert(0);
  return os;
}

void Postraitement_Forces_Interfaces_FT::set_param(Param& p)
{
  p.ajouter("calcul_forces", &calcul_forces_);
  p.ajouter("calcul_flux", &calcul_flux_);
  p.ajouter_non_std("methode_calcul_force_pression", (this), Param::REQUIRED);
  p.ajouter_non_std("methode_calcul_force_frottements", (this), Param::REQUIRED);
  p.ajouter_non_std("localisation_tenseur_contraintes", (this));
  p.ajouter("calcul_forces_theoriques_stokes", &calcul_forces_theoriques_stokes_);
  p.ajouter("distance_interpolation_pression_P1", &distance_interpolation_pression_P1_);
  p.ajouter("distance_interpolation_pression_P2", &distance_interpolation_pression_P2_);
  p.ajouter("distance_interpolation_temp_P1", &distance_interpolation_temp_P1_);
  p.ajouter("distance_interpolation_temp_P2", &distance_interpolation_temp_P2_);
  p.ajouter("distance_interpolation_gradU_P1", &distance_interpolation_gradU_P1_);
  p.ajouter("distance_interpolation_gradU_P2", &distance_interpolation_gradU_P2_);
  p.ajouter("pression_facettes", &flag_pression_facettes_);
  p.ajouter("force_pression_facettes", &flag_force_pression_facettes_);
  p.ajouter("force_frottements_facettes", &flag_force_frottements_facettes_);
  p.ajouter("tenseur_contraintes_facettes", &flag_tenseur_contraintes_facettes_);
}


int Postraitement_Forces_Interfaces_FT::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{

  if (mot=="methode_calcul_force_pression")
    {
      Motcles mots;
      mots.add("trilineaire_lineaire");
      mots.add("taylor_lagrande_ordre_3");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading methode_calcul_force_pression : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          methode_calcul_force_pression_ = Postraitement_Forces_Interfaces_FT::TRILINEAIRE_LINEAIRE;
          break;
        case 1:
          methode_calcul_force_pression_ = Postraitement_Forces_Interfaces_FT::TAYLOR_LAGRANGE_ORDRE3_PRESSION;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot=="methode_calcul_force_frottements")
    {
      Motcles mots;
      mots.add("trilineaire_lineaire_tenseur_complet");
      mots.add("trilineaire_tenseur_projete");
      mots.add("taylor_lagrande_ordre_3");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading methode_calcul_force_frottements : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          methode_calcul_force_frottements_ = Postraitement_Forces_Interfaces_FT::TRILINEAIRE_LINEAIRE_TENSEUR_COMPLET;
          break;
        case 1:
          methode_calcul_force_frottements_ = Postraitement_Forces_Interfaces_FT::TRILINENAIRE_TENSEUR_PROJETE;
          break;
        case 2:
          methode_calcul_force_frottements_ = Postraitement_Forces_Interfaces_FT::TAYLOR_LAGRANGE_ORDRE3_FROTTEMENTS;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot=="localisation_tenseur_contraintes")
    {
      Motcles mots;
      mots.add("faces_normale_x");
      mots.add("elements");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading methode_interpolation : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          localisation_tenseur_contraintes_ = Postraitement_Forces_Interfaces_FT::FACES_NORMALE_X;
          break;
        case 1:
          localisation_tenseur_contraintes_ = Postraitement_Forces_Interfaces_FT::ELEMENTS;
          break;
        case -1:
          localisation_tenseur_contraintes_ = Postraitement_Forces_Interfaces_FT::ELEMENTS;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else
    {
      Cerr << mot << " is not a keyword understood by " << que_suis_je() << " in lire_motcle_non_standard"<< finl;
      exit();
    }
  return -1;
}


int Postraitement_Forces_Interfaces_FT::postraiter_forces()
{
  if (calcul_forces_) return 1;
  return 0;
}
int Postraitement_Forces_Interfaces_FT::postraiter_flux()
{
  if (calcul_flux_) return 1;
  return 0;
}

int Postraitement_Forces_Interfaces_FT::postraiter_forces() const
{
  if (calcul_forces_) return 1;
  return 0;
}
int Postraitement_Forces_Interfaces_FT::postraiter_flux() const
{
  if (calcul_flux_) return 1;
  return 0;
}

double Postraitement_Forces_Interfaces_FT::get_distance_interpolation_pression_P1() const
{
  return distance_interpolation_pression_P1_;
}

double Postraitement_Forces_Interfaces_FT::get_distance_interpolation_pression_P2() const
{
  return distance_interpolation_pression_P2_;
}

double Postraitement_Forces_Interfaces_FT::get_distance_interpolation_temperature_P1() const
{
  return distance_interpolation_temp_P1_;
}

double Postraitement_Forces_Interfaces_FT::get_distance_interpolation_gradU_P1() const
{
  return distance_interpolation_gradU_P1_;
}

double Postraitement_Forces_Interfaces_FT::get_distance_interpolation_gradU_P2() const
{
  return distance_interpolation_gradU_P2_;
}

double Postraitement_Forces_Interfaces_FT::get_distance_interpolation_temperature_P2() const
{
  return distance_interpolation_temp_P2_;
}


