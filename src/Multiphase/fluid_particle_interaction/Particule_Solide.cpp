/*
 * Particule_Solide.cpp
 *
 *  Created on: Apr 1, 2022
 *      Author: edouard
 */


#include <Particule_Solide.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Particule_Solide,"Particule_Solide",Fluide_Incompressible);


Particule_Solide::Particule_Solide():
  e_dry_(0.),
  monodisperse_(1), // EB : on le met a 1 pour le moment car on ne fait que du monodisperse
  rayon_frottements_(-1.),
  rayon_collision_(-1.),
  diametre_(0.)
{}

Sortie& Particule_Solide::printOn(Sortie& os) const
{
  os << "{"                           << finl;
  os << "e_dry " << e_dry_            << finl;
  Milieu_base::ecrire(os);
  os << "}"					          << finl;
  return os;
}

Entree& Particule_Solide::readOn(Entree& is)
{
  Fluide_Incompressible::readOn(is);
  if (monodisperse_==1)
    {
      set_diametre(2*rayon_collision_);
      set_volume_compo(4 * M_PI * pow(rayon_collision_, 3) / 3);
      const double rho_solide = masse_volumique().valeurs()(0, 0);
      set_masse_compo((4 * M_PI * pow(rayon_collision_, 3) / 3)*rho_solide);
    }
  return is;
}

void Particule_Solide::set_param(Param& param)
{
  Fluide_Incompressible::set_param(param);
  //La lecture de rho n est pas specifiee obligatoire ici car ce
  //champ ne doit pas etre lu pour un fluide dilatable
  param.ajouter("e_dry",&e_dry_,Param::REQUIRED);
  param.ajouter("rayon_frottements", &rayon_frottements_,Param::REQUIRED); // XD_ADD_P radius of a spherical particle;
  param.ajouter("rayon_collision", &rayon_collision_,Param::REQUIRED); // XD_ADD_P radius of a spherical particle;
  param.ajouter("monodisperse", &monodisperse_);
}

const double& Particule_Solide::e_dry() const { return e_dry_; }
void Particule_Solide::set_rayon(const double rayon) const { rayon_collision_=rayon; }
void Particule_Solide::set_diametre(const double diametre) const { diametre_=diametre; }
void Particule_Solide::set_volume_compo(const double volume) const { volume_compo_=volume; }
void Particule_Solide::set_masse_compo(const double masse) const { masse_compo_=masse; }
const double& Particule_Solide::rayon_frottements() const { return rayon_frottements_; }
const double& Particule_Solide::rayon_collision() const { return rayon_collision_; }
const double& Particule_Solide::diametre() const { return diametre_; }
const double& Particule_Solide::volume_compo() const { return volume_compo_; }
const double& Particule_Solide::masse_compo() const { return masse_compo_; }
const int& Particule_Solide::monodisperse() const { return monodisperse_; }

