/*
 * Particule_Solide.h
 *
 *  Created on: Apr 1, 2022
 *      Author: edouard
 */

#ifndef PARTICULE_SOLIDE_included
#define PARTICULE_SOLIDE_included
#include <Fluide_Incompressible.h>

class Particule_Solide : public Fluide_Incompressible
{
  Declare_instanciable_sans_constructeur(Particule_Solide);

public:
  Particule_Solide();
  void set_param(Param& param) override;
  const double& e_dry() const;
  void set_rayon(const double rayon) const;
  void set_diametre(const double rayon) const;
  void set_volume_compo(const double volume) const;
  void set_masse_compo(const double masse) const;
  const double& rayon_frottements() const;
  const double& rayon_collision() const;
  const double& diametre() const;
  const double& volume_compo() const;
  const double& masse_compo() const;
  const int& monodisperse() const;
protected :
  double e_dry_;
  int monodisperse_;
  mutable double rayon_frottements_;
  mutable double rayon_collision_;
  mutable double diametre_;
  mutable double volume_compo_;
  mutable double masse_compo_;
};

#endif /* PARTICULE_SOLIDE_H_ */
