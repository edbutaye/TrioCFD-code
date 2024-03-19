/*
 * VDF_plus_discretisation.cpp
 *
 *  Created on: May 20, 2023
 *      Author: edouard
 */

#include <VDF_plus_discretisation.h>
#include <Domaine_Cl_VDF.h>
#include <Navier_Stokes_std.h>
#include <Domaine_VDF_plus.h>

Implemente_instanciable(VDF_plus_discretisation,"VDF+",VDF_discretisation);


Entree& VDF_plus_discretisation::readOn(Entree& s)
{
  return s ;
}

Sortie& VDF_plus_discretisation::printOn(Sortie& s) const
{
  return s ;
}

void VDF_plus_discretisation::domaine_Cl_dis(Domaine_dis& d,
                                             Domaine_Cl_dis& dcl) const
{
  Cerr << "Discretisation des conditions limites" << finl;
  Domaine_VDF_plus& domaine_vdf_plus=ref_cast(Domaine_VDF_plus, d.valeur());
  dcl.typer("Domaine_Cl_VDF");
  Domaine_Cl_VDF& domaine_cl_vdf=ref_cast(Domaine_Cl_VDF, dcl.valeur());
  domaine_cl_vdf.associer(domaine_vdf_plus);
}
