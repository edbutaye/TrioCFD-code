/*
 * VDF_plus_discretisation.h
 *
 *  Created on: May 20, 2023
 *      Author: edouard
 */

#ifndef VDF_PLUS_DISCRETISATION_included
#define VDF_PLUS_DISCRETISATION_included


#include <VDF_discretisation.h>

class VDF_plus_discretisation : public VDF_discretisation
{
  Declare_instanciable(VDF_plus_discretisation);

public:

  void domaine_Cl_dis(Domaine_dis& d, Domaine_Cl_dis& dcl) const override;
};

#endif
