/*
 * Domaine_VDF_plus.h
 *
 *  Created on: May 20, 2023
 *      Author: edouard
 */


#ifndef DOMAINE_VDF_PLUS_included
#define DOMAINE_VDF_PLUS_included

#include <Domaine_VDF.h>

//////////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION class Domaine_VDF_plus
// Classe fille de Domaine_VDF incluant une structure parallele aux aretes
//////////////////////////////////////////////////////////////////////////////////

class Domaine_VDF_plus : public Domaine_VDF
{

  Declare_instanciable(Domaine_VDF_plus);

public :
  void discretiser() override;
  void modifier_pour_Cl(const Conds_lim& cl) override;

private:
  void genere_et_cree_aretes(); // EB : pour construire la structure parallelle
  void calcul_xa(); // EB
  void remplir_volumes_aretes(); // EB
  void cree_aretes_virtuelles(IntTab& Aretes_som, IntTab& Elem_Aretes, Aretes& les_aretes); // EB

};


#endif
