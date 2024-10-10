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
// File:        Mod_echelle_cv_naturelle.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Lois_Paroi/TBLE/Models
//
//////////////////////////////////////////////////////////////////////////////

#include <Mod_echelle_cv_naturelle.h>
#include <EFichier.h>
#include <SFichier.h>
#include <Fluide_Incompressible.h>


Implemente_instanciable(Mod_echelle_cv_naturelle,"Mod_echelle_convection_naturelle",Mod_echelle_LRM_base);

//     printOn()
/////

Sortie& Mod_echelle_cv_naturelle::printOn(Sortie& os) const
{
  return os << que_suis_je() << " " << le_nom();
}

//// readOn
//

Entree& Mod_echelle_cv_naturelle::readOn(Entree& is)
{

  return is;


}


double Mod_echelle_cv_naturelle::calculer_vv_bar(double y, double k, double u, double visco_cin)
{

  //Attention : vv_bar est calcule en y pour ParoiVDF_LRM et en yc pour Diffu_k !
  //Dans Mod_echelle_LRM il est calcule par defaut en y et donc a interpoler pour Diffu_k.

  double  y_star,vv_bar;
  y_star=y*sqrt(k)/visco_cin;
  vv_bar=k*(7.19*1.e-3*y_star-4.33*1.e-5*y_star*y_star+8.8*1.e-8*y_star*y_star*y_star);//convection naturelle
  return vv_bar;

}


double Mod_echelle_cv_naturelle::calculer_l_eps(double y, double k, double u, double visco_cin)
{

  //l_eps est calcule en y

  double  y_nu, l_eps, vv_bar;

  vv_bar=calculer_vv_bar(y, k, u, visco_cin);
  y_nu=y*sqrt(vv_bar)/visco_cin;
  l_eps=(8.8*y)/(1+(10/y_nu)+5.15*1.e-2*y_nu);

  return l_eps;

}

double Mod_echelle_cv_naturelle::calculer_l_mu(double y, double k, double u, double visco_cin)
{

  //l_mu est calcule en yc

  double y_nu, l_mu, vv_bar;

  vv_bar=calculer_vv_bar(y, k, u, visco_cin);
  y_nu=y*sqrt(vv_bar)/visco_cin;
  l_mu=(0.544*y)/(1+5.025*1.e-4*pow(y_nu,1.65));

  return l_mu;

}


