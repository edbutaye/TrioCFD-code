/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Transport_K_Eps_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_Eps_base.h>
#include <Schema_Temps_base.h>
#include <Domaine_VF.h>
#include <Domaine_VEF.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Param.h>
#include <Debog.h>

Implemente_base(Transport_K_Eps_base, "Transport_K_Eps_base", Transport_2eq_base);
// X_D Transport_K_Eps_base Transport_2eq_base Transport_K_Eps_base 1 Base equation for RANS k-eps model. Should not be used directly

/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_Eps_base::printOn(Sortie& is) const
{ return is << que_suis_je() << "\n"; }

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_Eps_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}
void Transport_K_Eps_base::set_param(Param& param)
{
  param.ajouter_flag("exit_on_negative_k_eps", &exit_on_negative_k_eps_); // XD_ADD_P flag Flag to exit (with postprocessing of fields) if a negative value is found for k or epsilon
  param.ajouter_flag("exit_on_big_eps", &exit_on_big_eps_); // XD_ADD_P flag Flag to exit (with postprocessing of fields) if an excessively big values of epsilon is found
  Transport_2eq_base::set_param(param);
}


void Transport_K_Eps_base::discretiser()
{
  if (sub_type(Discret_Thyd,discretisation()))
    {
      Cerr << "K,Eps transport equation ("<< que_suis_je() <<") discretization" << finl;
      Cerr << "K_Eps field discretization" << finl;
      Noms noms(2);
      Noms unit(2);
      noms[0]="K";
      noms[1]="eps";
      unit[0]="m2/s2";
      unit[1]="m2/s3";

      discretisation().discretiser_champ("temperature",domaine_dis(),multi_scalaire,noms,unit,2,schema_temps().nb_valeurs_temporelles(),schema_temps().temps_courant(),le_champ_K_Eps);
      le_champ_K_Eps->nommer("K_Eps");
      champs_compris_.ajoute_champ(le_champ_K_Eps);
      if (modele_turbulence().equation().calculate_time_derivative())
        {
          set_calculate_time_derivative(1);
        }

      Equation_base::discretiser();
    }
  else
    {
      Cerr<<" Transport_K_Eps_base::discretiser "<<finl;
      Cerr<<"Discretization "<<discretisation().que_suis_je()<<" not recognized."<<finl;
      Process::exit();
    }
}

/*! @brief Controle le champ inconnue K-epsilon en forcant a zero les valeurs du champ
 *
 *     inferieurs a 1.e-10.
 *
 * @return (int) renvoie toujours 1 (c'est tres utile)
 */
int Transport_K_Eps_base::controler_K_Eps()
{
  DoubleTab& K_Eps = le_champ_K_Eps->valeurs();

  // size == nb_elem in VDF or nb_faces in VEF
  int size = K_Eps.dimension(0);

  // can this situation really happen ?
  if (size < 0)
    {
      if (sub_type(Champ_Inc_P0_base, le_champ_K_Eps.valeur()))
        {
          size = le_champ_K_Eps->equation().domaine_dis().domaine().nb_elem();
        }
      else
        {
          Cerr << "Unsupported K_Eps field in Transport_K_Eps_base::controler_K_Eps()" << finl;
          Process::exit(1);
        }
    }

  // neg will store the amount of problematic values of k or eps found
  // neg[0] : amount of negative k
  // neg[1] : amount of negative eps
  // neg[2] : amount of eps too big (> modele_turbulence().get_EPS_MAX())
  ArrOfInt neg(3);
  neg = 0;


  // On interdit K-Eps negatif pour le K-Eps seulement
  // Les autres modeles (2 couches, Launder, ne sont pas assez valides)

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF,domaine_dis());
  double LeEPS_MIN = modele_turbulence().get_EPS_MIN();
  double LeEPS_MAX = modele_turbulence().get_EPS_MAX();
  double LeK_MIN = modele_turbulence().get_K_MIN();
  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces = domaine_vf.elem_faces();


  // PL on ne fixe au seuil minimum que si negatifs
  // car la loi de paroi peut fixer a des valeurs tres petites
  // et le rapport K*K/eps est coherent
  // Changement: 13/12/07: en cas de valeurs negatives pour k OU eps
  // on fixe k ET eps a une valeur moyenne des 2 elements voisins

  Debog::verifier("Transport_K_Eps_base::controler_K_Eps K_Eps before", K_Eps);

  for (int n = 0; n < size; n++)
    {
      double& k   = K_Eps(n, 0);
      double& eps = K_Eps(n, 1);

      // correct big values of epsilon
      if (eps > LeEPS_MAX)
        {
          neg[2] += 1;
          eps = LeEPS_MAX;
        }

      // correct negative values of k or epsilon
      if ((k < 0 || eps < 0) )
        {
          neg[0] += (  k<0 ? 1 : 0);
          neg[1] += (eps<0 ? 1 : 0);

          // On impose une valeur plus physique (moyenne des elements voisins)
          k = 0;
          eps = 0;
          int nk = 0;
          int neps = 0;
          int nb_faces_elem = elem_faces.line_size();

          if (sub_type(Domaine_VEF, domaine_vf))
            {
              // Here we are in VEF discretization
              for (int i = 0; i < 2; i++)
                {
                  int elem = face_voisins(n, i);
                  if (elem != -1)
                    for (int j = 0; j < nb_faces_elem; j++)
                      if (j != n)
                        {
                          double k_face = K_Eps(elem_faces(elem, j), 0);
                          if (k_face > LeK_MIN)
                            {
                              k += k_face;
                              nk++;
                            }
                          double e_face = K_Eps(elem_faces(elem, j), 1);
                          if (e_face > LeEPS_MIN)
                            {
                              eps += e_face;
                              neps++;
                            }
                        }
                }
            }

          if (nk != 0) {k /= nk;}
          else {k = LeK_MIN;}
          if (neps != 0) { eps /= neps; }
          else { eps = LeEPS_MIN; }

        } // fin (k < 0 || eps < 0)
    }

  K_Eps.echange_espace_virtuel();


  Debog::verifier("Transport_K_Eps_base::controler_K_Eps K_Eps after", K_Eps);


  if (schema_temps().limpr() && !modele_turbulence().get_lquiet())
    {
      if (neg[0] > 0 || neg[1] > 0)
        {
          const double time = le_champ_K_Eps->temps();

          Journal() << "WARNING: Some values values forced for k and eps:" << finl;
          if (neg[0])
            Journal() << "Negative values of k found on   :" << neg[0] << "/" << size << " nodes at time " << time << finl;
          if (neg[1])
            Journal() << "Negative values of eps found on :" << neg[1] << "/" << size << " nodes at time " << time << finl;

          // Warning if more than 0.01% of nodes are values fixed
          double ratio_k = 100. * neg[0] / size;
          double ratio_eps = 100. * neg[1] / size;
          if ((ratio_k > 0.01 || ratio_eps > 0.01))
            {
              Cerr << "WARNING: Found high ratio of invalid values for k and/or epsilon (more that 0.01%) on process" << Process::me() << finl;
              Cerr << "Check journal log file for more information. These messages can be disabled with the flag 'quiet' in modele_turbulence." << finl;

              Journal() << "WARNING: Found high ratio of invalid values for k and/or epsilon (more that 0.01%)." << finl;
              Journal() << "It is possible your initial and/or boundary conditions on k and/or eps are wrong." << finl;
              Journal() << "Check the initial and boundary values for k and eps by using:" << finl;
              Journal() << "k~" << (dimension == 2 ? "" : "3/2*") << "(t*U)^2 (t turbulence rate, U mean velocity) ";
              Journal() << "and eps~Cmu^0.75 k^1.5/l with l turbulent length scale and Cmu a k-eps model parameter whose value is typically given as 0.09." << finl;
              Journal() << "Remark : by giving the velocity field (u) and the hydraulic diameter (d), by using boundary_field_uniform_keps_from_ud and field_uniform_keps_from_ud,  " << finl;
              Journal() << "respectively for boudnaries and initial conditions, TrioCFD will determine automatically values for k and eps." << finl;
              if (probleme().is_dilatable() == 1)
                {
                  Journal() << "Please, don't forget (sorry for this TrioCFD syntax weakness) that when using Quasi-Compressible module" << finl;
                  Journal() << "the unknowns for which you define initial and boundary conditions are rho*k and rho*eps." << finl;
                }
            }

          if (exit_on_negative_k_eps_)
            {
              // actually, this is writing the corrected field, so not very useful
              probleme().postraiter(1);
              Process::exit();
            };
        }

      if (neg[2]> 0)
        {
          const double time = le_champ_K_Eps->temps();
          Journal() << "WARNING: Some values values forced for k and eps:" << finl;
          Journal() << "Excessive values of eps found on " << neg[2] << "/" << size << " nodes at time " << time << finl;

          if (exit_on_big_eps_)
            {
              // actually, this is writing the corrected field, so not very useful
              probleme().postraiter(1);
              Process::exit();
            };
        }
    }

  return 1;
}

/*! @brief Method to correct the field K_Eps after an iteration
 *
 *  Calls Transport_K_Eps_base::controler_K_Eps() which does the work.
 *
 *  WARNING: The method controler_K_Eps is also called from Modele_turbulence_hyd_K_Eps_XXX::controler()
 *
 */
void Transport_K_Eps_base::valider_iteration()
{
  controler_K_Eps();
}
