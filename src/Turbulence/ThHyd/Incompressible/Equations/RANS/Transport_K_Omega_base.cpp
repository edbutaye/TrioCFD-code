/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Transport_K_Omega_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Domaine_VEF.h>
#include <Param.h>
#include <Debog.h>

Implemente_base(Transport_K_Omega_base, "Transport_K_Omega_base", Transport_2eq_base);

// X_D Transport_K_Omega_base Transport_2eq_base Transport_K_Omega_base 1 Base equation for RANS k-omega model. Should not be used directly

/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_Omega_base::printOn(Sortie& is) const
{ return is << que_suis_je() << "\n"; }

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_Omega_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}
void Transport_K_Omega_base::set_param(Param& param)
{
  Transport_2eq_base::set_param(param);
  param.ajouter_flag("exit_on_negative_k_omega", &exit_on_negative_k_omega_); // XD_ADD_P flag Flag to exit (with postprocessing of fields) if a negative value is found for k or omega
  param.ajouter_flag("exit_on_big_omega", &exit_on_big_omega_); // XD_ADD_P flag Flag to exit (with postprocessing of fields) if an excessively big values of omega are found
}

void Transport_K_Omega_base::discretiser()
{
  if (!sub_type(Discret_Thyd, discretisation()))
    {
      Cerr << " Transport_K_Omega_base::discretiser " << finl;
      Cerr << "Discretization " << discretisation().que_suis_je() << " not recognized." << finl;
      Process::exit();
    }

  Cerr << "K-Omega transport equation (" << que_suis_je() << ") discretization" << finl;
  Cerr << "K_Omega field discretization" << finl;
  Noms noms(2);
  Noms unit(2);
  noms[0] = "K";
  noms[1] = "omega";
  unit[0] = "m2/s2";
  unit[1] = "1/s1";

  // cAlan : possibilite de mutualiser ca dans Transport_RANS_2eq
  discretisation().discretiser_champ("temperature",  domaine_dis(), multi_scalaire,
                                     noms, unit, 2, schema_temps().nb_valeurs_temporelles(),
                                     schema_temps().temps_courant(), le_champ_K_Omega);
  le_champ_K_Omega->nommer("K_Omega");
  champs_compris_.ajoute_champ(le_champ_K_Omega);
  if (modele_turbulence().equation().calculate_time_derivative())
    set_calculate_time_derivative(1);

  Equation_base::discretiser();
}

/*! @brief Controle le champ inconnue K-Omega en forcant a zero les valeurs du champ
 *
 *     inferieurs a 1.e-10.
 *
 * @return (int) renvoie toujours 1
 */
int Transport_K_Omega_base::controler_K_Omega()
{
  DoubleTab& K_Omega = le_champ_K_Omega->valeurs();
  int size = K_Omega.dimension(0);
  if (size < 0)
    {
      if (!sub_type(Champ_Inc_P0_base, le_champ_K_Omega.valeur()))
        Process::exit("Unsupported K_Omega field in Transport_K_Omega_base::controler_K_Omega()");
      size = le_champ_K_Omega->equation().domaine_dis().domaine().nb_elem();
    }


  // neg will store the amount of problematic values of k or omega found
  // neg[0] : amount of negative k
  // neg[1] : amount of negative omega
  // neg[2] : amount of omega too big (> modele_turbulence().get_OMEGA_MAX())
  ArrOfInt neg(3);
  neg = 0;

  const int lquiet = modele_turbulence().get_lquiet(); // cAlan remonter ce lquiet dans modele_turbu

  // cAlan, le 20/01/2023 : on force les valeurs au min et max comme pour le K_Eps.
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
  const double OMEGA_MIN = modele_turbulence().get_OMEGA_MIN();
  const double OMEGA_MAX = modele_turbulence().get_OMEGA_MAX();
  const double K_MIN = modele_turbulence().get_K_MIN();
  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces = domaine_vf.elem_faces();


  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega before", K_Omega);

  for (int n = 0; n < size; n++)
    {
      double& enerK = K_Omega(n, 0);
      double& omega = K_Omega(n, 1);

      // correct big omega
      if (omega > OMEGA_MAX)
        {
          neg[2] += 1;
          omega = OMEGA_MAX;
        }

      // correct negative k or omega
      if ((enerK < 0 || omega < 0))
        {
          neg[0] += (enerK < 0 ? 1 : 0);
          neg[1] += (omega < 0 ? 1 : 0);


          // On impose une valeur plus physique (moyenne des elements voisins)
          enerK = 0;
          omega = 0;
          int nenerK = 0;
          int nomega = 0;
          const int nb_faces_elem = elem_faces.line_size();

          if (sub_type(Domaine_VEF, domaine_vf))
            {
              // cAlan : faire une fonction dans Transport_RANS_2eq qui fait la meme chose ?
              // K-Eps on faces (eg:VEF)
              for (int i = 0; i < 2; i++)
                {
                  int elem = face_voisins(n, i);
                  if (elem != -1)
                    for (int j = 0; j < nb_faces_elem; j++)
                      if (j != n)
                        {
                          double k_face = K_Omega(elem_faces(elem, j), 0);
                          if (k_face > K_MIN)
                            {
                              enerK += k_face;
                              nenerK++;
                            }
                          double o_face = K_Omega(elem_faces(elem, j), 1);
                          if (o_face > OMEGA_MIN)
                            {
                              omega += o_face;
                              nomega++;
                            }
                        }
                }

            }

          if (nenerK != 0)
            {enerK /= nenerK;}
          else
            {  enerK = K_MIN;}

          if (nomega != 0)
            { omega /= nomega;}
          else
            {omega = OMEGA_MIN;}


        }
    }
  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega middle", K_Omega);
  K_Omega.echange_espace_virtuel();
  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega after", K_Omega);


  if (schema_temps().limpr() && !modele_turbulence().get_lquiet())
    {
      if (neg[0] || neg[1])
        {


          const double time = le_champ_K_Omega->temps();
          Journal() << "Values forced for k and omega because:" << finl;
          if (neg[0])
            {
              Journal() << "Negative values found for k on "
                        << neg[0] << "/" << size << " nodes at time "
                        << time << finl;
            }
          if (neg[1])
            {
              Journal() << "Negative values found for omega on "
                        << neg[1] << "/" << size << " nodes at time "
                        << time << finl;
            }

          // Warning if more than 0.01% of nodes are values fixed
          // cAlan : mettre une variable "experte" dans le jdd pour ajuster ce seuil ?
          const double ratio_k = 100. * neg[0] / size;
          const double ratio_omega = 100. * neg[1] / size;
          if ((ratio_k > 0.01 || ratio_omega > 0.01) && !lquiet)
            {
              Cerr << "WARNING: Found high ratio of invalid values for k and/or omega (more that 0.01%) on process" << Process::me() << finl;
              Cerr << "Check journal log file for more information. These messages can be disabled with the flag 'quiet' in modele_turbulence." << finl;
              // cAlan : adapter le texte pour omega
              Journal() << "It is possible your initial and/or boundary conditions on k and/or omega are wrong." << finl;
            }

          if (exit_on_negative_k_omega_)
            {
              Cerr << "The problem is postprocessed in order to find the nodes where K or Omega values go below 0." << finl;
              probleme().postraiter(1);
              Process::exit();
            };
        }
      if (neg[2])
        {
          const double time = le_champ_K_Omega->temps();
          Journal() << "Values forced for omega because:" << finl;
          Journal() << "Maximum values found for omega on " << neg[2] << "/" << size << " nodes at time " << time << finl;

          if (exit_on_big_omega_)
            {
              Cerr << "The problem is postprocessed in order to find the nodes where Omega values too high." << finl;
              probleme().postraiter(1);
              Process::exit();
            };

        }
    }
  return 1;
}

/*! @brief on remet omega et K positifs
 *
 */
void Transport_K_Omega_base::valider_iteration()
{
  controler_K_Omega();
}
