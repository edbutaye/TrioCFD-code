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
  param.ajouter_flag("exit_on_negative_k_omega", &exit_on_negative_k_omega_); // X_D_ADD_P flag Flag to exit (with postprocessing of fields) if a negative value is found for k or omega
  param.ajouter_flag("exit_on_big_omega", &exit_on_big_omega_); // X_D_ADD_P flag Flag to exit (with postprocessing of fields) if an excessively big values of omega are found
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


  // these will store the amount of problematic values of k or omega found
  // for reporting at the end
  int count_negative_k = 0;
  int count_omega_under_threshold = 0;
  int count_omega_too_big = 0;

  const int lquiet = modele_turbulence().get_lquiet(); // cAlan remonter ce lquiet dans modele_turbu

  // cAlan, le 20/01/2023 : on force les valeurs au min et max comme pour le K_Eps.
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
  const double OMEGA_MIN = modele_turbulence().get_OMEGA_MIN();
  const double OMEGA_MAX = modele_turbulence().get_OMEGA_MAX();
  const double K_MIN = modele_turbulence().get_K_MIN();


  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega before", K_Omega);

  // Two-phase parallel algorithm following Transport_K_Eps_base pattern

  // Phase 1: Parallel detection and simple corrections
  struct CountValues
  {
    int count_negative_k;
    int count_omega_under_threshold;
    int count_omega_too_big;

    KOKKOS_INLINE_FUNCTION CountValues() : count_negative_k(0), count_omega_under_threshold(0), count_omega_too_big(0) {}
    KOKKOS_INLINE_FUNCTION CountValues(const CountValues& rhs) : count_negative_k(rhs.count_negative_k), count_omega_under_threshold(rhs.count_omega_under_threshold), count_omega_too_big(rhs.count_omega_too_big) {}
    KOKKOS_INLINE_FUNCTION CountValues& operator=(const CountValues&) = default;
    KOKKOS_INLINE_FUNCTION void operator+=(const CountValues& rhs)
    {
      count_negative_k += rhs.count_negative_k;
      count_omega_under_threshold += rhs.count_omega_under_threshold;
      count_omega_too_big += rhs.count_omega_too_big;
    }
  };

  DoubleTabView tab_K_Omega = K_Omega.view_rw();
  CountValues result;
  Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), size,
                          KOKKOS_LAMBDA(const int n, CountValues& local_count)
  {
    double& enerK = tab_K_Omega(n, 0);
    double& omega = tab_K_Omega(n, 1);

    // Count and correct big omega (simple correction, no dependencies)
    if (omega > OMEGA_MAX)
      {
        local_count.count_omega_too_big++;
        omega = OMEGA_MAX;
      }
    // Count and correct small omega (simple correction, no dependencies)
    if (omega < OMEGA_MIN)
      {
        local_count.count_omega_under_threshold++;
        omega = OMEGA_MIN;
      }

    // Count negative k (correction will be done in phase 2)
    if (enerK < 0) local_count.count_negative_k++;

  }, result);
  end_gpu_timer(__KERNEL_NAME__);

  count_negative_k = result.count_negative_k;
  count_omega_under_threshold = result.count_omega_under_threshold;
  count_omega_too_big = result.count_omega_too_big;

  // Phase 2: Parallel correction of negative k with neighbor averaging
  // Create a snapshot of original values to ensure consistent neighbor averaging
  DoubleTrav tab_K_Omega_original(K_Omega);
  tab_K_Omega_original = K_Omega;
  CDoubleTabView K_Omega_orig = tab_K_Omega_original.view_ro();
  const bool vef_algo = sub_type(Domaine_VEF, domaine_vf) && !disable_VEF_mean_value_corrections_;
  const int nb_faces_elem = domaine_vf.elem_faces().line_size();
  CIntTabView face_voisins = domaine_vf.face_voisins().view_ro();
  CIntTabView elem_faces = domaine_vf.elem_faces().view_ro();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size, KOKKOS_LAMBDA(const int n)
  {
    double& enerK = tab_K_Omega(n, 0);
    double& omega = tab_K_Omega(n, 1);

    // correct negative k
    if (K_Omega_orig(n, 0) < 0)
      {
        enerK = 0;
        omega = 0;
        int nenerK = 0;
        int nomega = 0;

        // in VEF disc, we compute the mean value of neighbours
        if (vef_algo)
          {
            for (int i = 0; i < 2; i++)
              {
                int elem = face_voisins(n, i);
                if (elem != -1)
                  for (int j = 0; j < nb_faces_elem; j++)
                    {
                      int face_j = elem_faces(elem, j);
                      if (face_j != n)
                        {
                          double k_face = K_Omega_orig(face_j, 0);
                          if (k_face > K_MIN)
                            {
                              enerK += k_face;
                              nenerK++;
                            }
                          double o_face = K_Omega_orig(face_j, 1);
                          if (o_face > OMEGA_MIN)
                            {
                              omega += o_face;
                              nomega++;
                            }
                        }
                    }
              }
          }

        if (nenerK != 0)
          enerK /= nenerK;
        else
          enerK = K_MIN;

        if (nomega != 0)
          omega /= nomega;
        else
          omega = OMEGA_MIN;
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega middle", K_Omega);
  K_Omega.echange_espace_virtuel();
  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega after", K_Omega);


  if (schema_temps().limpr() && !modele_turbulence().get_lquiet())
    {
      if (count_negative_k || count_omega_under_threshold)
        {
          const double time = le_champ_K_Omega->temps();
          Journal() << "Values forced for k and omega because:" << finl;
          if (count_negative_k)
            {
              Journal() << "Negative values found for k on "
                        << count_negative_k << "/" << size << " nodes at time "
                        << time << finl;
            }
          if (count_omega_under_threshold)
            {
              Journal() << "Negative values found for omega on "
                        << count_omega_under_threshold << "/" << size << " nodes at time "
                        << time << finl;
            }
          // Warning if more than 0.01% of nodes are values fixed
          // cAlan : mettre une variable "experte" dans le jdd pour ajuster ce seuil ?
          const double ratio_k = 100. * count_negative_k / size;
          const double ratio_omega = 100. * count_omega_under_threshold / size;
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
      if (count_omega_too_big)
        {
          const double time = le_champ_K_Omega->temps();
          Journal() << "Values forced for omega because:" << finl;
          Journal() << "Maximum values found for omega on " << count_omega_too_big << "/" << size << " nodes at time " << time << finl;

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
