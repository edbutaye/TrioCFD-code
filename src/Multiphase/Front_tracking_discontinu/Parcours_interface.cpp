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
// File:        Parcours_interface.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////
#include <Parcours_interface.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <TRUST_Deriv.h>
#include <Comm_Group.h>
#include <Connectivite_frontieres.h>
#include <Param.h>
#include <Scatter.h>
#include <Debog.h>
#include <Domaine_VDF.h> // EB
#include <Transport_Interfaces_FT_Disc.h>

//#define EXTENSION_TRIANGLE_POUR_CALCUL_INDIC_AVEC_CONSERVATION_FACETTE_COIN
Implemente_instanciable_sans_constructeur(Parcours_interface,"Parcours_interface",Objet_U);

static int flag_warning_code_missing=1;

const double Parcours_interface::Erreur_relative_maxi_ = 1.E-13;

Entree& Parcours_interface::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter_flag("correction_parcours_thomas", &correction_parcours_thomas_);
  param.lire_avec_accolades_depuis(is);
  return is;
}

Sortie& Parcours_interface::printOn(Sortie& os) const
{
  return os;
}

// ===================================================================
//  INITIALISATION DE L'OBJET
// ===================================================================

Parcours_interface::Parcours_interface()
{
  domaine_elem_ptr = 0;
  domaine_sommets_ptr = 0;
  Valeur_max_coordonnees_ = 0.;
  Erreur_max_coordonnees_ = 0.;
  correction_parcours_thomas_ = 0;
}

/*! @brief La construction par copie est interdite !
 *
 */
Parcours_interface::Parcours_interface(const Parcours_interface& a): Objet_U(a)
{
  exit();
}

/*! @brief L'operateur = est interdit !
 *
 */
const Parcours_interface& Parcours_interface::operator=(const Parcours_interface&)
{
  exit();
  return *this;
}

/*! @brief Remplissage des variables persistantes de la classe (refdomaine_vf_, nb_faces_elem_, nb_elements_reels_, type_element_,
 *
 *    equations_plans_faces_).
 *
 */
void Parcours_interface::associer_domaine_dis(const Domaine_dis& domaine_dis)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
  const Domaine& domaine = ref_cast(Domaine, domaine_dis.domaine()); // EB
  assert(! refdomaine_vf_.non_nul());
  refdomaine_vf_ = domaine_vf;
  nb_faces_elem_ = domaine_vf.domaine().nb_faces_elem();
  nb_elements_reels_ = domaine_vf.domaine().nb_elem();
  nb_sommets_par_face_ = domaine_vf.nb_som_face();
  nb_faces_reelles_ = domaine_vf.nb_faces(); // EB
  drapeaux_elements_parcourus_.resize_array(nb_elements_reels_);
  drapeaux_elements_parcourus_ = 0;
  domaine_elem_ptr = 0;
  domaine_sommets_ptr = 0;
  // debut EB
  drapeaux_faces_parcourues_x_.resize_array(nb_faces_reelles_);
  drapeaux_faces_parcourues_x_ = 0;
  drapeaux_faces_parcourues_y_.resize_array(nb_faces_reelles_);
  drapeaux_faces_parcourues_y_ = 0;
  drapeaux_faces_parcourues_z_.resize_array(nb_faces_reelles_);
  drapeaux_faces_parcourues_z_ = 0;

  nb_aretes_reelles_=domaine.aretes_som().dimension(0);
  drapeaux_aretes_parcourues_x_.resize_array(nb_aretes_reelles_);
  drapeaux_aretes_parcourues_x_ = 0;
  drapeaux_aretes_parcourues_y_.resize_array(nb_aretes_reelles_);
  drapeaux_aretes_parcourues_y_ = 0;
  drapeaux_aretes_parcourues_z_.resize_array(nb_aretes_reelles_);
  drapeaux_aretes_parcourues_z_ = 0;
  // fin EB
  // Calcul de la valeur maximale des coordonnees du maillage fixe:
  {
    const DoubleTab& coord_som = domaine_vf.domaine().coord_sommets();
    Valeur_max_coordonnees_ = mp_max_abs_vect(coord_som);
    Erreur_max_coordonnees_ = Valeur_max_coordonnees_ * Erreur_relative_maxi_;
  }

  // Remplissage du tableau equations_plans_faces_
  {
    const int nb_faces = domaine_vf.nb_faces();
    const int dimension3 = (Objet_U::dimension == 3);
    int i;
    double a, b, c = 0., d;
    double x, y = 0., z = 0.;
    double inverse_norme;
    equations_plans_faces_.resize(nb_faces, 4);
    //equations_plans_faces_face_.resize(nb_faces,2*Objet_U::dimension,4); // EB // TODO : il faut faire un tableau equations_plans_faces_face_(nb_faces_faces,4)
    for (i = 0; i < nb_faces; i++)
      {
        // Calcul de l'equation des faces des elements
        x = domaine_vf.xv(i, 0);
        if (Objet_U::bidim_axi && x == 0.)
          {
            // Cas particulier: la normale est nulle sur l'axe (surface nulle)
            // la normale pointe de gauche a droite en VDF.
            a = 1.;
            b = 0.;
          }
        else
          {
            a = domaine_vf.face_normales(i, 0);
            b = domaine_vf.face_normales(i, 1);
            y = domaine_vf.xv(i, 1);
            if (dimension3)
              {
                c = domaine_vf.face_normales(i, 2);
                z = domaine_vf.xv(i, 2);
              }
            inverse_norme = 1. / sqrt(a*a + b*b + c*c);
            a *= inverse_norme;
            b *= inverse_norme;
            c *= inverse_norme;
          }
        d = - (a * x + b * y + c * z);
        equations_plans_faces_(i, 0) = a;
        equations_plans_faces_(i, 1) = b;
        equations_plans_faces_(i, 2) = c;
        equations_plans_faces_(i, 3) = d;
      }
  }
  const Nom& nom_elem = domaine_vf.domaine().type_elem().valeur().que_suis_je();
  if (nom_elem == "Rectangle")
    type_element_ = RECTANGLE;
  else if (nom_elem == "Rectangle_2D_axi")
    type_element_ = RECTANGLE;
  else if (nom_elem == "Triangle")
    type_element_ = TRIANGLE;
  else if (nom_elem == "Tetraedre")
    type_element_ = TETRA;
  else if (nom_elem == "Hexaedre")
    type_element_ = HEXA;
  else
    {
      Cerr << "Parcours_interface::associer_domaine_vf\n";
      Cerr << " Le type d'element " << nom_elem;
      Cerr << " n'est pas reconnu !\n";
      exit();
    }
}

void Parcours_interface::associer_connectivite_frontieres(const Connectivite_frontieres& connect)
{
  refconnect_front_ = connect;
}

// ===================================================================
//   METHODE PUBLIQUE : parcours de l'interface
// ===================================================================

// Boucle sur les facettes :
//  Remplissage de intersections_elem_facettes_
//  Completion du maillage (facettes sur les "processeurs pauvres")
//   (ajout de noeuds virtuels, de facettes virtuelles, et completion
//    des espaces distants et virtuels)

void Parcours_interface::parcourir(Maillage_FT_Disc& maillage) const
{
  assert(refdomaine_vf_.non_nul());

  const Domaine_VF& domaine_vf = refdomaine_vf_.valeur();
  domaine_elem_ptr = & domaine_vf.domaine().les_elems();
  domaine_sommets_ptr = & domaine_vf.domaine().les_sommets();

  const Domaine_dis& domaine_dis = maillage.refdomaine_dis_.valeur();

  static int iteration_parcourir_faces=0;
  static int iteration_parcourir_aretes=0;
  int calcul_precis_indic_face=0,calcul_precis_indic_arete=0;
  if (maillage.refequation_transport_.non_nul())
    {
      calcul_precis_indic_face=maillage.equation_transport().calcul_precis_indic_faces();
      calcul_precis_indic_arete=maillage.equation_transport().calcul_precis_indic_aretes();
    }

  const double temps=maillage.temps();
  DoubleTab copie_sommets_maillage;
  if (correction_parcours_thomas_)
    {
      // Modification pour eliminer les erreurs de calcul de l'indicatrice
      //  Le calcul produit des resultats aleatoires si les sommets du maillage lagrangien
      //  se trouvent exactement sur une face d'un element eulerien (on risque de compter
      //  deux fois des contributions de volume). On deplace donc legerement les sommets
      //  pour les eloigner des faces.
      //
      // REMARQUE IMPORTANTE : a la fin de la procedure, on remettre
      // les sommets du maillage lagrangien dans leur etat initial
      copie_sommets_maillage = maillage.sommets();
      eloigner_sommets_des_faces(maillage);
    }



  // RAZ des intersections elements-facettes
  {
    const int nb_elem = domaine_vf.nb_elem();
    const int nb_facettes = maillage.facettes_.dimension(0);
    maillage.intersections_elem_facettes_.reset(nb_elem, nb_facettes);
  }

  // debut EB
  // RAZ des intersections face-facettes
  {
    const int nb_faces = domaine_vf.nb_faces();
    const int nb_facettes = maillage.facettes_.dimension(0);
    maillage.intersections_face_facettes_x_.reset(nb_faces, nb_facettes);
    maillage.intersections_face_facettes_y_.reset(nb_faces, nb_facettes);
    maillage.intersections_face_facettes_z_.reset(nb_faces, nb_facettes);
  }

  // RAZ des intersections arete-facettes
  {
    const int nb_aretes = domaine_vf.domaine().nb_aretes();
    const int nb_facettes = maillage.facettes_.dimension(0);
    maillage.intersections_arete_facettes_x_.reset(nb_aretes, nb_facettes);
    maillage.intersections_arete_facettes_y_.reset(nb_aretes, nb_facettes);
    maillage.intersections_arete_facettes_z_.reset(nb_aretes, nb_facettes);
  }
  // fin EB

  // Facettes et elements d'arrivee a envoyer aux processeurs voisins pour l'iteration
  // suivante:
  // ELEMENTS
  static ArrOfIntFT echange_facettes_numfacette;
  static ArrOfIntFT echange_facettes_numelement;
  // Facettes et elements d'arrivee de la liste de facette en cours de traitement
  // a l'iteration courante:
  static ArrOfIntFT facettes_a_traiter_numfacette;
  static ArrOfIntFT facettes_a_traiter_numelement;

  compteur_erreur_grossiere = 0;
  //int iteration = 0;
  // debut EB
  // FACES
  // Face d'arrivees a envoyer aux processeurs voisins pour l'iteration
  // suivante:
  static ArrOfIntFT echange_facettes_numfacette_x;
  static ArrOfIntFT echange_facettes_numfacette_y;
  static ArrOfIntFT echange_facettes_numfacette_z;

  static ArrOfIntFT echange_facettes_numface_x;
  static ArrOfIntFT echange_facettes_numface_y;
  static ArrOfIntFT echange_facettes_numface_z;
  // Faces d'arrivee de la liste de facette en cours de traitement
  // a l'iteration courante:
  static ArrOfIntFT facettes_a_traiter_numfacette_x;
  static ArrOfIntFT facettes_a_traiter_numfacette_y;
  static ArrOfIntFT facettes_a_traiter_numfacette_z;

  static ArrOfIntFT facettes_a_traiter_numface_x;
  static ArrOfIntFT facettes_a_traiter_numface_y;
  static ArrOfIntFT facettes_a_traiter_numface_z;

  // ARETES

  static ArrOfIntFT echange_facettes_numfacette_arete_x;
  static ArrOfIntFT echange_facettes_numfacette_arete_y;
  static ArrOfIntFT echange_facettes_numfacette_arete_z;

  static ArrOfIntFT echange_facettes_numarete_x;
  static ArrOfIntFT echange_facettes_numarete_y;
  static ArrOfIntFT echange_facettes_numarete_z;
  // Faces d'arrivee de la liste de facette en cours de traitement
  // a l'iteration courante:
  static ArrOfIntFT facettes_a_traiter_numfacette_arete_x;
  static ArrOfIntFT facettes_a_traiter_numfacette_arete_y;
  static ArrOfIntFT facettes_a_traiter_numfacette_arete_z;

  static ArrOfIntFT facettes_a_traiter_numarete_x;
  static ArrOfIntFT facettes_a_traiter_numarete_y;
  static ArrOfIntFT facettes_a_traiter_numarete_z;

  // fin EB


  int nb_facettes_echangees_elem = 0;

  // FACES
  int nb_facettes_echangees_x = 0;
  int nb_facettes_echangees_y = 0;
  int nb_facettes_echangees_z = 0;

  // ARETES
  int nb_facettes_echangees_arete_x = 0;
  int nb_facettes_echangees_arete_y = 0;
  int nb_facettes_echangees_arete_z = 0;

  int iteration_elem = 0, iteration_x=0, iteration_y=0, iteration_z=0;
  int iteration_arete_x=0, iteration_arete_y=0, iteration_arete_z=0;
  // On est oblige de faire une boucle pour les elements et une boucle par orientation
  // sinon on converge jamais en // a cause de la boucle while et de l'echange des numfacettes/numelem/numface
  // On commence par les elements

  do
    {
      echange_facettes_numfacette.resize_array(0);
      echange_facettes_numelement.resize_array(0);

      // Au premier passage, la liste des facettes a traiter est vide,
      // on traite toutes les facettes du domaine
      if (iteration_elem == 0)
        {
          int nb_facettes = maillage.facettes_.dimension(0);
          for (int i = 0; i < nb_facettes; i++)
            {
              // Le parcours commence par l'element qui contient le premier
              // sommet, si la facette m'appartient (element_depart >= 0)
              int premier_sommet = maillage.facettes_(i,0);
              int element_depart = maillage.sommet_elem_[premier_sommet];
              if (element_depart >= 0)
                parcours_facette(domaine_vf, maillage,
                                 echange_facettes_numfacette,
                                 echange_facettes_numelement,
                                 i, element_depart);
            }
        }
      else
        {
          // Aux passages suivants, on traite uniquement la liste
          int nb_facettes = facettes_a_traiter_numfacette.size_array();
          for (int i = 0; i < nb_facettes; i++)
            {
              int num_facette = facettes_a_traiter_numfacette[i];
              int element_depart = facettes_a_traiter_numelement[i];
              parcours_facette(domaine_vf, maillage,
                               echange_facettes_numfacette,
                               echange_facettes_numelement,
                               num_facette, element_depart);
            }
        }

      // Echange des facettes entre processeurs, on recupere
      // la liste des facettes a traiter et l'element d'arrivee.
      maillage.echanger_facettes(echange_facettes_numfacette,
                                 echange_facettes_numelement,
                                 facettes_a_traiter_numfacette,
                                 facettes_a_traiter_numelement);
      // Arret lorsque plus aucun processeur n'a de facettes a traiter.
      nb_facettes_echangees_elem =
        Process::mp_sum(facettes_a_traiter_numfacette.size_array());
      Process::Journal() << " nombre de faces echangees : ";
      Process::Journal() << facettes_a_traiter_numfacette.size_array() << " (total ";
      Process::Journal() << nb_facettes_echangees_elem << ")" << finl;

      iteration_elem++;

    }
  while (nb_facettes_echangees_elem > 0);

  if (calcul_precis_indic_face)
    {
      maillage.update_sommet_face();
      //maillage.update_sommet_arete();
      if (iteration_parcourir_faces>1 || temps>0) // lors des premieres iterations, bug parce qu'on a pas encore rempli sommet_face_ (fait uniquement apres
        // la premiere maj de Transport_Interfaces_FT_Disc
        {
          // face de normale x
          do
            {
              echange_facettes_numface_x.resize_array(0);
              echange_facettes_numfacette_x.resize_array(0);

              // EB
              // Au premier passage, la liste des facettes a traiter est vide,
              // on traite toutes les facettes du domaine
              if (iteration_x == 0)
                {

                  int nb_facettes = maillage.facettes_.dimension(0);
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      // On definit face_depart_x, face_depart_y et face_depart_z car un sommet est present dans le volume
                      // de controle de 3 faces, il faut donc appeler parcours_facette_face 3 fois
                      int premier_sommet = maillage.facettes_(i,0);
                      int face_depart_x = maillage.sommet_face_(premier_sommet,0);
                      if (face_depart_x >= 0)
                        {
                          parcours_facette_face_x(domaine_vf, maillage,
                                                  echange_facettes_numfacette_x,
                                                  echange_facettes_numface_x,
                                                  i, face_depart_x, maillage.intersections_face_facettes_x_, drapeaux_faces_parcourues_x_);
                        }
                    }
                }
              else
                {
                  // On definit face_depart_x, face_depart_y et face_depart_z car un sommet est present dans le volume
                  // de controle de 3 faces, il faut donc appeler parcours_facette_face 3 fois
                  int nb_facettes = facettes_a_traiter_numfacette_x.size_array();

                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int num_facette = facettes_a_traiter_numfacette_x[i];
                      int face_depart_x = facettes_a_traiter_numface_x[i];
                      if (face_depart_x >= 0 )
                        {
                          parcours_facette_face_x(domaine_vf, maillage,
                                                  echange_facettes_numfacette_x,
                                                  echange_facettes_numface_x,
                                                  num_facette, face_depart_x, maillage.intersections_face_facettes_x_, drapeaux_faces_parcourues_x_);
                        }
                    }
                }
              // Echange des facettes entre processeurs, on recupere
              // la liste des facettes a traiter et l'element d'arrivee.
              maillage.echanger_facettes_face_x(echange_facettes_numfacette_x,
                                                echange_facettes_numface_x,
                                                facettes_a_traiter_numfacette_x,
                                                facettes_a_traiter_numface_x);
              // maillage.update_sommet_face();
              // Arret lorsque plus aucun processeur n'a de facettes a traiter.
              nb_facettes_echangees_x =
                Process::mp_sum(facettes_a_traiter_numfacette_x.size_array());
              Process::Journal() << " nombre de faces echangees_x : ";
              Process::Journal() << facettes_a_traiter_numfacette_x.size_array() << " (total ";
              Process::Journal() << nb_facettes_echangees_x << ")" << finl;
              iteration_x++;
            }
          while (nb_facettes_echangees_x > 0);

          // faces de normale y
          do
            {
              echange_facettes_numface_y.resize_array(0);
              echange_facettes_numfacette_y.resize_array(0);

              if (iteration_y == 0)
                {
                  int nb_facettes = maillage.facettes_.dimension(0);
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int premier_sommet = maillage.facettes_(i,0);
                      int face_depart_y = maillage.sommet_face_(premier_sommet,1);

                      if (face_depart_y >= 0)
                        parcours_facette_face_y(domaine_vf, maillage,
                                                echange_facettes_numfacette_y,
                                                echange_facettes_numface_y,
                                                i, face_depart_y, maillage.intersections_face_facettes_y_, drapeaux_faces_parcourues_y_);
                    }

                }
              else
                {
                  // Aux passages suivants, on traite uniquement la liste
                  int nb_facettes = facettes_a_traiter_numfacette_y.size_array();
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int num_facette = facettes_a_traiter_numfacette_y[i];
                      int face_depart_y = facettes_a_traiter_numface_y[i];
                      if (face_depart_y >= 0 )
                        parcours_facette_face_y(domaine_vf, maillage,
                                                echange_facettes_numfacette_y,
                                                echange_facettes_numface_y,
                                                num_facette, face_depart_y, maillage.intersections_face_facettes_y_, drapeaux_faces_parcourues_y_);
                    }
                }
              maillage.echanger_facettes_face_y(echange_facettes_numfacette_y,
                                                echange_facettes_numface_y,
                                                facettes_a_traiter_numfacette_y,
                                                facettes_a_traiter_numface_y);
              nb_facettes_echangees_y =
                Process::mp_sum(facettes_a_traiter_numfacette_y.size_array());

              Process::Journal() << " nombre de faces echangees_y : ";
              Process::Journal() << facettes_a_traiter_numfacette_y.size_array() << " (total ";
              Process::Journal() << nb_facettes_echangees_y << ")" << finl;

              iteration_y++;
            }
          while (nb_facettes_echangees_y > 0);

          // faces de normale z
          do
            {
              echange_facettes_numface_z.resize_array(0);
              echange_facettes_numfacette_z.resize_array(0);

              if (iteration_z == 0)
                {
                  int nb_facettes = maillage.facettes_.dimension(0);
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int premier_sommet = maillage.facettes_(i,0);
                      int face_depart_z = maillage.sommet_face_(premier_sommet,2);
                      if (face_depart_z >= 0 )
                        parcours_facette_face_z(domaine_vf, maillage,
                                                echange_facettes_numfacette_z,
                                                echange_facettes_numface_z,
                                                i, face_depart_z, maillage.intersections_face_facettes_z_,drapeaux_faces_parcourues_z_);
                    }
                }
              else
                {
                  // Aux passages suivants, on traite uniquement la liste
                  int nb_facettes = facettes_a_traiter_numfacette_z.size_array();
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int num_facette = facettes_a_traiter_numfacette_z[i];
                      int face_depart_z = facettes_a_traiter_numface_z[i];
                      if (face_depart_z >= 0)
                        parcours_facette_face_z(domaine_vf, maillage,
                                                echange_facettes_numfacette_z,
                                                echange_facettes_numface_z,
                                                num_facette, face_depart_z, maillage.intersections_face_facettes_z_,drapeaux_faces_parcourues_z_);
                    }
                }

              maillage.echanger_facettes_face_z(echange_facettes_numfacette_z,
                                                echange_facettes_numface_z,
                                                facettes_a_traiter_numfacette_z,
                                                facettes_a_traiter_numface_z);
              nb_facettes_echangees_z =
                Process::mp_sum(facettes_a_traiter_numfacette_z.size_array());
              Process::Journal() << " nombre de faces echangees_z : ";
              Process::Journal() << facettes_a_traiter_numfacette_z.size_array() << " (total ";
              Process::Journal() << nb_facettes_echangees_z << ")" << finl;

              iteration_z++;
            }
          while (nb_facettes_echangees_z > 0);

        }

      iteration_parcourir_faces++;
    }

  /*
  Cerr << "AVANT PARCOURS ARETES\n{" << finl;
  Cerr << "\tmaillage.facettes_.dimension(0) "  << maillage.facettes_.dimension(0) << finl;
  Cerr << "\tmaillage.sommet_elem_.size_array() " << maillage.sommet_elem_.size_array() << finl;
  Cerr << "\tmaillage.sommet_face_.dimension(0) " << maillage.sommet_face_.dimension(0) << "\n}" << finl;
  */

  //Cerr << "avant calcul_precis_indic_arete" << finl;

  if (calcul_precis_indic_arete)
    {
      maillage.update_sommet_arete();
      const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF,domaine_dis.valeur());
      const IntVect& type_arete=domaine_vdf.type_arete();

      if (iteration_parcourir_aretes>1|| temps>0) // lors des premieres iterations, bug parce qu'on a pas encore rempli sommet_face_ (fait uniquement apres
        // la premiere maj de Transport_Interfaces_FT_Disc
        {
          // face de normale x
          do
            {
              echange_facettes_numarete_x.resize_array(0);
              echange_facettes_numfacette_arete_x.resize_array(0);

              // EB
              // Au premier passage, la liste des facettes a traiter est vide,
              // on traite toutes les facettes du domaine
              if (iteration_arete_x == 0)
                {

                  int nb_facettes = maillage.facettes_.dimension(0);
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      // On definit face_depart_x, face_depart_y et face_depart_z car un sommet est present dans le volume
                      // de controle de 3 faces, il faut donc appeler parcours_facette_face 3 fois
                      int premier_sommet = maillage.facettes_(i,0);
                      int arete_depart_x = maillage.sommet_arete_(premier_sommet,0);
                      if (type_arete(arete_depart_x)!=2) continue; // on ne calcule que pour les aretes_internes
                      if (arete_depart_x >= 0)
                        {
                          parcours_facette_arete_x(domaine_vf, maillage,
                                                   echange_facettes_numfacette_arete_x,
                                                   echange_facettes_numarete_x,
                                                   i, arete_depart_x, maillage.intersections_arete_facettes_x_, drapeaux_aretes_parcourues_x_);
                        }
                    }
                }
              else
                {
                  // On definit face_depart_x, face_depart_y et face_depart_z car un sommet est present dans le volume
                  // de controle de 3 faces, il faut donc appeler parcours_facette_face 3 fois
                  int nb_facettes = facettes_a_traiter_numfacette_arete_x.size_array();

                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int num_facette = facettes_a_traiter_numfacette_arete_x[i];
                      int arete_depart_x = facettes_a_traiter_numarete_x[i];
                      if (type_arete(arete_depart_x)!=2) continue; // on ne calcule que pour les aretes_internes
                      if (arete_depart_x >= 0 )
                        {
                          parcours_facette_arete_x(domaine_vf, maillage,
                                                   echange_facettes_numfacette_arete_x,
                                                   echange_facettes_numarete_x,
                                                   num_facette, arete_depart_x, maillage.intersections_arete_facettes_x_, drapeaux_aretes_parcourues_x_);
                        }
                    }
                }
              // Echange des facettes entre processeurs, on recupere
              // la liste des facettes a traiter et l'element d'arrivee.
              /*
              Cerr << "avant echanger_facettes_arete_x" << finl;
              Cerr << "echange_facettes_numarete_x.size_array() " << echange_facettes_numarete_x.size_array() << finl;
              for (int i=0; i<echange_facettes_numarete_x.size_array(); i++)
                {
                  int arete=echange_facettes_numarete_x(i);
                  Cerr << "arete " << arete;
                  Cerr << " cg " << domaine_vf.xa(arete,0) << " " << domaine_vf.xa(arete,1) << " " << domaine_vf.xa(arete,2) << finl;
                }
                */
              maillage.echanger_facettes_arete_x(echange_facettes_numfacette_arete_x,
                                                 echange_facettes_numarete_x,
                                                 facettes_a_traiter_numfacette_arete_x,
                                                 facettes_a_traiter_numarete_x);
              // Cerr << "apres echanger_facettes_arete_x" << finl;
              // Arret lorsque plus aucun processeur n'a de facettes a traiter.
              nb_facettes_echangees_arete_x =
                Process::mp_sum(facettes_a_traiter_numfacette_arete_x.size_array());

              Process::Journal() << " nombre de d'aretes echangees x : ";
              Process::Journal() << facettes_a_traiter_numfacette_arete_x.size_array() << " (total ";
              Process::Journal() << nb_facettes_echangees_arete_x << ")" << finl;

              iteration_arete_x++;
            }
          while (nb_facettes_echangees_arete_x > 0);

          //Cerr << "fin parcours arete x" << finl;
          /*
          Cerr << "APRES PARCOURS ARETES x \n{" << finl;
          Cerr << "\tmaillage.facettes_.dimension(0) "  << maillage.facettes_.dimension(0) << finl;
          Cerr << "\tmaillage.sommet_elem_.size_array() " << maillage.sommet_elem_.size_array() << finl;
          Cerr << "\tmaillage.sommet_face_.dimension(0) " << maillage.sommet_face_.dimension(0) << "\n}" << finl;
          Cerr << "debut parcours y "<< finl;
          */

          // faces de normale y
          do
            {
              echange_facettes_numarete_y.resize_array(0);
              echange_facettes_numfacette_arete_y.resize_array(0);

              if (iteration_arete_y == 0)
                {
                  int nb_facettes = maillage.facettes_.dimension(0);
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int premier_sommet = maillage.facettes_(i,0);
                      int arete_depart_y = maillage.sommet_arete_(premier_sommet,1);
                      if (type_arete(arete_depart_y)!=2) continue; // on ne calcule que pour les aretes_internes
                      if (arete_depart_y >= 0)
                        parcours_facette_arete_y(domaine_vf, maillage,
                                                 echange_facettes_numfacette_arete_y,
                                                 echange_facettes_numarete_y,
                                                 i, arete_depart_y, maillage.intersections_arete_facettes_y_, drapeaux_aretes_parcourues_y_);
                    }
                }
              else
                {
                  // Aux passages suivants, on traite uniquement la liste
                  int nb_facettes = facettes_a_traiter_numfacette_arete_y.size_array();
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int num_facette = facettes_a_traiter_numfacette_arete_y[i];
                      int arete_depart_y = facettes_a_traiter_numarete_y[i];
                      if (type_arete(arete_depart_y)!=2) continue; // on ne calcule que pour les aretes_internes
                      if (arete_depart_y >= 0 )
                        parcours_facette_arete_y(domaine_vf, maillage,
                                                 echange_facettes_numfacette_arete_y,
                                                 echange_facettes_numarete_y,
                                                 num_facette, arete_depart_y, maillage.intersections_arete_facettes_y_, drapeaux_aretes_parcourues_y_);
                    }
                }
              /*
              Cerr << "echange_facettes_numarete_y.size_array() " << echange_facettes_numarete_y.size_array() << finl;
              for (int arete=0; arete<echange_facettes_numarete_y.size_array(); arete++)
                {
                  int ind_arete=echange_facettes_numarete_y(arete);
                  Cerr << "arete " << ind_arete << "/" << domaine_vf.domaine().aretes_som().dimension(0) << "\t" << domaine_vf.xa(ind_arete,0) << " " << domaine_vf.xa(ind_arete,1) << " " << domaine_vf.xa(ind_arete,2) << finl;
                }
              */
              maillage.echanger_facettes_arete_y(echange_facettes_numfacette_arete_y,
                                                 echange_facettes_numarete_y,
                                                 facettes_a_traiter_numfacette_arete_y,
                                                 facettes_a_traiter_numarete_y);
              /*
                            Cerr << "APRES ECHANGE"<< finl;
                            for (int arete=0; arete<facettes_a_traiter_numarete_y.size_array(); arete++)
                              {
                                int ind_arete=facettes_a_traiter_numarete_y(arete);
                                Cerr << "arete " << ind_arete  << "/" << domaine_vf.domaine().aretes_som().dimension(0) << "\t" << domaine_vf.xa(ind_arete,0) << " " << domaine_vf.xa(ind_arete,1) << " " << domaine_vf.xa(ind_arete,2) << finl;
                              }

                            Cerr <<"echange_facettes_numfacette_arete_y\n " << echange_facettes_numfacette_arete_z << finl;
                            Cerr <<"facettes_a_traiter_numfacette_arete_y\n " << facettes_a_traiter_numfacette_arete_z << finl;
              */
              nb_facettes_echangees_arete_y =
                Process::mp_sum(facettes_a_traiter_numfacette_arete_y.size_array());
              //Cerr << "iteration_arete_x " << iteration_arete_x << finl;

              Process::Journal() << " nombre de d'aretes echangees y : ";
              Process::Journal() << facettes_a_traiter_numfacette_arete_y.size_array() << " (total ";
              Process::Journal() << nb_facettes_echangees_arete_y << ")" << finl;

              iteration_arete_y++;
            }
          while (nb_facettes_echangees_arete_y > 0);

          /*
          Cerr << "APRES PARCOURS ARETES y \n{" << finl;
          Cerr << "\tmaillage.facettes_.dimension(0) "  << maillage.facettes_.dimension(0) << finl;
          Cerr << "\tmaillage.sommet_elem_.size_array() " << maillage.sommet_elem_.size_array() << finl;
          Cerr << "\tmaillage.sommet_face_.dimension(0) " << maillage.sommet_face_.dimension(0) << "\n}" << finl;
          */

          // faces de normale z
          do
            {
              echange_facettes_numarete_z.resize_array(0);
              echange_facettes_numfacette_arete_z.resize_array(0);

              if (iteration_arete_z == 0)
                {
                  int nb_facettes = maillage.facettes_.dimension(0);
                  //Cerr << "ici" << finl;
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      //Cerr << "i = " << i  << "/" <<nb_facettes << finl;
                      int premier_sommet = maillage.facettes_(i,0);
                      //Cerr << "premier_sommet " << premier_sommet << finl;
                      int arete_depart_z = maillage.sommet_arete_(premier_sommet,2);
                      //Cerr << "arete_depart_z " << arete_depart_z << finl;
                      int type_arete_depart=type_arete(arete_depart_z);
                      //Cerr << "type_arete_depart " << type_arete_depart << finl;
                      if (type_arete_depart!=2) continue; // on ne calcule que pour les aretes_internes
                      //Cerr << "i = " << i << " arete_depart " << arete_depart_z << finl;
                      if (arete_depart_z >= 0 )
                        parcours_facette_arete_z(domaine_vf, maillage,
                                                 echange_facettes_numfacette_arete_z,
                                                 echange_facettes_numarete_z,
                                                 i, arete_depart_z, maillage.intersections_arete_facettes_z_,drapeaux_aretes_parcourues_z_);
                    }
                  //Cerr << "la" << finl;
                }
              else
                {

                  // Aux passages suivants, on traite uniquement la liste
                  int nb_facettes = facettes_a_traiter_numfacette_arete_z.size_array();
                  for (int i = 0; i < nb_facettes; i++)
                    {
                      int num_facette = facettes_a_traiter_numfacette_arete_z[i];
                      int arete_depart_z = facettes_a_traiter_numarete_z[i];
                      if (type_arete(arete_depart_z)!=2) continue; // on ne calcule que pour les aretes_internes
                      if (arete_depart_z >= 0)
                        parcours_facette_arete_z(domaine_vf, maillage,
                                                 echange_facettes_numfacette_arete_z,
                                                 echange_facettes_numarete_z,
                                                 num_facette, arete_depart_z, maillage.intersections_arete_facettes_z_,drapeaux_aretes_parcourues_z_);
                    }

                }
              //Cerr << "avant echanger fa7_arete z" << finl;
              maillage.echanger_facettes_arete_z(echange_facettes_numfacette_arete_z,
                                                 echange_facettes_numarete_z,
                                                 facettes_a_traiter_numfacette_arete_z,
                                                 facettes_a_traiter_numarete_z);
              //Cerr << "apres echanger fa7_arete z" << finl;
              nb_facettes_echangees_arete_z =
                Process::mp_sum(facettes_a_traiter_numfacette_arete_z.size_array());
              //Cerr << "nb_facettes_echangees_arete_z " << nb_facettes_echangees_arete_z << finl;
              Process::Journal() << " nombre d'aretes echangees z : ";
              Process::Journal() << facettes_a_traiter_numfacette_arete_z.size_array() << " (total ";
              Process::Journal() << nb_facettes_echangees_arete_z << ")" << finl;

              iteration_arete_z++;
            }
          while (nb_facettes_echangees_arete_z > 0);

          /*
          Cerr << "APRES PARCOURS ARETES z \n{" << finl;
          Cerr << "\tmaillage.facettes_.dimension(0) "  << maillage.facettes_.dimension(0) << finl;
          Cerr << "\tmaillage.sommet_elem_.size_array() " << maillage.sommet_elem_.size_array() << finl;
          Cerr << "\tmaillage.sommet_face_.dimension(0) " << maillage.sommet_face_.dimension(0) << "\n}" << finl;
          */
        }
      iteration_parcourir_aretes++;
    }
  //Cerr << "apres calcul_precis_indic_arete" << finl;


  if (correction_parcours_thomas_)
    {
      //Modification pour eliminer les erreurs de calcul de l'indicatrice
      //REMARQUE IMPORTANTE : on remet le maillage dans son etat initial
      //Ceci est possible car la classe Parcours_interface est une classe
      //amie de Maillage_FT_Disc
      // Il y a de nouveaux sommets virtuels dans le maillage (mais pas de changement de proprietaire)
      // il faut mettre a jour le tableau.
      const int ni = maillage.sommets_.dimension(0);
      const int nj = maillage.sommets_.dimension(1);
      copie_sommets_maillage.resize(ni, nj);
      maillage.desc_sommets().echange_espace_virtuel(copie_sommets_maillage);
      maillage.sommets_ = copie_sommets_maillage;
    }

  Process::Journal() << " comptage erreurs_grossieres=";
  Process::Journal() << compteur_erreur_grossiere << finl;
  Process::Journal() << "FIN Parcours_interface::parcourir"<<finl;

}

// ==========================================================================
// ==========================================================================
// ==========================================================================

static void effacer_drapeaux_elements_parcourus(
  ArrOfBit& drapeaux_elements_parcourus,
  const ArrOfInt& liste_elements_parcourus)
{
  // Remise a zero des drapeaux pour la prochaine fois
  // (efficace : on ne reinitialise que les drapeaux qui sont mis)
  int n = liste_elements_parcourus.size_array();
  for (int i = 0; i < n; i++)
    {
      const int element = liste_elements_parcourus[i];
      drapeaux_elements_parcourus.clearbit(element);
    }
}
// debut EB
static void effacer_drapeaux_faces_parcourues(
  ArrOfBit& drapeaux_faces_parcourues,
  const ArrOfInt& liste_faces_parcourues)
{
  // Remise a zero des drapeaux pour la prochaine fois
  // (efficace : on ne reinitialise que les drapeaux qui sont mis)
  int n = liste_faces_parcourues.size_array();
  for (int i = 0; i < n; i++)
    {
      const int face = liste_faces_parcourues[i];
      drapeaux_faces_parcourues.clearbit(face);
    }
}
// fin EB

// Calcul de l'equation de plan d'une face d'un element eulerien. Pour un point (x,y,z),
// la fonction (a*x+b*y+c*z+d) * signe est positive a l'interieur de l'element,
// negative a l'exterieur. La normale (a,b,c) est unitaire.
// Parametre:     domaine_vf
// Signification: reference au Domaine_VF. Attention, ce doit etre la meme domaine que celle qui
//                a ete associee au parcours !!!
// Parametre:     num_element
// Signification: numero d'un element reel du domaine (0 <= num_element < elem_faces().dimension(0))
// Parametre:     num_face_element
// Signification: numero d'une face de l'element (0 <= num_face_element < elem_faces().dimension(1))
// Parametre:     a,b,c,d
// Signification: on y met les coefficients du plan: normale (a,b,c) et constante d
// Valeur de retour : signe (1. ou -1.)
inline double Parcours_interface::calcul_eq_plan(const Domaine_VF& domaine_vf,
                                                 const int num_element,
                                                 const int num_face_element,
                                                 double& a, double& b, double& c, double& d) const
{
  assert((&domaine_vf) == (&(refdomaine_vf_.valeur())));
  // Numero de la face i
  const int num_face = domaine_vf.elem_faces(num_element, num_face_element);
  // Numero du premier element voisin de la face
  const int elem_voisin = domaine_vf.face_voisins(num_face, 0);
  // Equation du plan contenant la face
  a = equations_plans_faces_(num_face, 0);
  b = equations_plans_faces_(num_face, 1);
  c = equations_plans_faces_(num_face, 2);
  d = equations_plans_faces_(num_face, 3);

  if (elem_voisin == num_element)
    return -1.;
  else
    return 1.;
}

// debut EB
// beaucoup de doublons dans ce tableau... Il faudrait implementer un connectivite arete/num_face_arete au meme titre
// que la connectivite element/face
inline double Parcours_interface::calcul_eq_aretes(const Domaine_VDF& domaine_vdf,
                                                   const int num_arete, const int num_face_arete,
                                                   double& a, double& b, double& c, double& d) const
{
  // Equation du plan contenant la face
  const int premiere_arete_interne=domaine_vdf.premiere_arete_interne();
  const int arete=num_arete-premiere_arete_interne;

  a = equations_plans_faces_arete_(arete, num_face_arete, 0);
  b = equations_plans_faces_arete_(arete, num_face_arete, 1);
  c = equations_plans_faces_arete_(arete, num_face_arete,  2);
  d = equations_plans_faces_arete_(arete, num_face_arete, 3);

  return ((num_face_arete<dimension) ? 1 : -1);
}
// fin EB

// debut EB
// Calcul de l'equation de plan d'une face d'une face eulerienne. Pour un point (x,y,z),
// la fonction (a*x+b*y+c*z+d) * signe est positive a l'interieur de la face,
// negative a l'exterieur. La normale (a,b,c) est unitaire.
// Parametre:     domaine_VDF
// Signification: reference a la Domaine_VF. Attention, ce doit etre le meme domaine que celle qui
//                a ete associee au parcours !!!
// Parametre:     num_element
// Signification: numero d'une face reelle de la zone (0 <= num_face <face_voisins().dimension(0))
// Parametre:     num_face_face
// Signification: numero d'une face de la face (0 <= num_face_face < elem_faces().dimension(1))
// Parametre:     a,b,c,d
// Signification: on y met les coefficients du plan: normale (a,b,c) et constante d
// Valeur de retour : signe (1. ou -1.)
inline double Parcours_interface::calcul_eq_plan_face(const Domaine_VF& domaine_vf,
                                                      const int num_face,
                                                      const int num_face_face,
                                                      double& a, double& b, double& c, double& d) const
{
  a = (num_face_face%dimension==0);
  b = (num_face_face%dimension==1);
  c = (num_face_face%dimension==2);

  const int face_double=domaine_vf.faces_doubles()[num_face];

  int face_bord0 = (domaine_vf.face_voisins(num_face,0)<0 || domaine_vf.face_voisins(num_face,0)>=domaine_vf.nb_elem()) ; // pas d'acces a l'element de gauche, bas, arriere ou face_double
  int face_bord1 = (domaine_vf.face_voisins(num_face,1)<0 || domaine_vf.face_voisins(num_face,1)>=domaine_vf.nb_elem()); // pas d'acces a l'element de droite, haut, avant ou face double
  int direction=(num_face_face>=dimension);
  //const int nb_elem_reel=domaine_vf.nb_elem();

  if (!face_bord0 && !face_bord1 &&  !face_double)
    {
      int elem = domaine_vf.face_voisins(num_face,direction);
      if (num_face_face%dimension==domaine_vf.orientation(num_face)) d = - domaine_vf.xp(elem,num_face_face%dimension);
      else
        {
          int face_eulerienne=domaine_vf.elem_faces(elem, num_face_face);
          d = - domaine_vf.xv(face_eulerienne,num_face_face%dimension);
        }
    }
  else // on est sur le volume de controle d'une face de bord ou d'une face_double
    {
      int elem_eulerien=-1;
      const int elem0=domaine_vf.face_voisins(num_face,0);
      const int elem1=domaine_vf.face_voisins(num_face,1);
      int bord=-1;
      elem_eulerien = face_bord0 ? elem1 : elem0;
      bord = face_bord0 ? 0 : 1;
      //}

      if (num_face_face%dimension==domaine_vf.orientation(num_face))
        {
          if (bord==0 && (num_face_face<dimension)) d=-domaine_vf.xv(num_face,num_face_face%dimension);
          else if (bord==0 && (num_face_face>=dimension)) d = -domaine_vf.xp(elem_eulerien,num_face_face%dimension);
          else if (bord==1 && (num_face_face<dimension))  d = -domaine_vf.xp(elem_eulerien,num_face_face%dimension);
          else d=-domaine_vf.xv(num_face,num_face_face%dimension);
        }
      else
        {
          int face_eulerienne= domaine_vf.elem_faces(elem_eulerien, num_face_face);
          d = - domaine_vf.xv(face_eulerienne,num_face_face%dimension);
        }

    }

  return ((num_face_face<dimension) ? 1 : -1);
}

void Parcours_interface::remplir_equation_plan_faces_aretes_internes(const Domaine_dis& domaine_dis)
{
// const Domaine& domaine = domaine_dis.domaine();
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis.valeur());
  const IntVect& orientation_aretes=domaine_vdf.orientation_aretes();
  const IntVect& aretes_multiples=domaine_vdf.aretes_multiples();
  const IntVect& type_arete=domaine_vdf.type_arete();
  const IntTab& Qdm = domaine_vdf.Qdm();
  const IntTab& Aretes_Som=domaine_vdf.domaine().aretes_som();
  const DoubleTab& coord_som=domaine_vdf.domaine().coord_sommets();
  //const int nb_aretes=domaine.nb_aretes();
  const int nb_elems_reels=domaine_vdf.nb_elem();
  const int nb_aretes_internes=domaine_vdf.nb_aretes_internes();
  const int premiere_arete_interne=domaine_vdf.premiere_arete_interne();
  const int nb_faces_arete=2*dimension;
  equations_plans_faces_arete_.resize(nb_aretes_internes,nb_faces_arete,4);
  equations_plans_faces_arete_=-1;
  double a=0, b=0, c=0, d=0;
  int arete_interne=0;
  for (int arete=premiere_arete_interne; arete<premiere_arete_interne+nb_aretes_internes; arete++)
    {
      int ori_arete=dimension-orientation_aretes(arete)-1;

      int face1=Qdm(arete,0);
      int face2=Qdm(arete,1);
      int face3=Qdm(arete,2);
      int face4=Qdm(arete,3);
      int elem1,elem2,elem3,elem4;
      elem1=domaine_vdf.face_voisins(face1,0);
      elem2=domaine_vdf.face_voisins(face2,0);
      elem3=domaine_vdf.face_voisins(face1,1);
      elem4=domaine_vdf.face_voisins(face2,1);
      int s1=Aretes_Som(arete,0);
      int s2=Aretes_Som(arete,1);
      //Cerr << "cg_arete_interne " << domaine_vdf.xa(arete,0) << " " << domaine_vdf.xa(arete,1) << " " << domaine_vdf.xa(arete,2) << " ori_arete " << ori_arete << finl;
      int arete_multiple=aretes_multiples(arete);
      int le_type_arete=type_arete(arete);

      /*
      if (domaine_vdf.xa(arete,0)<-2e-5 && domaine_vdf.xa(arete,0)>-3e-5 &&  domaine_vdf.xa(arete,2)==0 && domaine_vdf.xa(arete,1)<1e-5 && domaine_vdf.xa(arete,1)>-1e-5 )
        {
          Cerr << "ori_arete " << ori_arete << " arete_multiple " << arete_multiple << " le_type_arete " << le_type_arete << finl;
        }
      */
      for (int face=0; face<nb_faces_arete; face++)
        {
          a = (face%dimension==0);
          b = (face%dimension==1);
          c = (face%dimension==2);
          int dir=(face>=dimension);
          int ori_face=face%dimension;

          if (le_type_arete==2 && arete_multiple==0)// arete interne pas de joint,
            {
              if (ori_arete==0)  // arete // a l'axe x (arete YZ)
                {
                  if (face%dimension==ori_arete)
                    {
                      face1=Qdm(arete,0);
                      elem1=domaine_vdf.face_voisins(face1,0);
                      int la_face=domaine_vdf.elem_faces(elem1,ori_arete+dir*dimension);
                      d= -domaine_vdf.xv(la_face,ori_arete);
                    }
                  else if (face%dimension==1)
                    {
                      int la_face= !dir ? Qdm(arete,0) : Qdm(arete,1);
                      d= -domaine_vdf.xv(la_face,ori_face);
                    }
                  else if (face%dimension==2)
                    {
                      int la_face=!dir ? Qdm(arete,2) : Qdm(arete,3);
                      d= -domaine_vdf.xv(la_face,ori_face);
                    }
                }

              else if (ori_arete==1) // arete // a l'axe y (arete XZ)
                {

                  if (face%dimension==0)
                    {
                      int la_face= !dir ? Qdm(arete,0) : Qdm(arete,1);
                      d= -domaine_vdf.xv(la_face,ori_face);
                    }
                  else if (face%dimension==ori_arete)
                    {
                      face1=Qdm(arete,0);
                      elem1=domaine_vdf.face_voisins(face1,0);
                      int la_face=domaine_vdf.elem_faces(elem1,ori_arete+dir*dimension);
                      d= -domaine_vdf.xv(la_face,ori_arete);
                    }
                  else if (face%dimension==2)
                    {
                      int la_face= !dir ? Qdm(arete,2) : Qdm(arete,3);
                      d= -domaine_vdf.xv(la_face,ori_face);
                    }
                }
              else if (ori_arete==2) // arete // a l'axe z (arete XY)
                {
                  if (face%dimension==0)
                    {
                      int la_face=(face<dimension) ? Qdm(arete,0) : Qdm(arete,1);
                      d= -domaine_vdf.xv(la_face,ori_face);
                    }
                  else if (face%dimension==1)
                    {
                      int la_face=(face<dimension) ? Qdm(arete,2) : Qdm(arete,3);
                      d= -domaine_vdf.xv(la_face,ori_face);
                    }
                  else if (face%dimension==ori_arete)
                    {
                      face1=Qdm(arete,0);
                      elem1=domaine_vdf.face_voisins(face1,0);
                      int la_face=domaine_vdf.elem_faces(elem1,ori_arete+dir*dimension);
                      d= -domaine_vdf.xv(la_face,ori_arete);
                    }
                }

            }
          else if (le_type_arete==2 && arete_multiple>0)
            {
              if (ori_arete==0)
                {
                  if (face%dimension==0)
                    {
                      int sommet=(face<dimension) ? s1 : s2;
                      d=-coord_som(sommet,0);
                      //continue;
                    }
                  if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem2<nb_elems_reels && elem4<nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,1);
                          else d=-domaine_vdf.xv(face2,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem2>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem3<nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,1);
                          else d=-domaine_vdf.xa(arete,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3<nb_elems_reels && elem3<nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,1);
                          else d=-domaine_vdf.xv(face2,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem3>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem2<nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,1);
                          else d=-domaine_vdf.xv(face2,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xa(arete,2);
                        }
                    }
                  else if (elem2>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,1);
                          else d=-domaine_vdf.xa(arete,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xa(arete,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,1);
                          else d=-domaine_vdf.xv(face2,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xa(arete,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3>=nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,1);
                          else d=-domaine_vdf.xv(face2,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,1);
                          else d=-domaine_vdf.xa(arete,1);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                }
              else if (ori_arete==1)
                {
                  if (face%dimension==1)
                    {
                      int sommet=(face<dimension) ? s1 : s2;
                      d=-coord_som(sommet,1);
                      //continue;
                    }
                  if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem2<nb_elems_reels && elem4<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem2>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem3<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xa(arete,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3<nb_elems_reels && elem3<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem3>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem2<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xa(arete,2);
                        }
                    }
                  else if (elem2>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xa(arete,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xa(arete,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,0);
                          else d=-domaine_vdf.xa(face2,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,2);
                          else d=-domaine_vdf.xa(arete,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xa(arete,0);
                        }
                      else if (face%dimension==2)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,2);
                          else d=-domaine_vdf.xv(face4,2);
                        }
                    }
                }
              else if (ori_arete==2)
                {
                  if (face%dimension==2)
                    {
                      int sommet=(face<dimension) ? s1 : s2;
                      d=-coord_som(sommet,2);
                      // continue;
                    }
                  if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem2<nb_elems_reels && elem4<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,1);
                          else d=-domaine_vdf.xv(face4,1);
                        }
                    }
                  else if (elem2>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem3<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xa(arete,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,1);
                          else d=-domaine_vdf.xv(face4,1);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3<nb_elems_reels && elem3<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,1);
                          else d=-domaine_vdf.xv(face4,1);
                        }
                    }
                  else if (elem3>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem2<nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,1);
                          else d=-domaine_vdf.xa(arete,1);
                        }
                    }
                  else if (elem2>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xa(arete,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,1);
                          else d=-domaine_vdf.xa(arete,1);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face3,1);
                          else d=-domaine_vdf.xa(arete,1);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,0);
                          else d=-domaine_vdf.xv(face2,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,1);
                          else d=-domaine_vdf.xv(face4,1);
                        }
                    }
                  else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem4>=nb_elems_reels)
                    {
                      if (face%dimension==0)
                        {
                          if (face<dimension) d=-domaine_vdf.xv(face1,0);
                          else d=-domaine_vdf.xa(arete,0);
                        }
                      else if (face%dimension==1)
                        {
                          if (face<dimension) d=-domaine_vdf.xa(arete,1);
                          else d=-domaine_vdf.xv(face4,1);
                        }
                    }
                }
            }

          equations_plans_faces_arete_(arete_interne,face,0)=a;
          equations_plans_faces_arete_(arete_interne,face,1)=b;
          equations_plans_faces_arete_(arete_interne,face,2)=c;
          equations_plans_faces_arete_(arete_interne,face,3)=d;

          // Cerr << "a b c d " << a << " " << b << " " << c << " " << d << finl;

        }
      arete_interne++;
    }
}

// fin EB
// Parcours des elements euleriens traverses par la facette specifiee,
// de proche en proche, en commencant par l'element_depart.

void Parcours_interface::parcours_facette(const Domaine_VF& domaine_vf,
                                          Maillage_FT_Disc& maillage,
                                          ArrOfInt& echange_facettes_numfacette,
                                          ArrOfInt& echange_facettes_numelement,
                                          int num_facette,
                                          int element_depart) const
{
  const int dimension3 = (Objet_U::dimension == 3);

  Intersections_Elem_Facettes& intersections =
    maillage.intersections_elem_facettes_;

  static ArrOfIntFT liste_elements_parcourus;
  static ArrOfIntFT elements_a_traiter;
  static ArrOfIntFT new_elements_a_traiter;

  // Initialisation des drapeaux des elements deja traites pour cette facette
  // Hypothese : on suppose que drapeaux_elements_parcourus_ est initialise
  //             a zero quand on arrive ici
  {
    intersections.get_liste_elements_traverses(num_facette,
                                               liste_elements_parcourus);
    int n = liste_elements_parcourus.size_array();
    for (int i = 0; i < n; i++)
      {
        const int element = liste_elements_parcourus[i];
        drapeaux_elements_parcourus_.setbit(element);
      }
  }

  if (drapeaux_elements_parcourus_[element_depart])
    {
      effacer_drapeaux_elements_parcourus(drapeaux_elements_parcourus_,
                                          liste_elements_parcourus);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  elements_a_traiter.resize_array(1);
  elements_a_traiter[0] = element_depart;
  liste_elements_parcourus.append_array(element_depart);
  drapeaux_elements_parcourus_.setbit(element_depart);

  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_elements_a_traiter.resize_array(0);
      const int n = elements_a_traiter.size_array();

      // Boucle sur les elements du front
      // ************************************************
      for (int i = 0; i < n; i++)
        {
          const int element = elements_a_traiter[i];

          // Calcul de l'intersection facette_element
          int code_retour;
          if (dimension3)
            code_retour = calcul_intersection_facelem_3D(domaine_vf, maillage,
                                                         num_facette, element);
          else
            code_retour = calcul_intersection_facelem_2D(domaine_vf, maillage,
                                                         num_facette, element);
          // Decryptage du code retour
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.
          if (code_retour < 0)
            {
              Journal() << "code_retour<0" << finl;
              continue;
            }
          for (int j = 0; j < nb_faces_elem_; j++)
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // Numero de l'element voisin
              const int num_face = domaine_vf.elem_faces(element, j);
              const int elem_voisin = domaine_vf.face_voisins(num_face, 0)
                                      + domaine_vf.face_voisins(num_face, 1)
                                      - element;

              // S'il n'y a pas de voisin, on est au bord du domaine
              if (elem_voisin < 0)
                continue;

              if (elem_voisin < nb_elements_reels_)
                {
                  // Si c'est un element reel, et qu'il n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des elements
                  // a parcourir dans le front suivant. On marque l'element
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_elements_parcourus_.testsetbit(elem_voisin))
                    {
                      liste_elements_parcourus.append_array(elem_voisin);
                      new_elements_a_traiter.append_array(elem_voisin);
                    }
                }
              else
                {
                  // Si le voisin est un element virtuel, il faut transmettre
                  // au processeur voisin. On stocke le numero de la facette
                  // et l'element voisin.
                  echange_facettes_numfacette.append_array(num_facette);
                  echange_facettes_numelement.append_array(elem_voisin);
                }
            }
        }
      // ************************************************
      // Mise a jour du "front" des elements a traiter
      elements_a_traiter = new_elements_a_traiter;
    }
  while (elements_a_traiter.size_array() > 0);
  // ********************************************************************

  effacer_drapeaux_elements_parcourus(drapeaux_elements_parcourus_,
                                      liste_elements_parcourus);
}

// Parcours des faces euleriens traverses par la facette specifiee,
// de proche en proche, en commencant par la face_depart.

void Parcours_interface::parcours_facette_face_x(const Domaine_VF& domaine_vf,
                                                 Maillage_FT_Disc& maillage,
                                                 ArrOfInt& echange_facettes_numfacette,
                                                 ArrOfInt& echange_facettes_numface,
                                                 int num_facette,
                                                 int face_depart, Intersections_Face_Facettes& intersections_face,
                                                 ArrOfBit& drapeaux_faces_parcourues) const
{

  const int ori=0;
  const int dimension3 = (Objet_U::dimension == 3);

  static ArrOfIntFT liste_faces_parcourues_x;
  static ArrOfIntFT faces_a_traiter_x;
  static ArrOfIntFT new_faces_a_traiter_x;

  const ArrOfInt& faces_doubles = domaine_vf.faces_doubles();
  const IntTab& faces_doubles_virt_pe_num = domaine_vf.faces_doubles_virt_pe_num();

  // Initialisation des drapeaux des faces deja traitees pour cette facette
  // Hypothese : on suppose que drapeaux_faces_parcourues_i_ est initialise
  //             a zero quand on arrive ici

  {
    intersections_face.get_liste_faces_traversees(num_facette,
                                                  liste_faces_parcourues_x);
    int n = liste_faces_parcourues_x.size_array();
    for (int i = 0; i < n; i++)
      {
        const int face = liste_faces_parcourues_x[i];
        drapeaux_faces_parcourues.setbit(face);
      }
  }

  if (drapeaux_faces_parcourues[face_depart])
    {
      effacer_drapeaux_faces_parcourues(drapeaux_faces_parcourues,
                                        liste_faces_parcourues_x);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  faces_a_traiter_x.resize_array(1);
  faces_a_traiter_x[0] = face_depart;
  liste_faces_parcourues_x.append_array(face_depart);
  drapeaux_faces_parcourues.setbit(face_depart);
  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_faces_a_traiter_x.resize_array(0);
      const int n = faces_a_traiter_x.size_array();

      // Boucle sur les faces du front
      // ************************************************
      for (int i = 0; i < n; i++)
        {
          const int face = faces_a_traiter_x[i];

          const int face_double = faces_doubles(face);
          // Calcul de l'intersection facette_face
          int code_retour=0;
          if (dimension3)
            code_retour = calcul_intersection_facette_face_3D(domaine_vf, maillage,
                                                              num_facette, face, intersections_face);
          else
            {
              exit();
            }
          // Decryptage du code retour
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.

          if (code_retour < 0)
            {
              //Journal() << "code_retour<0 parcours face x" << finl;
              continue;
            }
          for (int j = 0; j < nb_faces_elem_; j++) // on a aussi nb_faces_elem_ faces pour le volume de controle des faces
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // On cherche l'element voisin
              const int direction = (j>=dimension);
              const int elem0=domaine_vf.face_voisins(face,0);
              const int elem1=domaine_vf.face_voisins(face,1);

              int num_elem;
              int face_voisine;

              if (j%dimension==ori)
                {
                  num_elem = (direction==0) ? elem0 : elem1;
                  if (num_elem < 0)
                    continue;
                  face_voisine = domaine_vf.elem_faces(num_elem, ori+direction*dimension);
                }
              else
                {
                  const int elem_2 = (elem0>=0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem0,j),direction) : -1;
                  const int elem_3 = (elem1>=0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem1,j),direction) : -1;
                  num_elem = std::max(elem_2,elem_3);
                  // S'il n'y a pas de voisin, on est au bord du domaine
                  if (num_elem < 0)
                    continue;
                  face_voisine =  (elem_2>=0) ? domaine_vf.elem_faces(elem_2, ori+dimension) : domaine_vf.elem_faces(elem_3, ori);
                }


              // si face est une face_double
              if (face_double && elem0 >= domaine_vf.nb_elem() && j==ori)
                {
                  face_voisine=face;
                }
              if (face_double && elem1 >= domaine_vf.nb_elem() && j==ori+dimension)
                {
                  face_voisine=face;
                }

              const int face_voisine_double = (faces_doubles_virt_pe_num(face_voisine,0)>0);

              if (face_voisine < nb_faces_reelles_ )
                {
                  // Si c'est une face reelle, et qu'elle n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des faces
                  // a parcourir dans le front suivant. On marque la face
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_faces_parcourues.testsetbit(face_voisine))
                    {
                      liste_faces_parcourues_x.append_array(face_voisine);
                      new_faces_a_traiter_x.append_array(face_voisine);
                    }
                }

              const int face_a_echanger=(elem0 >= domaine_vf.nb_elem() && j==ori) || (elem1 >= domaine_vf.nb_elem() && j==ori+dimension);
              if (face_voisine >= nb_faces_reelles_ || (domaine_vf.faces_doubles()(face)&&face_voisine<nb_faces_reelles_&&face_a_echanger))
                {
                  // Si le voisin est une face virtuelle, il faut transmettre
                  // au processeur voisin. On stocke le numero de la facette
                  // et l'element voisin.
                  if (face_voisine_double && face_voisine >= nb_faces_reelles_) //face_voisine_double  on met 2 fois la fa7 et la face voisine pour les envoyer aux 2 procs dans Maillage_FT_Disc::echanger_facettes_face_x (y, z)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                    }
                  else
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                    }

                }
            }
        }
      // ************************************************
      // Mise a jour du "front" des faces a traiter
      faces_a_traiter_x = new_faces_a_traiter_x;
    }
  while (faces_a_traiter_x.size_array() > 0);
  // ********************************************************************
  effacer_drapeaux_faces_parcourues(drapeaux_faces_parcourues,
                                    liste_faces_parcourues_x);
}

void Parcours_interface::parcours_facette_face_y(const Domaine_VF& domaine_vf,
                                                 Maillage_FT_Disc& maillage,
                                                 ArrOfInt& echange_facettes_numfacette,
                                                 ArrOfInt& echange_facettes_numface,
                                                 int num_facette,
                                                 int face_depart, Intersections_Face_Facettes& intersections_face,
                                                 ArrOfBit& drapeaux_faces_parcourues) const
{

  const int ori=1;
  const int dimension3 = (Objet_U::dimension == 3);

  static ArrOfIntFT liste_faces_parcourues_y;
  static ArrOfIntFT faces_a_traiter_y;
  static ArrOfIntFT new_faces_a_traiter_y;

  const ArrOfInt faces_doubles = domaine_vf.faces_doubles();
  const IntTab& faces_doubles_virt_pe_num = domaine_vf.faces_doubles_virt_pe_num();

  // Initialisation des drapeaux des faces deja traitees pour cette facette
  // Hypothese : on suppose que drapeaux_faces_parcourues_i_ est initialise
  //             a zero quand on arrive ici
  {
    intersections_face.get_liste_faces_traversees(num_facette,
                                                  liste_faces_parcourues_y);
    int n = liste_faces_parcourues_y.size_array();
    for (int i = 0; i < n; i++)
      {
        const int face = liste_faces_parcourues_y[i];
        drapeaux_faces_parcourues.setbit(face);
      }
  }

  if (drapeaux_faces_parcourues[face_depart])
    {
      effacer_drapeaux_faces_parcourues(drapeaux_faces_parcourues,
                                        liste_faces_parcourues_y);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  faces_a_traiter_y.resize_array(1);
  faces_a_traiter_y[0] = face_depart;
  liste_faces_parcourues_y.append_array(face_depart);
  drapeaux_faces_parcourues.setbit(face_depart);
  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_faces_a_traiter_y.resize_array(0);
      const int n = faces_a_traiter_y.size_array();

      // Boucle sur les faces du front
      // ************************************************
      for (int i = 0; i < n; i++)
        {

          const int face = faces_a_traiter_y[i];

          const int face_double = faces_doubles(face);
          // Calcul de l'intersection facette_face
          int code_retour=0;
          if (dimension3)
            code_retour = calcul_intersection_facette_face_3D(domaine_vf, maillage,
                                                              num_facette, face, intersections_face);
          else
            {
              exit();
            }
          // Decryptage du code retourF
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.

          if (code_retour < 0)
            {
              //Journal() << "code_retour<0 parcours face y"<< finl;
              continue;
            }

          for (int j = 0; j < nb_faces_elem_; j++) // on a aussi nb_faces_elem_ faces pour le volume de controle des faces
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // On cherche l'element voisin
              const int direction = (j>=dimension);
              const int elem0=domaine_vf.face_voisins(face,0);
              const int elem1=domaine_vf.face_voisins(face,1);

              int num_elem;
              int face_voisine;

              if (j%dimension==ori)
                {
                  num_elem = (direction==0) ? elem0 : elem1;
                  if (num_elem < 0)
                    continue;
                  face_voisine = domaine_vf.elem_faces(num_elem, ori+direction*dimension);
                }
              else
                {
                  const int elem_2 = (elem0>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem0,j),direction) : -1;
                  const int elem_3 = (elem1>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem1,j),direction) : -1;
                  num_elem = std::max(elem_2,elem_3);
                  // S'il n'y a pas de voisin, on est au bord du domaine
                  if (num_elem < 0)
                    continue;
                  face_voisine =  (elem_2>0) ? domaine_vf.elem_faces(elem_2, ori+dimension) : domaine_vf.elem_faces(elem_3, ori);
                }

              if (face_double && elem0 >= domaine_vf.nb_elem() && j==ori)
                {
                  face_voisine=face;
                }
              if (face_double && elem1 >= domaine_vf.nb_elem() && j==ori+dimension)
                {
                  face_voisine=face;
                }
              const int face_voisine_double = (faces_doubles_virt_pe_num(face_voisine,0)>0);

              if (face_voisine < nb_faces_reelles_ )
                {
                  // Si c'est une face reelle, et qu'elle n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des faces
                  // a parcourir dans le front suivant. On marque la face
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_faces_parcourues.testsetbit(face_voisine))
                    {
                      liste_faces_parcourues_y.append_array(face_voisine);
                      new_faces_a_traiter_y.append_array(face_voisine);
                    }
                }
              const int face_a_echanger=(elem0 >= domaine_vf.nb_elem() && j==ori) || (elem1 >= domaine_vf.nb_elem() && j==ori+dimension);
              if (face_voisine >= nb_faces_reelles_ || (domaine_vf.faces_doubles()(face_voisine)&&face_voisine<nb_faces_reelles_&&face_a_echanger))
                {
                  // Si le voisin est une face virtuelle, il faut transmettre
                  // au processeur voisin. On stocke le numero de la facette
                  // et l'element voisin.
                  if (face_voisine_double && face_voisine >= nb_faces_reelles_ ) //  && on met 2 fois la fa7 et la face voisine pour les envoyer aux 2 procs dans Maillage_FT_Disc::echanger_facettes_face_x (y, z)
                    {

                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                    }
                  else
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                    }
                }
            }
        }


      // ************************************************
      // Mise a jour du "front" des faces a traiter
      faces_a_traiter_y = new_faces_a_traiter_y;
    }
  while (faces_a_traiter_y.size_array() > 0);
  // ********************************************************************
  effacer_drapeaux_faces_parcourues(drapeaux_faces_parcourues,
                                    liste_faces_parcourues_y);
}

void Parcours_interface::parcours_facette_face_z(const Domaine_VF& domaine_vf,
                                                 Maillage_FT_Disc& maillage,
                                                 ArrOfInt& echange_facettes_numfacette,
                                                 ArrOfInt& echange_facettes_numface,
                                                 int num_facette,
                                                 int face_depart, Intersections_Face_Facettes& intersections_face,
                                                 ArrOfBit& drapeaux_faces_parcourues) const
{
  const int ori=2;
  const int dimension3 = (Objet_U::dimension == 3);

  static ArrOfIntFT liste_faces_parcourues_z;
  static ArrOfIntFT faces_a_traiter_z;
  static ArrOfIntFT new_faces_a_traiter_z;

  const ArrOfInt faces_doubles = domaine_vf.faces_doubles();
  const IntTab& faces_doubles_virt_pe_num = domaine_vf.faces_doubles_virt_pe_num();

  // Initialisation des drapeaux des faces deja traitees pour cette facette
  // Hypothese : on suppose que drapeaux_faces_parcourues_i_ est initialise
  //             a zero quand on arrive ici
  {
    intersections_face.get_liste_faces_traversees(num_facette,
                                                  liste_faces_parcourues_z);
    int n = liste_faces_parcourues_z.size_array();
    for (int i = 0; i < n; i++)
      {
        const int face = liste_faces_parcourues_z[i];
        drapeaux_faces_parcourues.setbit(face);
      }
  }

  if (drapeaux_faces_parcourues[face_depart])
    {
      effacer_drapeaux_faces_parcourues(drapeaux_faces_parcourues,
                                        liste_faces_parcourues_z);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  faces_a_traiter_z.resize_array(1);
  faces_a_traiter_z[0] = face_depart;
  liste_faces_parcourues_z.append_array(face_depart);
  drapeaux_faces_parcourues.setbit(face_depart);
  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_faces_a_traiter_z.resize_array(0);
      const int n = faces_a_traiter_z.size_array();

      // Boucle sur les faces du front
      // ************************************************
      for (int i = 0; i < n; i++)
        {
          const int face = faces_a_traiter_z[i];

          const int face_double = faces_doubles(face) && face<domaine_vf.nb_faces();
          // Calcul de l'intersection facette_face
          int code_retour=0;
          if (dimension3)
            code_retour = calcul_intersection_facette_face_3D(domaine_vf, maillage,
                                                              num_facette, face, intersections_face);
          else
            {
              exit();
            }
          // Decryptage du code retour
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.
          if (code_retour < 0)
            {
              //Journal() << "code_retour<0 parcours face z" << finl;
              continue;
            }
          for (int j = 0; j < nb_faces_elem_; j++) // on a aussi nb_faces_elem_ faces pour le volume de controle des faces
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // On cherche l'element voisin
              const int direction = (j>=dimension);
              const int elem0=domaine_vf.face_voisins(face,0);
              const int elem1=domaine_vf.face_voisins(face,1);

              int num_elem;
              int face_voisine;

              if (j%dimension==ori)
                {
                  num_elem = (direction==0) ? elem0 : elem1;
                  if (num_elem < 0)
                    continue;
                  face_voisine = domaine_vf.elem_faces(num_elem, ori+direction*dimension);
                }
              else
                {
                  const int elem_2 = (elem0>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem0,j),direction) : -1;
                  const int elem_3 = (elem1>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem1,j),direction) : -1;
                  num_elem = std::max(elem_2,elem_3);
                  // S'il n'y a pas de voisin, on est au bord du domaine
                  if (num_elem < 0)
                    continue;
                  face_voisine =  (elem_2>0) ? domaine_vf.elem_faces(elem_2, ori+dimension) : domaine_vf.elem_faces(elem_3, ori);
                }

              if (face_double && elem0 >= domaine_vf.nb_elem() && j==ori)
                {
                  face_voisine=face;
                }
              else if (face_double && elem1 >= domaine_vf.nb_elem() && j==ori+dimension)
                {
                  face_voisine=face;
                }

              const int face_voisine_double = (faces_doubles_virt_pe_num(face_voisine,0)>0);

              if (face_voisine < nb_faces_reelles_ )
                {
                  // Si c'est une face reelle, et qu'elle n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des faces
                  // a parcourir dans le front suivant. On marque la face
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_faces_parcourues.testsetbit(face_voisine))
                    {
                      liste_faces_parcourues_z.append_array(face_voisine);
                      new_faces_a_traiter_z.append_array(face_voisine);
                    }
                }
              const int face_a_echanger=(elem0 >= domaine_vf.nb_elem() && j==ori) || (elem1 >= domaine_vf.nb_elem() && j==ori+dimension); // si face_voisine=face alors on l'envoie a l'autre proc proprietaire
              if (face_voisine >= nb_faces_reelles_ || (domaine_vf.faces_doubles()(face_voisine)&&face_voisine<nb_faces_reelles_&&face_a_echanger))
                {
                  // Si le voisin est une face virtuelle, il faut transmettre
                  // au processeur voisin. On stocke le numero de la facette
                  // et l'element voisin.
                  if (face_voisine_double && face_voisine >= nb_faces_reelles_ ) // face_voisine_double on met 2 fois la fa7 et la face voisine pour les envoyer aux 2 procs dans Maillage_FT_Disc::echanger_facettes_face_x (y, z)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                    }
                  else
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numface.append_array(face_voisine);
                    }
                }
            }
        }


      // ************************************************
      // Mise a jour du "front" des faces a traiter
      faces_a_traiter_z = new_faces_a_traiter_z;
    }
  while (faces_a_traiter_z.size_array() > 0);
  // ********************************************************************
  effacer_drapeaux_faces_parcourues(drapeaux_faces_parcourues,
                                    liste_faces_parcourues_z);
}

void Parcours_interface::parcours_facette_arete_x(const Domaine_VF& domaine_vf,
                                                  Maillage_FT_Disc& maillage,
                                                  ArrOfInt& echange_facettes_numfacette,
                                                  ArrOfInt& echange_facettes_numarete,
                                                  int num_facette,
                                                  int arete_depart, Intersections_Arete_Facettes& intersections_arete,
                                                  ArrOfBit& drapeaux_aretes_parcourues) const
{
  const int ori_arete=0;
  const int dimension3 = (Objet_U::dimension == 3);

  static ArrOfIntFT liste_aretes_parcourues_arete_x;
  static ArrOfIntFT arete_a_traiter_arete_x;
  static ArrOfIntFT new_aretes_a_traiter_arete_x;

  const Maillage_FT_Disc& maillage_const=maillage;
  const Domaine_dis& domaine_dis = maillage_const.refdomaine_dis_.valeur();
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis.valeur());

  const IntTab& Elem_Aretes=domaine_vf.domaine().elem_aretes();
  const IntTab& Qdm=domaine_vdf.Qdm();
  const IntVect& orientation_arete=domaine_vdf.orientation_aretes();
  const IntVect& type_arete=domaine_vdf.type_arete();
  const IntVect& aretes_multiples=domaine_vdf.aretes_multiples();

  // Initialisation des drapeaux des faces deja traitees pour cette facette
  // Hypothese : on suppose que drapeaux_faces_parcourues_i_ est initialise
  //             a zero quand on arrive ici

  {
    intersections_arete.get_liste_aretes_traversees(num_facette,
                                                    liste_aretes_parcourues_arete_x);
    int n = liste_aretes_parcourues_arete_x.size_array();
    for (int i = 0; i < n; i++)
      {
        const int face = liste_aretes_parcourues_arete_x[i];
        drapeaux_aretes_parcourues.setbit(face);
      }
  }

  if (drapeaux_aretes_parcourues[arete_depart])
    {
      effacer_drapeaux_faces_parcourues(drapeaux_aretes_parcourues,
                                        liste_aretes_parcourues_arete_x);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  arete_a_traiter_arete_x.resize_array(1);
  arete_a_traiter_arete_x[0] = arete_depart;
  liste_aretes_parcourues_arete_x.append_array(arete_depart);
  drapeaux_aretes_parcourues.setbit(arete_depart);

  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_aretes_a_traiter_arete_x.resize_array(0);
      const int n = arete_a_traiter_arete_x.size_array();

      // Boucle sur les aretes
      // ************************************************
      for (int i = 0; i < n; i++)
        {
          const int arete = arete_a_traiter_arete_x[i];
          if (type_arete(arete)!=2) continue; // on ne calcule que pour les aretes_internes

          const int arete_multiple=(aretes_multiples(arete)>0);

          // Calcul de l'intersection facette_face
          int code_retour=0;
          if (dimension3)
            code_retour = calcul_intersection_facette_arete_3D(domaine_vf, maillage,
                                                               num_facette, arete, intersections_arete); // ne fonctionne que si l'arete est interne
          else
            {
              exit();
            }
          // Decryptage du code retour
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.

          if (code_retour < 0)
            {
              //Journal() << "code_retour<0 parcours face x" << finl;
              continue;
            }

          int face1=Qdm(arete,0);
          int face2=Qdm(arete,1);

          int elem1=domaine_vf.face_voisins(face1,0);
          int elem2=domaine_vf.face_voisins(face2,0);
          int elem3=domaine_vf.face_voisins(face1,1);
          int elem4=domaine_vf.face_voisins(face2,1);

          int face1_av=domaine_vf.elem_faces(elem1, ori_arete+dimension);
          int face1_arr=domaine_vf.elem_faces(elem1, ori_arete);

          int elem_av=domaine_vf.face_voisins(face1_av,1);
          int elem_arr=domaine_vf.face_voisins(face1_arr,0);


          for (int j = 0; j < nb_faces_elem_; j++) // on a aussi nb_faces_elem_ faces pour le volume de controle des faces
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // On cherche l'element voisin
              //const int direction = (j>=dimension);

              int arete_voisine;

              if (j==ori_arete) arete_voisine= (elem_arr>=0) ? Elem_Aretes(elem_arr,orientation_arete(arete)) : -1; // voir numerotation Aretes dans Hexaedre.cpp
              else if (j==ori_arete+dimension) arete_voisine=(elem_av>=0) ? Elem_Aretes(elem_av,orientation_arete(arete)):-1;
              else if (j==1) arete_voisine= (elem1>=0) ? Elem_Aretes(elem1,11) : -1;
              else if (j==1+dimension) arete_voisine=(elem2>=0) ? Elem_Aretes(elem2,2) : -1;
              else if (j==2) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,8) : -1;
              else arete_voisine=(elem3>=0) ? Elem_Aretes(elem3,2) : -1;

              //if (arete_voisine>=nb_aretes_reelles_) Cerr << "coucou" << finl;
              int arete_a_echanger=0; // Si une arete voisine est multiple, alors on l'envoie a tous les procs
              if (arete_multiple)
                {
                  if ((elem2>=nb_elements_reels_ && elem4>=nb_elements_reels_) && j==1+dimension) arete_voisine=arete,arete_a_echanger=1;
                  else if ((elem1>=nb_elements_reels_ && elem3>=nb_elements_reels_) && j==1) arete_voisine=arete,arete_a_echanger=1;

                  else if ((elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_) &&  j==2) arete_voisine=arete,arete_a_echanger=1;
                  else if ((elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_) &&  j==2+dimension) arete_voisine=arete,arete_a_echanger=1;

                  else if (elem2>=nb_elements_reels_ && elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==1+dimension || j==2+dimension)) arete_voisine=arete,arete_a_echanger=1;
                  else if (elem1>=nb_elements_reels_ && elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==1 || j==2+dimension)) arete_voisine=arete,arete_a_echanger=1;

                  else if (elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_ && elem3>=nb_elements_reels_ && (j==1 || j==2)) arete_voisine=arete,arete_a_echanger=1;
                  else if (elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==1+dimension || j==2)) arete_voisine=arete,arete_a_echanger=1;
                }

              if (arete_voisine < nb_aretes_reelles_ )
                {
                  // Si c'est une face reelle, et qu'elle n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des faces
                  // a parcourir dans le front suivant. On marque la face
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_aretes_parcourues.testsetbit(arete_voisine))
                    {
                      liste_aretes_parcourues_arete_x.append_array(arete_voisine);
                      new_aretes_a_traiter_arete_x.append_array(arete_voisine);
                    }
                }
              if (arete_voisine >= nb_aretes_reelles_  || arete_a_echanger)
                {
                  int arete_mult=aretes_multiples(arete_voisine);
                  if (arete_mult<=1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  else if (arete_mult>1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  else if (arete_mult>2)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  if (1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                }
            }
        }
      //Cerr << "compteur " << compteur << finl;
      // ************************************************
      // Mise a jour du "front" des faces a traiter
      arete_a_traiter_arete_x = new_aretes_a_traiter_arete_x;
    }
  while (arete_a_traiter_arete_x.size_array() > 0);
  // ********************************************************************
  effacer_drapeaux_faces_parcourues(drapeaux_aretes_parcourues,
                                    liste_aretes_parcourues_arete_x);
  //Cerr << "FIN  Parcours_interface::parcours_facette_arete_x"  << finl;
}

void Parcours_interface::parcours_facette_arete_y(const Domaine_VF& domaine_vf,
                                                  Maillage_FT_Disc& maillage,
                                                  ArrOfInt& echange_facettes_numfacette,
                                                  ArrOfInt& echange_facettes_numarete,
                                                  int num_facette,
                                                  int arete_depart, Intersections_Arete_Facettes& intersections_arete,
                                                  ArrOfBit& drapeaux_aretes_parcourues) const
{

  const int ori_arete=1;
  const int dimension3 = (Objet_U::dimension == 3);

  static ArrOfIntFT liste_aretes_parcourues_arete_y;
  static ArrOfIntFT arete_a_traiter_arete_y;
  static ArrOfIntFT new_aretes_a_traiter_arete_y;

  const Maillage_FT_Disc& maillage_const=maillage;
  const Domaine_dis& domaine_dis = maillage_const.refdomaine_dis_.valeur();
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis.valeur());

  const IntTab& Elem_Aretes=domaine_vf.domaine().elem_aretes();
  const IntTab& Qdm=domaine_vdf.Qdm();
  const IntVect& orientation_arete=domaine_vdf.orientation_aretes();
  const IntVect& type_arete=domaine_vdf.type_arete();
  const IntVect& aretes_multiples=domaine_vdf.aretes_multiples();

  // Initialisation des drapeaux des faces deja traitees pour cette facette
  // Hypothese : on suppose que drapeaux_faces_parcourues_i_ est initialise
  //             a zero quand on arrive ici

  {
    intersections_arete.get_liste_aretes_traversees(num_facette,
                                                    liste_aretes_parcourues_arete_y);
    int n = liste_aretes_parcourues_arete_y.size_array();
    for (int i = 0; i < n; i++)
      {
        const int face = liste_aretes_parcourues_arete_y[i];
        drapeaux_aretes_parcourues.setbit(face);
      }
  }

  if (drapeaux_aretes_parcourues[arete_depart])
    {
      effacer_drapeaux_faces_parcourues(drapeaux_aretes_parcourues,
                                        liste_aretes_parcourues_arete_y);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  arete_a_traiter_arete_y.resize_array(1);
  arete_a_traiter_arete_y[0] = arete_depart;
  liste_aretes_parcourues_arete_y.append_array(arete_depart);
  drapeaux_aretes_parcourues.setbit(arete_depart);

  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_aretes_a_traiter_arete_y.resize_array(0);
      const int n = arete_a_traiter_arete_y.size_array();

      // Boucle sur les aretes
      // ************************************************
      for (int i = 0; i < n; i++)
        {
          const int arete = arete_a_traiter_arete_y[i];



          if (type_arete(arete)!=2) continue; // on ne calcule que pour les aretes_internes

          const int arete_multiple=(aretes_multiples(arete)>0);
          /* if (domaine_vf.xa(arete,0)<-2e-5 && domaine_vf.xa(arete,0)>-3e-5 &&  domaine_vf.xa(arete,2)==0 && domaine_vf.xa(arete,1)<1e-5 && domaine_vf.xa(arete,1)>-1e-5 ) Cerr << "arete_multiple " <<  arete_multiple << "cg " << domaine_vf.xa(arete,0)
                 << " " << domaine_vf.xa(arete,1) << " " << domaine_vf.xa(arete,2) << " facette " << maillage.cg_fa7_(num_facette,0) << " " <<  maillage.cg_fa7_(num_facette,1) << " " << maillage.cg_fa7_(num_facette,2) << finl;
           */
          // Calcul de l'intersection facette_face
          int code_retour=0;
          if (dimension3)
            code_retour = calcul_intersection_facette_arete_3D(domaine_vf, maillage,
                                                               num_facette, arete, intersections_arete); // ne fonctionne que si l'arete est interne
          else
            {
              exit();
            }
          // Decryptage du code retour
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.

          if (code_retour < 0)
            {
              //Journal() << "code_retour<0 parcours arete y" << finl;
              continue;
            }

          int face1=Qdm(arete,0);
          int face2=Qdm(arete,1);

          int elem1=domaine_vf.face_voisins(face1,0);
          int elem2=domaine_vf.face_voisins(face2,0);
          int elem3=domaine_vf.face_voisins(face1,1);
          int elem4=domaine_vf.face_voisins(face2,1);

          int face1_av=domaine_vf.elem_faces(elem1, ori_arete+dimension);
          int face1_arr=domaine_vf.elem_faces(elem1, ori_arete);

          int elem_av=domaine_vf.face_voisins(face1_av,1);
          int elem_arr=domaine_vf.face_voisins(face1_arr,0);


          for (int j = 0; j < nb_faces_elem_; j++) // on a aussi nb_faces_elem_ faces pour le volume de controle des faces
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // On cherche l'element voisin
              //const int direction = (j>=dimension);

              int arete_voisine=-1;

              if (j==ori_arete) arete_voisine= (elem_arr>=0) ? Elem_Aretes(elem_arr,orientation_arete(arete)) : -1; // voir numerotation Aretes dans Hexaedre.cpp
              else if (j==ori_arete+dimension) arete_voisine= (elem_av>=0) ? Elem_Aretes(elem_av,orientation_arete(arete)) : -1;

              else if (j==0) arete_voisine= (elem1>=0) ? Elem_Aretes(elem1,10) : -1; // voir numerotation Aretes dans Hexaedre.cpp
              else if (j==0+dimension) arete_voisine=(elem2>=0) ?Elem_Aretes(elem2,1):-1;
              else if (j==2) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,7):-1;
              else arete_voisine=(elem3>=0) ? Elem_Aretes(elem3,1):-1;

              int arete_a_echanger=0;
              if (arete_multiple)
                {
                  if ((elem2>=nb_elements_reels_ && elem4>=nb_elements_reels_) && j==0+dimension) arete_voisine=arete,arete_a_echanger=1;
                  else if ((elem1>=nb_elements_reels_ && elem3>=nb_elements_reels_) && j==0) arete_voisine=arete,arete_a_echanger=1;

                  else if ((elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_) &&  j==2) arete_voisine=arete,arete_a_echanger=1;
                  else if ((elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_) &&  j==2+dimension) arete_voisine=arete,arete_a_echanger=1;

                  else if (elem2>=nb_elements_reels_ && elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==0+dimension || j==2+dimension)) arete_voisine=arete,arete_a_echanger=1;
                  else if (elem1>=nb_elements_reels_ && elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==0 || j==2+dimension)) arete_voisine=arete,arete_a_echanger=1;

                  else if (elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_ && elem3>=nb_elements_reels_ && (j==0 || j==2)) arete_voisine=arete,arete_a_echanger=1;
                  else if (elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==0+dimension || j==2)) arete_voisine=arete,arete_a_echanger=1;
                }

              if (arete_voisine < nb_aretes_reelles_ )
                {
                  // Si c'est une face reelle, et qu'elle n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des faces
                  // a parcourir dans le front suivant. On marque la face
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_aretes_parcourues.testsetbit(arete_voisine))
                    {
                      liste_aretes_parcourues_arete_y.append_array(arete_voisine);
                      new_aretes_a_traiter_arete_y.append_array(arete_voisine);
                    }
                }
              if (arete_voisine >= nb_aretes_reelles_ || arete_a_echanger)
                {
                  int arete_mult=aretes_multiples(arete_voisine);

                  if (arete_mult<=1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  else if (arete_mult>1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  else if (arete_mult>2)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }

                  if (1) //(arete_voisine >= nb_aretes_reelles_)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                }
            }
        }
      // ************************************************
      // Mise a jour du "front" des faces a traiter
      arete_a_traiter_arete_y = new_aretes_a_traiter_arete_y;
    }
  while (arete_a_traiter_arete_y.size_array() > 0);
  // ********************************************************************
  effacer_drapeaux_faces_parcourues(drapeaux_aretes_parcourues,
                                    liste_aretes_parcourues_arete_y);
  //Cerr << "FIN  Parcours_interface::parcours_facette_arete_y"  << finl;
}

void Parcours_interface::parcours_facette_arete_z(const Domaine_VF& domaine_vf,
                                                  Maillage_FT_Disc& maillage,
                                                  ArrOfInt& echange_facettes_numfacette,
                                                  ArrOfInt& echange_facettes_numarete,
                                                  int num_facette,
                                                  int arete_depart, Intersections_Arete_Facettes& intersections_arete,
                                                  ArrOfBit& drapeaux_aretes_parcourues) const
{
  const int ori_arete=2;
  const int dimension3 = (Objet_U::dimension == 3);

  static ArrOfIntFT liste_aretes_parcourues_arete_z;
  static ArrOfIntFT arete_a_traiter_arete_z;
  static ArrOfIntFT new_aretes_a_traiter_arete_z;

  const Maillage_FT_Disc& maillage_const=maillage;
  const Domaine_dis& domaine_dis = maillage_const.refdomaine_dis_.valeur();
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis.valeur());

  const IntTab& Elem_Aretes=domaine_vf.domaine().elem_aretes();
  const IntTab& Qdm=domaine_vdf.Qdm();
  const IntVect& orientation_arete=domaine_vdf.orientation_aretes();
  const IntVect& type_arete=domaine_vdf.type_arete();
  const IntVect& aretes_multiples=domaine_vdf.aretes_multiples();


  // Initialisation des drapeaux des faces deja traitees pour cette facette
  // Hypothese : on suppose que drapeaux_faces_parcourues_i_ est initialise
  //             a zero quand on arrive ici

  {
    intersections_arete.get_liste_aretes_traversees(num_facette,
                                                    liste_aretes_parcourues_arete_z);
    int n = liste_aretes_parcourues_arete_z.size_array();
    for (int i = 0; i < n; i++)
      {
        const int face = liste_aretes_parcourues_arete_z[i];
        drapeaux_aretes_parcourues.setbit(face);
      }
  }

  if (drapeaux_aretes_parcourues[arete_depart])
    {
      effacer_drapeaux_faces_parcourues(drapeaux_aretes_parcourues,
                                        liste_aretes_parcourues_arete_z);
      return;
    }

  // Initialisation du "front" d'elements a traiter
  arete_a_traiter_arete_z.resize_array(1);
  arete_a_traiter_arete_z[0] = arete_depart;
  liste_aretes_parcourues_arete_z.append_array(arete_depart);
  drapeaux_aretes_parcourues.setbit(arete_depart);

  // Boucle sur les fronts successifs d'elements traverses par la facette
  // ********************************************************************
  do
    {
      new_aretes_a_traiter_arete_z.resize_array(0);
      const int n = arete_a_traiter_arete_z.size_array();

      // Boucle sur les aretes
      // ************************************************
      for (int i = 0; i < n; i++)
        {
          const int arete = arete_a_traiter_arete_z[i];
          if (type_arete(arete)!=2) continue; // on ne calcule que pour les aretes_internes

          const int arete_multiple=(aretes_multiples(arete)>0);

          // Calcul de l'intersection facette_face
          int code_retour=0;
          if (dimension3)
            code_retour = calcul_intersection_facette_arete_3D(domaine_vf, maillage,
                                                               num_facette, arete, intersections_arete); // ne fonctionne que si l'arete est interne
          else
            {
              exit();
            }
          // Decryptage du code retour
          // S'il est negatif, erreur grossiere dans l'algorithme...
          // En debug, on regarde de pres, en production, on ferme les yeux.

          if (code_retour < 0)
            {
              //Journal() << "code_retour<0 parcours arete z" << finl;
              continue;
            }

          int face1=Qdm(arete,0);
          int face2=Qdm(arete,1);

          int elem1=domaine_vf.face_voisins(face1,0);
          int elem2=domaine_vf.face_voisins(face2,0);
          int elem3=domaine_vf.face_voisins(face1,1);
          int elem4=domaine_vf.face_voisins(face2,1);

          int face1_av=domaine_vf.elem_faces(elem1, ori_arete+dimension);
          int face1_arr=domaine_vf.elem_faces(elem1, ori_arete);

          int elem_av=domaine_vf.face_voisins(face1_av,1);
          int elem_arr=domaine_vf.face_voisins(face1_arr,0);


          for (int j = 0; j < nb_faces_elem_; j++) // on a aussi nb_faces_elem_ faces pour le volume de controle des faces
            {
              const int masque = 1 << j;

              // La facette traverse-t-elle la face j ?
              if ((code_retour & masque) == 0)
                continue;

              // On cherche l'element voisin
              int arete_voisine;
              if (j==ori_arete) arete_voisine=  (elem_arr>=0) ? Elem_Aretes(elem_arr,orientation_arete(arete)) : -1; // voir numerotation Aretes dans Hexaedre.cpp
              else if (j==ori_arete+dimension) arete_voisine = (elem_av>=0) ? Elem_Aretes(elem_av,orientation_arete(arete)):-1;
              else if (j==0) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,9) : -1; // voir numerotation Aretes dans Hexaedre.cpp
              else if (j==0+dimension) arete_voisine=(elem2>=0) ? Elem_Aretes(elem2,0) : -1;
              else if (j==1) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,6) : -1;
              else arete_voisine=(elem3>=0) ? Elem_Aretes(elem3,0) : -1;

              int arete_a_echanger=0;
              if (arete_multiple)
                {
                  if ((elem2>=nb_elements_reels_ && elem4>=nb_elements_reels_) && j==0+dimension) arete_voisine=arete,arete_a_echanger=1;
                  else if ((elem1>=nb_elements_reels_ && elem3>=nb_elements_reels_) && j==0) arete_voisine=arete,arete_a_echanger=1;

                  else if ((elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_) &&  j==1) arete_voisine=arete,arete_a_echanger=1;
                  else if ((elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_) &&  j==1+dimension) arete_voisine=arete,arete_a_echanger=1;

                  else if (elem2>=nb_elements_reels_ && elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==0+dimension || j==1+dimension)) arete_voisine=arete,arete_a_echanger=1;
                  else if (elem1>=nb_elements_reels_ && elem3>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==0 || j==1+dimension)) arete_voisine=arete,arete_a_echanger=1;

                  else if (elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_ && elem3>=nb_elements_reels_ && (j==0 || j==1)) arete_voisine=arete,arete_a_echanger=1;
                  else if (elem1>=nb_elements_reels_ && elem2>=nb_elements_reels_ && elem4>=nb_elements_reels_ && (j==0+dimension || j==1)) arete_voisine=arete,arete_a_echanger=1;

                }
              if (arete_voisine < nb_aretes_reelles_ )
                {
                  // Si c'est une arete reelle, et qu'elle n'a pas encore ete
                  // parcouru, on l'ajoute a la liste des aretes
                  // a parcourir dans le front suivant. On marque l'arete
                  // pour ne pas l'ajouter plusieurs fois a la liste.
                  if ( ! drapeaux_aretes_parcourues.testsetbit(arete_voisine))
                    {
                      liste_aretes_parcourues_arete_z.append_array(arete_voisine);
                      new_aretes_a_traiter_arete_z.append_array(arete_voisine);
                    }
                }
              if (arete_voisine >= nb_aretes_reelles_ || arete_a_echanger)
                {
                  int arete_mult=aretes_multiples(arete_voisine);
                  if (arete_mult<=1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  else if (arete_mult>1)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  else if (arete_mult>2)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                  if (1) //(arete_voisine >= nb_aretes_reelles_)
                    {
                      echange_facettes_numfacette.append_array(num_facette);
                      echange_facettes_numarete.append_array(arete_voisine);
                    }
                }


            }
        }
      // ************************************************
      // Mise a jour du "front" des aretes a traiter
      arete_a_traiter_arete_z = new_aretes_a_traiter_arete_z;
    }
  while (arete_a_traiter_arete_z.size_array() > 0);
  // ********************************************************************
  effacer_drapeaux_faces_parcourues(drapeaux_aretes_parcourues,
                                    liste_aretes_parcourues_arete_z);
  //Cerr << "FIN  Parcours_interface::parcours_facette_arete_z"  << finl;
}

// ==========================================================================
// ==========================================================================

// Calcul de l'intersection entre une facette et un element eulerien.
// Code de retour :
//  -1 si l'intersection est vide (pas normal)
//  >=0 : les bits du code retour indiquent les plans de l'element
//        qui coupent la facette (bit 0 : premiere face de l'element,
//        bit 1 : deuxieme face, etc...)

int Parcours_interface::calcul_intersection_facelem_2D(
  const Domaine_VF& domaine_vf,
  Maillage_FT_Disc& maillage,
  int num_facette, int num_element) const
{
  assert(Objet_U::dimension == 2);

  const int sommet0 = maillage.facettes_(num_facette, 0);
  const int sommet1 = maillage.facettes_(num_facette, 1);
  // Extremites du segment.
  double x0 = maillage.sommets_(sommet0, 0);
  double y0 = maillage.sommets_(sommet0, 1);
  double x1 = maillage.sommets_(sommet1, 0);
  double y1 = maillage.sommets_(sommet1, 1);

  // Extremites du segment tronque (coordonnees barycentriques)
  double u0 = 1.;
  double u1 = 0.;

  // Estimation de l'incertitude sur l'evaluation de la fonction plan
  double Zero = - Erreur_max_coordonnees_;
  if (correction_parcours_thomas_)
    Zero = 0.;
  // Extension du segment si un sommet au bord
  {
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    double l = sqrt(dx * dx + dy * dy);
    if (l > 0.)
      l = 1. / l;
    const double ext = 4* Erreur_max_coordonnees_ * l;
    if (maillage.sommet_ligne_contact(sommet0))
      {
        x0 -= ext * dx;
        y0 -= ext * dy;
      }
    else if (maillage.sommet_ligne_contact(sommet1))
      {
        x1 += ext * dx;
        y1 += ext * dy;
      }
  }

  // Numero de la face qui genere l'extremite de l'intersection
  // (-1 -> pas d'intersection avec le bord de l'element)
  int plan_coupe0 = -1;
  int plan_coupe1 = -1;

  // *********************************************************************
  // *********************************************************************
  // Calcul de l'intersection de la facette avec les plans qui
  // definissent l'element
  for (int num_plan = 0; num_plan < nb_faces_elem_; num_plan++)
    {
      // Coefficients de la fonction qui definit le plan (c est nul en 2D).
      double a, b, c, d;
      double signe = calcul_eq_plan(domaine_vf, num_element, num_plan, a, b, c, d);
      // Calcul de la fonction f aux deux extremites du segment d'origine
      const double f_0 = (a * x0 + b * y0 + d) * signe;
      const double f_1 = (a * x1 + b * y1 + d) * signe;
      // Calcul de la fonction aux extremites du segment tronque.
      const double f0 = f_0 * u0 + f_1 * (1. - u0);
      const double f1 = f_0 * u1 + f_1 * (1. - u1);
      const int s0_dehors = inf_strict(f0,Zero) ? 1 : 0;
      const int s1_dehors = inf_strict(f1,Zero) ? 1 : 0;

      if (s0_dehors + s1_dehors == 1)
        {
          // Un point dedans et un dehors, on calcule l'intersection
          // Coordonnee barycentrique de l'intersection
          double t = f1 / (f1 - f0);
          if (t < 0.)
            t = 0.;
          else if (t > 1.)
            t = 1.;
          if (s0_dehors)
            {
              u0 = u0 * t + u1 * (1.-t);
              plan_coupe0 = num_plan;
            }
          else
            {
              u1 = u0 * t + u1 * (1.-t);
              plan_coupe1 = num_plan;
            }
        }
      else if (s0_dehors + s1_dehors == 2)
        {
          // Les deux points sont a l'exterieur => pas d'intersection,
          // c'est un probleme interne de l'algorithme.
          compteur_erreur_grossiere++;
          return -1;
        }
    }
  // *********************************************************************
  // *********************************************************************
  // Calcul du centre de gravite de l'intersection et de la contribution
  // de volume.

  // Coordonnee barycentrique du centre de gravite de l'intersection:
  double u_centre = (u0 + u1) * 0.5;
  // Surface d'intersection (fraction entre 0 et 1):
  double surface = u0 - u1;
  // En cas d'erreur d'arrondi ...
  if (inf_strict(surface,0.)) surface = 0.;
  // Contribution de volume :
  double volume = 0.;
  // Computation of barycenter in phase 1 :
  double barycentre_phase1[3] = {0.,0.,0.};

  // Si on passe par un coin, on ne veut pas faire de calcul de contribution,
  // il risque d'etre faux. Une erreur relative a chaque extremite. On ajoute
  // Un facteur 2.5 pour etre sur.
  if ((surface > 5. * Erreur_max_coordonnees_) || correction_parcours_thomas_)
    {
      // Calcul des coordonnees des extremites de l'intersection:
      double ix0 = x0 * u0 + x1 * (1.-u0);
      double iy0 = y0 * u0 + y1 * (1.-u0);
      double ix1 = x0 * u1 + x1 * (1.-u1);
      double iy1 = y0 * u1 + y1 * (1.-u1);
      switch(type_element_)
        {
        case RECTANGLE:
          volume = volume_barycentre_rectangle(domaine_vf, num_element,
                                               ix0, iy0, ix1, iy1,
                                               Erreur_max_coordonnees_,
                                               barycentre_phase1);
          break;
        case TRIANGLE:
          if (flag_warning_code_missing)
            {
              Cerr << "WARNING : barycentre_phase1 not filled properly!!" << finl;
              Cerr << "WARNING : Calculation of barycentre_phase1 not implemented yet." << finl;
            }
          volume = volume_triangle(domaine_vf, num_element,
                                   ix0, iy0, ix1, iy1,
                                   plan_coupe0, plan_coupe1);
          break;
        default:
          Process::exit();// qu'est-ce qu'on fout la ?
          break;
        }
    }

  // Stockage des donnees relatives a cette intersection
  maillage.intersections_elem_facettes_.ajoute_intersection(num_facette,
                                                            num_element,
                                                            surface,
                                                            volume,
                                                            barycentre_phase1,
                                                            u_centre,
                                                            1. - u_centre,
                                                            0.);
  // Calcul du code retour
  int code_retour = 0;
  if (plan_coupe0 >= 0)
    code_retour |= 1 << plan_coupe0;
  if (plan_coupe1 >= 0)
    code_retour |= 1 << plan_coupe1;
  return code_retour;
}


// Calcul de la contribution de volume d'une facette a la valeur de
// l'indicatrice dans un element. C'est une fraction du volume de l'element
// comprise entre Erreur_relative_maxi et 1.-Erreur_relative_maxi

// 4 cas possibles :
//     v=B-A            v=A+B          v=B-A         v=A+B
//    ----1------    ----0------    --0--------    --1--------  <- y_haut
//   |   /:      |  |   /:      |  |  :\       |  |  :\       |
//   |  0 :      |  |  1 :      |  |  : 1      |  |  : 0      |
//   |  : :      |  |  : :      |  |  : :      |  |  : :      |
//   |B :A:      |  |  :A:  B   |  |  :A:  B   |  | B:A:      |
//   |+ :-:      |  |  :+:  +   |  |  :-:  +   |  | +:+:      |
//   |  : :      |  |  : :      |  |  : :      |  |  : :      |
//    -----------    -----------    -----------    -----------  <- y_bas
//   ^           ^
// x_gauche  x_droite

// Parametres : x0, y0, x1, y1 : coordonnees de l'extremite de l'intersection
//                               segment/element
//              plan_coupe0, plan_coupe1 : numeros de la face de l'element sur
//                                           laquelle se trouve l'extremite.

inline double Parcours_interface::volume_rectangle(const Domaine_VF& domaine_vf,
                                                   int num_element,
                                                   double x0, double y0,
                                                   double x1, double y1,
                                                   double epsilon) const
{
  const double angle_bidim_axi = bidim_axi ? Maillage_FT_Disc::angle_bidim_axi() : 0.;
  const double un_tiers = 1. / 3.;

  // Conventions TRUST VDF :
  static const int NUM_FACE_GAUCHE = 0;
  static const int NUM_FACE_BAS = 1;
  static const int NUM_FACE_HAUT = 3;
  static const int NUM_FACE_DROITE = 2;
  int face_bas = domaine_vf.elem_faces(num_element, NUM_FACE_BAS);
  int face_gauche = domaine_vf.elem_faces(num_element, NUM_FACE_GAUCHE);
  int face_droite = domaine_vf.elem_faces(num_element, NUM_FACE_DROITE);
  int face_haut = domaine_vf.elem_faces(num_element, NUM_FACE_HAUT);
  double y_bas = domaine_vf.xv(face_bas, 1);
  double x_gauche = domaine_vf.xv(face_gauche, 0);
  double x_droite = domaine_vf.xv(face_droite, 0);
  double y_haut = domaine_vf.xv(face_haut, 1);

  // Volume de l'element
  double v_elem = (x_droite - x_gauche) * (y_haut - y_bas);
  if (bidim_axi)
    v_elem *= (x_droite + x_gauche) * 0.5 * angle_bidim_axi;

  // Volume de la partie A (avec signe)
  double v = (x0 - x1) * ((y0 + y1) * 0.5 - y_bas);
  if (bidim_axi)
    {
      double s0 = (x1 - x0) * (y0 - y_bas);
      double s1 = (x1 - x0) * (y1 - y_bas);
      if (s0 + s1 != 0)
        v *= ((x0 + x0 + x1) * s0 + (x0 + x1 + x1) * s1) / (s0 + s1);
      v *= un_tiers * angle_bidim_axi;
    }

  // Ajout du volume de B si on touche le bord haut
  double v2 = 0.;
  if (y0 > y_haut - epsilon)
    {
      v2 = (x_droite - x0) * (y_haut - y_bas);
      if (bidim_axi)
        v2 *= (x_droite + x0) * 0.5 * angle_bidim_axi;
    }
  else if (y1 > y_haut - epsilon)
    {
      v2 = (x1 - x_gauche) * (y_haut - y_bas);
      if (bidim_axi)
        v2 *= (x1 + x_gauche) * 0.5 * angle_bidim_axi;
    }
  else if (x1 > x0)
    {
      v2 = v_elem;
    }

  v = (v + v2) / v_elem;

  // On force la valeur entre 0 et 1 strictement.
  if (v < Erreur_relative_maxi_)
    v = Erreur_relative_maxi_;
  else if (v > 1. - Erreur_relative_maxi_)
    v = 1. - Erreur_relative_maxi_;

  return v;
}

/*! @brief Calcul de la matrice 2x2 de transformation pour passer d'une coordonnee dans le repere (x,y) global a une coordonnee (u,v,1-u-v) barycentrique
 *
 *  dans un element fini triangulaire.
 *  Le point 0 du triangle aura pour coordonnees (0,0,1)
 *  Le point 1                                   (1,0,0)
 *  Le point 2                                   (0,1,0)
 *
 */
// Comme volume_rectangle, cette methode calcul le barycentre
// It is identical to volume, but with multiplication by x or y when needed to obtain the barycenter...
double Parcours_interface::volume_barycentre_rectangle(const Domaine_VF& domaine_vf,
                                                       int num_element,
                                                       double x0, double y0,
                                                       double x1, double y1,
                                                       double epsilon,
                                                       double bary[3]) const
{
  const double angle_bidim_axi = bidim_axi ? Maillage_FT_Disc::angle_bidim_axi() : 0.;
  const double un_tiers = 1. / 3.;

  // Conventions TRUST VDF :
  static const int NUM_FACE_GAUCHE = 0;
  static const int NUM_FACE_BAS = 1;
  static const int NUM_FACE_HAUT = 3;
  static const int NUM_FACE_DROITE = 2;
  int face_bas = domaine_vf.elem_faces(num_element, NUM_FACE_BAS);
  int face_gauche = domaine_vf.elem_faces(num_element, NUM_FACE_GAUCHE);
  int face_droite = domaine_vf.elem_faces(num_element, NUM_FACE_DROITE);
  int face_haut = domaine_vf.elem_faces(num_element, NUM_FACE_HAUT);
  double y_bas = domaine_vf.xv(face_bas, 1);
  double x_gauche = domaine_vf.xv(face_gauche, 0);
  double x_droite = domaine_vf.xv(face_droite, 0);
  double y_haut = domaine_vf.xv(face_haut, 1);

  // Volume de l'element
  double v_elem = (x_droite - x_gauche) * (y_haut - y_bas);
  if (bidim_axi)
    v_elem *= (x_droite + x_gauche) * 0.5 * angle_bidim_axi;

  // Volume de la partie A (avec signe)
  double v = (x0 - x1) * ((y0 + y1) * 0.5 - y_bas);
  bary[0] = (x0 - x1) * ((y0 + y1) * 0.5 - y_bas);
  bary[1] = (x0 - x1) * ((y0 + y1) * 0.5 - y_bas);
  // Integral(2pi*r**2 dr) / Integral(2pi*r**2 dr) between r0 and r1
  if (bidim_axi)
    {
      double s0 = (x1 - x0) * (y0 - y_bas);
      double s1 = (x1 - x0) * (y1 - y_bas);
      if (s0 + s1 != 0)
        v *= ((x0 + x0 + x1) * s0 + (x0 + x1 + x1) * s1) / (s0 + s1);
      v *= un_tiers * angle_bidim_axi;
      double xx1 = x1*x1;
      double xx0 = x0*x0;
      if (fabs(xx1-xx0)>DMINFLOAT)
        bary[0] =  2.*un_tiers*(xx1*x1-xx0*x0)/(xx1-xx0);
      else
        bary[0] = x0;
    }
  else
    bary[0] =  (x0 + x1) *0.5;

  bary[1] = ((y0 + y1) * 0.5 + y_bas) * 0.5;

  // Ajout du volume de B si on touche le bord haut
  double v2 = 0.;
  double bary2[2] = {0.,0.};
  if (y0 > y_haut - epsilon)
    {
      v2 = (x_droite - x0) * (y_haut - y_bas);
      if (bidim_axi)
        v2 *= (x_droite + x0) * 0.5 * angle_bidim_axi;
    }
  else if (y1 > y_haut - epsilon)
    {
      v2 = (x1 - x_gauche) * (y_haut - y_bas);
      if (bidim_axi)
        v2 *= (x1 + x_gauche) * 0.5 * angle_bidim_axi;
    }
  else if (x1 > x0)
    {
      v2 = v_elem;
    }

  if (y0 > y_haut - epsilon)
    {
      if (bidim_axi)
        {
          double xx1 = x_droite*x_droite;
          double xx0 = x0*x0;
          if (fabs(xx1-xx0)>DMINFLOAT)
            bary2[0] =  2.*un_tiers*(xx1*x_droite-xx0*x0)/(xx1-xx0);
        }
      else
        bary2[0] = (x_droite + x0) * 0.5;

      bary2[1] = (y_haut + y_bas) * 0.5;
    }
  else if (y1 > y_haut - epsilon)
    {
      if (bidim_axi)
        {
          double xx1 = x1*x1;
          double xx0 = x_gauche*x_gauche;
          if (fabs(xx1-xx0)>DMINFLOAT)
            bary2[0] =  2.*un_tiers*(xx0*x_gauche-xx1*x1)/(xx0-xx1);
        }
      else
        bary2[0] = (x_gauche + x0) * 0.5;

      bary2[1] = (y_haut + y_bas) * 0.5;

    }
  else if (x1 > x0)
    {
      bary2[0] = (x_gauche + x_droite) * 0.5;
      if (bidim_axi)
        if (fabs(x_gauche-x_droite)>DMINFLOAT)
          bary2[0] = 2.*un_tiers*(x_droite*x_droite*x_droite-x_gauche*x_gauche*x_gauche)/(x_droite*x_droite-x_gauche*x_gauche);
      bary2[1] = (y_haut + y_bas) * 0.5;
    }

  if ((v+v2)>DMINFLOAT)
    for (int i=0; i<2; i++)
      bary[i] = (bary[i]*v + bary2[i]*v2) / (v+v2);

  bary[2] = 0.; // In 2D, nothing on z.
  v = (v + v2) / v_elem;

// On force la valeur entre 0 et 1 strictement.
  if (v < Erreur_relative_maxi_)
    v = Erreur_relative_maxi_;
  else if (v > 1. - Erreur_relative_maxi_)
    v = 1. - Erreur_relative_maxi_;

  return v;
}

// Description:
inline void Parcours_interface::matrice_triangle(int num_element,
                                                 FTd_vecteur2& origine,
                                                 FTd_matrice22& matrice,
                                                 double& surface) const
{
  const int s0 = (*domaine_elem_ptr)(num_element,0);
  const int s1 = (*domaine_elem_ptr)(num_element,1);
  const int s2 = (*domaine_elem_ptr)(num_element,2);
  const double x0 = (*domaine_sommets_ptr)(s0, 0);
  const double y0 = (*domaine_sommets_ptr)(s0, 1);
  origine[0] = x0;
  origine[1] = y0;
  const double dx1 = (*domaine_sommets_ptr)(s1, 0) - x0;
  const double dy1 = (*domaine_sommets_ptr)(s1, 1) - y0;
  const double dx2 = (*domaine_sommets_ptr)(s2, 0) - x0;
  const double dy2 = (*domaine_sommets_ptr)(s2, 1) - y0;

  //                         [ dx1  dx2 ](-1)   [ matrice[0][0]  matrice[0][1] ]
  // Inversion de la matrice:[          ]     = [                              ]
  //                         [ dy1  dy2 ]       [ matrice[1][0]  matrice[1][1] ]
  double det = dx1 * dy2 - dx2 * dy1;
  surface = det * 0.5;
  double inv_det = 1. / det;
  matrice[0][0] = dy2 * inv_det;
  matrice[0][1] = - dx2 * inv_det;
  matrice[1][0] = - dy1 * inv_det;
  matrice[1][1] = dx1 * inv_det;
}

/*! @brief Applique la transformation calculee par matrice_triangle a une coordonnee (x,y).
 *
 * Le resultat est stocke dans (u,v). La troisieme
 *  coordonnee barycentrique est implicitement egale a 1-u-v.
 *
 */
inline void Parcours_interface::transformation_2d(const FTd_vecteur2& origine,
                                                  const FTd_matrice22& matrice,
                                                  double x, double y,
                                                  double& u, double& v) const
{
  x -= origine[0];
  y -= origine[1];
  u = matrice[0][0] * x + matrice[0][1] * y;
  v = matrice[1][0] * x + matrice[1][1] * y;
}

/*! @brief Calcul de la contribution de volume d'une facette a la valeur de l'indicatrice dans un element.
 *
 * (x0,y0) et (x1,y1) sont les coordonnees
 *  de l'extremite du segment (obligatoirement a l'interieur du triangle
 *  ou sur un bord). plan_coupe* donne le numero de la face du triangle
 *  sur laquelle se trouve chacun des deux sommets du segment ou -1 si le
 *  sommet est strictement a l'interieur du triangle.
 *  C'est une fraction du volume de l'element comprise entre
 *  Erreur_relative_maxi et 1.-Erreur_relative_maxi
 *
 *  sommet2
 *   ^ ordonnee=v
 *   |\.
 *   | \.
 *   |  \.
 *   |   \. face 0 (opposee au sommet 0)
 *   |    \.
 *   | 0---1.
 *   | : A :\.
 *   | :   :B\.
 *    --------  -> abscisse=u  (coordonnees dans le triangle de reference : u,v)
 *  sommet0  sommet1
 *
 */
inline double Parcours_interface::volume_triangle(const Domaine_VF& domaine_vf,
                                                  int num_element,
                                                  double x0, double y0,
                                                  double x1, double y1,
                                                  int plan_coupe0,
                                                  int plan_coupe1) const
{
  assert(!bidim_axi); // Pas prevu pour l'instant

  static const int FACE_ZERO = 0; // Numero de la face opposee au sommet zero
  double origine[2];
  double matrice[2][2];
  double surface_triangle;
  double u0, v0, u1, v1;
  // Calcul des coordonnees de l'intersection dans un repere local au triangle
  matrice_triangle(num_element, origine, matrice, surface_triangle);
  transformation_2d(origine, matrice, x0, y0, u0, v0);
  transformation_2d(origine, matrice, x1, y1, u1, v1);

  // Volume de la partie A (avec signe) sur le triangle de reference
  double vol = (u0 - u1) * (v0 + v1) * 0.5;

  // Ajout du volume B si on touche la face 0
  if (plan_coupe0 == FACE_ZERO)
    vol += (1. - u0) * (1. - u0) * 0.5;
  else if (plan_coupe1 == FACE_ZERO)
    vol -= (1. - u1) * (1. - u1) * 0.5;

  if (vol < 0.)
    vol += 0.5;

  if (vol < Erreur_relative_maxi_)
    vol = Erreur_relative_maxi_;
  else if (vol > 0.5 - Erreur_relative_maxi_)
    vol = 0.5 - Erreur_relative_maxi_;

  // le triangle de reference a un volume de 0.5, on ramene vol a une fraction:
  vol = vol * 2.;
  return vol;
}

/*! @brief Cette methode permet de calculer l'intersection entre une facette et un element du maillage eulerien
 *
 * Precondition: dimension = 3
 *
 * @param (domaine_vf) domaine du calcul
 * @param (maillage) description du maillage de l'interface
 * @param (num_facette) indice de la facette intersectant
 * @param (num_element) indice de l'element intersecte
 * @return (int) 1 si ok, 0 si ??.
 */
int Parcours_interface::calcul_intersection_facelem_3D(
  const Domaine_VF& domaine_vf,
  Maillage_FT_Disc& maillage,
  int num_facette, int num_element) const
{
  assert(Objet_U::dimension == 3);
  int i,k;

  // Polygone d'intersection entre la facette et l'element.
  // On l'initialise avec la facette en entier.
  // Il contient des coordonnees barycentriques.
  static DoubleTabFT poly_(20,3);
  static DoubleTabFT new_poly_(20,3);
  // Polygone d'intersection entre la facette et l'element
  // contenant les coordonnees reelles de l'element
  static DoubleTabFT poly_reelles_(20,3);

  poly_.resize(3,3);
  poly_(0,0) = 1.;
  poly_(0,1) = 0.;
  poly_(0,2) = 0.;
  poly_(1,0) = 0.;
  poly_(1,1) = 1.;
  poly_(1,2) = 0.;
  poly_(2,0) = 0.;
  poly_(2,1) = 0.;
  poly_(2,2) = 1.;
  double coord_som[3][3];
  const int sommets[3] = { maillage.facettes_(num_facette, 0),
                           maillage.facettes_(num_facette, 1),
                           maillage.facettes_(num_facette, 2)
                         };
  {
    for (int ii = 0; ii < 3; ii++)
      {
        int s = sommets[ii];
        for (int j = 0; j < 3; j++)
          {
            coord_som[ii][j] = maillage.sommets_(s, j);
          }
      }
  }
  // Extension du triangle si deux sommets au bord
  {
    const double fact_mult = 1.;
    int isom,isom_s,isom_ss, kk;
    for (isom=0 ; isom<3 ; isom++)
      {
        isom_s = (isom+1)%3;
        if (maillage.sommet_ligne_contact(sommets[isom]) && maillage.sommet_ligne_contact(sommets[isom_s]))
          {
            isom_ss = (isom_s+1)%3;
            for (kk=0 ; kk<dimension ; kk++)
              {
                coord_som[isom][kk]   += fact_mult * (coord_som[isom][kk]  -coord_som[isom_ss][kk]);
                coord_som[isom_s][kk] += fact_mult * (coord_som[isom_s][kk]-coord_som[isom_ss][kk]);
              }
            break;//sortir au cas ou le 3e sommet soit aussi sur une ligne de contact (ce qui ne devrait normalement pas arriver...)
          }
      }
  }
  // Estimation de l'incertitude sur l'evaluation de la fonction plan
  double Zero = - Erreur_max_coordonnees_;
  if (correction_parcours_thomas_)
    Zero = 0.;

  // polygone_plan_coupe contient pour chaque segment du polygone
  // le numero du plan de l'element qui a genere ce segment.
  // polygone_plan_coupe_(1) est le plan qui genere le segment poly(0) - poly(1)
  // Il contient -1 si c'est un segment de la facette
  // d'origine.
  static ArrOfIntFT polygone_plan_coupe_(20);
  static ArrOfIntFT new_poly_plan_coupe_(20);
  polygone_plan_coupe_.resize_array(3);
  polygone_plan_coupe_ = -1;

  // *********************************************************************
  // *********************************************************************
  // Calcul de l'intersection du polygone avec les plans qui
  // definissent l'element
  for (int num_plan = 0; num_plan < nb_faces_elem_; num_plan++)
    {
      new_poly_.resize(0,3);
      new_poly_plan_coupe_.resize_array(0);
      int n = 0;
      // Coefficients de la fonction qui definit le plan
      // "fonction plan" (x,y,z) = a*x + b*y + c*d + d
      double a, b, c, d;
      double signe = calcul_eq_plan(domaine_vf, num_element, num_plan, a, b, c, d);
      // Calcul de la "fonction plan" aux sommets de la facette
      const double f0 = (a * coord_som[0][0] + b * coord_som[0][1] + c * coord_som[0][2] + d) * signe;
      const double f1 = (a * coord_som[1][0] + b * coord_som[1][1] + c * coord_som[1][2] + d) * signe;
      const double f2 = (a * coord_som[2][0] + b * coord_som[2][1] + c * coord_som[2][2] + d) * signe;

      // Recherche des points du polygone qui sont dans le demi-espace
      const int nb_sommets_poly = poly_.dimension(0);
      i = nb_sommets_poly - 1; // Le dernier sommet du polygone
      double u = poly_(i,0);
      double v = poly_(i,1);
      double w = poly_(i,2);
      // Calcul de la "fonction plan" au sommet du polygone
      double f = f0 * u + f1 * v + f2 * w;
      int sommet_dehors = inf_strict(f,Zero) ? 1 : 0;

      for (i = 0; i < nb_sommets_poly; i++)
        {
          int sommet_precedent_dehors = sommet_dehors;
          double u_prec = u;
          double v_prec = v;
          double f_prec = f;
          u = poly_(i,0);
          v = poly_(i,1);
          w = poly_(i,2);
          // Calcul de la "fonction plan" au sommet du polygone
          f = f0 * u + f1 * v + f2 * w;
          sommet_dehors = inf_strict(f,Zero) ? 1 : 0;
          if (sommet_dehors + sommet_precedent_dehors == 1)
            {
              // Le dernier segment du polygone coupe le plan.
              // On calcule le point d'intersection
              // Coordonnee barycentrique de l'intersection sur le segment du polygone
              double t = f / (f - f_prec);
              if (t < 0.)
                t = 0.;
              else if (t > 1.)
                t = 1.;
              new_poly_.resize(n+1,3);
              // Coordonnee barycentrique de l'intersection sur la facette
              double new_u = u_prec * t + u * (1.-t);
              double new_v = v_prec * t + v * (1.-t);
              double new_w = 1. - new_u - new_v;
              new_poly_(n,0) = new_u;
              new_poly_(n,1) = new_v;
              new_poly_(n,2) = new_w;
              n++;
              const int num_plan_coupe = sommet_dehors ? polygone_plan_coupe_[i] : num_plan;
              new_poly_plan_coupe_.append_array(num_plan_coupe);
            }
          if (! sommet_dehors)
            {
              // Le dernier sommet est dedans, je le garde
              new_poly_.resize(n+1,3);
              new_poly_(n,0) = u;
              new_poly_(n,1) = v;
              new_poly_(n,2) = w;
              n++;
              const int num_plan_coupe = polygone_plan_coupe_[i];
              new_poly_plan_coupe_.append_array(num_plan_coupe);
            }
        }
      if (new_poly_.dimension(0) == 0)
        {
          // L'intersection est vide: c'est anormal
          compteur_erreur_grossiere++;
          return -1;
        }
      poly_ = new_poly_;
      polygone_plan_coupe_ = new_poly_plan_coupe_;
    }
  // *********************************************************************
  const int nb_sommets_poly = poly_.dimension(0);
  // Calcul du centre de gravite et de la surface de l'intersection
  double volume = 0.;
  double surface = 0.;
  double u_centre = 0.;
  double v_centre = 0.;
  {
    i = nb_sommets_poly - 1; // Le dernier sommet du polygone
    double u = poly_(i,0);
    double v = poly_(i,1);
    for (i = 0; i < nb_sommets_poly; i++)
      {
        double u_prec = u;
        double v_prec = v;
        u = poly_(i,0);
        v = poly_(i,1);
        //calcul de la contribution de l'arete a l'aire :
        // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe u=0)
        double contrib_surface = (v - v_prec) * (u + u_prec) * 0.5;
        if (contrib_surface!=0.)
          {
            //calcul du centre de gravite X du trapeze :
            //    v    > |------------\.
            //           |            |\.
            //           |            | \.
            //    vX   > |       X    |  \.
            //           |            |   \.
            //H v_prec > |------------|----\.
            //           ^       ^    ^    ^
            //           0       uX   u    u_prec
            double B = std::max(u,u_prec);    //grande base du trapeze
            double b = std::min(u,u_prec);    //petite base du trapeze
            double h = std::fabs(v - v_prec); //hauteur du trapeze
            //le CdG du trapeze : le CdG des CdG du rectangle et du triangle
            double Sr = b * h;
            double CdGrX = b/2.;
            double CdGrY = (v+v_prec)/2.;
            double St = (B-b) * h/2.;
            double CdGtX = (2.*b+B)/3.;
            double CdGtY = v+v_prec;
            if (u<u_prec)
              {
                CdGtY += v_prec;
              }
            else
              {
                CdGtY += v;
              }
            CdGtY /= 3.;

            double S = Sr + St;
            double uX = (Sr * CdGrX + St * CdGtX) / S;
            double vX = (Sr * CdGrY + St * CdGtY) / S;

            double contrib_u_centre = uX * contrib_surface;
            double contrib_v_centre = vX * contrib_surface;

            surface += contrib_surface;
            u_centre += contrib_u_centre;
            v_centre += contrib_v_centre;
          }
      }
  }
  // En cas d'erreur d'arrondi ...
  if (surface < 0.) surface = 0.;

  // calcul des contributions.
  double barycentre_phase1[3] = {0.,0.,0.};
  if (correction_parcours_thomas_ || (sqrt(surface) > 5. * Erreur_max_coordonnees_))
    {
      //normalisation du centre de gravite
      if (sup_strict(surface,0.))
        {
          u_centre /= surface;
          v_centre /= surface;
        }
      else
        {
          u_centre = v_centre = 1./3.;
        }
      double w_centre = 1. - u_centre - v_centre;
      //calcul des coordonnees reelles du centre de gravite
      FTd_vecteur3 centre_de_gravite;
      {
        for (int ii = 0; ii < 3; ii++)
          {
            centre_de_gravite[ii] =
              u_centre * coord_som[0][ii]
              + v_centre * coord_som[1][ii]
              + w_centre * coord_som[2][ii];
          }
      }
      //calcul des coordonnees reelles du polygone d'intersection
      poly_reelles_.resize(nb_sommets_poly,3);
      for (i=0 ; i<nb_sommets_poly ; i++)
        {
          for (k=0 ; k<dimension ; k++)
            {
              poly_reelles_(i,k) =
                poly_(i,0) * coord_som[0][k]
                + poly_(i,1) * coord_som[1][k]
                + poly_(i,2) * coord_som[2][k];
            }
        }
      FTd_vecteur3 norme;
      //calcul des contributions de la facette pour l'indicatrice,
      //en fonction du type d'element eulerien
      switch(type_element_)
        {
        case TETRA:
          if (flag_warning_code_missing)
            {
              Cerr << "WARNING : barycentre_phase1 not filled properly!!" << finl;
              Cerr << "WARNING : Calculation of barycentre_phase1 not implemented yet for TETRA." << finl;
              flag_warning_code_missing=0;
            }
          volume = volume_tetraedre(domaine_vf, num_element,num_facette,
                                    maillage,
                                    poly_reelles_,
                                    centre_de_gravite,
                                    Erreur_max_coordonnees_);
          break;
        case HEXA:
          if (flag_warning_code_missing)
            {
              Cerr << "WARNING : barycentre_phase1 not filled properly!!" << finl;
              Cerr << "WARNING : Calculation of barycentre_phase1 not implemented yet for HEXA." << finl;
              flag_warning_code_missing=0;
            }
          maillage.calcul_normale_3D(num_facette,norme);
          volume = volume_hexaedre(domaine_vf, num_element,
                                   poly_reelles_,
                                   norme, centre_de_gravite,
                                   polygone_plan_coupe_,
                                   Erreur_max_coordonnees_);
          break;
        default:
          // qu'est-ce qu'on fout la ?
          Process::exit();
          break;
        }
      // "surface" contient la surface du polygone d'intersection dans le triangle
      // de reference. Or ce triangle est de surface 1/2. Pour avoir la fraction
      // de surface par rapport au triangle, on multiplie par 2.
      maillage.intersections_elem_facettes_.ajoute_intersection(num_facette,
                                                                num_element,
                                                                surface * 2.,
                                                                volume,
                                                                barycentre_phase1,
                                                                u_centre,
                                                                v_centre,
                                                                1. - u_centre - v_centre);
    }
  else
    {
      // Intersection de surface nulle. On l'enregistre pour ne pas traiter
      // cette facette a nouveau lors du parcours.
      maillage.intersections_elem_facettes_.ajoute_intersection(num_facette,
                                                                num_element,
                                                                0.,
                                                                0.,
                                                                barycentre_phase1,
                                                                0., 0., 0.);
    }

  int code_retour = 0;
  {
    for (i = 0; i < nb_sommets_poly; i++)
      {
        const int plan_coupe = polygone_plan_coupe_[i];
        if (plan_coupe >= 0)
          code_retour |= 1 << plan_coupe;
      }
  }
  return code_retour;
}

// debut EB
// Description:
//    Cette methode permet de calculer l'intersection entre une facette et une face (volume de controle) du maillage eulerien
// Precondition: dimension = 3
// Parametre: domaine_vf
//    Signification: zone du calcul
//    Valeurs par defaut: NA
//    Contraintes: NA
//    Acces: lecture
// Parametre: maillage
//    Signification: description du maillage de l'interface
//    Valeurs par defaut: NA
//    Contraintes: NA
//    Acces: lecture/ecriture
// Parametre: num_facette
//    Signification: indice de la facette intersectant
//    Valeurs par defaut: NA
//    Contraintes: NA
//    Acces: lecture
// Parametre: num_face
//    Signification: indice de la face intersecte
//    Valeurs par defaut: NA
//    Contraintes: NA
//    Acces: lecture
/// Retour: int
//    Signification: 1 si ok, 0 si ??.
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
int Parcours_interface::calcul_intersection_facette_face_3D(const Domaine_VF& domaine_vf,
                                                            Maillage_FT_Disc& maillage,
                                                            int num_facette,
                                                            int num_face, Intersections_Face_Facettes& intersections_face) const
{
  assert(Objet_U::dimension == 3);
  int i,k;

  // Polygone d'intersection entre la facette et la face
  // On l'initialise avec la facette en entier.
  // Il contient des coordonnees barycentriques.
  static DoubleTabFT poly_(20,3);
  static DoubleTabFT new_poly_(20,3);
  // Polygone d'intersection entre la facette et la face
  // contenant les coordonnees reelles de la face
  static DoubleTabFT poly_reelles_(20,3);

  poly_.resize(3,3);
  poly_(0,0) = 1.;
  poly_(0,1) = 0.;
  poly_(0,2) = 0.;
  poly_(1,0) = 0.;
  poly_(1,1) = 1.;
  poly_(1,2) = 0.;
  poly_(2,0) = 0.;
  poly_(2,1) = 0.;
  poly_(2,2) = 1.;
  double coord_som[3][3];
  const int sommets[3] = { maillage.facettes_(num_facette, 0),
                           maillage.facettes_(num_facette, 1),
                           maillage.facettes_(num_facette, 2)
                         };
  {
    for (int ii = 0; ii < 3; ii++)
      {
        int s = sommets[ii];
        for (int j = 0; j < 3; j++)
          {
            coord_som[ii][j] = maillage.sommets_(s, j);
          }
      }
  }
  // Extension du triangle si deux sommets au bord
  {
    const double fact_mult = 1.;
    int isom,isom_s,isom_ss, kk;
    for (isom=0 ; isom<3 ; isom++)
      {
        isom_s = (isom+1)%3;
        if (maillage.sommet_ligne_contact(sommets[isom]) && maillage.sommet_ligne_contact(sommets[isom_s]))
          {
            isom_ss = (isom_s+1)%3;
            for (kk=0 ; kk<dimension ; kk++)
              {
                coord_som[isom][kk]   += fact_mult * (coord_som[isom][kk]  -coord_som[isom_ss][kk]);
                coord_som[isom_s][kk] += fact_mult * (coord_som[isom_s][kk]-coord_som[isom_ss][kk]);
              }
            break;//sortir au cas ou le 3e sommet soit aussi sur une ligne de contact (ce qui ne devrait normalement pas arriver...)
          }
      }
  }

  // Estimation de l'incertitude sur l'evaluation de la fonction plan
  double Zero = - Erreur_max_coordonnees_;
  if (correction_parcours_thomas_)
    Zero = 0.;

  // polygone_plan_coupe contient pour chaque segment du polygone
  // le numero du plan de la face qui a genere ce segment.
  // polygone_plan_coupe_(1) est le plan qui genere le segment poly(0) - poly(1)
  // Il contient -1 si c'est un segment de la facette
  // d'origine.
  static ArrOfIntFT polygone_plan_coupe_(20);
  static ArrOfIntFT new_poly_plan_coupe_(20);
  polygone_plan_coupe_.resize_array(3);
  polygone_plan_coupe_ = -1;

  // *********************************************************************
  // *********************************************************************
  // Calcul de l'intersection du polygone avec les plans qui
  // definissent la face
  Domaine_dis& domaine_dis = maillage.refdomaine_dis_.valeur();
  Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF,domaine_dis.valeur());
  for (int num_plan=0; num_plan<nb_faces_elem_; num_plan++) // tout comme les elements, les faces euleriennes ont 6 faces decrivant leur volume de controle
    {
      new_poly_.resize(0,3);
      new_poly_plan_coupe_.resize_array(0);
      int n = 0;
      // Coefficients de la fonction qui definit le plan
      double a,b,c,d;
      double signe = calcul_eq_plan_face(domaine_vf, num_face, num_plan, a, b, c, d);

      // EB

      // Calcul de la "fonction plan" aux sommets de la facette
      const double f0 = (a * coord_som[0][0] + b * coord_som[0][1] + c * coord_som[0][2] + d) * signe;
      const double f1 = (a * coord_som[1][0] + b * coord_som[1][1] + c * coord_som[1][2] + d) * signe;
      const double f2 = (a * coord_som[2][0] + b * coord_som[2][1] + c * coord_som[2][2] + d) * signe;

      // Recherche des points du polygone qui sont dans le demi-espace
      const int nb_sommets_poly = poly_.dimension(0);
      i = nb_sommets_poly - 1; // Le dernier sommet du polygone
      double u = poly_(i,0);
      double v = poly_(i,1);
      double w = poly_(i,2);
      // Calcul de la "fonction plan" au sommet du polygone
      double f = f0 * u + f1 * v + f2 * w;
      int sommet_dehors = inf_strict(f,Zero) ? 1 : 0;

      for (i = 0; i < nb_sommets_poly; i++)
        {
          int sommet_precedent_dehors = sommet_dehors;
          double u_prec = u;
          double v_prec = v;
          double f_prec = f;
          u = poly_(i,0);
          v = poly_(i,1);
          w = poly_(i,2);
          // Calcul de la "fonction plan" au sommet du polygone
          f = f0 * u + f1 * v + f2 * w;
          sommet_dehors = inf_strict(f,Zero) ? 1 : 0;
          if (sommet_dehors + sommet_precedent_dehors == 1)
            {
              // Le dernier segment du polygone coupe le plan.
              // On calcule le point d'intersection
              // Coordonnee barycentrique de l'intersection sur le segment du polygone
              double t = f / (f - f_prec);
              if (t < 0.)
                t = 0.;
              else if (t > 1.)
                t = 1.;
              new_poly_.resize(n+1,3);
              // Coordonnee barycentrique de l'intersection sur la facette
              double new_u = u_prec * t + u * (1.-t);
              double new_v = v_prec * t + v * (1.-t);
              double new_w = 1. - new_u - new_v;
              new_poly_(n,0) = new_u;
              new_poly_(n,1) = new_v;
              new_poly_(n,2) = new_w;
              n++;
              const int num_plan_coupe = sommet_dehors ? polygone_plan_coupe_[i] : num_plan;
              new_poly_plan_coupe_.append_array(num_plan_coupe);
            }
          if (! sommet_dehors)
            {
              // Le dernier sommet est dedans, je le garde
              new_poly_.resize(n+1,3);
              new_poly_(n,0) = u;
              new_poly_(n,1) = v;
              new_poly_(n,2) = w;
              n++;
              const int num_plan_coupe = polygone_plan_coupe_[i];
              new_poly_plan_coupe_.append_array(num_plan_coupe);
            }
        }
      if (new_poly_.dimension(0) == 0)
        {
          // L'intersection est vide: c'est anormal
          compteur_erreur_grossiere++;
          return -1;
        }
      poly_ = new_poly_;
      polygone_plan_coupe_ = new_poly_plan_coupe_;
    }
  // *********************************************************************
  const int nb_sommets_poly = poly_.dimension(0);
  // Calcul du centre de gravite et de la surface de l'intersection
  double volume = 0.;
  double surface = 0.;
  double u_centre = 0.;
  double v_centre = 0.;
  {
    i = nb_sommets_poly - 1; // Le dernier sommet du polygone
    double u = poly_(i,0);
    double v = poly_(i,1);
    for (i = 0; i < nb_sommets_poly; i++)
      {
        double u_prec = u;
        double v_prec = v;
        u = poly_(i,0);
        v = poly_(i,1);
        //calcul de la contribution de l'arete a l'aire :
        // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe u=0)
        double contrib_surface = (v - v_prec) * (u + u_prec) * 0.5;
        if (contrib_surface!=0.)
          {
            //calcul du centre de gravite X du trapeze :
            //    v    > |------------\.
            //           |            |\.
            //           |            | \.
            //    vX   > |       X    |  \.
            //           |            |   \.
            //H v_prec > |------------|----\.
            //           ^       ^    ^    ^
            //           0       uX   u    u_prec
            double B = std::max(u,u_prec);    //grande base du trapeze
            double b = std::min(u,u_prec);    //petite base du trapeze
            double h = std::abs(v - v_prec); //hauteur du trapeze
            //le CdG du trapeze : le CdG des CdG du rectangle et du triangle
            double Sr = b * h;
            double CdGrX = b/2.;
            double CdGrY = (v+v_prec)/2.;
            double St = (B-b) * h/2.;
            double CdGtX = (2.*b+B)/3.;
            double CdGtY = v+v_prec;
            if (u<u_prec)
              {
                CdGtY += v_prec;
              }
            else
              {
                CdGtY += v;
              }
            CdGtY /= 3.;

            double S = Sr + St;
            double uX = (Sr * CdGrX + St * CdGtX) / S;
            double vX = (Sr * CdGrY + St * CdGtY) / S;

            double contrib_u_centre = uX * contrib_surface;
            double contrib_v_centre = vX * contrib_surface;

            surface += contrib_surface;
            u_centre += contrib_u_centre;
            v_centre += contrib_v_centre;
          }
      }
  }
  // En cas d'erreur d'arrondi ...
  if (surface < 0.) surface = 0.;
  // calcul des contributions.
  //  Modif BM 9/9/2004 : calcul d'une surface reelle et non d'une fraction!
  double normale[3];
  const double surface_facette = maillage.calcul_normale_3D(num_facette, normale);
  if (correction_parcours_thomas_ || (surface * surface_facette > 5. * Erreur_max_coordonnees_))
    {
      //normalisation du centre de gravite
      if (sup_strict(surface,0.))
        {
          u_centre /= surface;
          v_centre /= surface;
        }
      else
        {
          u_centre = v_centre = 1./3.;
        }
      double w_centre = 1. - u_centre - v_centre;
      //calcul des coordonnees reelles du centre de gravite
      FTd_vecteur3 centre_de_gravite;
      {
        for (int ii = 0; ii < 3; ii++)
          {
            centre_de_gravite[ii] =
              u_centre * coord_som[0][ii]
              + v_centre * coord_som[1][ii]
              + w_centre * coord_som[2][ii];
          }
      }
      //calcul des coordonnees reelles du polygone d'intersection
      poly_reelles_.resize(nb_sommets_poly,3);
      for (i=0 ; i<nb_sommets_poly ; i++)
        {
          for (k=0 ; k<dimension ; k++)
            {
              poly_reelles_(i,k) =
                poly_(i,0) * coord_som[0][k]
                + poly_(i,1) * coord_som[1][k]
                + poly_(i,2) * coord_som[2][k];
            }
        }
      FTd_vecteur3 norme;
      //calcul des contributions de la facette pour l'indicatrice,
      //en fonction du type d'element eulerien
      switch(type_element_)
        {
        case TETRA:
          Cerr << "Calcul de l'indicatrice aux faces non code pour les elements de types TETRA !" << finl;
          exit();
          break;
        case HEXA:
          maillage.calcul_normale_3D(num_facette,norme);
          volume = volume_hexaedre_face(domaine_vdf, num_face,
                                        poly_reelles_,
                                        norme, centre_de_gravite,
                                        polygone_plan_coupe_,
                                        Erreur_max_coordonnees_);
          break;
        default:
          // qu'est-ce qu'on fout la ?
          exit();
          break;
        }
      // "surface" contient la surface du polygone d'intersection dans le triangle
      // de reference. Or ce triangle est de surface 1/2. Pour avoir la fraction
      // de surface par rapport au triangle, on multiplie par 2.
      intersections_face.ajoute_intersection(num_facette,
                                             num_face,
                                             surface * 2.,
                                             volume,
                                             u_centre,
                                             v_centre,
                                             1. - u_centre - v_centre);
    }
  else
    {
      // Intersection de surface nulle. On l'enregistre pour ne pas traiter
      // cette facette a nouveau lors du parcours.
      intersections_face.ajoute_intersection(num_facette,
                                             num_face,
                                             0.,
                                             0.,
                                             0., 0., 0.);
    }

  int code_retour = 0;
  {
    for (i = 0; i < nb_sommets_poly; i++)
      {
        const int plan_coupe = polygone_plan_coupe_[i];
        if (plan_coupe >= 0)
          code_retour |= 1 << plan_coupe;
      }
  }

  return code_retour;
}

int Parcours_interface::calcul_intersection_facette_arete_3D(const Domaine_VF& domaine_vf,
                                                             Maillage_FT_Disc& maillage,
                                                             int num_facette,
                                                             int num_arete, Intersections_Arete_Facettes& intersections_arete) const
{
  assert(Objet_U::dimension == 3);
  int i,k;

  // Polygone d'intersection entre la facette et la face
  // On l'initialise avec la facette en entier.
  // Il contient des coordonnees barycentriques.
  static DoubleTabFT poly_(20,3);
  static DoubleTabFT new_poly_(20,3);
  // Polygone d'intersection entre la facette et la face
  // contenant les coordonnees reelles de la face
  static DoubleTabFT poly_reelles_(20,3);

  poly_.resize(3,3);
  poly_(0,0) = 1.;
  poly_(0,1) = 0.;
  poly_(0,2) = 0.;
  poly_(1,0) = 0.;
  poly_(1,1) = 1.;
  poly_(1,2) = 0.;
  poly_(2,0) = 0.;
  poly_(2,1) = 0.;
  poly_(2,2) = 1.;
  double coord_som[3][3];
  const int sommets[3] = { maillage.facettes_(num_facette, 0),
                           maillage.facettes_(num_facette, 1),
                           maillage.facettes_(num_facette, 2)
                         };
  {
    for (int ii = 0; ii < 3; ii++)
      {
        int s = sommets[ii];
        for (int j = 0; j < 3; j++)
          {
            coord_som[ii][j] = maillage.sommets_(s, j);
          }
      }
  }
  // Extension du triangle si deux sommets au bord
  {
    const double fact_mult = 1.;
    int isom,isom_s,isom_ss, kk;
    for (isom=0 ; isom<3 ; isom++)
      {
        isom_s = (isom+1)%3;
        if (maillage.sommet_ligne_contact(sommets[isom]) && maillage.sommet_ligne_contact(sommets[isom_s]))
          {
            isom_ss = (isom_s+1)%3;
            for (kk=0 ; kk<dimension ; kk++)
              {
                coord_som[isom][kk]   += fact_mult * (coord_som[isom][kk]  -coord_som[isom_ss][kk]);
                coord_som[isom_s][kk] += fact_mult * (coord_som[isom_s][kk]-coord_som[isom_ss][kk]);
              }
            break;//sortir au cas ou le 3e sommet soit aussi sur une ligne de contact (ce qui ne devrait normalement pas arriver...)
          }
      }
  }

  // Estimation de l'incertitude sur l'evaluation de la fonction plan
  double Zero = - Erreur_max_coordonnees_;
  if (correction_parcours_thomas_)
    Zero = 0.;

  // polygone_plan_coupe contient pour chaque segment du polygone
  // le numero du plan de la face qui a genere ce segment.
  // polygone_plan_coupe_(1) est le plan qui genere le segment poly(0) - poly(1)
  // Il contient -1 si c'est un segment de la facette
  // d'origine.
  static ArrOfIntFT polygone_plan_coupe_(20);
  static ArrOfIntFT new_poly_plan_coupe_(20);
  polygone_plan_coupe_.resize_array(3);
  polygone_plan_coupe_ = -1;

  // *********************************************************************
  // *********************************************************************
  // Calcul de l'intersection du polygone avec les plans qui
  // definissent la face
  Domaine_dis& domaine_dis = maillage.refdomaine_dis_.valeur();
  Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF,domaine_dis.valeur());

  //Cerr << "arete " << num_arete << " cg " << domaine_vf.xa(num_arete,0) << " " <<  domaine_vf.xa(num_arete,1) << " " << domaine_vf.xa(num_arete,2) << finl;
  for (int num_plan=0; num_plan<nb_faces_elem_; num_plan++) // tout comme les elements, les faces euleriennes ont 6 faces decrivant leur volume de controle
    {
      new_poly_.resize(0,3);
      new_poly_plan_coupe_.resize_array(0);
      int n = 0;
      // Coefficients de la fonction qui definit le plan
      double a,b,c,d;
      double signe = calcul_eq_aretes(domaine_vdf, num_arete, num_plan, a, b, c, d);

      // if (domaine_vf.xa(num_arete,0)<-2.8e-5 && domaine_vf.xa(num_arete,0)>-3e-5 &&  domaine_vf.xa(num_arete,2)==0 && domaine_vf.xa(num_arete,1)<1e-5 && domaine_vf.xa(num_arete,1)>-1e-5 ) Cerr << "num_plan "
      //       << num_plan << " a " << a << " b " << b << " c " << c << " d " << d << finl;

      // EB

      // Calcul de la "fonction plan" aux sommets de la facette
      const double f0 = (a * coord_som[0][0] + b * coord_som[0][1] + c * coord_som[0][2] + d) * signe;
      const double f1 = (a * coord_som[1][0] + b * coord_som[1][1] + c * coord_som[1][2] + d) * signe;
      const double f2 = (a * coord_som[2][0] + b * coord_som[2][1] + c * coord_som[2][2] + d) * signe;

      // Recherche des points du polygone qui sont dans le demi-espace
      const int nb_sommets_poly = poly_.dimension(0);
      i = nb_sommets_poly - 1; // Le dernier sommet du polygone
      double u = poly_(i,0);
      double v = poly_(i,1);
      double w = poly_(i,2);
      // Calcul de la "fonction plan" au sommet du polygone
      double f = f0 * u + f1 * v + f2 * w;
      int sommet_dehors = inf_strict(f,Zero) ? 1 : 0;

      for (i = 0; i < nb_sommets_poly; i++)
        {
          int sommet_precedent_dehors = sommet_dehors;
          double u_prec = u;
          double v_prec = v;
          double f_prec = f;
          u = poly_(i,0);
          v = poly_(i,1);
          w = poly_(i,2);
          // Calcul de la "fonction plan" au sommet du polygone
          f = f0 * u + f1 * v + f2 * w;
          sommet_dehors = inf_strict(f,Zero) ? 1 : 0;
          if (sommet_dehors + sommet_precedent_dehors == 1)
            {
              // Le dernier segment du polygone coupe le plan.
              // On calcule le point d'intersection
              // Coordonnee barycentrique de l'intersection sur le segment du polygone
              double t = f / (f - f_prec);
              if (t < 0.)
                t = 0.;
              else if (t > 1.)
                t = 1.;
              new_poly_.resize(n+1,3);
              // Coordonnee barycentrique de l'intersection sur la facette
              double new_u = u_prec * t + u * (1.-t);
              double new_v = v_prec * t + v * (1.-t);
              double new_w = 1. - new_u - new_v;
              new_poly_(n,0) = new_u;
              new_poly_(n,1) = new_v;
              new_poly_(n,2) = new_w;
              n++;
              const int num_plan_coupe = sommet_dehors ? polygone_plan_coupe_[i] : num_plan;
              new_poly_plan_coupe_.append_array(num_plan_coupe);
            }
          if (! sommet_dehors)
            {
              // Le dernier sommet est dedans, je le garde
              new_poly_.resize(n+1,3);
              new_poly_(n,0) = u;
              new_poly_(n,1) = v;
              new_poly_(n,2) = w;
              n++;
              const int num_plan_coupe = polygone_plan_coupe_[i];
              new_poly_plan_coupe_.append_array(num_plan_coupe);
            }
        }
      if (new_poly_.dimension(0) == 0)
        {
          // L'intersection est vide: c'est anormal
          compteur_erreur_grossiere++;
          return -1;
        }
      poly_ = new_poly_;
      polygone_plan_coupe_ = new_poly_plan_coupe_;
    }
  // *********************************************************************
  const int nb_sommets_poly = poly_.dimension(0);
  // Calcul du centre de gravite et de la surface de l'intersection
  double volume = 0.;
  double surface = 0.;
  double u_centre = 0.;
  double v_centre = 0.;
  {
    i = nb_sommets_poly - 1; // Le dernier sommet du polygone
    double u = poly_(i,0);
    double v = poly_(i,1);
    for (i = 0; i < nb_sommets_poly; i++)
      {
        double u_prec = u;
        double v_prec = v;
        u = poly_(i,0);
        v = poly_(i,1);
        //calcul de la contribution de l'arete a l'aire :
        // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe u=0)
        double contrib_surface = (v - v_prec) * (u + u_prec) * 0.5;
        if (contrib_surface!=0.)
          {
            //calcul du centre de gravite X du trapeze :
            //    v    > |------------\.
            //           |            |\.
            //           |            | \.
            //    vX   > |       X    |  \.
            //           |            |   \.
            //H v_prec > |------------|----\.
            //           ^       ^    ^    ^
            //           0       uX   u    u_prec
            double B = std::max(u,u_prec);    //grande base du trapeze
            double b = std::min(u,u_prec);    //petite base du trapeze
            double h = std::abs(v - v_prec); //hauteur du trapeze
            //le CdG du trapeze : le CdG des CdG du rectangle et du triangle
            double Sr = b * h;
            double CdGrX = b/2.;
            double CdGrY = (v+v_prec)/2.;
            double St = (B-b) * h/2.;
            double CdGtX = (2.*b+B)/3.;
            double CdGtY = v+v_prec;
            if (u<u_prec)
              {
                CdGtY += v_prec;
              }
            else
              {
                CdGtY += v;
              }
            CdGtY /= 3.;

            double S = Sr + St;
            double uX = (Sr * CdGrX + St * CdGtX) / S;
            double vX = (Sr * CdGrY + St * CdGtY) / S;

            double contrib_u_centre = uX * contrib_surface;
            double contrib_v_centre = vX * contrib_surface;

            surface += contrib_surface;
            u_centre += contrib_u_centre;
            v_centre += contrib_v_centre;
          }
      }
  }
  // En cas d'erreur d'arrondi ...
  if (surface < 0.) surface = 0.;
  // calcul des contributions.
  //  Modif BM 9/9/2004 : calcul d'une surface reelle et non d'une fraction!
  double normale[3];
  const double surface_facette = maillage.calcul_normale_3D(num_facette, normale);
  if (correction_parcours_thomas_ || (surface * surface_facette > 5. * Erreur_max_coordonnees_))
    {
      //normalisation du centre de gravite
      if (sup_strict(surface,0.))
        {
          u_centre /= surface;
          v_centre /= surface;
        }
      else
        {
          u_centre = v_centre = 1./3.;
        }
      double w_centre = 1. - u_centre - v_centre;
      //calcul des coordonnees reelles du centre de gravite
      FTd_vecteur3 centre_de_gravite;
      {
        for (int ii = 0; ii < 3; ii++)
          {
            centre_de_gravite[ii] =
              u_centre * coord_som[0][ii]
              + v_centre * coord_som[1][ii]
              + w_centre * coord_som[2][ii];
          }
      }
      //calcul des coordonnees reelles du polygone d'intersection
      poly_reelles_.resize(nb_sommets_poly,3);
      for (i=0 ; i<nb_sommets_poly ; i++)
        {
          for (k=0 ; k<dimension ; k++)
            {
              poly_reelles_(i,k) =
                poly_(i,0) * coord_som[0][k]
                + poly_(i,1) * coord_som[1][k]
                + poly_(i,2) * coord_som[2][k];
            }
        }
      FTd_vecteur3 norme;
      //calcul des contributions de la facette pour l'indicatrice,
      //en fonction du type d'element eulerien
      switch(type_element_)
        {
        case TETRA:
          Cerr << "Calcul de l'indicatrice aux faces non code pour les elements de types TETRA !" << finl;
          exit();
          break;
        case HEXA:
          maillage.calcul_normale_3D(num_facette,norme);
          volume = volume_hexaedre_arete_interne(domaine_vdf, num_arete,
                                                 poly_reelles_,
                                                 norme, centre_de_gravite,
                                                 polygone_plan_coupe_,
                                                 Erreur_max_coordonnees_);
          break;
        default:
          // qu'est-ce qu'on fout la ?
          exit();
          break;
        }
      // "surface" contient la surface du polygone d'intersection dans le triangle
      // de reference. Or ce triangle est de surface 1/2. Pour avoir la fraction
      // de surface par rapport au triangle, on multiplie par 2.
      intersections_arete.ajoute_intersection(num_facette,
                                              num_arete,
                                              surface * 2.,
                                              volume,
                                              u_centre,
                                              v_centre,
                                              1. - u_centre - v_centre);
    }
  else
    {
      // Intersection de surface nulle. On l'enregistre pour ne pas traiter
      // cette facette a nouveau lors du parcours.
      intersections_arete.ajoute_intersection(num_facette,
                                              num_arete,
                                              0.,
                                              0.,
                                              0., 0., 0.);
    }

  int code_retour = 0;
  {
    for (i = 0; i < nb_sommets_poly; i++)
      {
        const int plan_coupe = polygone_plan_coupe_[i];
        if (plan_coupe >= 0)
          code_retour |= 1 << plan_coupe;
      }
  }

  return code_retour;
}
// fin EB

// Description:
//   Cette methode (statique) permet d'inverser une matrice 3x3
// Precondition: des denominateurs non nuls
// Parametre: matrice
//    Signification: matrice 3x3 a inverser
//    Valeurs par defaut: NA
//    Contraintes: NA
//    Acces: lecture
// Parametre: matrice_inv
//    Signification: inverse de la matrice 3x3
//    Valeurs par defaut: NA
//    Contraintes: NA
//    Acces: ecriture
/// Retour: rien
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
/*! @brief Cette methode (statique) permet d'inverser une matrice 3x3
 *
 * Precondition: des denominateurs non nuls
 *
 * @param (matrice) matrice 3x3 a inverser
 * @param (matrice_inv) inverse de la matrice 3x3
 * @return (rien)
 */
void Parcours_interface::calcul_inverse_matrice33(const FTd_matrice33& matrice, FTd_matrice33& matrice_inv)
{
  const double a00 = matrice[0][0];
  const double a01 = matrice[0][1];
  const double a02 = matrice[0][2];
  const double a10 = matrice[1][0];
  const double a11 = matrice[1][1];
  const double a12 = matrice[1][2];
  const double a20 = matrice[2][0];
  const double a21 = matrice[2][1];
  const double a22 = matrice[2][2];
  //calcul de valeurs temporaires pour optimisation
  const double t4 = a00*a11;
  const double t6 = a00*a12;
  const double t8 = a01*a10;
  const double t10 = a02*a10;
  const double t12 = a01*a20;
  const double t14 = a02*a20;
  const double t = t4*a22-t6*a21-t8*a22+t10*a21+t12*a12-t14*a11;
  if (t==0.)
    {
      Cerr<<"PB : matrice non inversible : t= "<<t<<finl<<"  matrice= "<<finl;
      int i,j;
      for (i=0 ; i<3 ; i++)
        {
          for (j=0 ; j<3 ; j++)
            {
              Cerr<<"  "<<matrice[i][j];
            }
          Cerr<<finl;
        }
      exit();
    }
  const double t17 = 1/(t);

  //calcul de la matrice inverse
  matrice_inv[0][0] = (a11*a22-a12*a21)*t17;
  matrice_inv[0][1] = -(a01*a22-a02*a21)*t17;
  matrice_inv[0][2] = -(-a01*a12+a02*a11)*t17;
  matrice_inv[1][0] = (-a10*a22+a12*a20)*t17;
  matrice_inv[1][1] = (a00*a22-t14)*t17;
  matrice_inv[1][2] = -(t6-t10)*t17;
  matrice_inv[2][0] = -(-a10*a21+a11*a20)*t17;
  matrice_inv[2][1] = -(a00*a21-t12)*t17;
  matrice_inv[2][2] = (t4-t8)*t17;
}

/*! @brief Cette methode (statique) permet de calculer le produit d'une matrice 3x3 avec un vecteur 3
 *
 * @param (matrice) matrice 3x3 a inverser
 * @param (vect) vecteur 3 a multiplier
 * @param (vect) vecteur 3 resultat
 * @return (rien)
 */
void Parcours_interface::calcul_produit_matrice33_vecteur(const FTd_matrice33& matrice, const FTd_vecteur3& vect, FTd_vecteur3& res)
{
  int i,j;
  for (i=0 ; i<3 ; i++)
    {
      double x = 0.;
      for (j=0 ; j<3 ; j++)
        {
          x += matrice[i][j] * vect[j];
        }
      res[i] = x;
    }
}

/*! @brief Calcul de la contribution de volume d'une facette a la valeur de l'indicatrice dans un element.
 *
 * C'est une fraction du volume de l'element
 *    comprise entre epsilon et 1.-epsilon
 *
 * Precondition: dimension = 3
 *
 * @param (domaine_vf) domaine du calcul
 * @param (num_element) indice de l'element intersecte
 * @param (poly_reelles) coordonnees (reelles) des sommets definissant une surface contenue dans l'element (en pratique : surface d'intersection entre une facette d'interface et l'element)
 * @param (centre_de_gravite) centre de gravite de la surface
 * @param (epsilon) erreur relative
 * @return (double) contribution du volume engendre par la surface dans l'element
 */
double Parcours_interface::volume_tetraedre(const Domaine_VF& domaine_vf,
                                            int num_element,
                                            int num_facette,
                                            const Maillage_FT_Disc& maillage,
                                            const DoubleTab& poly_reelles,
                                            const FTd_vecteur3& centre_de_gravite,
                                            double epsilon) const
{
  static const int SOM_Z = 0;  // indice du sommet qui servira d'origine
  static const int SOM_A = 1;     // indice du sommet qui servira pour le premier vecteur
  static const int SOM_B = 2;     // indice du sommet qui servira pour le second vecteur
  static const int SOM_C = 3;     // indice du sommet qui servira pour le troisieme vecteur

  const int som_Z = (*domaine_elem_ptr)(num_element,SOM_Z);
  const int som_A = (*domaine_elem_ptr)(num_element,SOM_A);
  const int som_B = (*domaine_elem_ptr)(num_element,SOM_B);
  const int som_C = (*domaine_elem_ptr)(num_element,SOM_C);

  FTd_matrice33 matrice, matrice_inv;
  //on commence par calculer la matrice de transformation de la base (OX,OY,OZ) vers (ZA,ZB,ZC)
  //colonne 0 : vecteur OA
  matrice[0][0] = (*domaine_sommets_ptr)(som_A, 0) - (*domaine_sommets_ptr)(som_Z, 0); //OA(x)
  matrice[1][0] = (*domaine_sommets_ptr)(som_A, 1) - (*domaine_sommets_ptr)(som_Z, 1); //OA(y)
  matrice[2][0] = (*domaine_sommets_ptr)(som_A, 2) - (*domaine_sommets_ptr)(som_Z, 2); //OA(z)
  //colonne 1 : vecteur OB
  matrice[0][1] = (*domaine_sommets_ptr)(som_B, 0) - (*domaine_sommets_ptr)(som_Z, 0); //OB(x)
  matrice[1][1] = (*domaine_sommets_ptr)(som_B, 1) - (*domaine_sommets_ptr)(som_Z, 1); //OB(y)
  matrice[2][1] = (*domaine_sommets_ptr)(som_B, 2) - (*domaine_sommets_ptr)(som_Z, 2); //OB(z)
  //colonne 2 : vecteur OC
  matrice[0][2] = (*domaine_sommets_ptr)(som_C, 0) - (*domaine_sommets_ptr)(som_Z, 0); //OC(x)
  matrice[1][2] = (*domaine_sommets_ptr)(som_C, 1) - (*domaine_sommets_ptr)(som_Z, 1); //OC(y)
  matrice[2][2] = (*domaine_sommets_ptr)(som_C, 2) - (*domaine_sommets_ptr)(som_Z, 2); //OC(z)

  //on inverse ensuite la matrice de transformation
  calcul_inverse_matrice33(matrice,matrice_inv);

  int i,k;

  //calcul des coordonnes du polygone dans le repere de reference
  const int nb_sommets_poly = poly_reelles.dimension(0);
  DoubleTabFT poly_reelles_ref(nb_sommets_poly,3);
  FTd_vecteur3 poly, poly_ref;
  for (i=0 ; i<nb_sommets_poly ; i++)
    {
      for (k=0 ; k<3 ; k++)
        {
          poly[k] = poly_reelles(i,k) - (*domaine_sommets_ptr)(som_Z, k);
        }
      calcul_produit_matrice33_vecteur(matrice_inv,poly, poly_ref);
      for (k=0 ; k<3 ; k++)
        {
          poly_reelles_ref(i,k) = poly_ref[k];
        }
    }

  //calcul de la normale a la facette dans le repere de reference
  FTd_vecteur3 normale_ref;
  double AB[3], AC[3], AB_ref[3], AC_ref[3];
  const int somA = maillage.facettes_(num_facette, 0);
  const int somB = maillage.facettes_(num_facette, 1);
  const int somC = maillage.facettes_(num_facette, 2);
  for (k=0 ; k<dimension ; k++)
    {
      AB[k] = maillage.sommets_(somB,k) - maillage.sommets_(somA,k);
      AC[k] = maillage.sommets_(somC,k) - maillage.sommets_(somA,k);
    }
  calcul_produit_matrice33_vecteur(matrice_inv,AB, AB_ref);
  calcul_produit_matrice33_vecteur(matrice_inv,AC, AC_ref);
  //calcul pdt vect AB_ref x AC_ref
  normale_ref[0] = AB_ref[1]*AC_ref[2] - AB_ref[2]*AC_ref[1];
  normale_ref[1] = AB_ref[2]*AC_ref[0] - AB_ref[0]*AC_ref[2];
  normale_ref[2] = AB_ref[0]*AC_ref[1] - AB_ref[1]*AC_ref[0];

  //calcul de la contribution de volume d'une facette a la valeur de l'indicatrice dans l'element de reference
  //le calcul est fait en decomposant les polygones en triangles elementaires
  FTd_vecteur3 centre_de_gravite_triangle_ref;
  static DoubleTabFT triangle_ref(3,3);
  double contrib_volume = 0., vol;
  for (k=0 ; k<dimension ; k++)
    {
      triangle_ref(0,k) = poly_reelles_ref(0,k);
    }
  for (i=1 ; i<nb_sommets_poly-1 ; i++)
    {
      for (k=0 ; k<dimension ; k++)
        {
          triangle_ref(1,k) = poly_reelles_ref(i,k);
          triangle_ref(2,k) = poly_reelles_ref(i+1,k);
          centre_de_gravite_triangle_ref[k] = (triangle_ref(0,k)+triangle_ref(1,k)+triangle_ref(2,k)) /3.;
        }
      vol = volume_tetraedre_reference(triangle_ref,normale_ref,centre_de_gravite_triangle_ref, epsilon);
      contrib_volume += vol;
    }

  return contrib_volume;
}

/*! @brief Calcul de la contribution de volume d'une facette a la valeur de l'indicatrice dans le tetraedre de reference.
 *
 *    Dans ce tetraedre, on utilise le plan OXZ comme plan de projection (projection selon OY donc)
 *
 * Precondition: la surface doit etre un triangle (nb_sommets_poly==3)
 *
 * @param (poly_reelles_ref) coordonnees (reelles) des sommets l'element de reference
 * @param (norme) normale a la surface, dans l'element de reference
 * @param (centre_de_gravite, dans l'element de reference) centre de gravite de la surface
 * @param (epsilon) erreur relative
 * @return (double) contribution du volume engendre par la surface dans l'element
 */
double Parcours_interface::volume_tetraedre_reference(const DoubleTab& poly_reelles_ref,
                                                      const FTd_vecteur3& norme_ref,
                                                      const FTd_vecteur3& centre_de_gravite_ref,
                                                      double epsilon) const
{
  //determination des signes pour les differentes composantes :
  double signe_princ;
  if (norme_ref[1]>0.)
    {
      signe_princ = -1.;
    }
  else
    {
      signe_princ = 1.;
    }
  // Volume de l'element = base * hauteur /3
  //  avec base = 0.5, hauteur = 1
  double v_elem = 0.5 /3.;
  assert(v_elem>0.);

  const int nb_sommets_poly = poly_reelles_ref.dimension(0);
  int k,i, i_prec = nb_sommets_poly -1;
  //on entre dans cette fonction qu'avec des triangles
  assert(nb_sommets_poly==3);
  //calcul de la normale au triangle dans le repere de reference
  FTd_vecteur3 normale_tr;
  double AB[3], AC[3];
  assert(dimension==3);
  for (k=0 ; k<3 ; k++)
    {
      AB[k] = poly_reelles_ref(1,k) - poly_reelles_ref(0,k);
      AC[k] = poly_reelles_ref(2,k) - poly_reelles_ref(0,k);
    }
  //calcul pdt vect AB_ref x AC_ref
  normale_tr[0] = AB[1]*AC[2] - AB[2]*AC[1];
  normale_tr[1] = AB[2]*AC[0] - AB[0]*AC[2];
  normale_tr[2] = AB[0]*AC[1] - AB[1]*AC[0];
  const double aire_tr = normale_tr[0]*normale_tr[0]+normale_tr[1]*normale_tr[1]+normale_tr[2]*normale_tr[2];
  if (aire_tr<Erreur_relative_maxi_)
    {
      return 0;
    }
  int coupe_face_opp = -1;
  int coupe_face_opp_p1 = -1;
  //calcule l'aire de la surface projetee sur face_bas
  double aire_projetee = 0.;
  for (i=0 ; i<nb_sommets_poly ; i++)
    {
      //teste si l'arete [i_prec,i] est sur la face_opposee
      //l'equation de plan de la face opposee est : x+y+z-1=0
      double f_prec =  poly_reelles_ref(i_prec,0) +  poly_reelles_ref(i_prec,1) +  poly_reelles_ref(i_prec,2) - 1.;
      double f =  poly_reelles_ref(i,0) +  poly_reelles_ref(i,1) +  poly_reelles_ref(i,2) - 1.;
      if (f<Erreur_relative_maxi_ && f>-Erreur_relative_maxi_ && f_prec<Erreur_relative_maxi_ && f_prec>-Erreur_relative_maxi_)
        {
          //les somet i_prec et i sont dans le plan de la face opposee au sommet 0
          //comme le polygone est dans le tetraedre, les somets sont dans la face...
          coupe_face_opp = i_prec;
          coupe_face_opp_p1 = i;
        }

      //calcul de la contribution de l'arete a l'aire projetee dans le plan (OX,OZ) :
      // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe x=0)
      double contrib_aire_projetee = (poly_reelles_ref(i,2) - poly_reelles_ref(i_prec,2)) * (poly_reelles_ref(i,0) + poly_reelles_ref(i_prec,0)) * 0.5;

      aire_projetee += contrib_aire_projetee;
      i_prec = i;
    }
  double v = 0.;

  if (aire_projetee!=0.)
    {
      // Volume de la partie "projection sur la face face_bas (y_bas=0)
      v = signe_princ * std::fabs(aire_projetee) * (centre_de_gravite_ref[1]);
    }
  if (coupe_face_opp!=-1)
    {
      //la surface possede une arete sur la face opposee au sommet 0 (face 0)
      // -> il faut ajouter un volume complementaire
      //ce volume est defini par :
      // - les sommets coupe_face_opp et coupe_face_opp_p1 (notes som0 et som1)
      // - leur projection orthogonale sur la face_bas (ie ces sommets avec y=0) : proj_som0 et proj_som1
      // - la projection des sommets coupe_face_opp et coupe_face_opp_p1 sur face_bas, a partir du sommet (0,1,0) : proj2_som0 et proj2_som1
      //(cette derniere projection etant d'ailleurs equivalente a la projection de proj_som0 et proj_som1 sur l'arete (1,0,0)(0,0,1) a prtir du sommet (0,0,0)
      FTd_vecteur3 som0,som1, proj_som0,proj_som1, proj2_som0,proj2_som1;
      //recuperation des sommets sur la face 0
      for (k=0 ; k<dimension ; k++)
        {
          som0[k] = poly_reelles_ref(coupe_face_opp,k);
          som1[k] = poly_reelles_ref(coupe_face_opp_p1,k);
        }
      //calcul de leur projection sur y=0
      proj_som0[0] = som0[0];
      proj_som1[0] = som1[0];
      proj_som0[1] = 0.;
      proj_som1[1] = 0.;
      proj_som0[2] = som0[2];
      proj_som1[2] = som1[2];
      //calcul de leur projection sur l'arete (1,0,0)-(0,0,1)
      //    0                     1      X
      // 0  *---------------------------->
      //    |  ---              /
      //    |!    ---         /
      //    |!       -1-    /
      //    | !         --+
      //    | !         /  ---
      //    |  O      /
      //    |  !    /
      //    |   ! /
      //    |   +
      //    | /  !
      // 1  |    !
      //    |
      // Z  V
      //on va donc calculer
      //  -l'intersection des droites : z = z0/x0 . x  et z = 1 - x
      //  -l'intersection des droites : z = z1/x1 . x  et z = 1 - x
      if (som0[0]<epsilon)   //x0==0
        {
          proj2_som0[0] = 0.;
        }
      else
        {
          //x20 = 1 / (1 + z0/x0)
          assert((1. + som0[2]/som0[0]) !=0.);
          proj2_som0[0] = 1./ (1. + som0[2]/som0[0]);
        }
      if (som1[0]<epsilon)   //x1==0
        {
          proj2_som1[0] = 0.;
        }
      else
        {
          //x21 = 1 / (1 + z1/x1)
          assert((1. + som1[2]/som1[0]) !=0.);
          proj2_som1[0] = 1./ (1. + som1[2]/som1[0]);
        }
      proj2_som0[1] = 0.;
      proj2_som1[1] = 0.;
      proj2_som0[2] = 1. - proj2_som0[0];
      proj2_som1[2] = 1. - proj2_som1[0];

      //determination du signe des volumes complementaires
      FTd_vecteur3 normale_plan_proj;
      normale_plan_proj[0] = -(proj_som1[2]-proj_som0[2]);
      normale_plan_proj[1] = 0.;
      normale_plan_proj[2] = proj_som1[0]-proj_som0[0];
      double pdtscal = 0.;
      FTd_vecteur3 vect;
      for (i=0 ; i<dimension ; i++)
        {
          vect[i] = proj2_som0[i]- proj_som0[i];
          pdtscal += vect[i]*vect[i];
        }
      if (pdtscal<Erreur_relative_maxi_)
        {
          for (i=0 ; i<dimension ; i++)
            {
              vect[i] = proj2_som1[i]- proj_som1[i];
            }
        }
      pdtscal = 0.;
      for (i=0 ; i<dimension ; i++)
        {
          pdtscal += normale_plan_proj[i] * vect[i];
        }
      double signe_2 = 0.;
      if (pdtscal>0.)
        {
          signe_2 = -1.;
        }
      else
        {
          signe_2 = 1.;
        }
      //une fois les pojetes calcules, il faut calculer le volume correspondant
      //ce volume est constitue par les 4 projetes et les 2 sommets opposes
      //on decompose alors ce volume en :
      // - une pyramide a base quadrangulaire (som0, som1 et leur projete orthogonal qui definissent une base trapezoidale) et le sommet oppose proj2_som0
      // - un tetraedre defini par les sommets proj2_som0, som1 et les 2 projetes de som1
      //calcul du volume de la pyramide
      double l,L,base,h, vol;
      //  calcul de l'aire de la base de la pyramide
      l = sqrt((som1[0]-som0[0])*(som1[0]-som0[0]) + (som1[2]-som0[2])*(som1[2]-som0[2])); //hauteur du trapeze
      L = 0.5* (som0[1] + som1[1]); //moyenne des longueurs des bases du trapeze
      base = L*l;  //base = aire du trapeze (som0,som1,proj_som0,proj_som1)
      if (l>=Erreur_relative_maxi_)
        {
          //  calcul de la hauteur de la pyramide
          double L2 = (proj2_som0[0]-som0[0])*(proj2_som0[0]-som0[0]) + (proj2_som0[2]-som0[2])*(proj2_som0[2]-som0[2]);
          double l2 = ((som1[0]-som0[0])*(proj2_som0[0]-som0[0]) + (som1[2]-som0[2])*(proj2_som0[2]-som0[2])) / l;
          h = L2 - l2*l2;
          if (h>=Erreur_relative_maxi_*Erreur_relative_maxi_)
            {
              h = sqrt(h);
              //  calcul du volume de la pyramide : base * hauteur / 3.
              vol = base * h /3.;
              v += signe_2 * vol;
            }
        }

      //calcul du volume du tetraedre
      //  calcul de l'aire de la base du tetraedre : aire du triangle proj_som1,proj2_som0,proj2_som1
      base = (proj_som1[2]-proj2_som0[2]) * (proj2_som1[0]-proj2_som0[0]) - (proj_som1[0]-proj2_som0[0]) * (proj2_som1[2]-proj2_som0[2]);
      base = std::fabs(base * 0.5);
      //  calcul de la hauteur du tetraedre
      h = som1[1];
      //  calcul du volume du tetraedre : base * hauteur / 3.
      vol = base * h /3.;
      v += signe_2 * vol;
    }

  //normalisation par le volume de l'element
  v /= v_elem;

  // On force la valeur entre 0 et 1 strictement.
  if (std::fabs(v) < Erreur_relative_maxi_)
    v = Erreur_relative_maxi_;
  else if (v > 1. - Erreur_relative_maxi_)
    v = 1. - Erreur_relative_maxi_;
  else if (v < -1. + Erreur_relative_maxi_)
    v = -1. + Erreur_relative_maxi_;

  return v;
}

/*! @brief Calcul de la contribution de volume d'une facette a la valeur de l'indicatrice dans un element.
 *
 * C'est une fraction du volume de l'element
 *    comprise entre epsilon et 1.-epsilon
 *
 * Precondition: dimension = 3
 *
 * @param (domaine_vf) domaine du calcul
 * @param (num_element) indice de l'element intersecte
 * @param (poly_reelles) coordonnees (reelles) des sommets definissant une surface contenue dans l'element (en pratique : surface d'intersection entre une facette d'interface et l'element)
 * @param (norme) normale a la facette
 * @param (centre_de_gravite) centre de gravite de la surface
 * @param (epsilon) erreur relative
 * @return (double) contribution du volume engendre par la surface dans l'element
 */
double Parcours_interface::volume_hexaedre(const Domaine_VF& domaine_vf,
                                           int num_element,
                                           const DoubleTab& poly_reelles,
                                           const FTd_vecteur3& norme,
                                           const FTd_vecteur3& centre_de_gravite,
                                           const ArrOfInt& polygone_plan_coupe,
                                           double epsilon) const
{
  // Conventions TRUST VDF :
  static const int NUM_FACE_GAUCHE = 0;
  static const int NUM_FACE_DROITE = 3;
  static const int NUM_FACE_BAS = 1;
  static const int NUM_FACE_HAUT = 4;
  static const int NUM_FACE_ARRIERE = 2;
  static const int NUM_FACE_AVANT = 5;
  int face_bas = domaine_vf.elem_faces(num_element, NUM_FACE_BAS);
  int face_haut = domaine_vf.elem_faces(num_element, NUM_FACE_HAUT);
  int face_gauche = domaine_vf.elem_faces(num_element, NUM_FACE_GAUCHE);
  int face_droite = domaine_vf.elem_faces(num_element, NUM_FACE_DROITE);
  int face_arriere = domaine_vf.elem_faces(num_element, NUM_FACE_ARRIERE);
  int face_avant = domaine_vf.elem_faces(num_element, NUM_FACE_AVANT);
  double y_bas = domaine_vf.xv(face_bas, 1);
  double y_haut = domaine_vf.xv(face_haut, 1);
  double x_gauche = domaine_vf.xv(face_gauche, 0);
  double x_droite = domaine_vf.xv(face_droite, 0);
  double z_arriere = domaine_vf.xv(face_arriere, 2);
  double z_avant = domaine_vf.xv(face_avant, 2);

  //determination des signes pour les differentes composantes :
  double signe_princ, signe_compl0,signe_compl1;
  if (norme[1]>0.)
    {
      signe_princ = -1;
    }
  else if (norme[1]<0.)
    {
      signe_princ = 1;
    }
  else
    {
      signe_princ = 0;
    }
  if (norme[0]>0.)
    {
      signe_compl0 = -1;
    }
  else if (norme[0]<0.)
    {
      signe_compl0 = +1;
    }
  else
    {
      signe_compl0 = 0;
    }
  if (norme[2]>0.)
    {
      signe_compl1 = -1;
    }
  else if (norme[2]<0.)
    {
      signe_compl1 = 1;
    }
  else
    {
      signe_compl1 = 0;
    }

  // Volume de l'element
  double v_elem = (x_droite - x_gauche) * (y_haut - y_bas) * (z_avant - z_arriere);
  assert(v_elem>0.);

  const int nb_sommets_poly = poly_reelles.dimension(0);
  int i, i_prec = nb_sommets_poly -1;
  //calcule l'aire de la surface projetee sur face_bas
  double aire_projetee = 0.;
  for (i=0 ; i<nb_sommets_poly ; i++)
    {
      //calcul de la contribution de l'arete a l'aire projetee dans le plan (X,Z) :
      // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe x=0)
      double contrib_aire_projetee = (poly_reelles(i,2) - poly_reelles(i_prec,2)) * ((poly_reelles(i,0) + poly_reelles(i_prec,0)) * 0.5);

      aire_projetee += contrib_aire_projetee;
      i_prec = i;
    }
  double v = 0.;

  // Volume de la partie "projection sur la face face_bas"
  v = signe_princ * std::fabs(aire_projetee) * (centre_de_gravite[1] - y_bas);

  // Recherche d'un segment coupant la face du haut
  for (i = 0; i < nb_sommets_poly; i++)
    {
      if (polygone_plan_coupe[i] == NUM_FACE_HAUT)
        {
          const int i_precedent = (i-1) < 0 ? nb_sommets_poly-1 : (i-1);
          const int i_suivant   = (i+1) >= nb_sommets_poly ? 0 : (i+1);
          // Le segment [i-1,i] et sur la face du haut
          const int coupe_face_haut = i_precedent;
          const int coupe_face_haut_p1 = i;
          // Un segment voisin est-il sur la face de droite ?
          int coupe_face_droite = -1;
          if (polygone_plan_coupe[i_precedent] == NUM_FACE_DROITE)
            coupe_face_droite = i_precedent;
          else if (polygone_plan_coupe[i_suivant] == NUM_FACE_DROITE)
            coupe_face_droite = i;//_suivant;

          // Les sommets coupe_haut et coupe_haut_p1 coupent la face face_haut
          // -> ajoute la composante de volume associee = moyenne_x * (y_haut-y_bas) * diff_z
          double vol_compl0 =
            (0.5*(poly_reelles(coupe_face_haut,0) + poly_reelles(coupe_face_haut_p1,0)) - x_gauche)
            * (y_haut - y_bas)
            * (poly_reelles(coupe_face_haut,2) - poly_reelles(coupe_face_haut_p1,2));

          v += signe_compl0 * std::fabs(vol_compl0);

          if (coupe_face_droite >= 0)
            {
              //cette arete coupe aussi la face xmax
              //ajoute la composante de volume associee = (x_droite-x_gauche) * (y_haut-y_bas) * (z-z_arriere)
              double vol_compl1 =
                (x_droite - x_gauche)
                * (y_haut - y_bas)
                * (poly_reelles(coupe_face_droite,2) - z_arriere);

              v += signe_compl1 * std::fabs(vol_compl1);
            }
        }
    }

  //normalisation par le volume de l'element
  v /= v_elem;

  // On force la valeur entre 0 et 1 strictement.
#if 0
  // B.Math 10/09/2004: il faut corriger le calcul pour que les
  // contributions soient entre 0 et 1, ensuite on pourra faire ca:
  if (v < Erreur_relative_maxi_)
    v = Erreur_relative_maxi_;
  else if (v > 1. - Erreur_relative_maxi_)
    v = 1. - Erreur_relative_maxi_;
#endif

  return v;
}

double Parcours_interface::volume_hexaedre_arete_interne(const Domaine_VDF& domaine_vdf, int num_arete,
                                                         const DoubleTab& poly_reelles,
                                                         const FTd_vecteur3& norme,
                                                         const FTd_vecteur3& centre_de_gravite,
                                                         const ArrOfInt& polygone_plan_coupe,
                                                         double epsilon) const
{
  // Conventions TRUST VDF :
  static const int NUM_FACE_DROITE_ARETE = 3;
  static const int NUM_FACE_HAUT_ARETE = 4;

  const int nb_elems_reels=domaine_vdf.nb_elem();

  const IntTab& Qdm=domaine_vdf.Qdm();
  const IntTab& Aretes_Som=domaine_vdf.domaine().aretes_som();
  const DoubleTab& coord_som=domaine_vdf.domaine().coord_sommets();
  const IntVect& orientation_aretes=domaine_vdf.orientation_aretes();

  const int arete_multiple=(domaine_vdf.aretes_multiples()(num_arete)>0); // peu importe si l'arete est multiple 1, 2 ou 3, on decoupe en 4 dans tous les cas

  double y_bas=0,y_haut=0,x_gauche=0,x_droite=0,z_arriere=0,z_avant=0;



  int ori_arete=dimension-orientation_aretes(num_arete)-1;
  int face1=Qdm(num_arete,0);
  int face2=Qdm(num_arete,1);
  int face3=Qdm(num_arete,2);
  int face4=Qdm(num_arete,3);

  int s1=Aretes_Som(num_arete,0);
  int s2=Aretes_Som(num_arete,1);

  assert(face1>=0 && face2>=0 && face3>=0 && face4>=0);

  if (!arete_multiple)
    {
      if (ori_arete==0)
        {
          x_gauche=coord_som(s1,0);
          x_droite=coord_som(s2,0);
          y_haut=domaine_vdf.xv(face2,1);
          y_bas=domaine_vdf.xv(face1,1);
          z_avant=domaine_vdf.xv(face4,2);
          z_arriere=domaine_vdf.xv(face3,2);
        }
      else if (ori_arete==1)
        {
          y_bas=coord_som(s1,ori_arete);
          y_haut=coord_som(s2,ori_arete);
          x_gauche=domaine_vdf.xv(face1,0);
          x_droite=domaine_vdf.xv(face2,0);
          z_avant=domaine_vdf.xv(face4,2);
          z_arriere=domaine_vdf.xv(face3,2);
        }
      else if (ori_arete==2)
        {
          z_arriere=coord_som(s1,ori_arete);
          z_avant=coord_som(s2,ori_arete);
          x_gauche=domaine_vdf.xv(face1,0);
          x_droite=domaine_vdf.xv(face2,0);
          y_haut=domaine_vdf.xv(face4,1);
          y_bas=domaine_vdf.xv(face3,1);
        }
    }
  else
    {
      int elem1,elem2,elem3,elem4;
      elem1=domaine_vdf.face_voisins(face1,0);
      elem2=domaine_vdf.face_voisins(face2,0);
      elem3=domaine_vdf.face_voisins(face1,1);
      elem4=domaine_vdf.face_voisins(face2,1);

      if (ori_arete==0)
        {
          x_gauche=coord_som(s1,0);
          x_droite=coord_som(s2,0);
          if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem2<nb_elems_reels && elem4<nb_elems_reels)
            {
              y_bas=domaine_vdf.xa(num_arete,1);
              y_haut=domaine_vdf.xv(face2,1);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem2>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem3<nb_elems_reels)
            {
              y_haut=domaine_vdf.xa(num_arete,1);
              y_bas=domaine_vdf.xv(face1,1);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3<nb_elems_reels && elem3<nb_elems_reels)
            {
              y_bas=domaine_vdf.xv(face1,1);
              y_haut=domaine_vdf.xv(face2,1);
              z_arriere=domaine_vdf.xa(num_arete,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem3>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem2<nb_elems_reels)
            {
              y_bas=domaine_vdf.xv(face1,1);
              y_haut=domaine_vdf.xv(face2,1);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xa(num_arete,2);
            }
          else if (elem2>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              y_bas=domaine_vdf.xv(face1,1);
              y_haut=domaine_vdf.xa(num_arete,1);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xa(num_arete,2);
            }
          else if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              y_bas=domaine_vdf.xa(num_arete,1);
              y_haut=domaine_vdf.xv(face2,1);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xa(num_arete,2);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3>=nb_elems_reels)
            {
              y_bas=domaine_vdf.xa(num_arete,1);
              y_haut=domaine_vdf.xv(face2,1);
              z_arriere=domaine_vdf.xa(num_arete,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              y_bas=domaine_vdf.xv(face1,1);
              y_haut=domaine_vdf.xa(num_arete,1);
              z_arriere=domaine_vdf.xa(num_arete,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
        }

      else if (ori_arete==1)
        {
          y_bas=coord_som(s1,ori_arete);
          y_haut=coord_som(s2,ori_arete);
          if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem2<nb_elems_reels && elem4<nb_elems_reels)
            {
              x_gauche=domaine_vdf.xa(num_arete,0);
              x_droite=domaine_vdf.xv(face2,0);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem2>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem3<nb_elems_reels)
            {
              x_droite=domaine_vdf.xa(num_arete,0);
              x_gauche=domaine_vdf.xv(face1,0);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3<nb_elems_reels && elem3<nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xv(face2,0);
              z_arriere=domaine_vdf.xa(num_arete,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem3>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem2<nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xv(face2,0);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xa(num_arete,2);
            }
          else if (elem2>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xa(num_arete,0);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xa(num_arete,2);
            }
          else if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xa(num_arete,0);
              x_droite=domaine_vdf.xv(face2,0);
              z_arriere=domaine_vdf.xv(face3,2);
              z_avant=domaine_vdf.xa(num_arete,2);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xa(num_arete,0);
              x_droite=domaine_vdf.xv(face2,0);
              z_arriere=domaine_vdf.xa(num_arete,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xa(num_arete,0);
              z_arriere=domaine_vdf.xa(num_arete,2);
              z_avant=domaine_vdf.xv(face4,2);
            }
        }

      else if (ori_arete==2)
        {
          z_arriere=coord_som(s1,ori_arete);
          z_avant=coord_som(s2,ori_arete);
          if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem2<nb_elems_reels && elem4<nb_elems_reels)
            {
              x_gauche=domaine_vdf.xa(num_arete,0);
              x_droite=domaine_vdf.xv(face2,0);
              y_bas=domaine_vdf.xv(face3,1);
              y_haut=domaine_vdf.xv(face4,1);
            }
          else if (elem2>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem3<nb_elems_reels)
            {
              x_droite=domaine_vdf.xa(num_arete,0);
              x_gauche=domaine_vdf.xv(face1,0);
              y_bas=domaine_vdf.xv(face3,1);
              y_haut=domaine_vdf.xv(face4,1);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3<nb_elems_reels && elem3<nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xv(face2,0);
              y_bas=domaine_vdf.xa(num_arete,1);
              y_haut=domaine_vdf.xv(face4,1);
            }
          else if (elem3>=nb_elems_reels && elem4>=nb_elems_reels && elem1<nb_elems_reels && elem2<nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xv(face2,0);
              y_bas=domaine_vdf.xv(face3,1);
              y_haut=domaine_vdf.xa(num_arete,1);
            }
          else if (elem2>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xa(num_arete,0);
              y_bas=domaine_vdf.xv(face3,1);
              y_haut=domaine_vdf.xa(num_arete,1);
            }
          else if (elem1>=nb_elems_reels && elem3>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xa(num_arete,0);
              x_droite=domaine_vdf.xv(face2,0);
              y_bas=domaine_vdf.xv(face3,1);
              y_haut=domaine_vdf.xa(num_arete,1);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem3>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xa(num_arete,0);
              x_droite=domaine_vdf.xv(face2,0);
              y_bas=domaine_vdf.xa(num_arete,1);
              y_haut=domaine_vdf.xv(face4,1);
            }
          else if (elem1>=nb_elems_reels && elem2>=nb_elems_reels && elem4>=nb_elems_reels)
            {
              x_gauche=domaine_vdf.xv(face1,0);
              x_droite=domaine_vdf.xa(num_arete,0);
              y_bas=domaine_vdf.xa(num_arete,1);
              y_haut=domaine_vdf.xv(face4,1);
            }
        }
    }
  double signe_princ, signe_compl0,signe_compl1;
  if (norme[1]>0.)
    {
      signe_princ = -1;
    }
  else if (norme[1]<0.)
    {
      signe_princ = 1;
    }
  else
    {
      signe_princ = 0;
    }
  if (norme[0]>0.)
    {
      signe_compl0 = -1;
    }
  else if (norme[0]<0.)
    {
      signe_compl0 = +1;
    }
  else
    {
      signe_compl0 = 0;
    }
  if (norme[2]>0.)
    {
      signe_compl1 = -1;
    }
  else if (norme[2]<0.)
    {
      signe_compl1 = 1;
    }
  else
    {
      signe_compl1 = 0;
    }

  // Volume de l'element
  double v_elem = (x_droite - x_gauche) * (y_haut - y_bas) * (z_avant - z_arriere);
  assert(v_elem>0.);

  const int nb_sommets_poly = poly_reelles.dimension(0);
  int i, i_prec = nb_sommets_poly -1;
  //calcule l'aire de la surface projetee sur face_bas
  double aire_projetee = 0.;
  for (i=0 ; i<nb_sommets_poly ; i++)
    {
      //calcul de la contribution de l'arete a l'aire projetee dans le plan (X,Z) :
      // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe x=0)
      double contrib_aire_projetee = (poly_reelles(i,2) - poly_reelles(i_prec,2)) * ((poly_reelles(i,0) + poly_reelles(i_prec,0)) * 0.5);

      aire_projetee += contrib_aire_projetee;
      i_prec = i;
    }
  double v = 0.;

  // Volume de la partie "projection sur la face face_bas"
  v = signe_princ * fabs(aire_projetee) * (centre_de_gravite[1] - y_bas);

  // Recherche d'un segment coupant la face du haut
  for (i = 0; i < nb_sommets_poly; i++)
    {
      if (polygone_plan_coupe[i] == NUM_FACE_HAUT_ARETE)
        {
          const int i_precedent = (i-1) < 0 ? nb_sommets_poly-1 : (i-1);
          const int i_suivant   = (i+1) >= nb_sommets_poly ? 0 : (i+1);
          // Le segment [i-1,i] et sur la face du haut
          const int coupe_face_haut = i_precedent;
          const int coupe_face_haut_p1 = i;
          // Un segment voisin est-il sur la face de droite ?
          int coupe_face_droite = -1;
          if (polygone_plan_coupe[i_precedent] == NUM_FACE_DROITE_ARETE)
            coupe_face_droite = i_precedent;
          else if (polygone_plan_coupe[i_suivant] == NUM_FACE_DROITE_ARETE)
            coupe_face_droite = i;//_suivant;

          // Les sommets coupe_haut et coupe_haut_p1 coupent la face face_haut
          // -> ajoute la composante de volume associee = moyenne_x * (y_haut-y_bas) * diff_z
          double vol_compl0 =
            (0.5*(poly_reelles(coupe_face_haut,0) + poly_reelles(coupe_face_haut_p1,0)) - x_gauche)
            * (y_haut - y_bas)
            * (poly_reelles(coupe_face_haut,2) - poly_reelles(coupe_face_haut_p1,2));

          v += signe_compl0 * std::fabs(vol_compl0);

          if (coupe_face_droite >= 0)
            {
              //cette arete coupe aussi la face xmax
              //ajoute la composante de volume associee = (x_droite-x_gauche) * (y_haut-y_bas) * (z-z_arriere)
              double vol_compl1 =
                (x_droite - x_gauche)
                * (y_haut - y_bas)
                * (poly_reelles(coupe_face_droite,2) - z_arriere);

              v += signe_compl1 * std::fabs(vol_compl1);
            }
        }
    }

  //normalisation par le volume de l'element
  v /= v_elem;

  // On force la valeur entre 0 et 1 strictement.
#if 0
  // B.Math 10/09/2004: il faut corriger le calcul pour que les
  // contributions soient entre 0 et 1, ensuite on pourra faire ca:
  if (v < Erreur_relative_maxi_)
    v = Erreur_relative_maxi_;
  else if (v > 1. - Erreur_relative_maxi_)
    v = 1. - Erreur_relative_maxi_;
#endif

  return v;

}


// debut EB : copie de la fonction volume_hexaedre mais en adaptant la position des sommets pour les faces
double Parcours_interface::volume_hexaedre_face(const Domaine_VDF& domaine_vdf, int num_face,
                                                const DoubleTab& poly_reelles,
                                                const FTd_vecteur3& norme,
                                                const FTd_vecteur3& centre_de_gravite,
                                                const ArrOfInt& polygone_plan_coupe,
                                                double epsilon) const
{
  // Conventions TRUST VDF :
  static const int NUM_FACE_DROITE_FACE = 3;
  static const int NUM_FACE_HAUT_FACE = 4;
  /*
  static int affichage_x=0;
  static int affichage_y=0;
  static int affichage_z=0;
  */
  int face_gauche=0, face_droite=0, face_bas=0, face_haut=0, face_avant=0, face_arriere=0;
  double y_bas=0,y_haut=0,x_gauche=0,x_droite=0,z_arriere=0,z_avant=0;
  const int elem0 = domaine_vdf.face_voisins(num_face,0);
  const int elem1 = domaine_vdf.face_voisins(num_face,1);
  const int nb_elem_reel=domaine_vdf.nb_elem();
  // On verifie si on est sur une face de bord ou non
  const int face_bord0 = (elem0==-1 || (elem0>=nb_elem_reel)); // pas d'acces a l'element de gauche, bas, arriere
  const int face_bord1 = (elem1==-1 || (elem1>=nb_elem_reel)); // pas d'acces a l'element de droite, haut, avant
  int orientation = domaine_vdf.orientation(num_face);

  const int face_double=domaine_vdf.faces_doubles()[num_face];
  if (!face_bord0 && !face_bord1 && !face_double) // au debut on utilise elem0 ou elem1
    {
      face_gauche=domaine_vdf.elem_faces(elem0,0);
      face_droite=domaine_vdf.elem_faces(elem0,0+dimension);
      face_bas=domaine_vdf.elem_faces(elem0,1);
      face_haut=domaine_vdf.elem_faces(elem0,1+dimension);
      face_arriere=domaine_vdf.elem_faces(elem0,2);
      face_avant=domaine_vdf.elem_faces(elem0,2+dimension);
      x_gauche=domaine_vdf.xv(face_gauche,0);
      x_droite=domaine_vdf.xv(face_droite,0);
      y_bas=domaine_vdf.xv(face_bas,1);
      y_haut=domaine_vdf.xv(face_haut,1);
      z_avant=domaine_vdf.xv(face_avant,2);
      z_arriere=domaine_vdf.xv(face_arriere,2);
      // on corrige ensuite pour les faces du volume de controle de meme orientation que la face num_face
      if (orientation==0)
        {
          x_gauche=domaine_vdf.xp(elem0,0);
          x_droite=domaine_vdf.xp(elem1,0);
        }
      if (orientation==1)
        {
          y_bas=domaine_vdf.xp(elem0,1);
          y_haut=domaine_vdf.xp(elem1,1);
        }
      if (orientation==2)
        {
          z_arriere=domaine_vdf.xp(elem0,2);
          z_avant=domaine_vdf.xp(elem1,2);
        }
    }
  else
    {
      int elem=-1;// = face_bord0 ? elem1 : elem0;
      if (face_double)
        {
          if (elem0>=nb_elem_reel) elem=elem1;
          else if (elem1>=nb_elem_reel) elem=elem0;
          else
            {
              Cerr << "Parcours_interface::"
                   " : face_double detectee mais elem0 et elem1 sont reels" << finl;
              exit();
            }
        }
      else
        {
          elem = face_bord0 ? elem1 : elem0;
        }
      face_gauche=domaine_vdf.elem_faces(elem,0);
      face_droite=domaine_vdf.elem_faces(elem,0+dimension);
      face_bas=domaine_vdf.elem_faces(elem,1);
      face_haut=domaine_vdf.elem_faces(elem,1+dimension);
      face_arriere=domaine_vdf.elem_faces(elem,2);
      face_avant=domaine_vdf.elem_faces(elem,2+dimension);
      x_gauche=domaine_vdf.xv(face_gauche,0);
      x_droite=domaine_vdf.xv(face_droite,0);
      y_bas=domaine_vdf.xv(face_bas,1);
      y_haut=domaine_vdf.xv(face_haut,1);
      z_avant=domaine_vdf.xv(face_avant,2);
      z_arriere=domaine_vdf.xv(face_arriere,2);
      // on corrige ensuite pour les faces du volume de controle de meme orientation que la face num_face
      //int bord = face_bord0==1 ? 0 : 1;
      if (orientation==0)
        {
          if (face_bord0)
            {
              x_gauche=domaine_vdf.xv(num_face,0);
              x_droite=domaine_vdf.xp(elem,0);
            }
          else
            {
              x_gauche=domaine_vdf.xp(elem,0);
              x_droite=domaine_vdf.xv(num_face,0);
            }
        }
      if (orientation==1)
        {
          if (face_bord0)
            {
              y_bas=domaine_vdf.xv(num_face,1);
              y_haut=domaine_vdf.xp(elem,1);
            }
          else
            {
              y_bas=domaine_vdf.xp(elem,1);
              y_haut=domaine_vdf.xv(num_face,1);
            }
        }
      if (orientation==2)
        {
          if (face_bord0)
            {
              z_arriere=domaine_vdf.xv(num_face,2);
              z_avant=domaine_vdf.xp(elem,2);
            }
          else
            {
              z_arriere=domaine_vdf.xp(elem,2);
              z_avant=domaine_vdf.xv(num_face,2);
            }
        }
    }
  //if (face_double && orientation==2 && fabs(domaine_vf.xv(num_face,0))<0.2e-3 && domaine_vf.xv(num_face,1)<9e-3 && domaine_vf.xv(num_face,1)>8.6e-3 &&  domaine_vf.xv(num_face,2)<0 && domaine_vf.xv(num_face,2)>-0.5e-3 )
  //Cerr << "x_gauche : " << x_gauche << " x_droite : " << x_droite << " y_bas : " << y_bas << " y_haut : " << y_haut << " z_arriere : " << z_arriere << " z_avant : " << z_avant << finl;
  //determination des signes pour les differentes composantes :
  double signe_princ, signe_compl0,signe_compl1;
  if (norme[1]>0.)
    {
      signe_princ = -1;
    }
  else if (norme[1]<0.)
    {
      signe_princ = 1;
    }
  else
    {
      signe_princ = 0;
    }
  if (norme[0]>0.)
    {
      signe_compl0 = -1;
    }
  else if (norme[0]<0.)
    {
      signe_compl0 = +1;
    }
  else
    {
      signe_compl0 = 0;
    }
  if (norme[2]>0.)
    {
      signe_compl1 = -1;
    }
  else if (norme[2]<0.)
    {
      signe_compl1 = 1;
    }
  else
    {
      signe_compl1 = 0;
    }

  // Volume de l'element
  double v_elem = (x_droite - x_gauche) * (y_haut - y_bas) * (z_avant - z_arriere);
  assert(v_elem>0.);

  const int nb_sommets_poly = poly_reelles.dimension(0);
  int i, i_prec = nb_sommets_poly -1;
  //calcule l'aire de la surface projetee sur face_bas
  double aire_projetee = 0.;
  for (i=0 ; i<nb_sommets_poly ; i++)
    {
      //calcul de la contribution de l'arete a l'aire projetee dans le plan (X,Z) :
      // par la somme algebrique des aires des trapezes (generes par la projection de l'arete sur l'axe x=0)
      double contrib_aire_projetee = (poly_reelles(i,2) - poly_reelles(i_prec,2)) * ((poly_reelles(i,0) + poly_reelles(i_prec,0)) * 0.5);

      aire_projetee += contrib_aire_projetee;
      i_prec = i;
    }
  double v = 0.;

  // Volume de la partie "projection sur la face face_bas"
  v = signe_princ * fabs(aire_projetee) * (centre_de_gravite[1] - y_bas);

  // Recherche d'un segment coupant la face du haut
  for (i = 0; i < nb_sommets_poly; i++)
    {
      if (polygone_plan_coupe[i] == NUM_FACE_HAUT_FACE)
        {
          const int i_precedent = (i-1) < 0 ? nb_sommets_poly-1 : (i-1);
          const int i_suivant   = (i+1) >= nb_sommets_poly ? 0 : (i+1);
          // Le segment [i-1,i] et sur la face du haut
          const int coupe_face_haut = i_precedent;
          const int coupe_face_haut_p1 = i;
          // Un segment voisin est-il sur la face de droite ?
          int coupe_face_droite = -1;
          if (polygone_plan_coupe[i_precedent] == NUM_FACE_DROITE_FACE)
            coupe_face_droite = i_precedent;
          else if (polygone_plan_coupe[i_suivant] == NUM_FACE_DROITE_FACE)
            coupe_face_droite = i;//_suivant;

          // Les sommets coupe_haut et coupe_haut_p1 coupent la face face_haut
          // -> ajoute la composante de volume associee = moyenne_x * (y_haut-y_bas) * diff_z
          double vol_compl0 =
            (0.5*(poly_reelles(coupe_face_haut,0) + poly_reelles(coupe_face_haut_p1,0)) - x_gauche)
            * (y_haut - y_bas)
            * (poly_reelles(coupe_face_haut,2) - poly_reelles(coupe_face_haut_p1,2));

          v += signe_compl0 * std::fabs(vol_compl0);

          if (coupe_face_droite >= 0)
            {
              //cette arete coupe aussi la face xmax
              //ajoute la composante de volume associee = (x_droite-x_gauche) * (y_haut-y_bas) * (z-z_arriere)
              double vol_compl1 =
                (x_droite - x_gauche)
                * (y_haut - y_bas)
                * (poly_reelles(coupe_face_droite,2) - z_arriere);

              v += signe_compl1 * std::fabs(vol_compl1);
            }
        }
    }

  //normalisation par le volume de l'element
  v /= v_elem;

  // On force la valeur entre 0 et 1 strictement.
#if 0
  // B.Math 10/09/2004: il faut corriger le calcul pour que les
  // contributions soient entre 0 et 1, ensuite on pourra faire ca:
  if (v < Erreur_relative_maxi_)
    v = Erreur_relative_maxi_;
  else if (v > 1. - Erreur_relative_maxi_)
    v = 1. - Erreur_relative_maxi_;
#endif

  return v;
}

// fin EB

// ==========================================================================
// ==========================================================================

// Description:
//  Pour un point P0 (x0, y0, z0) a l'INTERIEUR de l'element num_element
//  et un autre point P1 (x1, y1, z1), calcule l'intersection du segment (P0,P1)
//  avec les bords de l'element.
//  Si le point P1 est sur un bord de l'element (a epsilon pres), on considere
//  qu'il est a l'interieur et on ne reporte aucune intersection.
//  Si on trouve une intersection I, on met dans pos_intersection la coordonnee
//  barycentrique de l'intersection definie par I = (1-pos) * P0 + pos * P1
//  Si on ne trouve pas d'intersection, pos_intersection est inchange.
// Valeur de retour:
//  Si une intersection a ete trouvee, numero de la face de sortie dans la zone_vf
//  (peut servir d'index dans face_voisins par exemple).
//  Sinon, renvoie -1.
/*! @brief Pour un point P0 (x0, y0, z0) a l'INTERIEUR de l'element num_element et un autre point P1 (x1, y1, z1), calcule l'intersection du segment (P0,P1)
 *
 *   avec les bords de l'element.
 *   Si le point P1 est sur un bord de l'element (a epsilon pres), on considere
 *   qu'il est a l'interieur et on ne reporte aucune intersection.
 *   Si on trouve une intersection I, on met dans pos_intersection la coordonnee
 *   barycentrique de l'intersection definie par I = (1-pos) * P0 + pos * P1
 *   Si on ne trouve pas d'intersection, pos_intersection est inchange.
 *  Valeur de retour:
 *   Si une intersection a ete trouvee, numero de la face de sortie dans le domaine_vf
 *   (peut servir d'index dans face_voisins par exemple).
 *   Sinon, renvoie -1.
 *
 */
int Parcours_interface::calculer_face_sortie_element(const Domaine_VF& domaine_vf,
                                                     const int num_element,
                                                     double x0, double y0, double z0,
                                                     double x1, double y1, double z1,
                                                     double& pos_intersection) const
{
  const int n_faces = nb_faces_elem_;

  double t_sortie = 2.;
  // Numero de la face de sortie (numero entre 0 et n_faces)
  int face_sortie = -1;

  for (int face = 0; face < n_faces; face++)
    {
      double a, b, c, d;
      double signe = calcul_eq_plan(domaine_vf, num_element, face, a, b, c, d);

      // f est toujours positif si le point est dans l'element,
      // sinon il est negatif pour au moins une face de l'element.
      double f0 = (a * x0 + b * y0 + c * z0 + d) * signe;
      double f1 = (a * x1 + b * y1 + c * z1 + d) * signe;
      // Le point x0, y0, z0 doit etre a l'interieur de l'element
      assert(f0 > - Erreur_max_coordonnees_ * 2.);
      if (f0 < 0.) f0 = 0.;
      // On considere que le point d'arrivee est dans l'element s'il est a une
      // distance inferieure a Erreur_max_coordonnees.
      if (f1 < - Erreur_max_coordonnees_)
        {
          // Le point d'arrivee est a l'exterieur de l'element.
          // Calcul de la coordonnee barycentrique de l'intersection
          // telle que x|y|z_intersection = x0|y0|z0 * (1.-t) + x1|y1|z1 * t
          // double t = f0 / (f0-f1); // Fixed bug: Arithmetic exception
          double t = t_sortie;
          if (std::fabs(f0-f1)>=DMINFLOAT) t = f0 / (f0-f1);
          if (t < t_sortie)
            {
              t_sortie = t;
              face_sortie = face;
            }
        }
    }
  // Conversion du numero de la face de l'element en numero global de face
  if (face_sortie >= 0)
    {
      face_sortie = domaine_vf.elem_faces() (num_element, face_sortie);
      pos_intersection = t_sortie;
    }

  return face_sortie;
}

/*! @brief Methode outil de Maillage_FT_Disc::deplacer_un_point dans le cas d'un marqueur de la ligne de contact.
 *
 *  P0 est la position initiale du marqueur en contact (sur la face_bord)
 *  P1 est la position finale visee apres le deplacement
 *   Pour un point P0(x0,y0,z0) sur la face "face_bord" et un point P1(x1,y1,z1),
 *   on determine la projection orthogonale p(P1) de P1 sur le plan contenant la
 *   face, et on calcule l'intersection (x,y,z) du segment [P0,p(P1)] avec
 *   les bords de la face. S'il n'y a pas d'intersection, (x,y,z)=p(P1) et
 *   la valeur de retour est -1, sinon la valeur de retour est le numero de
 *   l'arete de la face qui est coupee par le segment [P0,p(P1)] (aretes
 *   telles qu'elles sont definies dans la classe Connectivite_frontieres).
 *
 * @param (domaine_vf) le Domaine_VF a laquelle se rapportent les indices de face et d'element
 * @param (face_0) le numero de la face reelle sur laquelle se trouve le point P0
 * @param (num_element) un numero d'element voisin de la face face_0
 * @param (x0, y0, z0) Coordonnees du point P0 (en 2D, z0 est ignore)
 * @param (x1, y1, z1) Coordonnees du point P1 (en 2D, z1 est ignore)
 * @param (x, y, z) On met dans (x,y,z) les coordonnees de l'intersection s'il y en a une, sinon on y met les coordonnees de p(P1) qui est sur la face "face_0" Valeur de retour: -1 s'il n'y a pas d'intersection sinon le numero de la face de bord ou passe le sommet
 */
int Parcours_interface::calculer_sortie_face_bord(const int face_0,
                                                  const int num_element,
                                                  double x0, double y0, double z0,
                                                  double x1, double y1, double z1,
                                                  double& x, double& y, double& z) const
{
  // Principe de l'algo :
  // * On projete le point d'arrivee sur la face_0
  // * On calcule l'intersection du segment avec l'element adjacent
  // * On cherche le numero de l'arete correspondant a la face de sortie
  {
    // Coefficients du plan contenant la face_0 :
    const double a = equations_plans_faces_(face_0, 0);
    const double b = equations_plans_faces_(face_0, 1);
    const double c = equations_plans_faces_(face_0, 2);
    const double d = equations_plans_faces_(face_0, 3);

    // On projete les points 0 et 1 sur le plan
    // En 2D, le coefficient c est nul, donc (z0,z1) peut etre quelconque.
    const double f0 = a * x0 + b * y0 + c * z0 + d;
    x0 -= f0 * a;
    y0 -= f0 * b;
    z0 -= f0 * c;
    const double f1 = a * x1 + b * y1 + c * z1 + d;
    x1 -= f1 * a;
    y1 -= f1 * b;
    z1 -= f1 * c;
  }

  // Calcul de l'intersection avec l'element
  const Domaine_VF& domaine_vf = refdomaine_vf_.valeur();
  // Verifie que num_element est un voisin de face_0:
  assert(domaine_vf.face_voisins()(face_0,0) == num_element
         || domaine_vf.face_voisins()(face_0,1) == num_element);
  double pos_intersection = 1.;
  int face_sortie = calculer_face_sortie_element(domaine_vf, num_element,
                                                 x0, y0, z0,
                                                 x1, y1, z1,
                                                 pos_intersection);

  // Coordonnees de l'intersection
  x = (1. - pos_intersection) * x0 + pos_intersection * x1;
  y = (1. - pos_intersection) * y0 + pos_intersection * y1;
  z = (1. - pos_intersection) * z0 + pos_intersection * z1;

  int face_bord_sortie;

  // Le test face_sortie != face_0 permet de gerer des erreurs d'arrondi
  // qui font qu'on sort par la mauvaise face (deja vu)
  if (face_sortie >= 0 && face_sortie != face_0)
    {
      const IntTab& face_sommets = domaine_vf.face_sommets();
      const int dimension2 = (Objet_U::dimension == 2);
      //
      // Calcul du numero de l'arete de bord coupee :
      // c'est l'arete dont les deux sommets sont des sommets
      // de la face "face".
      // ----------------------------------------------------
      assert(nb_sommets_par_face_ <= 4);
      // On met dans face_sortie_sommets les numeros des sommets de la face
      // de sortie (indices de sommets dans domaines_vf.domaine().domaine().les_sommets()).
      int face_sortie_sommets[4];
      int i;
      for (i = 0; i < nb_sommets_par_face_; i++)
        face_sortie_sommets[i] = face_sommets(face_sortie, i);

      // Tableau des numeros des sommets communs a face_0 et a face_sortie,
      // tries dans l'ordre croissant (si 2 sommets communs)
      // Ce sont des numeros locaux de sommets sur la face_0 :
      int sommet_commun[2] = { -1, -1 };
      for (i = 0; i < nb_sommets_par_face_; i++)
        {
          const int face_bord_sommet = face_sommets(face_0, i);
          int s = -1;
          int j;
          for (j = 0; j < nb_sommets_par_face_; j++)
            {
              s = face_sortie_sommets[j];
              if (s == face_bord_sommet)
                break;
            }
          if (s == face_bord_sommet)
            {
              if (sommet_commun[0] < 0)
                {
                  sommet_commun[0] = i;
                  if (dimension2)
                    break; // En dimension 2, il n'y a qu'un sommet commun
                }
              else
                {
                  sommet_commun[1] = i;
                  // tri dans l'ordre croissant
                  if (sommet_commun[0] > sommet_commun[1])
                    {
                      sommet_commun[1] = sommet_commun[0];
                      sommet_commun[0] = i;
                    }
                  break; // On a trouve les deux sommets communs
                }
            }
        }

      // On doit avoir trouve un sommet commun en 2d, 2 en 3d.
      // Test aussi en optimise, sinon ca calcule n'importe-quoi
      {
        int n;
        n = (sommet_commun[0] >= 0) ? 1 : 0;
        n += (sommet_commun[1] >= 0) ? 1 : 0;
        if (n != Objet_U::dimension - 1)
          {
            Cerr << "Erreur sur PE " << Process::me()
                 << " face de sortie non trouvee" << finl;
            exit();
          }
      }
      // On cherche a quelle arete correspondent les deux sommets
      // communs que l'on vient de trouver :
      const Connectivite_frontieres& connect_front = refconnect_front_.valeur();
      const IntTab& def_face_arete = connect_front.def_face_aretes();
      const int nb_aretes_face = def_face_arete.dimension(0);
      for (i = 0; i < nb_aretes_face; i++)
        {
          // Les deux sommets de l'arete a tester (numeros locaux sur la face_0)
          int sommet0 = def_face_arete(i, 0);
          int sommet1 = def_face_arete(i, 1);
          // ... tries dans l'ordre croissant
          if (sommet1 >= 0 && sommet1 < sommet0)
            {
              int s = sommet1;
              sommet1 = sommet0;
              sommet0 = s;
            }
          // Est-ce que c'est la bonne arete ?
          if (sommet0 == sommet_commun[0] && sommet1 == sommet_commun[1])
            break;
        }
      // On verifie qu'on a trouve l'arete commune
      if (i >= nb_aretes_face)
        {
          Cerr << "Erreur sur PE " << Process::me()
               << " face de sortie non trouvee (2)" << finl;
          exit();
        }
      // Calcul du numero de la face voisine de la face_0 par l'arete i :
      face_bord_sortie = connect_front.faces_voisins()(face_0, i);
    }
  else
    {
      // On n'a pas trouve d'intersection, (x1,y1,z1) est sur la face_0
      face_bord_sortie = -1;
    }
  return face_bord_sortie;
}

// Description :
//  Renvoie la plus grande distance (signee) entre le sommet (x,y,z) et
//  les plans qui definissent l'element REEL num_element.
//  Si cette distance est positive, le sommet est a l'exterieur de l'element,
//  sinon il est a l'interieur. Une majoration de l'erreur sur ce resultat
//  est donnee par get_erreur_geometrique.
double Parcours_interface::distance_sommet_faces(const Domaine_VF& domaine_vf,
                                                 const int num_element,
                                                 double x, double y, double z) const
{
  const int n_faces = nb_faces_elem_;
  double f_min = DMAXFLOAT;
  for (int face = 0; face < n_faces; face++)
    {
      double a, b, c, d;
      double signe = calcul_eq_plan(domaine_vf, num_element, face, a, b, c, d);
      const double f = (a * x + b * y + c * z + d) * signe;
      if (f < f_min)
        f_min = f;
    }
  return - f_min;
}

/*! @brief Renvoie une estimation de l'erreur geometrique (valeur homogene a une distance).
 *
 * Par exemple si le domaine a une dimension caracteristique de 1e-3
 *   et si la precision relative des calculs Erreur_relative_maxi_ est fixee a 1E-14
 *   (valeur fixee dans le code source) alors get_erreur_geometrique renvoie 1E-17.
 *   Cette valeur est identique sur tous les processeurs et est calculee a partir
 *   de la taille du domaine.
 *
 */
double Parcours_interface::get_erreur_geometrique() const
{
  assert(Erreur_max_coordonnees_ > 0.);
  return Erreur_max_coordonnees_;
}

/*! @brief Methode outil utilisee pour le traitement des lignes de contact.
 *
 * Projection du vecteur x_, y_, z_ sur le plan parallele a la face num_face passant
 *   par l'origine (permet d'obtenir la direction de deplacement d'un sommet sur une
 *   ligne de contact).
 *
 */
void Parcours_interface::projeter_vecteur_sur_face(const int num_face,
                                                   double& x_,
                                                   double& y_,
                                                   double& z_) const
{
  const double x = x_;
  const double y = y_;
  const double z = z_;
  const double nx = equations_plans_faces_(num_face, 0);
  const double ny = equations_plans_faces_(num_face, 1);
  const double nz = equations_plans_faces_(num_face, 2);
  const double a  = x * nx + y * ny + z * nz;
  x_ = x - a * nx;
  y_ = y - a * ny;
  z_ = z - a * nz;
}

// Normale unitaire a la face de bord num_face, au point x,y,z, dirigee vers
// l'interieur du domaine.
void Parcours_interface::calculer_normale_face_bord(int num_face, double x, double y, double z,
                                                    double& nx_, double& ny_, double& nz_) const
{
  const IntTab& face_voisins = refdomaine_vf_.valeur().face_voisins();
  double signe;
  if (face_voisins(num_face, 0) < 0)
    signe = 1.;
  else
    signe = -1.;
  nx_ = equations_plans_faces_(num_face, 0) * signe;
  ny_ = equations_plans_faces_(num_face, 1) * signe;
  nz_ = equations_plans_faces_(num_face, 2) * signe;
}

/*! @brief Algorithme base sur une version initiale de Thomas (recode par BM) Ramene le point (x,y,z) a l'interieur de l'element elem du domaine_vf a une
 *
 *   distance >= Erreur_max_coordonnees_ par un algorithme d'Uzawa.
 *  Valeur de retour: distance finale du sommet aux faces de l'element
 *   (positive si le sommet est a l'interieur)
 *
 */
double Parcours_interface::uzawa2(const Domaine_VF& domaine_vf,
                                  const int elem,
                                  double& x, double& y, double& z, const int is_face) const
{
  const int nb_faces = domaine_vf.elem_faces().dimension(1);
  ArrOfDouble coef_lagrange(nb_faces);
  int i;
  double somme_m_x = 0., somme_m_y = 0., somme_m_z = 0.;
  double norme_infty = 0.;
  for (i = 0; i < nb_faces; i++)
    {
      double a = 0., b = 0., c = 0., d = 0.;
      if (is_face) calcul_eq_plan_face(domaine_vf, elem, i, a, b, c, d);
      else calcul_eq_plan(domaine_vf, elem, i, a, b, c, d);
      a = std::fabs(a);
      b = std::fabs(b);
      c = std::fabs(c);
      somme_m_x += a;
      somme_m_y += b;
      somme_m_z += c;
      double somme = a + b + c;
      norme_infty = (somme > norme_infty) ? somme : norme_infty;
    }
  double norme1 = (somme_m_y > somme_m_x) ? somme_m_y : somme_m_x;
  norme1 = (somme_m_z > norme1) ? somme_m_z : norme1;
  const double rho = 1. / (norme1 * norme_infty);
  double solution_x = x;
  double solution_y = y;
  double solution_z = z;
  double dist_min;
  const int max_iter = 10;
  int niter;
  for (niter = 0; niter < max_iter; niter++)
    {
      dist_min = 1e37;
      int finished = 1;
      double nsol_x = 0.;
      double nsol_y = 0.;
      double nsol_z = 0.;
      for (i = 0; i < nb_faces; i++)
        {
          double a = 0., b = 0., c = 0., d = 0.;
          double s = is_face ? calcul_eq_plan_face(domaine_vf, elem, i, a, b, c, d) : calcul_eq_plan(domaine_vf, elem, i, a, b, c, d);
          double dist = (a * solution_x + b * solution_y + c * solution_z + d) * s;
          // Le nombre suivant est un peu au pif: il ne faut pas prendre une valeur
          // qui risque de refaire des erreurs (si on prend 1., ca plante en vdf car
          // on a des angles de 45 degres qui font qu'on retombe sur des cas particuliers)
          // Il faut un facteur significativement superieur a 1 pour que l'algorithme converge
          // en moins de 10 iterations.
          double r = 2.0389787 * std::fabs(a) + 3.07687 * std::fabs(b) + 2.528764 * std::fabs(c);
          double coef = coef_lagrange[i] - rho * (dist - r * Erreur_max_coordonnees_);
          coef = (coef > 0.) ? coef : 0.;
          nsol_x += a * s * coef;
          nsol_y += b * s * coef;
          nsol_z += c * s * coef;
          coef_lagrange[i] = coef;
          if (dist < Erreur_max_coordonnees_) // Le point est trop pres de parois, on n'a pas converge
            finished = 0;
          dist_min = (dist < dist_min) ? dist : dist_min;
        }
      if (finished)
        break;
      solution_x = nsol_x * 0.5 + x;
      solution_y = nsol_y * 0.5 + y;
      solution_z = nsol_z * 0.5 + z;
    }
  if (niter == max_iter)
    {
      Journal(8) << "Parcours_interface::uzawa2("
                 << x << " " << y << " " << z << ") non converge. dist_min=" << dist_min << finl;
    }
  x = solution_x;
  y = solution_y;
  z = solution_z;
  return dist_min;
}

/*! @brief Pour chaque sommet, s'il est trop pres d'une face eulerienne, deplace le sommet pour l'en eloigner.
 *
 * Mise a jour de l'espace virtuel des sommets
 *
 */
int Parcours_interface::eloigner_sommets_des_faces(Maillage_FT_Disc& maillage) const
{
  //maillage.update_sommet_face(); // EB
  const int nb_sommets = maillage.nb_sommets();
  DoubleTab& sommets = maillage.sommets_;
  const ArrOfInt& som_elem = maillage.sommet_elem_;
  //const IntTab& som_face = maillage.sommet_face_; // EB
  const int dim3 = (sommets.dimension(1) == 3);
  const Domaine_VF& domaine_vf = refdomaine_vf_.valeur();
  int count = 0;
  for (int i = 0; i < nb_sommets; i++)
    {
      if (maillage.sommet_virtuel(i))
        continue;
      // Pour les sommets sur les bords du domaine, on les rentre a l'interieur. Avec
      //  le traitement supplementaire (extension des triangles au dela du domaine,
      //  cela semble donner un resultat correct. A voir...
      //if (maillage.sommet_ligne_contact(i))
      //  continue;
      const int elem = som_elem[i];
      double x = sommets(i, 0);
      double y = sommets(i, 1);
      double z = dim3 ? sommets(i, 2) : 0.;
      const double d = distance_sommet_faces(domaine_vf, elem, x, y, z);
      if (d < Erreur_max_coordonnees_)
        {
          // On deplace un peu ce sommet:
          uzawa2(domaine_vf, elem, x, y, z);
          sommets(i, 0) = x;
          sommets(i, 1) = y;
          if (dim3)
            sommets(i, 2) = z;
          count++;
        }
      /*
      // debut EB
      for (int dim=0; dim <dimension; dim++)
        {
          const int face = som_face(i,dim);
          const double d_face = distance_sommet_faces_face(zone_vf, face, x, y, z);
          if (d_face < Erreur_max_coordonnees_)
            {
              // On deplace un peu ce sommet:
              uzawa2(zone_vf, face, x, y, z,1); // 1 : pour specifier que l'on deplace une face
              sommets(i, 0) = x;
              sommets(i, 1) = y;
              if (dim3)
                sommets(i, 2) = z;
              count++;
            }
        }
       */
      // fin EB

    }
  const int nb_som_tot_deplaces = mp_sum(count);
  if (nb_som_tot_deplaces > 0)
    {
      if (je_suis_maitre())
        Journal(5) << "Parcours_interface::eloigner_sommets_des_faces deplacement " << count << " sommets" << finl;
      Cerr<< "Parcours_interface::eloigner_sommets_des_faces deplacement " << count << " sommets" << finl;
      maillage.desc_sommets().echange_espace_virtuel(sommets);
    }
  return nb_som_tot_deplaces;
}
