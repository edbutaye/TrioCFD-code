/*
 * Domaine_VDF_plus.cpp
 *
 *  Created on: May 20, 2023
 *      Author: edouard
 */


#include <Domaine_VDF_plus.h>

#include <Domaine_Cl_VDF.h>
#include <Periodique.h>
#include <Dirichlet_entree_fluide_leaves.h>
#include <Neumann_sortie_libre.h>
#include <Scatter.h>
#include <Rectangle.h>
#include <Hexaedre.h>

Implemente_instanciable(Domaine_VDF_plus, "Domaine_VDF+", Domaine_VDF);

Sortie& Domaine_VDF_plus::printOn(Sortie& os) const
{
  Domaine_VDF::printOn(os);
  return os ;
}

//// readOn
//

Entree& Domaine_VDF_plus::readOn(Entree& is)
{
  Domaine_VDF::readOn(is);
  return is ;
}

void Domaine_VDF_plus::discretiser()
{
  Domaine_VF::discretiser();

  Domaine& domaine_geom=domaine();

  // Verification de la coherence entre l'element geometrique et la discretisation

  const Elem_geom_base& elem_geom = domaine_geom.type_elem().valeur();

  if (dimension == 2)
    {
      if (!sub_type(Rectangle,elem_geom))
        {
          Cerr << " Le type de l'element geometrique " << elem_geom.que_suis_je() << " est incorrect" << finl;
          Cerr << " Seul le type Rectangle et Rectangle_Axi sont compatibles avec la discretisation VDF en dimension 2" << finl;
          exit();
        }
    }
  else if (dimension == 3)
    {
      if (!sub_type(Hexaedre,elem_geom))
        {
          Cerr << " Le type de l'element geometrique " << elem_geom.que_suis_je() << " est incorrect" << finl;
          Cerr << " Seul le type Hexaedre ou Hexadre_Axi sont compatibles avec la discretisation VDF en dimension 3" << finl;
          exit();
        }
    }
  // Calcul de l'orientation des faces virtuelles
  // Note BM: pour compatibilite avec la version 1.5, orientation_
  //  n'a pas d'espace virtuel.
  MD_Vector md_nul;
  creer_tableau_faces(orientation_);
  orientation_.echange_espace_virtuel();
  orientation_.set_md_vector(md_nul); // Detache la structure parallele

  // Application de la convention VDF sur face_voisin:
  // L'element face_voisin(i,0) doit avoir une coordonnee "ori" plus petite que la face
  // et l'element face_voisin(i,1) doit avoir une coordonnee plus grande.
  {
    const int nbr_faces_tot = face_voisins_.dimension_tot(0);
    for (int i_face = 0; i_face < nbr_faces_tot; i_face++)
      {
        const int elem0 = face_voisins_(i_face, 0);
        const int elem1 = face_voisins_(i_face, 1);
        const int ori = orientation_[i_face];
        const double x_face = xv(i_face, ori);
        double delta = 0.;
        // L'element 0 doit avoir une coordonnee "ori" plus petite que la face
        // et l'element 1 doit avoir une coordonnee plus grande.
        if (elem0 >= 0)
          {
            delta = x_face - xp(elem0, ori);
          }
        else
          {
            assert(elem1 >= 0); // Sinon, grosse erreur: elements voisins non renseignes
            delta = xp(elem1, ori) - x_face;
          }
        if (delta < 0)
          {
            // On inverse les deux elements
            face_voisins_(i_face, 0) = elem1;
            face_voisins_(i_face, 1) = elem0;
          }
      }
  }

  calculer_volumes_entrelaces();
  genere_et_cree_aretes(); // EB
  calcul_xa(); // EB
  remplir_volumes_aretes(); // EB


  calcul_h();
  Cerr << "L'objet de type Domaine_VDF a ete rempli avec succes " << finl;

}

void Domaine_VDF_plus::modifier_pour_Cl(const Conds_lim& conds_lim)
{
  // Cerr << "Domaine_VDF::Modifier_pour_Cl" << finl;
  Domaine_VF::modifier_pour_Cl(conds_lim);

  int fac3;
  int nb_cond_lim=conds_lim.size();
  IntTab aretes_coin_traitees(nb_aretes_coin());
  aretes_coin_traitees = -1;

  for (int num_cond_lim=0; num_cond_lim<nb_cond_lim; num_cond_lim++)
    {
      // for cl ...
      //Journal() << "On traite la cl num : " << num_cond_lim << finl;
      const Cond_lim_base& cl = conds_lim[num_cond_lim].valeur();
      if (sub_type(Periodique, cl))
        {
          // if cl perio

          // Modification du tableau Qdm_ pour les aretes de type periodicite

          const Domaine_Cl_VDF& domaine_Cl_VDF = ref_cast(Domaine_Cl_VDF,cl.domaine_Cl_dis());

          int ndeb_arete = premiere_arete_bord();
          int fac1,fac2,sign,num_face,elem1,n_type;
          for (int n_arete=ndeb_arete; n_arete<ndeb_arete+nb_aretes_bord(); n_arete++)
            {
              // for n_arete
              n_type=domaine_Cl_VDF.type_arete_bord(n_arete-ndeb_arete);

              if (n_type == TypeAreteBordVDF::PERIO_PERIO) // arete de type periodicite
                {
                  //if arete perio
                  fac1 = Qdm_(n_arete,0);
                  fac2 = Qdm_(n_arete,1);
                  fac3 = Qdm_(n_arete,2);
                  sign = Qdm_(n_arete,3);

                  const Front_VF& la_frontiere_dis = ref_cast(Front_VF,cl.frontiere_dis());
                  int ndeb = la_frontiere_dis.num_premiere_face();
                  int nfin = ndeb + la_frontiere_dis.nb_faces();
                  if ( ( ( ndeb <= fac1) && (fac1 < nfin) ) ||
                       ( ( ndeb <= fac2) && (fac2 < nfin) )
                     )
                    {
                      if (sign == 1)
                        elem1 = face_voisins(fac1,1);
                      else
                        elem1 = face_voisins(fac1,0);
                      num_face = dimension+orientation_[fac3];

                      if (sign == -1)
                        {
                          Qdm_(n_arete,3) = fac3;
                          Qdm_(n_arete,2) = elem_faces(elem1,num_face);
                        }
                      else
                        Qdm_(n_arete,3) = elem_faces(elem1,num_face);
                    }
                }
            }
          // prise  en compte des aretes coin
          int ndeb_arete_coin = premiere_arete_coin();
          Cerr << "Modifications de Qdm pour les aretes coins  num_cond_lim=" << num_cond_lim <<  finl;
          int fac4;

          for (int n_arete=ndeb_arete_coin; n_arete<ndeb_arete_coin+nb_aretes_coin(); n_arete++)
            {
              if (aretes_coin_traitees[n_arete] != 1)
                {
                  fac1 = Qdm_(n_arete,0);
                  fac2 = Qdm_(n_arete,1);
                  fac3 = Qdm_(n_arete,2);
                  fac4 = Qdm_(n_arete,3);

                  // On recupere le numero des faces qui ne sont pas egales a -1
                  IntVect f(2);
                  f = -2;
                  int i=0;
                  if (fac1 != -1)
                    {
                      f(i) = fac1;
                      i++;
                    }
                  if (fac2 != -1)
                    {
                      f(i) = fac2;
                      i++;
                    }
                  if (fac3 != -1)
                    {
                      f(i) = fac3;
                      i++;
                    }
                  if (fac4 != -1)
                    {
                      f(i) = fac4;
                      i++;
                    }

                  // On regarde si la face est une face de periodicite

                  const Front_VF& la_frontiere_dis = ref_cast(Front_VF,cl.frontiere_dis());
                  int ndeb = la_frontiere_dis.num_premiere_face();
                  int nfin = ndeb + la_frontiere_dis.nb_faces();
                  int  elem, fac,dim0,dim1,indic_f0=-100,indic_f1=-100;

                  n_type=domaine_Cl_VDF.type_arete_coin(n_arete-ndeb_arete_coin);
                  dim1 = orientation_[f(1)];
                  dim0 = orientation_[f(0)];

                  for (int j=0; j<2; j++)
                    for (int k=0; k<2; k++)
                      if ((face_voisins(f(0),j) == face_voisins(f(1),k))  && (face_voisins(f(0),j)!=-1) )
                        {
                          indic_f0 = j;
                          indic_f1 = k;
                        }

                  if ((n_type == TypeAreteCoinVDF::PERIO_PAROI) || (n_type == TypeAreteCoinVDF::PERIO_FLUIDE))// arete coin perio-paroi
                    {
                      if ((f(0) >= ndeb)&&(f(0) < nfin))
                        {
                          Qdm_(n_arete,2)=f(0);
                          Qdm_(n_arete,indic_f0)=f(1);

                          elem = face_voisins(f(0),1-indic_f0);
                          fac = elem_faces(elem,dim1+(1-indic_f1)*dimension);
                          Qdm_(n_arete,1-indic_f0)=fac;

                          Qdm_(n_arete,3)=1-2*indic_f1;

                          //                          Cerr << "n_arete=" << n_arete << "  OK!!" << finl;
                          aretes_coin_traitees[n_arete] = 1;
                        }
                      else
                        {
                          if ((f(1) >= ndeb)&&(f(1) < nfin))
                            {
                              Qdm_(n_arete,2)=f(1);
                              Qdm_(n_arete,indic_f1)=f(0);

                              elem = face_voisins(f(1),1-indic_f1);
                              fac = elem_faces(elem,dim0+(1-indic_f0)*dimension);
                              Qdm_(n_arete,1-indic_f1)=fac;

                              Qdm_(n_arete,3)=1-2*indic_f0;

                              //                              Cerr << "n_arete=" << n_arete << "  OK!!" << finl;
                              aretes_coin_traitees[n_arete] = 1;
                            }
                          else
                            {
                              // Cerr<<"Attention cas non traite lors du remplissage du tableau Qdm"<<finl;
                              // exit();
                              ;
                            }
                        }
                    }
                  else if (n_type == TypeAreteCoinVDF::PERIO_PERIO) // arete coin perio-perio
                    {
                      if ((f(0) >= ndeb)&&(f(0) < nfin))
                        {
                          Qdm_(n_arete,2+indic_f1)=f(0);
                          Qdm_(n_arete,indic_f0)=f(1);

                          elem = face_voisins(f(0),1-indic_f0);
                          fac = elem_faces(elem,dim1+(1-indic_f1)*dimension);
                          Qdm_(n_arete,1-indic_f0)=fac;

                          elem = face_voisins(f(1),1-indic_f1);
                          fac = elem_faces(elem,dim0+(1-indic_f0)*dimension);
                          Qdm_(n_arete,3-indic_f1)=fac;

                          //                            Cerr << "n_arete=" << n_arete << "  OK!!" << finl;
                          aretes_coin_traitees[n_arete] = 1;
                        }
                      else
                        {
                          if ((f(1) >= ndeb)&&(f(1) < nfin))
                            {
                              Qdm_(n_arete,2+indic_f0)=f(1);
                              Qdm_(n_arete,indic_f1)=f(0);

                              elem = face_voisins(f(1),1-indic_f1);
                              fac = elem_faces(elem,dim0+(1-indic_f0)*dimension);
                              Qdm_(n_arete,1-indic_f1)=fac;

                              elem = face_voisins(f(0),1-indic_f0);
                              fac = elem_faces(elem,dim1+(1-indic_f1)*dimension);
                              Qdm_(n_arete,3-indic_f0)=fac;

                              //                                Cerr << "n_arete=" << n_arete << "  OK!!" << finl;
                              aretes_coin_traitees[n_arete] = 1;
                            }
                        }
                    }

                  else
                    {
                      Cerr << "Attention : les cas concernant d autres types d aretes coin "  << finl;
                      Cerr << "que perio-perio ou perio-paroi ne sont pas traites!!" << finl;
                    }
                }
            }
        }

      // Modif OC 01/2005 pour traiter les coins de type paroi/fluide
      else if (sub_type(Dirichlet_entree_fluide, cl) || sub_type(Neumann_sortie_libre, cl) )
        {
          // if cl de type paroi

          const Domaine_Cl_VDF& domaine_Cl_VDF = ref_cast(Domaine_Cl_VDF,cl.domaine_Cl_dis());

          int fac1,fac2,n_type;
          int ndeb_arete_coin = premiere_arete_coin();
          Cerr << "Modifications de Qdm pour les aretes coins touchant une paroi num_cond_lim=" << num_cond_lim <<  finl;
          int fac4;

          for (int n_arete=ndeb_arete_coin; n_arete<ndeb_arete_coin+nb_aretes_coin(); n_arete++)
            {
              if (aretes_coin_traitees[n_arete] != 1)
                {
                  fac1 = Qdm_(n_arete,0);
                  fac2 = Qdm_(n_arete,1);
                  fac3 = Qdm_(n_arete,2);
                  fac4 = Qdm_(n_arete,3);

                  // On recupere le numero des faces qui ne sont pas egales a -1
                  IntVect f(2);
                  f = -2;
                  int i=0;
                  if (fac1 != -1)
                    {
                      f(i) = fac1;
                      i++;
                    }
                  if (fac2 != -1)
                    {
                      f(i) = fac2;
                      i++;
                    }
                  if (fac3 != -1)
                    {
                      f(i) = fac3;
                      i++;
                    }
                  if (fac4 != -1)
                    {
                      f(i) = fac4;
                      i++;
                    }

                  // On regarde si la face est une face de type paroi

                  const Front_VF& la_frontiere_dis = ref_cast(Front_VF,cl.frontiere_dis());
                  int ndeb = la_frontiere_dis.num_premiere_face();
                  int nfin = ndeb + la_frontiere_dis.nb_faces();
                  int indic_f0=-100,indic_f1=-100;

                  n_type=domaine_Cl_VDF.type_arete_coin(n_arete-ndeb_arete_coin);

                  for (int j=0; j<2; j++)
                    for (int k=0; k<2; k++)
                      if ((face_voisins(f(0),j) == face_voisins(f(1),k)) && (face_voisins(f(0),j)!=-1) )
                        {
                          indic_f0 = j;
                          indic_f1 = k;
                        }

                  if ((n_type == TypeAreteCoinVDF::PAROI_FLUIDE) || (n_type == TypeAreteCoinVDF::FLUIDE_PAROI) || (n_type == TypeAreteCoinVDF::FLUIDE_FLUIDE))// arete coin paroi-fluide
                    {
                      if ((f(0) >= ndeb)&&(f(0) < nfin))
                        {
                          Qdm_(n_arete,0)=f(1);
                          Qdm_(n_arete,1)=f(1);
                          Qdm_(n_arete,2)=f(0);
                          Qdm_(n_arete,3)=1-2*indic_f1;
                          aretes_coin_traitees[n_arete] = 1;
                        }
                      else
                        {
                          if ((f(1) >= ndeb)&&(f(1) < nfin))
                            {
                              Qdm_(n_arete,0)=f(0);
                              Qdm_(n_arete,1)=f(0);
                              Qdm_(n_arete,2)=f(1);
                              Qdm_(n_arete,3)=1-2*indic_f0;
                              aretes_coin_traitees[n_arete] = 1;
                            }
                          else
                            {
                              // Cerr<<"Attention cas non traite lors du remplissage du tableau Qdm"<<finl;
                              // exit();
                              ;
                            }
                        }
                    }
                }
            }
        }

      //   Cerr << "Qdm : " << Qdm_ << finl;
      //   Cerr << "aretes_coin_traitees : " << aretes_coin_traitees << finl;
    }
  //Journal() << "le domaine apres modif pour Cl : " << *this << finl;


  // PQ : 10/10/05 : les faces periodiques etant a double contribution
  //                      l'appel a marquer_faces_double_contrib s'effectue dans cette methode
  //                      afin de pouvoir beneficier de conds_lim.

  Domaine_VF::marquer_faces_double_contrib(conds_lim);
// debut EB
  Domaine_VF& zvf = ref_cast(Domaine_VF,*this); // EB
  Joints& joints=zvf.domaine().faces_joint();
  const int nb_voisins = joints.size();
  //const DoubleTab& coord=domaine().domaine().coord_sommets();
  DoubleTab cg_aretes(1,dimension);
  ArrOfInt aretes(1);
  const IntTab& Aretes_Som=domaine().aretes_som();

  for (int voisin=0; voisin<nb_voisins; voisin++)
    {
      Joint& joint     = joints[voisin];
      //const int pe = joint.PEvoisin();
      ArrOfInt& items_communs = joint.set_joint_item(Joint::ARETE).set_items_communs(); // arete multiple
      items_communs.set_smart_resize(1);
      items_communs.resize(0);
      int compteur_aretes_joint=0;

      const ArrOfInt& indices_sommets =
        joint.joint_item(Joint::SOMMET).items_communs();

      for (int ind_S1=0; ind_S1< indices_sommets.size_array(); ind_S1++)
        {
          int S1 = indices_sommets(ind_S1);
          for (int ind_S2=ind_S1; ind_S2< indices_sommets.size_array(); ind_S2++)
            {
              int S2 = indices_sommets(ind_S2);
              compteur_aretes_joint=items_communs.size_array();
              for (int k=0; k<Aretes_Som.dimension(0); k++)
                if (Aretes_Som(k,0)==S1 && Aretes_Som(k,1)==S2)
                  {
                    items_communs.append_array(k);
                    compteur_aretes_joint++;
                  }
            }
        }

      items_communs.append_array(-1);
      items_communs.set_smart_resize(0);
      items_communs.resize(compteur_aretes_joint);
    }
  Scatter::calculer_renum_items_communs(joints,Joint::ARETE); // EB

  Domaine_VF::marquer_aretes_multiple_contrib(conds_lim); // EB

  IntTab tmp;
  tmp=Domaine_VF::set_arete_virt_pe_num(nb_aretes_reelles_, nb_aretes_tot_);
  Domaine_VF::construire_arete_virt_pe_num(nb_aretes_reelles_, nb_aretes_tot_,tmp); // EB

  // fin EB
}

static inline int face_vois_plus(const Domaine_VDF& dvdf, const Domaine& domaine, int face, int i)
{
  int elem=dvdf.face_voisins(face,i);
  if(elem==-1)
    {
      face=domaine.face_bords_interne_conjuguee(face);
      if(face!=-1)
        {
          elem=dvdf.face_voisins(face,i);
          assert(elem!=-1);
        }
    }
  return elem;
}

// EB
void Domaine_VDF_plus::genere_et_cree_aretes()
{
  Domaine& mon_domaine=domaine();

  const IntTab& elem_som = mon_domaine.les_elems();
  IntTab& Aretes_som = mon_domaine.set_aretes_som();
  IntTab& Elem_Aretes = mon_domaine.set_elem_aretes();
  Aretes_som.set_smart_resize(1);
  Aretes_som.resize(0, 2);

  const int nb_aretes_elem=12;
  Elem_Aretes.resize(0, nb_aretes_elem);
  mon_domaine.creer_tableau_elements(Elem_Aretes, Array_base::NOCOPY_NOINIT);
  Elem_Aretes=-1;
  //  int nb_poly = mon_domaine.nb_elem();
  int nb_poly_tot = mon_domaine.nb_elem_tot();
  // Estimation avec majoration du nombre d'aretes : nb_aretes_plus
  nb_aretes_=-1; // cf Aretes::affecter
  int nb_aretes_plus=-1;
  if (dimension == 2)
    nb_aretes_plus = mon_domaine.nb_som_tot();
  else if (dimension == 3)
    nb_aretes_plus = mon_domaine.nb_som_tot()*3;
  Cerr << "Creation des aretes et de la structure parallele" << finl;
  Aretes les_aretes(nb_aretes_plus);

  // On balaie les elements :
  int el1, el2, el3, el4;
  int face12, face13, face34, face24;
  int s1,s2;
  int nb_dir;
  if(dimension==2) nb_dir=1;
  else nb_dir=3;
  IntVect gauche(nb_dir);
  gauche(0)=0;
  if(dimension==3)
    {
      gauche(1)=0;
      gauche(2)=1;
    }
  IntVect droite(nb_dir);
  droite(0)=dimension;
  if(dimension==3)
    {
      droite(1)=dimension;
      droite(2)=dimension+1;
    }
  IntVect haut(nb_dir);
  haut(0)=dimension+1;
  if(dimension==3)
    {
      haut(1)=dimension+2;
      haut(2)=dimension+2;
    }
  IntVect bas(nb_dir);
  bas(0)=1;
  if(dimension==3)
    {
      bas(1)=2;
      bas(2)=2;
    }
  IntVect ind_som1(nb_dir);
  int ind_som2=-1;

  IntVect ind_som1_bas_droite(nb_dir);
  IntVect ind_som2_bas_droite(nb_dir);
  IntVect ind_som1_haut_gauche(nb_dir);
  IntVect ind_som2_haut_gauche(nb_dir);
  IntVect ind_som1_bas_gauche(nb_dir);
  IntVect ind_som2_bas_gauche(nb_dir);
  {
    assert(nb_dir==3);
    ind_som1(0)=3;
    ind_som1(1)=5;
    ind_som1(2)=6;
    ind_som2=7;

    ind_som1_bas_droite(0)=1;
    ind_som1_bas_droite(1)=1;
    ind_som1_bas_droite(2)=2;
    ind_som2_bas_droite(0)=5;
    ind_som2_bas_droite(1)=3;
    ind_som2_bas_droite(2)=3;

    ind_som1_haut_gauche(0)=2;
    ind_som1_haut_gauche(1)=4;
    ind_som1_haut_gauche(2)=4;
    ind_som2_haut_gauche(0)=6;
    ind_som2_haut_gauche(1)=6;
    ind_som2_haut_gauche(2)=5;

    ind_som1_bas_gauche(0)=0;
    ind_som1_bas_gauche(1)=0;
    ind_som1_bas_gauche(2)=0;
    ind_som2_bas_gauche(0)=4;
    ind_som2_bas_gauche(1)=2;
    ind_som2_bas_gauche(2)=1;
  }

  int coin=-1;
  int bord=0;
  int mixte=1;
  int interne=2;
  // Detection des plaques (2 faces frontieres se superposent):
  IntVect est_une_plaque(nb_faces());
  creer_tableau_faces(est_une_plaque);
  est_une_plaque=0;
  // Boucles sur les faces frontieres
  ArrOfDouble P1(3), P2(3);
  P1=0;
  P2=0;
  int ok=0;
  double eps = 10.0*Objet_U::precision_geom;
  for (int face=0; face<premiere_face_int(); face++)
    {
      int ori = orientation(face);
      for (int i = 0; i < dimension; i++)
        {
          P1[i] = xv(face, i) + (ori == i ? eps : 0);
          P2[i] = xv(face, i) - (ori == i ? eps : 0);
        }
      int elem1 = domaine().chercher_elements(P1[0], P1[1], P1[2]);
      int elem2 = domaine().chercher_elements(P2[0], P2[1], P2[2]);
      if (elem1 >= 0 && elem2 >= 0 && elem1 != elem2)
        {
          // Find 2 points P1 and P2 in cells so:
          est_une_plaque(face) = 1;
          Journal() << "We detect an internal boundary on face " << face << " between elements " << elem1 << " and "
                    << elem2 << finl;
        }
    }
  est_une_plaque.echange_espace_virtuel();

  // Dans cette premiere boucle, on remplit Aretes_som et les aretes 0-2 de Elem_Aretes
  for(int dir=0; dir<nb_dir; dir++)
    for (el1=0; el1<nb_poly_tot; el1++)
      {
        // On doit generer l'arete en haut a droite de el1
        face12=elem_faces(el1,droite(dir));
        face13=elem_faces(el1,haut(dir));

        // sommets de cette aretes
        s1 = elem_som(el1,ind_som1(dir));
        s2 = elem_som(el1,ind_som2);

        el2=face_vois_plus(*this, mon_domaine, face12, 1);

        if(el2>-1)
          face24=elem_faces(el2,haut(dir));
        else
          face24=-1;

        el3=face_vois_plus(*this, mon_domaine, face13, 1);

        if(el3>-1)
          face34=elem_faces(el3,droite(dir));
        else
          face34=-1;

        if( (el2>-1) && (el3>-1))
          el4=face_vois_plus(*this, mon_domaine, face24, 1);
        else if((el2>-1))
          el4=face_vois_plus(*this, mon_domaine, face24, 1);
        else if((el3>-1))
          el4=face_vois_plus(*this, mon_domaine, face34, 1);
        else
          el4=-1;

        if ( (el2==-1) && (el4>=0) )
          face24=elem_faces(el4,bas(dir));
        if ( (el3==-1) && (el4>=0) )
          face34=elem_faces(el4,gauche(dir));

        const int nb_f = nb_faces();
        if (el2 > -1 && el3 > -1 && el4 > -1) // arete interne
          ok=les_aretes.affecter_aretes(nb_aretes_, dir, interne, nb_f, face13, face24, face12, face34, est_une_plaque);
        else if ( (el3 > -1 && el4 > -1) ||
                  (el2 > -1 && el4 > -1) ||
                  (el2 > -1 && el3 > -1) ) // arete mixte
          ok=les_aretes.affecter_aretes(nb_aretes_, dir, mixte, nb_f, face13, face24, face12, face34, est_une_plaque);
        else if (el2 > -1) // arete bord
          ok=les_aretes.affecter_aretes(nb_aretes_, dir, bord, nb_f, face13, face24, face12, 1, est_une_plaque);
        else if (el3 > -1) // arete bord
          ok=les_aretes.affecter_aretes(nb_aretes_, dir, bord, nb_f, face12, face34, face13, 1, est_une_plaque);
        else // arete coin
          ok=les_aretes.affecter_aretes(nb_aretes_, dir, coin, nb_f, face13, face24, face12, face34, est_une_plaque);

        if (ok)
          {
            Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
            Elem_Aretes(el1,dir)=nb_aretes_;
            //Cerr << "elements_haut_droit " << el1 << " " << el2 << " " << el3 << " " << el4 << finl;
            //Journal() << "Provisoire faces arete: " << face12 << " " << face13 << " " << face24 << " " << face34 << finl;
          }

        ok=0;
        // Pour les coins ou bords :
        face13 = elem_faces(el1, bas(dir));
        el3 = face_vois_plus(*this, mon_domaine, face13, 0);
        if (el3 < 0) // On doit generer l'arete en bas a droite de el1
          // si la maille en bas de el1 n'existe pas
          {
            if (el2 >= 0)
              {
                face24 = elem_faces(el2, bas(dir));
                el4 = face_vois_plus(*this, mon_domaine, face24, 0);
                if (el4 < 0) // arete bord
                  ok=les_aretes.affecter_aretes(nb_aretes_, dir, bord, nb_f, face13, face24, face12, -1, est_une_plaque);
                else // arete mixte
                  {
                    face34 = elem_faces(el4, gauche(dir));
                    if (face_vois_plus(*this, mon_domaine, face34, 0) == -1)
                      ok=les_aretes.affecter_aretes(nb_aretes_, dir, mixte, nb_f, face13, face24, face34, face12, est_une_plaque);
                  }
              }
            else   // arete coin
              ok=les_aretes.affecter_aretes(nb_aretes_, dir, coin, nb_f, face13, -1, -1, face12, est_une_plaque);


            if (ok)
              {
                s1 = elem_som(el1,ind_som1_bas_droite(dir));
                s2 = elem_som(el1,ind_som2_bas_droite(dir));
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,6+dir)=nb_aretes_; // dir=0 -> arete 6, dir=1 -> arete 7, dir=2 -> arete 8
              }
            if (el2 >= 0)
              {
                if (dir==0) Elem_Aretes(el2,3)=nb_aretes_;
                if (dir==1) Elem_Aretes(el2,4)=nb_aretes_;
                if (dir==2) Elem_Aretes(el2,5)=nb_aretes_;
              }
            if (el4 >= 0)
              {
                if (dir==0) Elem_Aretes(el4,9)=nb_aretes_;
                if (dir==1) Elem_Aretes(el4,10)=nb_aretes_;
                if (dir==2) Elem_Aretes(el4,11)=nb_aretes_;
              }
          }

        ok=0;
        face13 = elem_faces(el1, gauche(dir));
        el3 = face_vois_plus(*this, mon_domaine, face13, 0);
        if (el3 < 0) // On doit generer l'arete en haut a gauche de el1
          // si la maille en haut a gauche n'existe pas
          {
            face12 = elem_faces(el1, haut(dir));
            el2 = face_vois_plus(*this, mon_domaine, face12, 1);
            if (el2 >= 0)
              {
                face24 = elem_faces(el2, gauche(dir));
                el4 = face_vois_plus(*this, mon_domaine, face24, 0);
                if (el4 < 0) // arete bord
                  ok=les_aretes.affecter_aretes(nb_aretes_, dir, bord, nb_f, face13, face24, face12, -1, est_une_plaque);
              }
            else   // arete coin
              ok=les_aretes.affecter_aretes(nb_aretes_, dir, coin, nb_f, -1, face12, -1, face13, est_une_plaque);

            if (ok)
              {
                s1 = elem_som(el1,ind_som1_haut_gauche(dir));
                s2 = elem_som(el1,ind_som2_haut_gauche(dir));
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,9+dir)=nb_aretes_; // dir=0 -> arete 9, dir=1 -> arete 10, dir=2 -> arete 11
              }

            if (el2 >= 0)
              {
                if (dir==0) Elem_Aretes(el2,3)=nb_aretes_;
                if (dir==1) Elem_Aretes(el2,4)=nb_aretes_;
                if (dir==2) Elem_Aretes(el2,5)=nb_aretes_;
              }
            if (el4 >= 0)
              {
                if (dir==0) Elem_Aretes(el4,6)=nb_aretes_;
                if (dir==1) Elem_Aretes(el4,7)=nb_aretes_;
                if (dir==2) Elem_Aretes(el4,8)=nb_aretes_;
              }
          }

        ok=0;
        // On doit generer l'arete en bas a gauche de el1
        face12 = elem_faces(el1, gauche(dir));
        el2 = face_vois_plus(*this, mon_domaine, face12, 0);
        face13 = elem_faces(el1, bas(dir));
        el3 = face_vois_plus(*this, mon_domaine, face13, 0);
        if ((el2 < 0) && (el3 < 0)) // arete coin
          {
            ok=les_aretes.affecter_aretes(nb_aretes_, dir, coin, nb_f, face13, -1, face12, -1, est_une_plaque);
            if (ok)
              {
                s1 = elem_som(el1,ind_som1_bas_gauche(dir));
                s2 = elem_som(el1,ind_som2_bas_gauche(dir));
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,3+dir)=nb_aretes_; // dir=0 -> arete 9, dir=1 -> arete 10, dir=2 -> arete 11
              }
          }
      }

  // On complete Elem_Aretes

  for (int elem1=0; elem1<nb_poly_tot; elem1++)
    {
      int face_bas=elem_faces(elem1,2);
      int elem_bas=face_vois_plus(*this, mon_domaine, face_bas, 0);
      if (elem_bas>=0)
        {
          if (Elem_Aretes(elem1,8)<0) Elem_Aretes(elem1,8)=Elem_Aretes(elem_bas,2);
          if (Elem_Aretes(elem1,7)<0) Elem_Aretes(elem1,7)=Elem_Aretes(elem_bas,1);
        }

      int face_gauche= elem_faces(elem1,0);
      int elem_gauche=face_vois_plus(*this,mon_domaine, face_gauche, 0);
      if (elem_gauche>=0)
        {
          if (Elem_Aretes(elem1,9)<0) Elem_Aretes(elem1,9)=Elem_Aretes(elem_gauche,0);
          if (Elem_Aretes(elem1,10)<0) Elem_Aretes(elem1,10)=Elem_Aretes(elem_gauche,1);

          int face_avant_elem_gauche= elem_faces(elem_gauche,1);
          int elem_gauche_avant=face_vois_plus(*this, mon_domaine, face_avant_elem_gauche, 0);
          if (elem_gauche_avant>=0 && Elem_Aretes(elem1,3)<0) Elem_Aretes(elem1,3)=Elem_Aretes(elem_gauche_avant,0);

          int face_bas_elem_gauche= elem_faces(elem_gauche,2);
          int elem_gauche_bas=face_vois_plus(*this, mon_domaine, face_bas_elem_gauche, 0);
          if (elem_gauche_bas>=0 && Elem_Aretes(elem1,4)<0) Elem_Aretes(elem1,4)=Elem_Aretes(elem_gauche_bas,1);
        }

      int face_avant= elem_faces(elem1,1);
      int elem_avant=face_vois_plus(*this,mon_domaine, face_avant, 0);
      if (elem_avant>=0)
        {
          if (Elem_Aretes(elem1,6)<0) Elem_Aretes(elem1,6)=Elem_Aretes(elem_avant,0);
          if (Elem_Aretes(elem1,11)<0) Elem_Aretes(elem1,11)=Elem_Aretes(elem_avant,2);

          int face_bas_elem_avant=elem_faces(elem_avant,2);
          int elem_avant_bas=face_vois_plus(*this,mon_domaine,face_bas_elem_avant,0);
          if (elem_avant_bas>=0 && Elem_Aretes(elem1,5)<0) Elem_Aretes(elem1,5)=Elem_Aretes(elem_avant_bas,2);
        }

    }

  nb_aretes_reelles_=nb_aretes_+1;

  nb_aretes_tot_=nb_aretes_;
  cree_aretes_virtuelles(Aretes_som,Elem_Aretes,les_aretes);

  //Elem_Aretes.echange_espace_virtuel();
  nb_aretes_tot_++;
  nb_aretes_++;
  les_aretes.dimensionner(nb_aretes_tot_);
  Aretes_som.resize(nb_aretes_tot_,2);
  // Ajuste la taille du tableau Aretes_som
  const int n_aretes_tot = Aretes_som.dimension(0);
  Aretes_som.append_line(-1, -1);
  Aretes_som.set_smart_resize(0);
  Aretes_som.resize(n_aretes_tot, 2);

  // Cerr << "Tri des aretes " << finl;
  les_aretes.trier(nb_aretes_coin_,
                   nb_aretes_bord_,
                   nb_aretes_mixtes_,
                   nb_aretes_internes_, nb_aretes_reelles_,
                   nb_elem(), Aretes_som,Elem_Aretes);

  domaine().creer_structure_parallelle_aretes(nb_aretes_reelles_,Aretes_som,Elem_Aretes);
  md_vector_aretes_ = domaine().aretes_som().get_md_vector();

#ifdef SORT_POUR_DEBOG
  les_aretes.trier_pour_debog(nb_aretes_coin_,
                              nb_aretes_bord_,
                              nb_aretes_mixtes_,
                              nb_aretes_internes_,xv());
#endif
  nb_aretes_joint_=0;
  assert(nb_aretes_==nb_aretes_coin_+nb_aretes_bord_+nb_aretes_mixtes_
         +nb_aretes_internes_+nb_aretes_joint_);
  Qdm_.ref(les_aretes.faces());
  orientation_aretes_.ref(les_aretes.type1());
  type_arete_.ref(les_aretes.type2());


  //Elem_Aretes.echange_espace_virtuel();

}

// debut EB
void Domaine_VDF_plus::remplir_volumes_aretes()
{
  const int nb_aretes=domaine().nb_aretes();
  double y_bas=0,y_haut=0,x_gauche=0,x_droite=0,z_arriere=0,z_avant=0;
  int ori_arete,face1,face2,face3,face4;
  int s1,s2;
  const IntTab& Aretes_Som=domaine().aretes_som();
  const DoubleTab& coord_som=domaine().coord_sommets();
  const IntVect& types_aretes=type_arete();
  volumes_aretes_.resize(nb_aretes);
  int type_arete;
  for (int arete=0; arete<nb_aretes; arete++)
    {
      type_arete=types_aretes(arete);
      ori_arete=dimension-orientation_aretes()(arete)-1;
      face1=Qdm(arete,0);
      face2=Qdm(arete,1);
      face3=Qdm(arete,2);
      face4=Qdm(arete,3);
      if (type_arete==2)
        {
          s1=Aretes_Som(arete,0);
          s2=Aretes_Som(arete,1);

          assert(face1>=0 && face2>=0 && face3>=0 && face4>=0);

          if (ori_arete==0)
            {
              x_gauche=coord_som(s1,0);
              x_droite=coord_som(s2,0);
              y_haut=xv(face2,1);
              y_bas=xv(face1,1);
              z_avant=xv(face4,2);
              z_arriere=xv(face3,2);
            }
          else if (ori_arete==1)
            {
              y_bas=coord_som(s1,ori_arete);
              y_haut=coord_som(s2,ori_arete);
              x_gauche=xv(face1,0);
              x_droite=xv(face2,0);
              z_avant=xv(face4,2);
              z_arriere=xv(face3,2);
            }
          if (ori_arete==2)
            {
              z_arriere=coord_som(s1,ori_arete);
              z_avant=coord_som(s2,ori_arete);
              x_gauche=xv(face1,0);
              x_droite=xv(face2,0);
              y_bas=xv(face3,1);
              y_haut=xv(face4,1);
            }

          volumes_aretes_(arete)= (x_droite - x_gauche) * (y_haut - y_bas) * (z_avant - z_arriere);
        }
      else if (type_arete==1) volumes_aretes_(arete)=volumes_entrelaces(face3)/2; // bord

      else // coin, mixte (rip pour les aretes mixtes)
        {
          int elem1=(face4>=0) ? face_voisins(face4,0) : -1;
          int elem2=(face4>=0) ? face_voisins(face4,1) : -1;
          int elem3=(face3>=0) ? face_voisins(face3,0) : -1;
          int elem4=(face3>=0) ? face_voisins(face3,1) : -1;

          if (elem1>=0) volumes_aretes_(arete)=volumes(elem1)/4;
          else if (elem2>=0) volumes_aretes_(arete)=volumes(elem2)/4;
          else if (elem3>=0) volumes_aretes_(arete)=volumes(elem3)/4;
          else if (elem4>=0) volumes_aretes_(arete)=volumes(elem4)/4;
          else
            exit();
        }
    }



}
// fin EB
// Description:
// calcul du centre de gravite des aretes
void Domaine_VDF_plus::calcul_xa()
{
  const Domaine& dom = domaine();
  // Calcul des centres de gravite des aretes xa_ stockes dans la Domaine_VF
  const IntTab& aretes_som = domaine().aretes_som();
  const int nb_aretes_tot = aretes_som.dimension_tot(0);
  const DoubleTab& coord = dom.les_sommets();
  const int dim = coord.dimension(1);
  xa_.resize(nb_aretes_tot, dim);
  //creer_tableau_aretes(xa_, ArrOfDouble::NOCOPY_NOINIT);
  for (int i = 0; i < nb_aretes_tot; i++)
    {
      const int s0 = aretes_som(i, 0);
      const int s1 = aretes_som(i, 1);
      for (int j = 0; j < dim; j++)
        xa_(i, j) = (coord(s0, j) + coord(s1, j)) * 0.5;
    }
  //xa_.echange_espace_virtuel();
}
// fin EB


void Domaine_VDF_plus::cree_aretes_virtuelles(IntTab& Aretes_som, IntTab& Elem_Aretes, Aretes& les_aretes)
{
  Domaine& mon_domaine=domaine();
  const IntTab& elem_som = mon_domaine.les_elems();
  //const int nb_aretes_elem=12;

  int nb_poly = mon_domaine.nb_elem();
  int nb_poly_tot = mon_domaine.nb_elem_tot();
  // On balaie les elements :
  int el1, el2, el3, el4;
  int face12, face13, face34, face24;
  int s1,s2;
  int nb_dir;
  if(dimension==2) nb_dir=1;
  else nb_dir=3;
  IntVect gauche(nb_dir);
  gauche(0)=0;
  if(dimension==3)
    {
      gauche(1)=0;
      gauche(2)=1;
    }
  IntVect droite(nb_dir);
  droite(0)=dimension;
  if(dimension==3)
    {
      droite(1)=dimension;
      droite(2)=dimension+1;
    }
  IntVect haut(nb_dir);
  haut(0)=dimension+1;
  if(dimension==3)
    {
      haut(1)=dimension+2;
      haut(2)=dimension+2;
    }
  IntVect bas(nb_dir);
  bas(0)=1;
  if(dimension==3)
    {
      bas(1)=2;
      bas(2)=2;
    }
  IntVect ind_som1(nb_dir);
  int ind_som2=-1;

  IntVect ind_som1_bas_droite(nb_dir);
  IntVect ind_som2_bas_droite(nb_dir);
  IntVect ind_som1_haut_gauche(nb_dir);
  IntVect ind_som2_haut_gauche(nb_dir);
  IntVect ind_som1_bas_gauche(nb_dir);
  IntVect ind_som2_bas_gauche(nb_dir);
  {
    assert(nb_dir==3);
    ind_som1(0)=3;
    ind_som1(1)=5;
    ind_som1(2)=6;
    ind_som2=7;

    ind_som1_bas_droite(0)=1;
    ind_som1_bas_droite(1)=1;
    ind_som1_bas_droite(2)=2;
    ind_som2_bas_droite(0)=5;
    ind_som2_bas_droite(1)=3;
    ind_som2_bas_droite(2)=3;

    ind_som1_haut_gauche(0)=2;
    ind_som1_haut_gauche(1)=4;
    ind_som1_haut_gauche(2)=4;
    ind_som2_haut_gauche(0)=6;
    ind_som2_haut_gauche(1)=6;
    ind_som2_haut_gauche(2)=5;

    ind_som1_bas_gauche(0)=0;
    ind_som1_bas_gauche(1)=0;
    ind_som1_bas_gauche(2)=0;
    ind_som2_bas_gauche(0)=4;
    ind_som2_bas_gauche(1)=2;
    ind_som2_bas_gauche(2)=1;
  }

  int coin=-1;
  int bord=0;
  int mixte=1;
  int interne=2;
  int ok=0;
  // Dans cette premiere boucle, on remplit Aretes_som et les aretes 0-2 de Elem_Aretes
  for(int dir=0; dir<nb_dir; dir++)
    for (el1=nb_poly; el1<nb_poly_tot; el1++)
      {

        // On doit generer l'arete en haut a droite de el1
        face12=elem_faces(el1,droite(dir));
        face13=elem_faces(el1,haut(dir));

        // sommets de cette aretes
        s1 = elem_som(el1,ind_som1(dir));
        s2 = elem_som(el1,ind_som2);

        el2=face_vois_plus(*this, mon_domaine, face12, 1);

        if(el2>-1)
          face24=elem_faces(el2,haut(dir));
        else
          face24=-1;

        el3=face_vois_plus(*this, mon_domaine, face13, 1);

        if(el3>-1)
          face34=elem_faces(el3,droite(dir));
        else
          face34=-1;

        if( (el2>-1) && (el3>-1))
          el4=face_vois_plus(*this, mon_domaine, face24, 1);
        else if((el2>-1))
          el4=face_vois_plus(*this, mon_domaine, face24, 1);
        else if((el3>-1))
          el4=face_vois_plus(*this, mon_domaine, face34, 1);
        else
          el4=-1;

        if ( (el2==-1) && (el4>=0) )
          face24=elem_faces(el4,bas(dir));
        if ( (el3==-1) && (el4>=0) )
          face34=elem_faces(el4,gauche(dir));


        const int nb_f = nb_faces();
        if (Elem_Aretes(el1,dir)<0)
          {
            if (el2 > -1 && el3 > -1 && el4 > -1) // arete interne
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, interne, nb_f, face13, face24, face12, face34);
            else if ( (el3 > -1 && el4 > -1) ||
                      (el2 > -1 && el4 > -1) ||
                      (el2 > -1 && el3 > -1) ) // arete mixte
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, mixte, nb_f, face13, face24, face12, face34);
            else if (el2 > -1) // arete bord
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, bord, nb_f, face13, face24, face12, 1);
            else if (el3 > -1) // arete bord
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, bord, nb_f, face12, face34, face13, 1);
            else // arete coin
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, coin, nb_f, face13, face24, face12, face34);

            if (ok)
              {
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,dir)=nb_aretes_tot_;
                //Cerr << "elements_haut_droit " << el1 << " " << el2 << " " << el3 << " " << el4 << finl;
                //Journal() << "Provisoire faces arete: " << face12 << " " << face13 << " " << face24 << " " << face34 << finl;
              }
          }
        ok=0;
        // Pour les coins ou bords :
        face13 = elem_faces(el1, bas(dir));
        el3 = face_vois_plus(*this, mon_domaine, face13, 0);
        if (el3 < 0 && Elem_Aretes(el1,6+dir)<0) // On doit generer l'arete en bas a droite de el1
          // si la maille en bas de el1 n'existe pas
          {
            if (el2 >= 0)
              {
                face24 = elem_faces(el2, bas(dir));
                el4 = face_vois_plus(*this, mon_domaine, face24, 0);
                if (el4 < 0) // arete bord
                  ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, bord, nb_f, face13, face24, face12, -1);
                else // arete mixte
                  {
                    face34 = elem_faces(el4, gauche(dir));
                    if (face_vois_plus(*this, mon_domaine, face34, 0) == -1)
                      ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, mixte, nb_f, face13, face24, face34, face12);
                  }
              }
            else   // arete coin
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, coin, nb_f, face13, -1, -1, face12);


            if (ok)
              {
                s1 = elem_som(el1,ind_som1_bas_droite(dir));
                s2 = elem_som(el1,ind_som2_bas_droite(dir));
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,6+dir)=nb_aretes_tot_; // dir=0 -> arete 6, dir=1 -> arete 7, dir=2 -> arete 8
              }

            if (el2 >= 0)
              {
                if (dir==0) Elem_Aretes(el2,3)=nb_aretes_tot_;
                if (dir==1) Elem_Aretes(el2,4)=nb_aretes_tot_;
                if (dir==2) Elem_Aretes(el2,5)=nb_aretes_tot_;
              }
            if (el4 >= 0)
              {
                if (dir==0) Elem_Aretes(el4,9)=nb_aretes_tot_;
                if (dir==1) Elem_Aretes(el4,10)=nb_aretes_tot_;
                if (dir==2) Elem_Aretes(el4,11)=nb_aretes_tot_;
              }

          }

        ok=0;
        face13 = elem_faces(el1, gauche(dir));
        el3 = face_vois_plus(*this, mon_domaine, face13, 0);
        if (el3 < 0 && Elem_Aretes(el1,9+dir)<0) // On doit generer l'arete en haut a gauche de el1
          // si la maille en haut a gauche n'existe pas
          {
            face12 = elem_faces(el1, haut(dir));
            el2 = face_vois_plus(*this, mon_domaine, face12, 1);
            if (el2 >= 0)
              {
                face24 = elem_faces(el2, gauche(dir));
                el4 = face_vois_plus(*this, mon_domaine, face24, 0);
                if (el4 < 0) // arete bord
                  ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, bord, nb_f, face13, face24, face12, -1);
              }
            else   // arete coin
              ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, coin, nb_f, -1, face12, -1, face13);

            if (ok)
              {
                s1 = elem_som(el1,ind_som1_haut_gauche(dir));
                s2 = elem_som(el1,ind_som2_haut_gauche(dir));
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,9+dir)=nb_aretes_tot_; // dir=0 -> arete 9, dir=1 -> arete 10, dir=2 -> arete 11
              }

            if (el2 >= 0)
              {
                if (dir==0) Elem_Aretes(el2,3)=nb_aretes_tot_;
                if (dir==1) Elem_Aretes(el2,4)=nb_aretes_tot_;
                if (dir==2) Elem_Aretes(el2,5)=nb_aretes_tot_;
              }
            if (el4 >= 0)
              {
                if (dir==0) Elem_Aretes(el4,6)=nb_aretes_tot_;
                if (dir==1) Elem_Aretes(el4,7)=nb_aretes_tot_;
                if (dir==2) Elem_Aretes(el4,8)=nb_aretes_tot_;
              }

          }

        ok=0;
        // On doit generer l'arete en bas a gauche de el1
        face12 = elem_faces(el1, gauche(dir));
        el2 = face_vois_plus(*this, mon_domaine, face12, 0);
        face13 = elem_faces(el1, bas(dir));
        el3 = face_vois_plus(*this, mon_domaine, face13, 0);
        if ((el2 < 0) && (el3 < 0) &&  Elem_Aretes(el1,3+dir)<0) // arete coin
          {
            ok=les_aretes.affecter_aretes_virtuelle(nb_aretes_tot_, dir, coin, nb_f, face13, -1, face12, -1);
            if (ok)
              {
                s1 = elem_som(el1,ind_som1_bas_gauche(dir));
                s2 = elem_som(el1,ind_som2_bas_gauche(dir));
                Aretes_som.append_line(std::min(s1,s2),std::max(s1,s2)); // On ajoute l'arete
                Elem_Aretes(el1,3+dir)=nb_aretes_tot_; // dir=0 -> arete 9, dir=1 -> arete 10, dir=2 -> arete 11
              }
          }

      }

  // On complete Elem_Aretes
  for (int elem1=nb_elem(); elem1<nb_poly_tot; elem1++)
    {
      int face_bas=elem_faces(elem1,2);
      int elem_bas=face_vois_plus(*this, mon_domaine, face_bas, 0);
      if (elem_bas>=0)
        {
          if (Elem_Aretes(elem1,8)<0) Elem_Aretes(elem1,8)=Elem_Aretes(elem_bas,2);
          if (Elem_Aretes(elem1,7)<0) Elem_Aretes(elem1,7)=Elem_Aretes(elem_bas,1);
        }

      int face_gauche= elem_faces(elem1,0);
      int elem_gauche=face_vois_plus(*this,mon_domaine, face_gauche, 0);
      if (elem_gauche>=0)
        {
          if (Elem_Aretes(elem1,9)<0) Elem_Aretes(elem1,9)=Elem_Aretes(elem_gauche,0);
          if (Elem_Aretes(elem1,10)<0) Elem_Aretes(elem1,10)=Elem_Aretes(elem_gauche,1);

          int face_avant_elem_gauche= elem_faces(elem_gauche,1);
          int elem_gauche_avant=face_vois_plus(*this,mon_domaine, face_avant_elem_gauche, 0);
          if (elem_gauche_avant>=0 && Elem_Aretes(elem1,3)<0) Elem_Aretes(elem1,3)=Elem_Aretes(elem_gauche_avant,0);

          int face_bas_elem_gauche= elem_faces(elem_gauche,2);
          int elem_gauche_bas=face_vois_plus(*this,mon_domaine, face_bas_elem_gauche, 0);
          if (elem_gauche_bas>=0 && Elem_Aretes(elem1,4)<0) Elem_Aretes(elem1,4)=Elem_Aretes(elem_gauche_bas,1);
        }

      int face_avant= elem_faces(elem1,1);
      int elem_avant=face_vois_plus(*this,mon_domaine, face_avant, 0);
      if (elem_avant>=0)
        {
          if (Elem_Aretes(elem1,6)<0) Elem_Aretes(elem1,6)=Elem_Aretes(elem_avant,0);
          if (Elem_Aretes(elem1,11)<0) Elem_Aretes(elem1,11)=Elem_Aretes(elem_avant,2);

          int face_bas_elem_avant=elem_faces(elem_avant,2);
          int elem_avant_bas=face_vois_plus(*this,mon_domaine,face_bas_elem_avant,0);
          if (elem_avant_bas>=0 && Elem_Aretes(elem1,5)<0) Elem_Aretes(elem1,5)=Elem_Aretes(elem_avant_bas,2);
        }

    }

}

