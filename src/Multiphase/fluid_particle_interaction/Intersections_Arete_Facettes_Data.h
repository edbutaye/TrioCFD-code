/*
 * Intersections_Arete_Facettes_Data.h
 *
 *  Created on: Mar 23, 2023
 *      Author: edouard
 */

#ifndef INTERSECTIONS_ARETE_FACETTES_DATA_included
#define INTERSECTIONS_ARETE_FACETTES_DATA_included
#include <TRUSTTabFT.h>
// EB : Cette classe reprend a l'identique (ou presque) ce qui est fait dans Intersections_Elem_Facettes_Data.
// Cette classe permet de construire une indicatrice aux faces de la meme maniere que l'indicatrice aux elements
// Cela est requis pour les ecoulements fluide-particules avec de forts sauts de propriete a l'interface (viscosite, masse volumique)
// Dans ces cas, avoir une indicatrice aux faces calculee par une moyenne des indicatrices des elements de part et d'autre ne suffit plus
// idem aux aretes
// /!\ Une face est vue comme volume, le travail est donc similaire a ce qui est fait pour les elements, il faut "simplment" decaler les volumes de
//controle et gerer les bords

class Intersections_Arete_Facettes_Data
{
public: // On copie ce qui est fait dans Intersections_Elem_Facettes_Data (parce que je ne comprends pas tout)
#ifdef AVEC_BUG_SURFACES
  // Surface de l'intersection facette/element
  // !!!!! gros bug : en fait, quand on remplit la structure, on y met la fraction
  //  de surface.
  double surface_intersection_;
#else
  // Fraction surface_intersection_facette_arete / surface_totale_facette:
  double fraction_surface_intersection_;
#endif
  double contrib_volume_phase1_;
  // Coordonnees barycentriques du centre de gravite de l'intersection
  // par rapport aux trois sommets de la facette. On a toujours
  // barycentre[0] + barycentre[1] + barycentre[2] = 1.
  // En 2D, on a barycentre_[2] = 0.;
  double barycentre_[3];

  int index_arete_suivante_;  // -1 si derniere arete de la liste
  int index_facette_suivante_; // idem.
  int numero_facette_;
  int numero_arete_;
};

// ====================================================================
// .DESCRIPTION        : class Intersections_Arete_Facettes
//
// Cette classe contient les donnees des intersections entre les facettes
// de l'interface et les aretes euleriennes sous la forme d'une liste
// doublement chainee.
// Pour parcourir les facettes qui coupent une face "arete", on fait:
//
//   int index=index_arete()[arete];
//   while (index >= 0) {
//     const Intersections_Arete_Facettes_Data & data = data_intersection(index);
//     ... // faire quelque chose avec data
//     index = data.index_facette_suivante_;
//   }
//
// Pour parcourir les faces qui sont coupes par une facette "facette":
//
//   int index=index_facette()[facette];
//   while (index >= 0) {
//     const Intersections_Arete_Facettes_Data & data = data_intersection(index);
//     ... // faire quelque chose avec data
//     index = data.index_arete_suivante_;
//   }

class Intersections_Arete_Facettes
{
public:
  Intersections_Arete_Facettes();
  ~Intersections_Arete_Facettes();
  void get_liste_aretes_traversees(int num_facette,
                                   ArrOfInt& liste_aretes) const;
  void get_liste_facettes_traversantes(int num_arete,
                                       ArrOfInt& liste_facettes) const;
  void ajoute_intersection(int num_facette,
                           int num_arete,
                           double surface_intersection,
                           double contrib_volume_phase1,
                           double barycentre_u,
                           double barycentre_v,
                           double barycentre_w);

  // A faire avant la premiere utilisation et pour remettre a zero
  // nb_aretes_euleriennes doit etre correct.
  // nb_facettes peut etre une estimation seulement (en pratique lors du
  // parcours on ne peut pas prevoir le nombre final de facettes)
  void reset(int nb_aretes_euleriennes=0, int nb_facettes=0);

  const ArrOfInt& index_arete() const;
  const ArrOfInt& index_facette() const;
  inline const Intersections_Arete_Facettes_Data& data_intersection(int index) const;
  inline Intersections_Arete_Facettes_Data& get_set_data_intersection(int index);

  //operateur de copie
  const Intersections_Arete_Facettes& operator=(const Intersections_Arete_Facettes& iaf);

private:
  // Tableau de taille "nombre d'aretes euleriennes"
  // Pour chaque face, index de l'entree correspondant a la premiere facette
  // dans data.
  //  -1 si l'arete ne contient pas de facette
  ArrOfIntFT index_arete_facette_;
  // pour chaque facette, index de l'entree correspondant a la premiere arete
  //  -1 si la facette ne traverse aucune face (facette qui devrait etre effacee)
  ArrOfIntFT index_facette_arete_;

  //
  // Les donnees : une arete par couple facette/face traverse.
  // (politique d'allocation : on agrandit d'un facteur 2, on ne diminue jamais).
  int data_allocated_size; // Le nombre d'aretes allouees
  int data_real_size;      // Le nombre d'aretes reellement utilisees
  Intersections_Arete_Facettes_Data * data;
};

// Description:
//  Renvoie les donnees de l'intersection stockee a l'indice "index"
//  dans le tableau "data" ( 0 <= index < data_real_size_ )
inline const Intersections_Arete_Facettes_Data& Intersections_Arete_Facettes::data_intersection(int index) const
{
  //Cerr << "data_real_size " << data_real_size << finl;
  assert(index >= 0 && index < data_real_size);
  return data[index];
}
// Description:
//  Renvoie les donnees de l'intersection stockee a l'indice "index"
//  dans le tableau "data" ( 0 <= index < data_real_size_ )
// ATTENTION A SON UTILISATION !!!
inline Intersections_Arete_Facettes_Data& Intersections_Arete_Facettes::get_set_data_intersection(int index)
{
  assert(index >= 0 && index < data_real_size);
  return data[index];
}





#endif
