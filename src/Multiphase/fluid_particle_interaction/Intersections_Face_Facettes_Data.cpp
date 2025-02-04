/*
 * Intersections_Face_Facettes_Data.cpp
 *
 *  Created on: Feb 9, 2023
 *      Author: edouard
 */


#include <Intersections_Face_Facettes_Data.h>

// ===========================================================================
// Description:
// Ajoute une entree a la liste doublement chainee d'intersections entre
// la facette d'interface num_facette et la face eulerienne num_face.
// Le numero de la face doit verifier 0 <= num_face < zone.nb_faces()
// Si le numero de facette est superieur a la taille de l'index des facettes,
// on agrandit l'index.
void Intersections_Face_Facettes::ajoute_intersection(int num_facette,
                                                      int num_face,
                                                      double fraction_surface_intersection,
                                                      double contrib_volume_phase1,
                                                      double barycentre_u,
                                                      double barycentre_v,
                                                      double barycentre_w)
{
  // Verification de taille de l'index des facettes
  int nb_facettes = index_facette_face_.size_array();
  if (num_facette >= nb_facettes)
    {
      index_facette_face_.resize_array(num_facette+1);
      for (int i = nb_facettes; i <= num_facette; i++)
        index_facette_face_[i] = -1;
    }
  // Verification de taille du tableau de donnees des intersections
  if (data_allocated_size < data_real_size+1)
    {
      data_allocated_size = (data_real_size+1)*2;
      Intersections_Face_Facettes_Data *new_data =
        new Intersections_Face_Facettes_Data[data_allocated_size];
      if (data_real_size > 0)
        for (int i=0; i<data_real_size; i++)
          new_data[i] = data[i];
      if (data)
        delete[] data;
      data = new_data;
    }
  // Ajout de la nouvelle entree en debut de liste (pour eviter d'avoir a
  // chercher la fin...)
  {
    // Indice de l'intersection entre la face et la premiere facette
    int& lindex_face = index_face_facette_[num_face];
    // Indice de l'intersection entre la facette et la premiere face
    int& lindex_facette = index_facette_face_[num_facette];
    Intersections_Face_Facettes_Data& new_entry = data[data_real_size];
    // Indice de l'intersection de la meme facette avec la face suivante
    new_entry.index_face_suivante_ = lindex_facette;
    // Indice de l'intersection de la meme face avec la facette suivante
    new_entry.index_facette_suivante_ = lindex_face;
    new_entry.numero_facette_ = num_facette;
    new_entry.numero_face_ = num_face;
#ifdef AVEC_BUG_SURFACES
    new_entry.surface_intersection_ = fraction_surface_intersection;
#else
    new_entry.fraction_surface_intersection_ = fraction_surface_intersection;
#endif
    new_entry.contrib_volume_phase1_ = contrib_volume_phase1;
    new_entry.barycentre_[0] = barycentre_u;
    new_entry.barycentre_[1] = barycentre_v;
    new_entry.barycentre_[2] = barycentre_w;
    lindex_face = data_real_size;
    lindex_facette = data_real_size;
  }
  data_real_size++;
}


Intersections_Face_Facettes::Intersections_Face_Facettes()
  : data_allocated_size(0), data_real_size(0), data(0)
{
}

Intersections_Face_Facettes::~Intersections_Face_Facettes()
{
  if (data)
    delete[] data;
  data = 0;
}

void Intersections_Face_Facettes::reset(int nb_faces_euleriennes,
                                        int nb_facettes)
{
  data_real_size = 0;
  index_face_facette_.resize_array(nb_faces_euleriennes);
  index_face_facette_ = -1;
  index_facette_face_.resize_array(nb_facettes);
  index_facette_face_ = -1;
}

void Intersections_Face_Facettes::get_liste_faces_traversees(int num_facette,
                                                             ArrOfInt& liste_faces) const
{
  liste_faces.resize_array(0);
  if (num_facette >= index_facette_face_.size_array())
    // Aucune intersection n'a ete ajoutee pour cette facette,
    // elle ne figure meme pas dans l'index
    return;

  int index = index_facette_face_[num_facette];
  for (; index >= 0; index = data[index].index_face_suivante_)
    {
      const int face = data[index].numero_face_;
      liste_faces.append_array(face);
    }
}

void Intersections_Face_Facettes::get_liste_facettes_traversantes(int num_face,
                                                                  ArrOfInt& liste_facettes) const
{
  liste_facettes.resize_array(0);
  int index = index_face_facette_[num_face];
  for (; index >= 0; index = data[index].index_facette_suivante_)
    {
      const int facette = data[index].numero_facette_;
      liste_facettes.append_array(facette);
    }
}

// Description:
// Renvoie un tableau de taille zone.nb_faces():
//  pour une face 0 <= face < zone.nb_faces(),
//  index_face()[face] est l'indice de la premiere intersection entre la face
//  et les facettes du maillage lagrangien (voir description de la classe)
const ArrOfInt& Intersections_Face_Facettes::index_face() const
{
  return index_face_facette_;
}

// Description:
// Renvoie un tableau de taille "nombre de facettes de l'interface"
//  pour une face 0 <= facette < nb_facettes,
//  index_facette()[facette] est l'indice de la premiere intersection entre
//  la facette lagrangienne et les faces du maillage eulerien
//  (voir description de la classe)
const ArrOfInt& Intersections_Face_Facettes::index_facette() const
{
  return index_facette_face_;
}

// Description:
//operateur de copie
const Intersections_Face_Facettes& Intersections_Face_Facettes::operator=(const Intersections_Face_Facettes& iff)
{
  if (&iff != this)
    {
      index_face_facette_ = iff.index_face();
      assert(index_face_facette_.size_array()>0);
      index_facette_face_ = iff.index_facette();
      data_allocated_size = iff.data_allocated_size;
      data_real_size = iff.data_real_size;

      if (data)
        delete data;
      data = new Intersections_Face_Facettes_Data[data_allocated_size];
      int i;
      for (i=0 ; i<data_real_size; i++)
        {
          data[i] = iff.data[i];
        }
    }
  return *this;
}
