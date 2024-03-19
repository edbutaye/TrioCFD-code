/////////////////////////////////////////////////////////////////////////////
//
// File      : Modele_Collision_FT.cpp
// Directory : $FPI184_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Modele_Collision_FT.h>
#include <TRUST_Deriv.h>
#include <Objet_U.h>
#include <Param.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Statistiques.h>
#include <Connex_components_FT.h>
#include <communications.h>
#include <Connex_components.h>
#include <fstream>
#include <iomanip>
#include <SFichier.h>
#include <EFichier.h>
#include <MD_Vector_tools.h>
#include <EcritureLectureSpecial.h>

//#include <Ref_Navier_Stokes_FT_Disc.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Schema_Temps_base.h>
#include <Probleme_FT_Disc_gen.h>

DoubleVect Modele_Collision_FT::Longueurs_modele_collision=0;
IntVect Modele_Collision_FT::Nb_Noeuds_modele_collision=0;
DoubleVect Modele_Collision_FT::Origine_modele_collision=0;

Implemente_instanciable_sans_constructeur(Modele_Collision_FT,"Modele_Collision_FT",Objet_U);

Modele_Collision_FT::Modele_Collision_FT() :
  tau_coll_(0),
  decalage_bords_(1),
  f_elem_diph_(1),
  modele_lubrification_(0),
  sigma_(0),
  delta_n_(0),
  raideur_cst_(0),
  amortissement_cst_(0),
  d_act_lub_(0.),
  d_sat_lub_(0.),
  nb_compo_tot_(0),
  //sauvegarde_reprise_fichier_unique_(0),
  raideur_(0),
  e_eff_(0),
  F_old_(0),
  F_now_(0),
  nb_dt_Verlet_(0),
  dt_compute_Verlet_(-1),
  activate_linked_cell_(0),
  valeurs_decalage(0),
  compteur_collisions_(0)
{
  positions_bords_.resize(0,6);
  valeurs_decalage.resize(6);
}


Entree& Modele_Collision_FT::readOn (Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  //Cerr << "Modele_Collision_FT::set_param sauvegarde_reprise_fichier_unique_ " << sauvegarde_reprise_fichier_unique_ << finl;

  return is;
}

Sortie& Modele_Collision_FT::printOn(Sortie& os) const
{
  Cerr << "Erreur : ::printOn n'est pas code." << finl;
  assert(0);
  return os;
}

void Modele_Collision_FT::set_param(Param& p)
{
  p.ajouter_non_std("modele_collision", (this),Param::REQUIRED);
  p.ajouter_non_std("detection_collision_Verlet", (this),Param::REQUIRED);
  p.ajouter("tau_coll", &tau_coll_, Param::REQUIRED); //en secondes
  p.ajouter("decalage_bords", &decalage_bords_, Param::REQUIRED); //en % du diametre de la particule
  p.ajouter("force_sur_elem_diphasiques", &f_elem_diph_, Param::REQUIRED); // 1 (True), 0 (False)
  p.ajouter("modele_lubrification", &modele_lubrification_, Param::REQUIRED); // 1 (True), 0 (False)
  p.ajouter("sigma", &sigma_);
  p.ajouter("delta_n", &delta_n_);
  p.ajouter("raideur_cst", &raideur_cst_);
  p.ajouter("amortissement_cst", &amortissement_cst_);
}

int Modele_Collision_FT::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{

  if (mot=="modele_collision")
    {
      Motcles mots;
      mots.add("ressort_amorti_esi");
      mots.add("ressort_amorti_vva");
      mots.add("ressort_amorti_ee");
      mots.add("mohaghegh");
      mots.add("hybrid_esi");
      mots.add("hybrid_ee");
      mots.add("hybrid_vva");
      mots.add("breugem");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading modele_collisions_part : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          modele_collision_ = Modele_Collision_FT::RESSORT_AMORTI_ESI;
          break;
        case 1:
          modele_collision_ = Modele_Collision_FT::RESSORT_AMORTI_VVA;
          break;
        case 2:
          modele_collision_ = Modele_Collision_FT::RESSORT_AMORTI_EE;
          break;
        case 3:
          modele_collision_ = Modele_Collision_FT::MOHAGHEGH;
          break;
        case 4:
          modele_collision_ = Modele_Collision_FT::HYBRID_ESI;
          break;
        case 5:
          modele_collision_ = Modele_Collision_FT::HYBRID_EE;
          break;
        case 6:
          modele_collision_ = Modele_Collision_FT::HYBRID_VVA;
          break;
        case 7:
          modele_collision_ = Modele_Collision_FT::BREUGEM;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot=="detection_collision_Verlet")
    {
      Motcles mots;
      mots.add("is_active");
      mots.add("s_Verlet");
      mots.add("nb_pas_dt_max");
      mots.add("activate_linked_cell");
      mots.add("partitionnement");
      Motcle accouverte = "{" , accfermee = "}" ;
      Motcle motbis;
      is >> motbis;
      if (motbis==accouverte)
        {
          is >> motbis;
          while (motbis != accfermee)
            {
              int rang2 = mots.search(motbis);
              switch(rang2)
                {
                case 0:
                  {
                    is >> is_detection_Verlet_;
                    break;
                  }
                case 1:
                  {
                    is >> s_Verlet_;
                    break;
                  }
                case 2:
                  {
                    is >> nb_pas_dt_max_Verlet_;
                    break;
                  }
                case 3:
                  {
                    is >> activate_linked_cell_;
                    break;
                  }
                case 4:
                  {
                    is >> Px_;
                    is >> Py_;
                    is >> Pz_;
                    break;
                  }
                default:
                  Cerr << "Transport_Interfaces_FT_Disc::lire\n"
                       << " les options de detection_collision_Verlet sont :\n"
                       << mots;
                  exit();

                }
              is >> motbis;
            }
        }
      return 1;
    }
  else
    {
      Cerr << mot << " is not a keyword understood by " << que_suis_je() << " in lire_motcle_non_standard"<< finl;
      exit();
    }
  return -1;
}

// copie de la methode ecrire_tableau(Sortie& os, const DoubleTab& tab) de Sauvegarde_Reprise_Maillage_FT.cpp
void ecrire_tableau_donnee_modele_collision(Sortie& os, const DoubleTab& tab)
{
  const int dim0 = tab.dimension(0);
  if (Process::je_suis_maitre())
    os << dim0 << space << tab.dimension(1) << finl;
  os.put(tab.addr(), tab.size_array());
  os.syncfile();
}
//
void lire_tableau_donnee_modele_collision(Entree& is, DoubleTab& tab, Entree * fichier)
{
  int dim0;
  int dim1;
  (*fichier)  >> dim0  >> dim1;
  DoubleTab tmp;
  tmp.resize(dim0,dim1);
  fichier->get(tmp.addr(), tmp.size_array());
  tab=tmp;
}

void Modele_Collision_FT::reset()
{
  /*const Transport_Interfaces_FT_Disc& eq_transport =refequation_transport_.valeur();
  const Maillage_FT_Disc& maillage_interface=eq_transport.maillage_interface();
  const int nb_facettes = maillage_interface.nb_facettes();
  ArrOfInt compo_connexes_facettes(nb_facettes); // Init a zero
  int n = search_connex_components_local_FT(maillage_interface, compo_connexes_facettes); // */
  int nb_compo_tot = nb_compo_tot_;//compute_global_connex_components_FT(maillage_interface, compo_connexes_facettes, n); //
  int nb_bords=6;
  F_old_.resize(nb_compo_tot,nb_compo_tot+nb_bords);
  raideur_.resize(nb_compo_tot,nb_compo_tot+nb_bords);
  e_eff_.resize(nb_compo_tot,nb_compo_tot+nb_bords);
}


void ouvrir_fichier_collision(SFichier& os, const int& flag, const Transport_Interfaces_FT_Disc& equation, Nom fichier_sauvegarde)
{
  // flag nul on n'ouvre pas le fichier
  if (flag==0)
    return ;
  Nom fichier=fichier_sauvegarde;

  const Schema_Temps_base& sch=equation.probleme().schema_temps();
  const int& precision=sch.precision_impr();
  // On cree le fichier a la premiere impression avec l'en tete ou si le fichier n'existe pas

  os.ouvrir(fichier,std::_S_out); // std::_S_out pour ne conserver que la derniere position
  os.precision(precision);
  os.setf(ios::scientific);
}

void ouvrir_fichier_collision(EFichier& os, const int& flag, const Transport_Interfaces_FT_Disc& equation, Nom fichier_reprise)
{
  // flag nul on n'ouvre pas le fichier
  if (flag==0)
    return ;
  Nom fichier=fichier_reprise;//Objet_U::nom_du_cas()+"_modele_collision.sauv";

  const Schema_Temps_base& sch=equation.probleme().schema_temps();
  const int& precision=sch.precision_impr();
  // On cree le fichier a la premiere impression avec l'en tete ou si le fichier n'existe pas

  os.ouvrir(fichier,ios::app);

  os.precision(precision);
  os.setf(ios::scientific);
}

int Modele_Collision_FT::reprendre(Entree& is)
{
  Nom motlu;
  const int format_xyz = EcritureLectureSpecial::is_lecture_special();

  reset();
  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          EFichier fich_modele_collision;
          Cerr << "fichier_reprise_FT_ " << fichier_reprise_FT_ << finl;
          ouvrir_fichier_collision(fich_modele_collision,1,refequation_transport_.valeur(),fichier_reprise_FT_);
          fich_modele_collision >> motlu;
          if (motlu != que_suis_je())
            {
              Cerr << "Erreur dans Modele_Collision_FT::reprendre\n";
              Cerr << " On attendait le motcle " << que_suis_je();
              Cerr << "\n On a trouve " << motlu << finl;
              Process::exit();
            }
          fich_modele_collision >> motlu;
          F_old_.lit(fich_modele_collision);
          fich_modele_collision >> motlu;
          raideur_.lit(fich_modele_collision);
          fich_modele_collision >> motlu;
          e_eff_.lit(fich_modele_collision);
          fich_modele_collision.close();
        }
      envoyer_broadcast(F_old_,0);
      envoyer_broadcast(raideur_,0);
      envoyer_broadcast(e_eff_,0);
      barrier();

      return 1;
    }
  else
    {
      is >> motlu;
      if (motlu != que_suis_je())
        {
          Cerr << "Erreur dans Modele_Collision_FT::reprendre\n";
          Cerr << " On attendait le motcle " << que_suis_je();
          Cerr << "\n On a trouve " << motlu << finl;
          Process::exit();
        }
      is >> motlu;
      F_old_.lit(is);
      is >> motlu;
      raideur_.lit(is);
      is >> motlu;
      e_eff_.lit(is);
    }


  return 1;
}

int Modele_Collision_FT::sauvegarder(Sortie& os) const
{
  int special, afaire;
  const int format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);

  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          Nom fichier_sauvegarde_FT="donnees_particules_FT.sauv";
          SFichier fich_modele_collision;
          ouvrir_fichier_collision(fich_modele_collision,1,refequation_transport_.valeur(),fichier_sauvegarde_FT);
          fich_modele_collision << Nom(que_suis_je()) << finl;
          fich_modele_collision << "F_old" << finl;
          F_old_.ecrit(fich_modele_collision);
          fich_modele_collision << "raideur" << finl;
          raideur_.ecrit(fich_modele_collision);
          fich_modele_collision << "e_eff" << finl;
          e_eff_.ecrit(fich_modele_collision);
          fich_modele_collision.close();
        }
      return 0;
    }
  else
    {
      int bytes = 0;
      os << que_suis_je() << finl;
      os << "F_old" << finl;
      F_old_.ecrit(os);
      bytes += 8 * F_old_.size_array();
      os << "raideur" << finl;
      raideur_.ecrit(os);
      bytes += 8 * raideur_.size_array();
      os << "e_eff" << finl;
      e_eff_.ecrit(os);
      bytes += 8 * e_eff_.size_array();
      return bytes;
    }


  return 0;
}


void  Modele_Collision_FT::calculer_force_contact(DoubleTab& force_contact, int& isFirstStepOfCollision, double& dist_int, double& next_dist_int, DoubleTab& norm, DoubleTab& dUn, double& masse_eff, int& compo, int& voisin, double& Stb, double& ed, double& vitesseRelNorm, double& dt, double& prod_scal)
{
  switch (modele_collision_)
    {
    case Modele_Collision_FT::RESSORT_AMORTI_ESI:
      {
        //Cerr << "Modele ressort_amorti. \n" << finl;
        double raideur = get_raideur_cst();
        double amortisseur = get_amortissement_cst();

        for (int d = 0; d < dimension; d++)
          force_contact(d)=-1 * raideur * next_dist_int * norm(d) -1*amortisseur*dUn(d);
      }
      break;
    case Modele_Collision_FT::RESSORT_AMORTI_EE:
      {
        //Cerr << "Modele ressort_amorti. \n" << finl;
        double raideur = get_raideur_cst();
        double amortisseur = get_amortissement_cst();

        for (int d = 0; d < dimension; d++)
          force_contact(d)=-1 * raideur * dist_int * norm(d) -1*amortisseur*dUn(d);
      }
      break;

    case Modele_Collision_FT::RESSORT_AMORTI_VVA:
      {
        //Cerr << "Modele ressort_amorti. \n" << finl;

        double raideur = get_raideur_cst();
        double amortisseur = get_amortissement_cst();

        for (int d = 0; d < dimension; d++)
          {

            double F_n = -1 * raideur * dist_int * norm(d) -1*amortisseur*dUn(d);
            double y_np1 = dist_int + dUn(d) * dt +0.5 * F_n * dt * dt / masse_eff ;
            double u_np12 = dUn(d) + 0.5 * F_n * dt / masse_eff ;
            double F_np1 = -1 * raideur * y_np1 * norm(d) -1*amortisseur*u_np12;
            force_contact(d)=0.5*(F_np1+F_n);
          }
      }
      break;

    case Modele_Collision_FT::MOHAGHEGH:
      {
        //Cerr << "Modele Mohagheg. \n" << finl;
        DoubleTab& raideur=get_raideur();
        DoubleTab& e_eff=get_e_eff();
        //int isFirstStepOfCollision = Tool::F_now(compo, voisin) > Tool::F_old(compo, voisin);
        if (isFirstStepOfCollision)
          {
            raideur(compo,voisin)= masse_eff * pow(vitesseRelNorm / sigma_, 2);
            e_eff(compo,voisin)=Stb > 18 ? ed * (1 - 8.65 * pow(Stb, -0.75)) : 0;
          }
        int isPhasePenetration = prod_scal <= 0;
        double le_e_eff = isPhasePenetration ? 1 : e_eff(compo, voisin);
        double la_raideur = raideur(compo, voisin);

        for (int d = 0; d < dimension; d++)
          force_contact(d)=-1 * le_e_eff * le_e_eff * la_raideur * next_dist_int * norm(d);
      }
      break;
    case Modele_Collision_FT::HYBRID_EE:
      {
        Cerr << "Modele hybrid. \n" << finl;
        DoubleTab& raideur=get_raideur();
        DoubleTab& e_eff=get_e_eff();
        if (isFirstStepOfCollision)
          {
            //double  tau_c = Nc * dt ;
            e_eff(compo,voisin)=ed * exp(-35 / (Stb + 1e-6));
            raideur(compo,voisin)=(masse_eff * (pow(M_PI,2) + pow(log(ed), 2))) / pow(tau_coll_, 2);
            //Tool::e_eff(compo, voisin) = ed * exp(-35 / (Stb + 1e-6));
            //Tool::raideur(compo, voisin) = (masse_eff * (pow(M_PI,2) + pow(log(ed), 2))) / pow(tau_coll, 2);
          }
        int isPhaseCompression = prod_scal <= 0;

        double le_e_eff = isPhaseCompression ? 1 : e_eff(compo, voisin);
        double la_raideur = raideur(compo, voisin);
        for (int d = 0; d < dimension; d++)
          {
            force_contact(d)=-1 * le_e_eff * le_e_eff * la_raideur * dist_int * norm(d);
          }
      }
      break;
    case Modele_Collision_FT::HYBRID_ESI:
      {
        DoubleTab& raideur=get_raideur();
        DoubleTab& e_eff=get_e_eff();
        if (isFirstStepOfCollision)
          {
            int Nc=8;
            double  tau_c = Nc * dt ;
            raideur(compo,voisin)=(masse_eff * (pow(M_PI,2)+ pow(log(ed), 2))) / pow(tau_c, 2);
            e_eff(compo,voisin)=ed * exp(-35 / (Stb + 1e-6));
          }

        int isPhaseCompression = prod_scal <= 0;
        double le_e_eff = isPhaseCompression ? 1 : e_eff(compo, voisin);
        double la_raideur = raideur(compo, voisin);
        //double amortisseur = -1*(masse_eff * log(ed)) / (tau_coll); //EB

        for (int d = 0; d < dimension; d++)
          {
            force_contact(d)=-1 * le_e_eff * le_e_eff * la_raideur * next_dist_int * norm(d);
          }

      }
      break;
    case Modele_Collision_FT::HYBRID_VVA:
      {
        DoubleTab& raideur=get_raideur();
        DoubleTab& e_eff=get_e_eff();
        if (isFirstStepOfCollision)
          {
            raideur(compo,voisin)=(masse_eff * (pow(M_PI,2) + pow(log(ed), 2))) / pow(tau_coll_, 2);
            e_eff(compo,voisin)=ed * exp(-35 / (Stb + 1e-6));

          }

        int isPhaseCompression = prod_scal <= 0;
        double le_e_eff = isPhaseCompression ? 1 : e_eff(compo, voisin);
        double la_raideur = raideur(compo, voisin);

        for (int d = 0; d < dimension; d++)
          {
            double F_n = -1 * la_raideur * dist_int * norm(d) ;
            double y_np1 = dist_int + dUn(d) * dt +0.5 * F_n * dt * dt / masse_eff ;
            double F_np1 = -1 * le_e_eff * le_e_eff*la_raideur * y_np1 * norm(d) ;
            force_contact(d)=F_np1;
          }
      }

      break;


    case Modele_Collision_FT::BREUGEM:
      {
        double raideur = (masse_eff * (pow(M_PI,2) + pow(log(ed), 2))) / pow(tau_coll_, 2); // tau_coll tel que definit dans le jdd
        double amortisseur = -2*(masse_eff * log(ed)) / (tau_coll_);

        for (int d = 0; d < dimension; d++)
          force_contact(d)=-1 * raideur * next_dist_int * norm(d) -1*amortisseur*dUn(d);
      }
      break;


    default:
      Cerr << "The method specified for modele_collision in not recognized. \n" << finl;
      Process::exit();
    }
}

void Modele_Collision_FT::calculer_positions_bords(const DoubleVect& rayon_particule)
{

  int nb_compo_tot=rayon_particule.size();
  positions_bords_.resize(nb_compo_tot,6);
  valeurs_decalage.resize(6);
  for (int compo=0; compo<nb_compo_tot; compo++)

    {
      switch(decalage_bords_)
        {
        case 0:
          valeurs_decalage=0;
          break;
        case 1:
          if (delta_n_>0)
            {
              valeurs_decalage=delta_n_*(2*rayon_particule(compo))/100;
            }
          else
            {
              valeurs_decalage=(Longueurs_modele_collision(0)/(Nb_Noeuds_modele_collision(0)-1))/4;
            }
          break;

        default:
          Cerr << "error value of decalage_bords"  <<finl;
          Process::exit();
          break;
        }

      // les bords sont traites comme des particules fictives placer en dehors du domaine, en mirroir.
      positions_bords_(compo,0) = Origine_modele_collision(0) - rayon_particule(compo) + valeurs_decalage(0);
      positions_bords_(compo,1) = Origine_modele_collision(1) - rayon_particule(compo) + valeurs_decalage(1);
      positions_bords_(compo,2) = Origine_modele_collision(2) - rayon_particule(compo) + valeurs_decalage(2);
      positions_bords_(compo,3) = Origine_modele_collision(0) + Longueurs_modele_collision(0) + rayon_particule(compo) - valeurs_decalage(3);
      positions_bords_(compo,4) = Origine_modele_collision(1) + Longueurs_modele_collision(1) + rayon_particule(compo) - valeurs_decalage(4);
      positions_bords_(compo,5) = Origine_modele_collision(2) + Longueurs_modele_collision(2) + rayon_particule(compo) - valeurs_decalage(5);
    }
}

void Modele_Collision_FT::set_resize_parametres_geometriques()
{
  Longueurs_modele_collision.resize(dimension);
  Nb_Noeuds_modele_collision.resize(dimension);
  Origine_modele_collision.resize(dimension);
}
void Modele_Collision_FT::set_longueur(DoubleVect& Longueurs)
{
  Longueurs_modele_collision=Longueurs;
}
void Modele_Collision_FT::set_longueur(double Lx, double Ly, double Lz)
{
  Longueurs_modele_collision(0)=Lx;
  Longueurs_modele_collision(1)=Ly;
  if (dimension==3) Longueurs_modele_collision(2)=Lz;
}

void Modele_Collision_FT::set_origin(DoubleVect& Origin)
{
  Origine_modele_collision=Origin;
}
void Modele_Collision_FT::set_origin(double Ox, double Oy, double Oz)
{
  Origine_modele_collision(0)=Ox;
  Origine_modele_collision(1)=Oy;
  if (dimension==3) Origine_modele_collision(2)=Oz;
}

DoubleVect& Modele_Collision_FT::get_origin()
{
  return Origine_modele_collision;
}
DoubleVect& Modele_Collision_FT::get_longueurs()
{
  return Longueurs_modele_collision;
}
IntVect& Modele_Collision_FT::get_nb_noeuds()
{
  return Nb_Noeuds_modele_collision;
}

void Modele_Collision_FT::set_s_Verlet(double s_Verlet)
{
  s_Verlet_=s_Verlet;
}

const int& Modele_Collision_FT::is_detection_Verlet() const
{
  return is_detection_Verlet_;
}
const int& Modele_Collision_FT::is_LC_activated() const
{
  return activate_linked_cell_;
}
double& Modele_Collision_FT::get_s_Verlet()
{
  return s_Verlet_;
}
const int& Modele_Collision_FT::get_Px() const
{
  return Px_;
}
const int& Modele_Collision_FT::get_Py() const
{
  return Py_;
}
const int& Modele_Collision_FT::get_Pz() const
{
  return Pz_;
}

int& Modele_Collision_FT::get_nb_dt_Verlet()
{
  return nb_dt_Verlet_;
}
int& Modele_Collision_FT::get_dt_compute_Verlet()
{
  return dt_compute_Verlet_;
}
int& Modele_Collision_FT::get_nb_pas_dt_max_Verlet()
{
  return nb_pas_dt_max_Verlet_;
}
const double& Modele_Collision_FT::get_d_act_lub() const
{
  return d_act_lub_;
}
const double& Modele_Collision_FT::get_d_sat_lub() const
{
  return d_sat_lub_;
}
ArrOfIntFT& Modele_Collision_FT::get_liste_zone_sup()
{
  return liste_zone_sup_;
}

ArrOfIntFT& Modele_Collision_FT::get_liste_zone_inf()
{
  return liste_zone_inf_;
}

const ArrOfIntFT& Modele_Collision_FT::get_liste_zone_sup() const
{
  return liste_zone_sup_;
}

const ArrOfIntFT& Modele_Collision_FT::get_liste_zone_inf() const
{
  return liste_zone_inf_;
}

const double& Modele_Collision_FT::sigma() const
{
  return sigma_;
}
const double& Modele_Collision_FT::delta_n() const
{
  return delta_n_;
}
const double& Modele_Collision_FT::get_raideur_cst() const
{
  return raideur_cst_;
}
const double& Modele_Collision_FT::get_amortissement_cst() const
{
  return amortissement_cst_;
}
DoubleTab& Modele_Collision_FT::get_raideur()
{
  return raideur_;
}
DoubleTab& Modele_Collision_FT::get_e_eff()
{
  return e_eff_;
}
DoubleTab& Modele_Collision_FT::get_F_old()
{
  return F_old_;
}
DoubleTab& Modele_Collision_FT::get_F_now()
{
  return F_now_;
}
DoubleTab& Modele_Collision_FT::get_forces_solide()
{
  return forces_solide_;
}
DoubleVect& Modele_Collision_FT::get_collisions_detected()
{
  return collision_detected_;
}
const DoubleTab& Modele_Collision_FT::position_bords() const
{
  return positions_bords_;
}

int& Modele_Collision_FT::compteur_collisions()
{
  return compteur_collisions_;
}

const int& Modele_Collision_FT::modele_lubrification() const
{
  return modele_lubrification_;
}

const int& Modele_Collision_FT::force_elem_diphasique() const
{
  return f_elem_diph_;
}

const double& Modele_Collision_FT::tau_coll() const
{
  return tau_coll_;
}

int Modele_Collision_FT::checkForDuplicates(ArrOfInt& vector)
{
  int flag =0;
  ArrOfInt copy_vector(vector);
  const int size = copy_vector.size_array();
  copy_vector.ordonne_array();
  for (int i = 0; i < size-1; i++)
    {
      if (copy_vector(i)==copy_vector(i+1))
        {
          flag = 1;
          Cerr << copy_vector(i) << " is duplicate !!" << finl ;
        }
    }

  return flag;
}
void Modele_Collision_FT::set_param_geom(Domaine_VDF& domaine_vdf)
{
  const Domaine& domaine = domaine_vdf.domaine();
  DoubleTab BB=domaine.getBoundingBox();
  const Bords& bords=domaine.faces_bord();
  DoubleVect NiNj(dimension); // EB NiNj=(NyNz NxNz NxNy )
  Modele_Collision_FT::set_resize_parametres_geometriques();
  // 1. Nombre de noeuds par direction
  NiNj=0;
  for (int i=0; i<bords.nb_bords(); i++)
    {
      int nb_boundary_faces = mp_sum(ref_cast(Frontiere,bords(i)).nb_faces());
      int nb_boundary_faces_local=ref_cast(Frontiere,bords(i)).nb_faces();
      if (nb_boundary_faces_local>0)
        {
          int face1=ref_cast(Frontiere,bords(i)).num_premiere_face();
          int orientation_face1=domaine_vdf.orientation(face1);
          NiNj(orientation_face1)=nb_boundary_faces;
        }
    }

  // EB : les valeurs seront bidons si les conditions citees precedemment ne sont pas remplies
  long long NxNy= static_cast<int>(mp_max(NiNj(2))) ;
  long long NxNz= static_cast<int>(mp_max(NiNj(1)));
  long long NyNz= static_cast<int>(mp_max(NiNj(0)));

  int Nx,Ny,Nz;

  Ny= NxNz>0 ? static_cast<int>(sqrt(NxNy*NyNz/NxNz)) : 0; // nb elements dans la direction y
  Nz= NxNy>0 ? static_cast<int>((NxNz*Ny/NxNy)) : 0; // idem z, attention a l'ordre des operations car Nz est un entier
  Nx= Ny>0 ? static_cast<int>(NxNy/Ny) : 0; // idem x
  Nx++; // nb noeuds dans la direction x
  Ny++; // idem y
  Nz++; // idem z
  Nb_Noeuds_modele_collision(0)=Nx;
  Nb_Noeuds_modele_collision(1)=Ny;
  Nb_Noeuds_modele_collision(2)=Nz;

  // 2. Origine et Longueurs
  double Ox=0,Oy=0,Oz=0;
  double Lx=0,Ly=0,Lz=0;

  for (int j=0; j<dimension; j++)
    {
      double min_ = mp_min(BB(j,0));
      double max_ = mp_max(BB(j,1));
      if (j==0)
        {
          Ox=min_;
          Lx=max_-min_;
        }
      if (j==1)
        {
          Oy=min_;
          Ly=max_-min_;
        }
      if (j==2)
        {
          Oz=min_;
          Lz=max_-min_;
        }
    }

  set_origin(Ox,Oy,Oz);
  set_longueur(Lx,Ly,Lz);

  Cerr << "Ox Oy Oz " << Ox << " " << Oy << " " << Oz << finl;
  Cerr << "Lx Ly Lz " << Lx << " " << Ly << " " << Lz << finl;
  Cerr << "Nx Ny Nz " << Nx << " " << Ny << " " << Nz << finl;
}
void Modele_Collision_FT::associer_equation_transport(const Equation_base& equation)
{
  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc,equation);
  refequation_transport_ = eq;
}

void Modele_Collision_FT::set_nb_compo_tot(int nb_compo_tot)
{
  nb_compo_tot_=nb_compo_tot;
}

void Modele_Collision_FT::set_nom_fichier_reprise_FT(Nom fichier_reprise_FT)
{
  fichier_reprise_FT_=fichier_reprise_FT;
}

void Modele_Collision_FT::set_d_act_lub(double d_act_lub)
{
  d_act_lub_=d_act_lub;
}
void Modele_Collision_FT::set_d_sat_lub(double d_sat_lub)
{
  d_sat_lub_=d_sat_lub;
}
