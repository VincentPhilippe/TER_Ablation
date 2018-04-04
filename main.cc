#include <iostream>
#include "diffusion.h"
#include "maillage.h"
#include "read_data.h"

using namespace std;

int main(int argc, char** argv)
{

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  // --------------------------- Fichier de données --------------------------------/
  read_data* data_file = new read_data(data_file_name);
  data_file->read_datafile();

  double dx = data_file->Get_dx();
  double dz = data_file->Get_dz(); //une petite fonction pour get dy serait pas mal dans le cas 3D :3
  double Lx = data_file->Get_Lx();
  double Lz = data_file->Get_Lz();
  double Nx = data_file->Get_Nx();
  double dt = data_file->Get_dt();
  double tfinal = data_file->Get_tfinal();
  string file_name= data_file->Get_file_name();

  cout << "------------------- LECTURE DES DONNEES -------------------" << endl;
  cout << "Simulation -- Recul d'une surface composite"                 << endl;
  cout << "Fichier source : " << file_name                              << endl;
  cout << "Pas de temps dt = " << dt                                    << endl;
  cout << "Temps final Tmax = " << tfinal                               << endl;

  // -------------------------- Création du maillage -------------------------------/
  Cartesien2D maillage("Maillage 2D", dx, dz, Lx, Lz);
  cout << "----------------------- MAILLAGE -----------------------" << endl;
  cout << maillage.GetMeshName() << dx                               << endl;
  cout << "Pas d'espace dx =" << dx                                  << endl;
  cout << "Pas d'espace dz = " << dz                                 << endl;
  cout << "Nombre de noeuds = " << maillage.GetNx()*maillage.GetNz() << endl;

  // ----------------------------- Initialisation ----------------------------------/
  diffusion plateau(*data_file, maillage);

// Boucle en temps

  /*
  Résolution équation de diffusion
  */

  /*
  Calcul de variation de hauteur de la surface
  */

  /*
  Reconstruction de la surface
  */

  /*
  Sauvegarde de la concentration et de la surface
  */

// Fin de la boucle

  delete data_file;
  return 0;
}
