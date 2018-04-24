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
  double dt = data_file->Get_dt();
  double tfinal = data_file->Get_tfinal();
  string file_name= data_file->Get_file_name();

  int nb_iterations = tfinal/dt;

  cout << "------------------- LECTURE DES DONNEES -------------------" << endl;
  cout << "Simulation -- Recul d'une surface composite"                 << endl;
  cout << "Fichier source : " << file_name                              << endl;
  cout << "Pas de temps dt = " << dt                                    << endl;
  cout << "Temps final Tmax = " << tfinal                               << endl;
  cout << "Nombre d'itérations = " << nb_iterations                     << endl;

  // -------------------------- Création du maillage -------------------------------/
  Cartesien2D maillage("Maillage 2D", dx, dz, Lx, Lz);
  cout << "----------------------- MAILLAGE -----------------------" << endl;
  cout << maillage.GetMeshName() << dx                               << endl;
  cout << "Pas d'espace dx =" << dx                                  << endl;
  cout << "Pas d'espace dz = " << dz                                 << endl;
  cout << "Nombre de noeuds = " << maillage.GetNx()*maillage.GetNz() << endl;

  // ----------------------------- Initialisation ----------------------------------/
  diffusion* plateau = new diffusion(*data_file, maillage);
  cout << "diffusion" << endl;

  recul* precul = new recul(*data_file);
  cout << "recul" << endl;
  precul->cpositive();
  cout << "cpositive" << endl;

  plic* pplic = new plic();
  cout << "plic" << endl;
  //construction de la première interface
  pplic->update(precul);
  cout << "update" << endl;
  pplic->interf();
  cout << "interf" << endl;

  // Boucle en temps
  for (int i = 1; i < nb_iterations; i++)
  {
    cout << "i " << i << endl;
    //sauvergarde de la solution
    cout << "SaveSol" << endl;
    pplic->SaveSol(i);
    //calcul des nouvelles valeurs de concentrations et de la vitesse de recul
    cout << "resolution" << endl;
    plateau->resolution();
    //plateau->vitesse();
    //mise à jour des données à l'interface, recul de la surface
    cout << "update recul" << endl;
    precul->update(pplic, plateau);
    cout << "recul_surface" << endl;
    precul->recul_surface();
    //reconstruction de l'interface
    cout << "update plic" << endl;
    pplic->update(precul);
    cout << "interf" << endl;
    pplic->interf();
  }
  // Fin de la boucle
  
  // Sauvegarde de la dernière solution
  pplic->SaveSol(nb_iterations);

  delete data_file;
  delete precul;
  delete pplic;
  delete plateau;
  return 0;
}
