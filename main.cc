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
  cout << "diffusion" << endl;
  diffusion* plateau = new diffusion(*data_file, maillage);

  cout << "recul" << endl;
  recul* precul = new recul(*data_file);
  cout << "cpositive" << endl;
  precul->cpositive();

  cout << "plic" << endl;
  plic* pplic = new plic(*data_file);
  //construction de la première interface
  cout << "update" << endl;
  pplic->update(precul);
  cout << "interf" << endl;
  pplic->interf();

  // Boucle en temps
  for (int i = 1; i < nb_iterations; i++)
  {
    cout << endl << "--------------I= " << i << "----------------" <<endl << endl;
    //sauvergarde de la solution
    cout << "SaveSol" << endl;
    pplic->SaveSol(i);
    //calcul des nouvelles valeurs de concentrations et de la vitesse de recul
    //cout << "update resolution" << endl;
    plateau->update(pplic);
    cout << "resolution" << endl;
    plateau->resolution();
    //cout << "vitesse" << endl << plateau->GetVitesse() << endl;
    //mise à jour des données à l'interface, recul de la surface
    //cout << "update recul" << endl;
    precul->update(pplic, plateau);
    cout << "recul_surface" << endl;
    precul->recul_surface();
    //cout << "C_solide" << endl << precul->Get_C_solide() << endl;
    //cout << "ninterf" << endl << precul->Get_ninterf() << endl;
    //reconstruction de l'interface
    //cout << "update plic" << endl;
    pplic->update(precul);
    cout << "interf" << endl;
    pplic->interf();
    cout<<"N_INTERFACE APRES:"<<endl<<precul->Get_ninterf()<<endl;
    cout<<"C_SOLIDE APRES;"<<endl<<precul->Get_C_solide()<<endl;
    //cout << "interface" << endl << pplic->Get_interface() << endl;

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
