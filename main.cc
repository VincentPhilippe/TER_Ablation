#include <iostream>
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

  // ----------------------- Fichier de données --------------------------------
  read_data* data_file = new read_data(data_file_name);
  data_file->read_datafile();
  //------- /

  /* code */

  /*
  Lire paramètres .dat
  */

  /*
  Afficher menu
  */

  /*
  Création maillage
  */

  /*
  Initialisation concentration et surface
  +Sauvegarde
  */

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
