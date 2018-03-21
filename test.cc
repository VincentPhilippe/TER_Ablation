//fichier de test (si vous aviez pas compris)
#include "maillage.h"
#include <iostream>
using namespace std;
using namespace Eigen;

int main()
{
  Cartesien3D Cartesien("test", 1., 1., 1., 5., 5., 5.);
  cout << Cartesien.GetIndices() << endl;
  cout << Cartesien.GetCoordX() << endl;
  cout << Cartesien.GetCoordY() << endl;
  cout << Cartesien.GetCoordZ() << endl;
  cout << Cartesien.GetNx() << endl;
  cout << Cartesien.GetNy() << endl;
  cout << Cartesien.GetNz() << endl;
}
