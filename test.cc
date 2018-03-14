//fichier de test (si vous aviez pas compris)
#include "maillage.h"
#include <iostream>
using namespace std;
using namespace Eigen;

int main()
{
  Cartesien Cartesien("test", 1., 1., 10., 10.);
  cout << Cartesien.GetIndices() << endl;
  cout << Cartesien.GetCoordX() << endl;
  cout << Cartesien.GetCoordZ() << endl;
  cout << Cartesien.GetNx() << endl;
  cout << Cartesien.GetNz() << endl;
}
