//fichier de test (si vous aviez pas compris)
#include "maillage.h"
#include <iostream>
using namespace std;
using namespace Eigen;

int main()
{
  Cartesien3D Cartesien("test", 1., 1., 1., 2., 2., 2.);
  cout << Cartesien.GetIndices() << endl;
  cout << "----------------------     X     ----------------------" << endl;
  cout << Cartesien.GetCoordX() << endl;
  cout << "----------------------     Y     ----------------------" << endl;
  cout << Cartesien.GetCoordY() << endl;
  cout << "----------------------     Z     ----------------------" << endl;
  cout << Cartesien.GetCoordZ() << endl;
  cout << Cartesien.GetNx() << endl;
  cout << Cartesien.GetNz() << endl;
  cout << Cartesien.GetNz() << endl;
}
