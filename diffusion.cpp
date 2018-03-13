#include <iostream>

#include "dense"
#include "maillage.h"
#include "read_data.h"
#include "rebuild_surface.h"

using namespace std;
using namespace Eigen;

diffusion::diffusion(Maillage2DCarre& maillage, ReadData& data)
{
  _maillage = maillage;
  _data = data;

 _ninterf = data.Get
}
