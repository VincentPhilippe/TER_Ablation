#ifndef RECUL3D_H
#define RECUL3D_H

//#include "objet.h"
#include "Dense"
#include <cmath>
#include <iostream>
//#include "plic.h"
//#include "diffusion.h"
//#include "read_data.h"

class recul3D
{
public:
  //recul3D(read_data &read_data, std::vector<std::vector<std::vector<double> > > C_solide);
  recul3D(double dt, double dx, double dy, double dz, int nx, int ny, int nz, std::vector<std::vector<std::vector<double> > > C_solide);
  ~recul3D();
  //void update(plic *pplic, diffusion *pdiffusion){_plic=pplic; _diff=pdiffusion;};

  Eigen::VectorXd eqplan(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc);
  double volume_pyramide(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd);//Eigen::MatrixXd coord);
  double surface_triangle(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc);
  void recul3D_1(Eigen::MatrixXd repere, Eigen::MatrixXd coord, double vrdt);
  Eigen::MatrixXd repereglobal(int i, int j, int k);
  Eigen::MatrixXd rotationz(Eigen::MatrixXd repere_prec);
  Eigen::MatrixXd reductionrepere(Eigen::MatrixXd repere_prec);
  Eigen::MatrixXd rotationcoin(Eigen::MatrixXd repere_prec);

  inline std::vector<std::vector<std::vector<double> > > Get_C_solide(){return _C_solide;};
  inline std::vector<std::vector<std::vector<double> > > Get_ninterf(){return _ninterf;};

private:
  std::vector<std::vector<std::vector<double> > > _ninterf;
  Eigen::MatrixXd _interface;
  Eigen::VectorXd _vitesse;
  std::vector<std::vector<std::vector<double> > > _C_solide;
  double _dt;
  double _dtmax;
  double _dx;
  double _dy;
  double _dz;
  int _nx;
  int _ny;
  int _nz;
  //read_data& _read_data;
  //diffusion *_diff;
  //plic *_plic;
};

#endif
