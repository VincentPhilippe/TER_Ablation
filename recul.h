#ifndef RECUL_H
#define RECUL_H

#include "Dense"
#include <cmath>
#include "plic.h"
#include "diffusion.h"
#include "read_data.h"

class recul
{
public:
  recul(read_data &read_data, Eigen::MatrixXd C_solide);
  //recul(double dt, double dx, double dz, Eigen::MatrixXd C_solide);
  ~recul();
void update(plic *plic, diffusion *diffusion){_plic=plic, _diff=diffusion;};/////////////////////////////////////////////////////////////////

  void recul_surface();
void recul1(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul2(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul3(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul4(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul5(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul6(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul7(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul8(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul9(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul10(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul11(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul12(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul13(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul14(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul15(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul16(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul17(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void recul18(int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord);
void cpositive();

inline Eigen::MatrixXd Get_C_solide(){return _C_solide;};
inline Eigen::MatrixXd Get_ninterf(){return _ninterf;};

private:
  Eigen::MatrixXd _ninterf;
  Eigen::MatrixXd _interface;
  Eigen::VectorXd _vitesse;
  Eigen::MatrixXd _C_solide;
  double _dt;
  double _dtmax;
  double _dx;
  double _dz;
  int _nx;
  int _nz;
  read_data& _read_data;
  diffusion *_diff;
  plic *_plic;
};

#endif