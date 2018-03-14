#include "Dense"

class recul
{
public:
  //recul(read_data& data, diffusion& diffusion, plic& plic);
  recul(double dt, double dx, double dz, Eigen::MatrixXd C_solide);
  ~recul();
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

private:
  Eigen::MatrixXd _ninterf;
  Eigen::MatrixXd _interface;
  Eigen::VectorXd _vitesse;
  Eigen::MatrixXd _C_solide;
  double _dt;
  double _dx;
  double _dz;
  int _nx;
  int _nz;
  //readdata& _data;
  //diffusion& _diff;
  //plic& _plic;
};
