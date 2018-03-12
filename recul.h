#include "Dense"

void recul_surface(Eigen::MatrixXd Surface, Eigen::MatrixXd C_solide, double dt, double dx, double dy);
Eigen::MatrixXd recul1(Eigen::MatrixXd C_solide, int i, int j, double alpha, double l, double vrdt, double dx, double dy);
Eigen::MatrixXd recul2(Eigen::MatrixXd C_solide, int i, int j, double alpha, double l, double vrdt, double dx, double dy, double xc);
Eigen::MatrixXd recul3(Eigen::MatrixXd C_solide, int i, int j, double alpha, double l, double vrdt, double dx, double dy, double yd);
Eigen::MatrixXd recul4(Eigen::MatrixXd C_solide, int i, int j, double alpha, double l, double vrdt, double dx, double dy, double xc, double yd, int nx);
Eigen::MatrixXd recul5(Eigen::MatrixXd C_solide, int i, int j, double alpha, double l, double vrdt, double dx, double dy, double xc, double yc, double yd, int nx);
