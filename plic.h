#include <vector>
#include <string>
#include "Sparse"
#include "Dense"

class plic {
  public:
  plic(read_data* read_data);
  double grad_x(const int i, const int j);
  double grad_y(const int i, const int j);
  void interf(const int lon,const int lar);


  private:
  Eigen::Matrix<double, Eigen::Dynamic, 2> _phi; 
  std::vector<std::vector<int> > _phi;
  double p,nx,nxx,ny;
  int k,_pointsupl
  //std::vector<std::vector<int> > _interface;   //_inteface[i][j,ax,ay,bx,by]
  Eigen::Matrix<int, Eigen::Dynamic, 2> _ninterf, typinterf;   //_ninterf[i][j]=k   //numérotes des cases où il y a présence d'interface
  Eigen::Matrix<double, Eigen::Dynamic, 2> _interface; //_interface[k][ax,ay,bx,by], où k est le numéro de la case

  // Sauvegarde la solution
	//void SaveSol( int n);


/////////////////////RESTE A FAIRE///////////////////
//resize _interface
//SaveSol
}
