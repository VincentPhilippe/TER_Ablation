#ifndef PLIC_H
#define PLIC_H

#include "objet.h"
#include <vector>
#include <string>
#include "Sparse"
#include "Dense"
#include "recul.h"

class plic {
  public:
  plic();
  ~plic();
  void update(recul *recul){_recul=recul;};
  double grad_x(const int i, const int j);
  double grad_y(const int i, const int j);
  void interf();
  const Eigen::MatrixXd Get_interface() const {return _interface;};

  private:
  recul *_recul;
  read_data *_read_data;
  Eigen::MatrixXd _phi;
  double p,nx,nxx,ny;
  int _pointsupl,lon,lar,num,k;
  double dx,dz;
  //std::vector<std::vector<int> > _interface;   //_inteface[i][j,ax,ay,bx,by]
  Eigen::MatrixXd _ninterf, typinterf; //int  //_ninterf(i,j)=k   //numérotes des cases où il y a présence d'interface
  Eigen::MatrixXd _interface; //_interface[k][ax,ay,bx,by], où k est le numéro de la case, les coordonnées sont en locales, 0 en bas à gauche
  Eigen::MatrixXd tri,quad,penta;  //contient les somments pour chaque polygones
  Eigen::MatrixXd pttri,ptquad,ptpenta; //contient les coord des sommets
  Eigen::MatrixXd tri1,quad1,penta1; //outils
  // Sauvegarde la solution
	void SaveSol( int n);


/////////////////////RESTE A FAIRE///////////////////
//resize _interface
//SaveSol
};

#endif
