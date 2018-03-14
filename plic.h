#include <vector>
#include <string>
#include "Sparse"
#include "Dense"
#ifndef PLIC_H
#define PLIC_H

class plic {
  public:
  plic(recul& _recul);
  void update(recul& recul){_recul=recul};
  double grad_x(const int i, const int j);
  double grad_y(const int i, const int j);
  void interf(const int lon,const int lar);
  const Eigen::MatrixXd Get_interface() const {return _interface;};
  const Eigen::MatrixXd Get_ninterf() const {return _ninterf;};

  private:
  recul& _recul;
  Eigen::MatrixXd _phi;
  double p,nx,nxx,ny;
  int k,_pointsupl,lon,lar
  //std::vector<std::vector<int> > _interface;   //_inteface[i][j,ax,ay,bx,by]
  Eigen::MatrixXd _ninterf, typinterf; //int  //_ninterf[i][j]=k   //numérotes des cases où il y a présence d'interface
  Eigen::MatrixXd _interface; //_interface[k][ax,ay,bx,by], où k est le numéro de la case

  // Sauvegarde la solution
	//void SaveSol( int n);


/////////////////////RESTE A FAIRE///////////////////
//resize _interface
//SaveSol
}

#ENDIF
