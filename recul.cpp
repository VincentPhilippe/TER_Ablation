#include "recul.h"
#include "Dense"

using namespace std;
using namespace Eigen;

recul::recul(readdata& data, diffusion& diffusion, plic& plic)
{
  _data=data;
  _diff=diff;
  _plic=plic;
  _dt=_data.getdt();
  _dx=_data.getdx();
  _dz=_data.getdz();
  _ninterf=MatrixXd::Zero(1,1);
  _interface=MatrixXd::Zero(1,1);
  _vitesse=VectorXd::Zero(1);
  _C_solide=MatrixXd::Zero(1,1);
  _nx=1
  _nz=1
}

recul::~recul()
{}

void recul::recul_surface()
{
  //récuperer _ninterf _interface _vitesse avec des get
  //N_surface, Nx, Ny
  _nx = _ninterf.cols();
  _nz = _ninterf.rows();
  //boucle sur les surfaces
  for(int i=0; i<_nx; i++)
  {
    for(int j=0; j<_nz; j++)
    {
      double k;

      k = _ninterf(i,j);

      if (k>0) {

        //calcul de l'angle alpha
        double xa, ya, xb, yb, t_alpha, alpha;

        xa = _interface(1,k-1);
        ya = _interface(2,k-1);
        xb = _interface(3,k-1);
        yb = _interface(4,k-1);
        t_alpha = (yb-ya)/(xb-xa);
        alpha = atan(t_alpha);

        //calcul des coordonnées des points C et D
        double xc, yc, xd, yd, vr, vrdt;

        vr = _vitesse(k);
        vrdt = vr*dt;
        xc = xa + vrdt*sin(alpha);
        yc = ya - vrdt*cos(alpha);
        xd = xb + vrdt*sin(alpha);
        yd = yb - vrdt*cos(alpha);

        MatrixXd coord;
        coord.resize(4,2);
        coord(1,1)=xa;
        coord(1,2)=ya;
        coord(2,1)=xb;
        coord(2,2)=yb;
        coord(3,1)=xc;
        coord(3,2)=yc;
        coord(4,1)=xd;
        coord(4,2)=yd;

        //identification du cas et modification du tableau des concentrations en solide
        if (xc<0) {
          if (xd<0) {
            if (yc<0) {//yd forcément négatif
              recul3(i, j, alpha, vrdt, coord);
            } else {//yc>0
              if (yd<0) {
                recul2(i, j, alpha, vrdt, coord);
              } else {//yd>0
                recul1(i, j, alpha, vrdt, coord);
              }
            }
          } else {//0<xd<dx
            if (yc<0) {//yd forcémentnégatif
              recul9(i, j, alpha, vrdt, coord);
            } else {//yc>0
              if (yd<0) {
                /* code */ //////////////////////////////////////////////////////attention 2 cas possibles cas 7 et 8
              } else {//yd>0
                recul4(i, j, alpha, vrdt, coord);
              }
            }
          }
        } else if (xc<dx) {
          if (xd<dx) {
            if (yc<0) {
              if (yd<0) {
                recul11(i, j, alpha, vrdt, coord);
              } else {
                recul6(i, j, alpha, vrdt, coord);
              }
            } else {
              if (yd<0) {
                recul10(i, j, alpha, vrdt, coord);
              } else {
                recul5(i, j, alpha, vrdt, coord);
              }
            }
          } else {
            if (yc<0) {
              if (yd<0) {
                recul17(i, j, alpha, vrdt, coord);
              } else {
                /* code *///////////////////////////////////////////////////////////////////////cas 13 et 14
              }
            } else {
              recul12(i, j, alpha, vrdt, coord);
            }
          }
        } else {
          if (yc<0) {
            if (yd<0) {
              recul18(i, j, alpha, vrdt, coord);
            } else {
              recul16(i, j, alpha, vrdt, coord);
            }
          } else {
            recul15(i, j, alpha, vrdt, coord);
          }
        }

      }
    }
  }
}


void recul::recul1(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double Stot, S1, S2, l1, l2;
  Stot=vrdt*l;
  l1 = -xc*sin(alpha);
  l2 = -xd*sin(alpha);
  S1 = l*(l1+l2)/2;
  S2 = Stot-S1;
  C_solide(i,(j-1+nx)%nx)=C_solide(i,(j-1+nx%nx))-S1/(dx*dy);
  C_solide(i,j)=C_solide(i,j)-(S2)/(dx*dy);

  return C_solide;
}


void recul::recul2(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S3 = xb*xb/(2*tan(alpha));
  S1 = -xc*yc-yc*yc/(2*tan(alpha))+xc*xc/(2*tan(alpha));
  S2 = xd*yd-xd*xd/(2*tan(alpha))+yd*yd/(2*tan(alpha));
  S4 = Stot-(S1+S2+S3);
  C_solide(i,(j-1+nx)%nx)=C_solide(i,(j-1+nx%nx))-S1/(dx*dy);
  C_solide(i,j)=C_solide(i,j)-(S4)/(dx*dy);
  if (i-1>=0) {
    C_solide(i-1,(j-1+nx)%nx)=C_solide(i,(j-1+nx)%nx)-S2/(dx*dy);
    C_solide(i-1,j)=C_solide(i-1,j)-S3/(dx*dy);
  }

  return C_solide;
}


void recul::recul3(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double Stot, S1, S2, S3, S4;
  Stot=vrdt*l;
  S1 = ya*ya*tan(alpha)/2;
  S4 = ya*xb/2;
  S3 = xb*xb/(2*tan(alpha));
  S2 = Stot-(S1+S3+S4);
  C_solide(i,(j-1+nx)%nx)=C_solide(i,(j-1+nx%nx))-S1/(dx*dy);
  C_solide(i,j)=C_solide(i,j)-(S4)/(dx*dy);
  if (i-1>=0) {
    C_solide(i-1,(j-1+nx)%nx)=C_solide(i,(j-1+nx)%nx)-S2/(dx*dy);
    C_solide(i-1,j)=C_solide(i-1,j)-S3/(dx*dy);
  }

  return C_solide;
}

//C à gauche cas 4
void recul::recul4(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double S1;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  C_solide(i,(j-1+nx)%nx)=C_solide(i,(j-1+nx%nx))-S1/(dx*dy);
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S1)/(dx*dy);

  return C_solide;
}

//C et D dans Cij cas 5
void recul::recul5(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  C_solide(i,j)=C_solide(i,j)-l*vrdt/(dx*dy);

  return C_solide;
}

////////////////////////////////////////////////////////////////à remplir
void recul::recul6(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double S2;

  S2=yd*yd*(1/tan(alpha)+tan(alpha))/2;
  if (i-1>=0) {
    C_solide(i-1,j)=C_solide(i-1,j)-S2/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S2)/(dx*dy);

  return C_solide;
}

//C à gauche et D en dessous cas 7
void recul::recul7(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double S1, S2;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  S2=yd*yd*(1/tan(alpha)+tan(alpha))/2;
  C_solide(i,(j-1+nx)%nx)=C_solide(i,(j-1+nx)%nx)-S1/(dx*dy);
  if (i-1>=0) {
    C_solide(i-1,j)=C_solide(i-1,j)-S2/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S1-S2)/(dx*dy);

  return C_solide;
}

//C à gauche et D en dessous, depasse dans le coin cas 8
void recul::recul8(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double S1, S2, S3;

  S1=xc*xc*(1/tan(alpha)+tan(alpha))/2;
  S2=yd*yd*(1/tan(alpha)+tan(alpha))/2;
  S3=(abs(xc)*tan(alpha)-yc)*(abs(xc)*tan(alpha)-yc)/(2*tan(alpha));

  C_solide(i,(j-1+nx)%nx)=C_solide(i,(j-1+nx)%nx)-(S1-S3)/(dx*dy);
  if (i-1>=0) {
    C_solide(i-1,j)=C_solide(i-1,j)-(S2-S3)/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S1-S2+S3)/(dx*dy);

  return C_solide;
}

////////////////////////////////////////////////////////////////////////à remplir
void recul::recul9(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

//D en dessous cas 10
void recul::recul10(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  double S2;

  S2=yd*yd*(1/tan(alpha)+tan(alpha))/2;
  if (i-1>=0) {
    C_solide(i-1,j)=C_solide(i-1,j)-S2/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S2)/(dx*dy);

  return C_solide;
}

///////////////////////////////////////////////////////////////////////////////////à remplir
void recul::recul11(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

/////////////////////////////////////////////////////////////////////////////////////////à remplir
void recul::recul12(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

/////////////////////////////////////////////////////////////à remplir
void recul::recul13(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

//////////////////////////////////////////////////////////////////à remplir
void recul::recul14(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

//////////////////////////////////////////////////////////////////////à remplir
void recul::recul15(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

///////////////////////////////////////////////////////////////////////////////à remplir
void recul::recul16(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

////////////////////////////////////////////////////////////////////////////////à remplir
void recul::recul17(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}

/////////////////////////////////////////////////////////////////////////////////////à remplir
void recul::recul18(int i, int j, double alpha, double vrdt, MatrixXd coord)
{
  double xa,ya,xb,yb,xc,yc,xd,yd,l;

  xa=coord(1,1);
  ya=coord(1,2);
  xb=coord(2,1);
  yb=coord(2,2);
  xc=coord(3,1);
  yc=coord(3,2);
  xd=coord(4,1);
  yd=coord(4,2);

  l=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));

  //à compléter

  return C_solide;
}
