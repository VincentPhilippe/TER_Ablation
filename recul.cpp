#include "recul.h"
#include "Dense"

using namespace std;
using namespace Eigen;

void recul_surface(Eigen::MatrixXd Surface, Eigen::MatrixXd C_solide, double dt, double dx, double dy)
{
  //N_surface, Nx, Ny
  double N_surface, Nx, Ny;

  N_surface = Surface.rows();
  Nx = C_solide.cols();
  Ny = C_solide.rows();

  //boucle sur les surfaces
  for(int k=0; k<N_surface; k++)
  {
    double i, j;

    i = Surface(k,0);
    j = Surface(k,1);

    //calcul de l'angle alpha
    double xa, ya, xb, yb, t_alpha, alpha;

    xa = Surface(k,2);
    ya = Surface(k,3);
    xb = Surface(k,4);
    yb = Surface(k,5);
    t_alpha = (yb-ya)/(xb-xa);
    alpha = atan(t_alpha);

    //calcul des coordonnées des points C et D
    double xc, yc, xd, yd, vr, vrdt;

    vr = Surface(k,6);
    vrdt = vr*dt;
    xc = xa + vrdt*sin(alpha);
    yc = ya - vrdt*cos(alpha);
    xd = xb + vrdt*sin(alpha);
    yd = yb - vrdt*cos(alpha);

    //identification du cas et modification du tableau des concentrations en solide
    if (xc<0) {
      if (xd<0) {
        if (yc<0) {//yd forcément négatif
          /* code *///cas 3
        } else {//yc>0
          if (yd<0) {
            /* code *///cas 2
          } else {//yd>0
            /* code *///cas 1
          }
        }
      } else {//0<xd<dx
        if (yc<0) {//yd forcémentnégatif
          /* code *///cas 9
        } else {//yc>0
          if (yd<0) {
            /* code */ //attention 2 cas possibles cas 7 et 8
          } else {//yd>0
            /* code *///cas 4
          }
        }
      }
    } else if (xc<dx) {
      if (xd<dx) {
        if (yc<0) {
          if (yd<0) {
            /* code *///cas 11
          } else {
            /* code *///cas 6
          }
        } else {
          if (yd<0) {
            /* code *///cas 10
          } else {
            /* code *///cas 5
          }
        }
      } else {
        if (yc<0) {
          if (yd<0) {
            /* code *///cas 17
          } else {
            /* code *///cas 13 et 14
          }
        } else {
          /* code *///cas 12
        }
      }
    } else {
      if (yc<0) {
        if (yd<0) {
          /* code *///cas 18
        } else {
          /* code *///cas 16
        }
      } else {
        /* code *///cas 15
      }
    }

  }
}

//////////////////////////////////////////////////à remplir
Eigen::MatrixXd recul1(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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

/////////////////////////////////////////////////////////à remplir
Eigen::MatrixXd recul2(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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

//////////////////////////////////////////////////////////////à remplir
Eigen::MatrixXd recul3(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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

//C à gauche cas 4
Eigen::MatrixXd recul4(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul5(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul6(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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

//C à gauche et D en dessous cas 7
Eigen::MatrixXd recul7(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
    C_solide(i-1,j)=C_solide(i,j-1)-S2/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S1-S2)/(dx*dy);

  return C_solide;
}

//C à gauche et D en dessous, depasse dans le coin cas 8
Eigen::MatrixXd recul8(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
    C_solide(i-1,j)=C_solide(i,j-1)-(S2-S3)/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S1-S2+S3)/(dx*dy);

  return C_solide;
}

////////////////////////////////////////////////////////////////////////à remplir
Eigen::MatrixXd recul9(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul10(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
    C_solide(i-1,j)=C_solide(i,j-1)-S2/(dx*dy);
  }
  C_solide(i,j)=C_solide(i,j)-(l*vrdt-S2)/(dx*dy);

  return C_solide;
}

///////////////////////////////////////////////////////////////////////////////////à remplir
Eigen::MatrixXd recul11(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul12(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul13(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul14(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul15(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul16(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul17(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
Eigen::MatrixXd recul18(Eigen::MatrixXd C_solide, int i, int j, double alpha, double vrdt, Eigen::MatrixXd coord, double dx, double dy, int nx)
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
