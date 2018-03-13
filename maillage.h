#include "Dense"

class Maillage2DCarre
{
  protected:
    // nom du maillage
    std::string _name;
    // taille du maillage et raffinement
    int _deltaX, _deltaZ, _Lz, _Lx;
    // vecteurs coordonnées
    Eigen::MatrixXd _coordonnees;

    // C'est Thomas tu peux rajouter Nx et Nz dans les attributs et mettre une méthode GetNx / GetNz

  public: // UN CONSTRUCTEUR N'A PAS DE VOID DEVANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // constructeur par défaut
    Maillage2DCarre();
    // constructeur 1
    Maillage2DCarre(std::string name, double deltaX, double deltaZ, double _Lx, double _Lz);
    // destructeur
    ~Maillage2DCarre();
}
