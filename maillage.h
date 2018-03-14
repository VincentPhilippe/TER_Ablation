#include "Dense"

class Maillage2DCarre
{
  protected:
    // nom du maillage
    std::string _name;
    // taille du maillage et raffinement
    double _deltaX, _deltaZ, _Lz, _Lx;
    // vecteurs coordonnées
    Eigen::VectorXd _coordX, _coordZ;

    // C'est Thomas tu peux rajouter Nx et Nz dans les attributs et mettre une méthode GetNx / GetNz

  public:
    // constructeur par défaut
    Maillage2DCarre();
    // constructeur 1
    Maillage2DCarre(std::string name, double deltaX, double deltaZ, double _Lx, double _Lz);
    // destructeur
    ~Maillage2DCarre();

    //fonctions d'accès aux variables de l'objet
    Eigen::VectorXd GetCoordX(){return _coordX;};
    Eigen::VectorXd GetCoordZ(){return _coordZ;};

};
