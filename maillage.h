#include "Dense"

class Maillage2DCarre
{
  protected:
    // nom du maillage
    std::string _name;
    // taille du maillage et raffinement
    int _N, _nbPtsX, _nbPtsY;
    // pas d'espace
    double deltaX, deltaY;
    // matrice des concentrations
    Eigen::MatrixXd _Concentrations;

  public:
    // constructeur par défaut
    Maillage2DCarre();
    // constructeur 1
    Maillage2DCarre(std::string name, int nbPtsX, int nbPtsY, int N);
    // destructeur
    ~Maillage2DCarre();
}
