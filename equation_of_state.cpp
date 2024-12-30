#include <cmath>

double density(double P, std::string material) {
    double rho0, c, n;

    if (material == "iron") { // IRON CORE 
        rho0 = 8300; // kg/m^3
        c = 1.0e-12; // Pa^(-n)
        n = 0.528;
    } else if (material == "silicate") { // SILICATE MANTLE
        rho0 = 4000; // kg/m^3
        c = 5.0e-12; // Pa^(-n)
        n = 0.544;
    } else if (material == "crust") { // CRUST
        rho0 = 2700; // kg/m^3
        c = 6.0e-12; // Pa^(-n)
        n = 0.544;
    } else {
        throw std::invalid_argument("Unknown material type");
    }

    return rho0 + c * pow(P, n);
}
