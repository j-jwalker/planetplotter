#include <iostream>
#include <cmath>
#include "cpgplot.h"
#include "final_project_header.h"

// Values for the polytropic equation of state based on material in each layer
const float rho_core = 8300.00; // Values for Fe (core)
const float c_core = 0.00349;
const float n_core = 0.528;

const float rho_mantleL = 4500.00; // Values for MgSiO3 (mantle)
const float c_mantleL = 0.0018;
const float n_mantleL = 0.5;

const float rho_crust = 3220.00; // Values for SiC (crust)
const float c_crust = 0.00172;
const float n_crust = 0.537;

// Equation of State
float density(float P, float P_core_threshold, float P_crust_threshold) {
    float rho_0, c, n;

    if (P > P_core_threshold) {  // Core
        rho_0 = rho_core;
        c = c_core;
        n = n_core;
    } else if (P > P_crust_threshold) {  // Mantle
        rho_0 = rho_mantleL;
        c = c_mantleL;
        n = n_mantleL;
    } else if (P > 0) {  // Crust
        rho_0 = rho_crust;
        c = c_crust;
        n = n_crust;
    } else {
        std::cerr << "Pressure = 0. Surface reached.\n";
        return -1;  // integrate_and_plot loop termination
    }

    return rho_0 + c * pow(P, n);
}

void integrate_and_plot(float P_c, float dr, float r_max, float P_core_frac) {
    // Initial conditions
    float r = 1e-6; // Initial radius (small but not 0) in meters
    float M = 0.0; // Initial mass in kg
    float P = P_c; // Initial central pressure in Pa

    // Threshold pressures
    float P_core_threshold = P_c * P_core_frac; // Core-mantle boundary
    float P_crust_threshold = P_c * 0.01; // Crust threshold

    // Arrays for storing values
    const int MAX_POINTS = 10000;
    float r_values[MAX_POINTS], M_values[MAX_POINTS], P_values[MAX_POINTS], rho_values[MAX_POINTS];
    int point_count = 0;

    // Integration loop
    while (r < r_max && P > 0 && point_count < MAX_POINTS) {
        float rho = density(P, P_core_threshold, P_crust_threshold);

        if (rho == -1) break; // Stop the loop if surface is reached

        // Store values
        r_values[point_count] = r / 1000.0; // Convert to km
        M_values[point_count] = M / EARTH_MASS; // Convert to Earth masses
        P_values[point_count] = P / GPa_to_Pa; // Convert to GPa
        rho_values[point_count] = rho;
        point_count++;

        // Runge-Kutta (RK4) for integration
        float k1_M = 4.0 * M_PI * r * r * rho;
        float k1_P = -G * M * rho / (r * r);

        float r_k2 = r + 0.5 * dr;
        float M_k2 = M + 0.5 * dr * k1_M;
        float P_k2 = P + 0.5 * dr * k1_P;
        float rho_k2 = density(P_k2, P_core_threshold, P_crust_threshold);
        float k2_M = 4.0 * M_PI * r_k2 * r_k2 * rho_k2;
        float k2_P = -G * M_k2 * rho_k2 / (r_k2 * r_k2);

        float r_k3 = r + 0.5 * dr;
        float M_k3 = M + 0.5 * dr * k2_M;
        float P_k3 = P + 0.5 * dr * k2_P;
        float rho_k3 = density(P_k3, P_core_threshold, P_crust_threshold);
        float k3_M = 4.0 * M_PI * r_k3 * r_k3 * rho_k3;
        float k3_P = -G * M_k3 * rho_k3 / (r_k3 * r_k3);

        float r_k4 = r + dr;
        float M_k4 = M + dr * k3_M;
        float P_k4 = P + dr * k3_P;
        float rho_k4 = density(P_k4, P_core_threshold, P_crust_threshold);
        float k4_M = 4.0 * M_PI * r_k4 * r_k4 * rho_k4;
        float k4_P = -G * M_k4 * rho_k4 / (r_k4 * r_k4);

        M += (k1_M + 2 * k2_M + 2 * k3_M + k4_M) * dr / 6.0;
        P += (k1_P + 2 * k2_P + 2 * k3_P + k4_P) * dr / 6.0;
        r += dr;
    }

    // PGPLOT Visualization
    if (cpgopen("/xwindow") > 0) {
        cpgscr(0, 1.0, 1.0, 1.0);
        cpgscr(1, 0.0, 0.0, 0.0);

        // Plot Mass vs Radius
        cpgenv(0, r_values[point_count - 1], 0, M_values[point_count - 1], 0, 0);
        cpglab("Radius (km)", "Mass (Earth Masses)", "Mass vs Radius");
        cpgline(point_count, r_values, M_values);


        // Plot Pressure vs Radius
        cpgenv(0, r_values[point_count - 1], 0, P_values[0], 0, 0);
        cpglab("Radius (km)", "Pressure (GPa)", "Pressure vs Radius");
        cpgline(point_count, r_values, P_values);


        // Plot Density vs Radius
        cpgenv(0, r_values[point_count - 1], 0, rho_values[0], 0, 0);
        cpglab("Radius (km)", "Density (kg/m^3)", "Density vs Radius");
        cpgline(point_count, r_values, rho_values);

        cpgclos();
    } else {
        std::cerr << "PGPLOT failed to open.\n";
    }
}

int main() {
    // Inform the user about the process
    std::cout << "We will begin by calculating the central pressure of the planet using observed values of the exoplanet through the bisection method.\n";
    std::cout << "Please provide the total radius, total mass, and core fraction to perform this calculation.\n";

    // User inputs for central pressure bisection calculation
    float total_radius, total_mass, P_core_frac;

    std::cout << "Enter the total radius of the planet (km) (Earth ~6300km): ";
    std::cin >> total_radius;
    total_radius *= 1000;  // Convert to meters

    std::cout << "Enter the total mass of the planet (Earth masses): ";
    std::cin >> total_mass;
    total_mass *= EARTH_MASS;  // Convert to kg

    std::cout << "Enter the core fraction of central pressure (A core fraction of .40 means we switch from a core-density regime to a mantle-density regime when central pressure has decreased to .40 of the initial central pressure. Essentially, an analogue for the core/mantle pressure ratio. In earth models, this valus is ~0.40): ";
    std::cin >> P_core_frac;

    // Tolerance for bisection method in Pa
    float tolerance = 1e6;
    // Derive central pressure using bisection method
    float P_c = bisection_central_pressure(total_radius, total_mass, tolerance, P_core_frac);

    // Display the derived central pressure and get user approval
    std::cout << "Derived central pressure: " << P_c / 1e9 << " GPa.\n";
    std::cout << "Do you approve this central pressure? (y/n): ";
    char approval;
    std::cin >> approval;

    if (approval != 'y' && approval != 'Y') {
        std::cout << "Program terminated.\n";
        return 0;  // Exit the program
    }

    // Proceed with step size and maximum radius inputs
    std::cout << "Enter the step size (km): ";
    float dr;
    std::cin >> dr;
    dr *= 1000;

    std::cout << "Enter the maximum radius to integrate (km): ";
    float r_max;
    std::cin >> r_max;
    r_max *= 1000;

    // Throw over to integrate_and_plot function
    integrate_and_plot(P_c, dr, r_max, P_core_frac);

    return 0;
}
