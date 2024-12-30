#include <cmath>
#include "final_project_header.h"

float bisection_central_pressure(float total_radius, float total_mass, float tolerance, float P_core_frac) {
    // Define bounds for central pressure (reasonable starting guesses)
    float P_low = 1e9; // Lower bound (Pa)
    float P_high = 1e12; // Upper bound (Pa)
    float P_mid = (P_low + P_high) / 2.0;

    // Variables to hold calculated radius and mass
    float calculated_radius, calculated_mass;

    // Perform bisection
    while ((P_high - P_low) > tolerance) {
        P_mid = (P_low + P_high) / 2.0;  // Midpoint

        // Integrate outward from the guessed central pressure
        integrate_and_get_surface(P_mid, P_core_frac, calculated_radius, calculated_mass);

        // Compare calculated values to target values
        if (calculated_radius > total_radius || calculated_mass > total_mass) {
            // Guessed pressure is too high
            P_high = P_mid;
        } else {
            // Guessed pressure is too low
            P_low = P_mid;
        }
    }

    // Return the midpoint as the estimated central pressure
    return P_mid;
}

void integrate_and_get_surface(float P_c, float P_core_frac, float &calculated_radius, float &calculated_mass) {
    // Initial conditions
    float r = 1e-6;
    float M = 0.0;
    float P = P_c;

    float P_core_threshold = P_c * P_core_frac; // Core-mantle boundary
    float P_crust_threshold = P_c * 0.01; // Crust threshold
    float dr = 1000;

    // Integrate outward until pressure reaches 0
    while (P > 0) {
        float rho = density(P, P_core_threshold, P_crust_threshold);
        float dMdr = 4.0 * M_PI * r * r * rho;
        float dPdr = -G * M * rho / (r * r);

        // Update values
        M += dMdr * dr;
        P += dPdr * dr;
        r += dr;
    }

    // Set the outputs
    calculated_radius = r; // Final radius where pressure reaches 0
    calculated_mass = M; // Total enclosed mass
}
