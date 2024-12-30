#ifndef FINAL_PROJECT_HEADER_H
#define FINAL_PROJECT_HEADER_H

#include <iostream>
#include <cmath>

// Constants
const float G = 6.67430e-11;
const float GPa_to_Pa = 1e9;
const float EARTH_MASS = 5.972e24;
const float m_to_km = 0.001;

// Function declarations
float density(float P, float P_core_threshold, float P_crust_threshold);

float bisection_central_pressure(float total_radius, float total_mass, float tolerance, float P_core_frac);

void integrate_and_get_surface(float P_c, float P_core_frac, float &calculated_radius, float &calculated_mass);

#endif
