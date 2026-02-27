# planetplotter

A C++ simulation of the internal structure of terrestrial exoplanets. Given observable 
properties of a planet (mass and radius), the program estimates central pressure using 
the bisection method, then integrates outward from the core using fourth-order 
Runge-Kutta to model how density, pressure, and enclosed mass vary with radius.

## Background

Exoplanet interiors are impossible to observe directly. Mass and radius can only be 
inferred from transit data. This model bridges that gap using a modified polytropic 
equation of state derived from Seager et al. (2007), with material constants for an 
iron core, perovskite mantle, and SiC crust. Validated against Earth-analogue inputs, 
producing a central pressure of ~363 GPa and core radius of ~3400 km.

## How It Works

1. User provides total radius, mass (in Earth units), and a core fraction parameter
2. Bisection method estimates central pressure to within 1 MPa tolerance
3. RK4 integration propagates mass and pressure outward from center to surface
4. PGPLOT renders density, pressure, and mass profiles vs. radius

## Build & Run

**Dependencies:** PGPLOT, X11
```bash
make
./final_project
```

## Key Parameters

| Parameter | Description | Earth value |
|---|---|---|
| Total radius (km) | Observed from transit | ~6300 |
| Total mass (Earth masses) | Observed from transit | 1.0 |
| Core fraction | Core/mantle pressure boundary | ~0.40 |

## References

Hales, A. L., and J. L. Roberts. 1970. "Shear Velocities in the Lower Mantle and the Radius 
of the Core." *Bulletin of the Seismological Society of America* 60 (5): 1427–36. 
https://doi.org/10.1785/BSSA0600051427

Seager, S., M. Kuchner, C. A. Hier-Majumder, and B. Militzer. 2007. "Mass-Radius 
Relationships for Solid Exoplanets." *The Astrophysical Journal* 669 (2): 1279. 
https://doi.org/10.1086/521346
