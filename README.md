# A new coefficient to measure agreement between continuous variables

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![fastmatrix](https://img.shields.io/badge/fastmatrix-0.6-orange)](https://cran.r-project.org/package=fastmatrix)
[![MVT](https://img.shields.io/badge/MVT-0.3--81-orange)](https://cran.r-project.org/package=MVT)
[![L1pack](https://img.shields.io/badge/L1pack-0.60-orange)](https://cran.r-project.org/package=L1pack)
[![plot3D](https://img.shields.io/badge/plot3D-1.4.1-orange)](https://cran.r-project.org/package=plot3D)

Supplementary material to **A new coefficient to measure agreement between continuous variables** by Ronny Vallejos, Felipe Osorio and Clemente Ferrer

Code tested on R under development (2025-02-20 r87772), running Linux Mint 22.1 (64 bits)

Attached packages: fastmatrix 0.6, MVT 0.3-81, L1pack, plot3D

### Instructions: 
To create the Dynamically Loaded (DL) library, using the console prompt move to `/RNG` directory and enter:

`R CMD SHLIB -o RNG.so *.c -lblas -llapack -lm`

Next, copy `RNG.so` file to the working directory (in our case to `/simulation`), and execute the commands in `simul_L1ccc.R` file.
