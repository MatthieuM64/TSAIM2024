# Two-species active Ising model

<a href="https://doi.org/10.5281/zenodo.14257074" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14257074.svg" alt="DOI"></a>

Codes used in the scientific publication: M. Mangeat, S. Chatterjee, J. D. Noh, and H. Rieger, <i>Emergent complex phases in a discrete flocking model with reciprocal and non-reciprocal interactions</i>, <a href='https://www.nature.com/articles/s42005-025-02098-x'>Commun. Phys. <b>8</b>, 186 (2025)</a>. A preprint is available on <a href='https://arxiv.org/abs/2412.02501'>arXiv</a>.</br></br>

For each model, a C++ code to compute the numerical simulations of the microscopic model and a C++ code to compute the numerical solutions of the hydrodynamic equations are available in this repository. Some additional Python codes are also available to generate movies of the system dynamics.</br></br>
<b>Exportations:</b> density snapshots and profiles shown in the different figures of the paper.</br>
<b>Compile:</b> g++ filename.cpp -fopenmp -lgsl -lgslcblas -lm -O3 -s -o filename.out.</br>
<b>Run:</b> ./filename.out -parameter=value.</br>
<b>Generate the movie:</b> python filename.py -parameter=value.</br>

## TSAIM

Codes for the reciprocal two-species active Ising model without species flip.</br></br>
<b>List of parameters for the numerical simulations</b> (<i>TSAIM_omp.cpp</i>): beta, D, v, theta, rho0, mag0, LX, LY, init, tmax, ran, threads (details as comments in the code).</br>
<b>Codes to generate the movies for the numerical simulations:</b> <i>figure_TSAIM_dynamics1d.py</i> and <i>figure_TSAIM_dynamics2d_rect.py</i>.</br>
<b>List of parameters for the numerical solutions of hydrodynamic equations</b> (<i>TSAIM_hydro_omp.cpp</i>): beta, epsilon, rho0, mag0, LX, init, dt, tmax, dx, threads (details as comments in the code).</br>
<b>Code to generate the movies for numerical solutions of hydrodynamic equations:</b> <i>figure_TSAIM_hydro_dynamics1d.py</i>.

## TSAIM_species

Codes for the reciprocal two-species active Ising model with species flip.</br></br>
<b>List of parameters for the numerical simulations</b> (<i>TSAIM_species_omp.cpp</i>): beta1, beta2, D, v, theta, rho0, mag0, gamma, LX, LY, init, tmax, ran, threads (details as comments in the code).</br>
<b>Codes to generate the movies for the numerical simulations:</b> <i>figure_TSAIM_species_dynamics1d.py</i> and <i>figure_TSAIM_species_dynamics2d_rect.py</i>.</br>
<b>List of parameters for the numerical solutions of hydrodynamic equations</b> (<i>TSAIM_species_hydro_omp.cpp</i>): beta1, beta2, epsilon, rho0, mag0, gamma, LX, init, dt, tmax, dx, threads (details as comments in the code).</br>
<b>Code to generate the movies for numerical solutions of hydrodynamic equations:</b> <i>figure_TSAIM_species_hydro_dynamics1d.py</i>.

## NRTSAIM

Codes for the non-reciprocal two-species active Ising model.</br></br>
<b>List of parameters for the numerical simulations</b> (<i>NRTSAIM_omp.cpp</i>): beta, D, v, theta, rho0, mag0, JAB, JBA, LX, LY, init, tmax, ran, threads (details as comments in the code).</br>
<b>Codes to generate the movies for the numerical simulations:</b> <i>figure_NRTSAIM_dynamics1d.py</i> and <i>figure_NRTSAIM_dynamics2d_rect.py</i>.</br>
<b>List of parameters for the numerical solutions of hydrodynamic equations</b> (<i>NRTSAIM_hydro_omp.cpp</i>): beta, epsilon, rho0, mag0, JAB, JBA, LX, init, dt, tmax, dx, threads (details as comments in the code).</br>
<b>Code to generate the movies for numerical solutions of hydrodynamic equations:</b> <i>figure_NRTSAIM_hydro_dynamics1d.py</i>.
