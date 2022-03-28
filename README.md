# Higher-order-Kuramoto-model

Code for simulating the higher-order Kuramoto model on simplicial complexes and on networks. 
Accompaining code for "Higher-order simplicial synchronization of coupled topological signal"

If you use any of these codes, please cite the following two  papers:

1. Ghorbanchian, R., Restrepo, J.G., Torres, J.J. and Bianconi, G., 2021. Higher-order simplicial synchronization of coupled topological signals. Communications Physics, 4(1), pp.1-13.  https://doi.org/10.1038/s42005-021-00605-4

2.Millán, A.P., Torres, J.J. and Bianconi, G., 2020. Explosive higher-order Kuramoto dynamics on simplicial complexes. Physical Review Letters, 124(21), p.218301. https://doi.org/10.1103/PhysRevLett.124.218301

This repository contains the following two folders:

Models: contains MATLAB and C codes to generate networks and simplicial complexes on which to run the higher-order Kuramoto model 
Kuramoto: Contains MATLB code of the  Euler algorithm to simulate the higher-order Kuramoto model on networks and simplicial complexes.

The folder Models contains two subfolders:

  simplicial complexes: generating simplicial complexes from NGFs and Configuration model of simplicial complexes.
                      The main (MATLAB) code generating all the parameters necessary to run the Kuramoto simulation is Parameters_simplices.m

  random graphs: generates networks and clique complex of Poisson networks and scale-free networks.
                The main (MATLAB) code generating all the parameters necessary to run the Kuramoto simulation is Parameters_configuration.m
               
The folder Kuramoto contains the MATLAB code to run the higher-order Kuramoto model on networks and simplicial complexes

Code for the higher-order Kuramoto model on networks:globalcode_SC.m
Code for the higher-order Kuramoto model on simplicial complexes:globalcode_network.m


Code for the projected dynamics: projected_sim_Poisson.m (on Poisson network) projected_sim_SF.m (on scale-free network) projected_sim_NGF.m (on NGF).



