# Constrained wind simulation in a 2D plane
[![View condWindSim on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/128423-condwindsim)
<a href="https://www.buymeacoffee.com/echeynet" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 25px !important;width: 120px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>

## Summary
 
Matlab code for a conditional/constrained wind simulation in a 2D plane. The coherence function is modelled using the Davenport model. Only one velocity component is simulated at a time. It is possible to generate non-stationary event with time-varying coherence, but the simulation time will be significantly increased.
 
 
 <img src="./illustration.jpg" width="50%" height="50%" />
 
 
## Contents

The repository contains

  - coherence.m: Calculates the coherence for two signals.
  - condWindSim2D.m: Generates a 2D conditional/constrained wind field.
  - The Matlab Live Script Documentation.mlx, which documents how to use the function in 2D
  - The Matlab Live Script Documentation2.mlx, which documents how to use the function in 1D fr a vertical lines
  - The Matlab Live Script Documentation2.mlx, which documents how to use the function in 1D fr a horizontal lines
  


This is the first version of the repository, there may be many bugs. Do not hesitate to contact me if you find any.
