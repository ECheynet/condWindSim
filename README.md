# Constrained wind simulation in a 2D plane
[![View condWindSim on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/128423-condwindsim)
[![Donation](https://camo.githubusercontent.com/a37ab2f2f19af23730565736fb8621eea275aad02f649c8f96959f78388edf45/68747470733a2f2f77617265686f7573652d63616d6f2e636d68312e707366686f737465642e6f72672f316339333962613132323739393662383762623033636630323963313438323165616239616439312f3638373437343730373333613266326636393664363732653733363836393635366336343733326536393666326636323631363436373635326634343666366536313734363532643432373537393235333233303664363532353332333036313235333233303633366636363636363536353264373936353663366336663737363737323635363536653265373337363637)](https://www.buymeacoffee.com/echeynet)

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
