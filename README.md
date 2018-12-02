# Image Stabilization Method based on Variational Mode Decomposition and Relative Entropy (Entropy 2017)  
106-2 Digital Image Processing Final Project

## Requirement
 - MATLAB R2017a
 - [sepspyr](https://www.mathworks.com/matlabcentral/fileexchange/43909-separable-steerable-pyramid-toolbox?focused=3800382&tab=function)  

## Files
 - `src/README.m`, setup for `sepspyr` package and run main.m, main2.m 
 - `src/main.m`, implement image stabilization with noisy video 
 - `src/main2.m`, implement image stabilization with stable video
 - `report.pdf`, the report
 - `presentation.pptx`, the ppt for the final results
 - `src/data/`, input video
 - `src/output/`, output video
 - `src/frames/`, store frame data


## Usage

open MATLAB > set path to `src/` > run `README.m`

## Demo Video

[![2018 Spring DIP FinalDemo](https://i.imgur.com/DCue79D.png)](https://youtu.be/clO5tRXlWKA "2018 Spring DIP FinalDemo")


## References
 - [Hao, D., Li, Q., & Li, C. (2017). Digital Image Stabilization Method Based on Variational Mode Decomposition and Relative Entropy. Entropy, 19(11), 623.](https://www.mdpi.com/1099-4300/19/11/623)
 - [Dragomiretskiy, K., & Zosso, D. (2014). Variational mode decomposition. IEEE transactions on signal processing, 62(3), 531-544.](https://ieeexplore.ieee.org/document/6655981/)
 - [Mathematical statistics, the Kullback–Leibler divergence (also called relative entropy)](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)
 - [Scale-invariant feature transform](https://en.wikipedia.org/wiki/Scale-invariant_feature_transform)
