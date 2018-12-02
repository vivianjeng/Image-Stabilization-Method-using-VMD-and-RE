% DIP Final Project 
% June 26, 2018
% Team: 28
% Name: Liang Chien-wei, Jeng Ya-wen
% ID #: R06521521 R06922097
%###############################################################################% 
% Title: Digital Image Stabilization Method Based on 
%        Variational Mode Decomposition and Relative Entropy. (Entropy 2017)
% Setup for SIFT:
%       Download 'sepspyr' from
%       https://www.mathworks.com/matlabcentral/fileexchange/43909-separable-steerable-pyramid-toolbox?focused=3800382&tab=function
disp('Setup for "SIFT"...'); 
cd sepspyr/sepspyr/;
set_paths;
cd ../..;
disp('Done setup');
% 
% Implementation:
%       SIFT, feature matching, Calculating GMV, VMD, Relative Entropy
%
% main.m:
%       description: input a noisy video, and implement image stabilization
%       input: "data/IMG_3714_resize.mp4"
%       output: "output/3714_out1.avi"
%
disp('run stabilization without noising step'); 
main;
disp('Done sataibilization');
% main2.m:
%       description: input a stable video, output a noisy video, then 
%       implement image stabilization
%       input: "data/data1.mp4"
%       output: "output/data1_out.avi"
%
disp('run stabilization with noising step'); 
main2;
disp('Done sataibilization');
%###############################################################################%
 