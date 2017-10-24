% Phong Shading Demo: 
%  This demo highlights the interaction of light with 2 
%  spherical objects with different chromatic and
%  material properties. The sphere that appears
%  in the front has RGB = [1.00 0.63 0.4]
%  and the one in the back has RGB = [1.0 0 0]
%
%  Variables used here follow the formula given
%  on page 13 of the lecture slides. 
%            
% Calls:
%  phongShade.m
% Authors: ADJ CSC
% Fall 2001

% clean up workspace and close all figures
clear all;
close all;

% white light shines on 2 spheres
lightColor = [1 1 1]; 

%surface material = 'metal'
surfaceType = 'metal';
ka  = 0.1; % ambient  reflection coefficient
kd  = 0.1; % diffuse  reflection coefficient
ks  = 1.0; % specular reflection coefficient
ke  = 5.0; % spectral exponent
scr = 0.5; % reflected light a combination of 
           % illuminant and surface color
phongShade(surfaceType, lightColor, ka, kd, ks, ke, scr);
  
% surfaceType = shiny
surfaceType = 'shiny';
ka  = 0.1; 
kd  = 0.6; 
ks  = 0.7; 
ke  = 5.0; 
scr = 1.0; % reflected light pure illuminant color
phongShade(surfaceType, lightColor, ka, kd, ks, ke, scr);

% surfaceType = diffuse
surfaceType = 'diffuse';
ka  = 0.1; 
kd  = 0.7; 
ks  = 0.0; 
ke  = 1.0; % since ks = 0, the exact values for
scr = 1.0; % ke and scr do not matter
phongShade(surfaceType, lightColor, ka, kd, ks, ke, scr);

% surfaceType = ambient, observe complete lack of
% 3D information from the spheres. why is that?
surfaceType = 'ambient';
ka  = 1.0; 
kd  = 0.0; 
ks  = 0.0; 
ke  = 1.0; % since ks = 0, the exact values for
scr = 1.0; % ke and scr do not matter
phongShade(surfaceType, lightColor, ka, kd, ks, ke, scr);

% QUESTIONS:
%  Why is there no inter-reflection?
%  How are r(\lambda), I(\lambda) computed?
%  Check default values for reflection coefficients with 
%  the following commands in matlab:
%
%   help material
%   material('metal')
%   material('shiny')
%   material('dull')

