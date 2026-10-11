
% To create a new "Domain.mat" file, follow these steps:
% 1: Run the first section (pdetool) and create the 2D domain in the app. 
% 2: Inside the PDE Modeler app, draw the domain and display the boundary.
%    Then, export the geometry (Boundary -> Export Decomposed Geometry, ...).
% 3: Exit the pdetool app (it is not necessary to save the file).
% 4: Choose a name for your "Domain.mat" file and run the second section.  

%% First section
clc , close all, clear all

pdetool;

%% Second section

Name_domain = 'Domain3.mat';
save(Name_domain,'b','g')


%%
