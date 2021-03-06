%Copyright (C) 2022 by Frida Viset

current_folder=pwd;
if ispc
    addpath([current_folder,'\GPs']);
    addpath([current_folder,'\InertialNavigation']);
    addpath([current_folder,'\QuaternionAlgebra']);
    addpath([current_folder,'\RandomTrajectories']);
    addpath([current_folder,'\ReducedRankGPs']);
    addpath([current_folder,'\SimulateMeasurements']);
else
    addpath([current_folder,'/GPs']);
    addpath([current_folder,'/InertialNavigation']);
    addpath([current_folder,'/QuaternionAlgebra']);
    addpath([current_folder,'/RandomTrajectories']);
    addpath([current_folder,'/ReducedRankGPs']);
    addpath([current_folder,'/SimulateMeasurements']);
end