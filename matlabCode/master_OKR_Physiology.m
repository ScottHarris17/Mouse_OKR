%Analysis for single cell recordings for OKR project

%CHANGE the project root if necessary (change for each machine)
physiologyRootPath = 'C:/Users/mrsco/Box/DataAndAnalysis/labData/OKR/physiology';
addpath(genpath(physiologyRootPath)); %add all subfolders to matlab path

%CHANGE the location of the list of all cells
cellListPath = 'C:/Users/mrsco/Box/DataAndAnalysis/labData/OKR/batch/';


%% Load Data
% User selects cells and center data to analyze. .mat files chosen should
% all live in the same directory within .../labData/OKR/physiology/
% User can also select a .h5 file that will cause the program to parse and
% export as a .mat file for each cell in the .h5 file. After that user will
% have the opportunity to select the desired .mat files again.
% Running load_OKR_Physiology will return a structure called data that
% contains a subfield for each cell that will be analyzed, along with the
% centerData.
load_OKR_Physiology

%% Build Epoch Log
% Determine which epochs were played for each cell and the order in which
% they occured. Add this to the data structure as a field for each cell
% called epochLog that contains a cell array. The first row of epochLog
% contains the protocol name and the second row contains the indexes of
% epochs that occured during that protocol run. You can change the
% parameters that will result in a split by editing the
% splitOnTheseParameters variable within the
% determineEpochsPerCell_OKR_Physiology.m file.
determineEpochsPerCell_OKR_Physiology

%% Fill Center Data
% Produce polar coordinates for each cell based on the center data file. If
% there is no centerData file selected then this function will get skipped
% over if coordinates already exist for this cell, otherwise coordinates
% will be listed as 'unknown'.
computeCoordinates_OKR_Physiology

%% Process Raw Traces
% Open dataViewerApp_OKR_Physiology.mlapp to view and process data for each
% cell at a time.

% THERE ARE TWO VERSIONS OF THE DATAVIEWER APP: One does an object oriented
% analysis, the other just uses a series of functions. Analyses are
% largely the same, but are stored in different files.
processData_OKR_Physiology


