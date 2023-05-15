clear all
close all
clc
%% Mission parameter
disp('loading mission parameter...');
data_current;
disp('done');

%% Mission supervisor
disp('loading supervisor parameter...');
mission;
disp('done');

%% Controller
disp('loading controller parameter...');
controller;
disp('done');

%% Model
disp('loading model parameter...');
model;
disp('done');

%% Sensors
disp('loading sensors parameter...');
sensors;
disp('done');

%% Navigation
disp('loading navigation parameter...');
navigation;
disp('done');
disp('Ready for the simulation...');