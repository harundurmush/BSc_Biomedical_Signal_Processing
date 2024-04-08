clc;
clear all;
close all;
%% uploading images
%%
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');
%% step a) Applying a 10-point moving average filter
%%
size_window = 10; % Defined window size of the filter
nom = ones(1,size_window); % nominator coefficient of the filter
denom = size_window; % denominator coefficient of the filter
ecg_hfn_filtered = filter(nom, denom, ecg_hfn); % ecg_hfn is from ecg_2.mat file
figure (1)
plot(ecg_hfn);
hold on;
plot(ecg_hfn_filtered);