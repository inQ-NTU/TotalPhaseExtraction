close all; 
clear all;

load('g4_try.mat')

%Calculating M4 measure for each g4 function
for i = 1:1
    M4(i) = sum(abs(g4_full{i}-g4_disconnected{i}), 'all')/sum(abs(g4_disconnected{i}),'all');
end