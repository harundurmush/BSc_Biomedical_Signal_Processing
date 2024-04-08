clc;
clear all;
close all;
%%
n=-50:50;
tic;
for i=0:0.01:10
    x=cos(pi*i*n);
    plot(n,x);
    title([" i = " num2str(i)])
    drawnow;
end
toc;