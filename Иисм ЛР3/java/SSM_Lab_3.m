clc;
clear all;
X = fscanf(fopen('out.txt','r'), '%f');
subplot(3,1,1);
histogram(X)

Y = fscanf(fopen('out1.txt','r'), '%f');
subplot(3,1,2);
histogram(Y)

Z = fscanf(fopen('out2.txt','r'), '%f');
subplot(3,1,3);
histogram(Z)