clc;clear;close all;

A = importdata('l1_zncc.txt');
l1 = A(:,1);
zncc = A(:,2);
negZncc = -zncc;

norNegZncc = (negZncc - max(negZncc)) / (max(negZncc) - min(negZncc));
norL1 = (l1 - max(l1))/(max(l1) - min(l1));

figure;
plot(norNegZncc,'b');
hold on;
plot(norL1,'r')
legend('zncc','l1')