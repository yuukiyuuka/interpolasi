clc
clear all

%data percobaan

x0 = load('data_plot_x0.txt');
y0 = load('data_plot_y0.txt');

%data interpolasi lagrange
xl = load('data_plot_xl.txt');
yl = load('data_plot_yl.txt');

%data interpolasi lagrange kubik
xlk = load('data_plot_xlk.txt');
ylk = load('data_plot_ylk.txt');

%data interpolasi hermite kubik
xhk = load('data_plot_xhk.txt');
yhk = load('data_plot_yhk.txt');

plot(x0,y0,'.',xl,yl,xlk,ylk,xhk,yhk,MarkerSize=20,LineWidth=1.5)
title  ('Plot Data dan Interpolasinya','FontSize',15)
grid on
legend ('Data Percobaan','interpolasi lagrange','interpolasi lagrange kubik','interpolasi hermite kubik','location','northwest','FontSize',10)