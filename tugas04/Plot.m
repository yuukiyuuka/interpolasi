clc
clear all

%data percobaan
dp = [10            0.01794775
15             0.03808997
20             0.05516225
20             0.05516225
25             0.05598281
30             0.04795629
35             0.04807485
40             0.06273566
45             0.07853982
50             0.07395442
55             0.04201338];

nn = size(dp);
n = nn(1,1);

for i=1:n
 x(i) = dp(i,1);
 y(i) = dp(i,2);
end

%data x untuk plot

max = 60;
min = 0;

n0 = 100; %jumlah data
x0 = min:(max-min)/(n0-1):max;

%%fungsi polinomial%%

%orde 3

a3 = [8.9837134848484990E-003
   1.3244783881118901E-003
   3.4878087878787730E-005
  -8.0858097902097699E-007];

y3 = zeros([n0 1]);

for i=1:n0
    for j=1:4
        y3(i) = y3(i) + a3(j)*x0(i)^(j-1);
    end
end

%orde 5

a5 = [-3.5578006393798829E-002
   1.0619961883994427E-003
   8.8021860020276108E-004
  -5.8623326993851915E-005
   1.3620461940544838E-006
  -1.0639517550537521E-008];

y5 = zeros([n0 1]);

for i=1:n0
    for j=1:6
        y5(i) = y5(i) + a5(j)*x0(i)^(j-1);
    end
end

%orde 7

a7 = [0.18647490708191194     
  -4.6318235667165010E-002
   4.0076386811539466E-003
  -8.9855952489531894E-005
  -3.2305315366267318E-006
   1.9128144895704346E-007
  -3.2528728395514216E-009
   1.8761882791405654E-011];

y7 = zeros([n0 1]);

for i=1:n0
    for j=1:8
        y7(i) = y7(i) + a7(j)*x0(i)^(j-1);
    end
end

%orde9

a9 = [-1.7400076958852716E-002
   1.5902090578251694E-002
  -3.3937470904475552E-003
   3.3515641672599142E-004
  -1.3651779377898548E-005
   1.1203196447739404E-007
   8.2873938863914494E-009
  -2.7449233577996040E-010
   3.3197095109260723E-012
  -1.4602161647973160E-014];

y9 = zeros([n0 1]);

for i=1:n0
    for j=1:10
        y9(i) = y9(i) + a9(j)*x0(i)^(j-1);
    end
end

plot(x,y,'.',x0,y3,x0,y5,x0,y7,x0,y9,MarkerSize=20,LineWidth=1)
title  ('Plot Data dan Interpolasinya','FontSize',15)
legend ('data percobaan','polinomial derajat 3','polinomial derajat 5','polinomial derajat 7','polinomial derajat 9','location','northwest','FontSize',10)
grid on