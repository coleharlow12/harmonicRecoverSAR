
a=100;
b=5;
repeat = 10^4;

a = rand(repeat,1);
b = rand(repeat,1);
tic
%a.*b;
%a+b;
exp(a);
toc
