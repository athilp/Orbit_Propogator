function [output] = f(x)
b = 1;      
f1 = b/x(1)^2;
f2 = b/x(2)^4;
f3 = b/x(1)^2;
f4 = b/x(2)^4;
f5 = b/x(1)^2;
f6 = b/x(2)^4;
output = [f1;f2;f3;f4;f5;f6];
end