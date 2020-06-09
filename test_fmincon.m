function [W] = test_fmincon(A,B,x)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
W = A(1)*B(1)*x(1) + A(2)*B(2)*x(2);
W = W^2;
end

