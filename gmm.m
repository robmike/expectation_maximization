% Output a random vector selected from a bi-modal Gaussian mixing model.
% Parameters are the means, variances, and probabilities of each mode and m
% is the number of samples to be generated.
% Vectorized.

function y=gmm(P,u,var,m)
p1=P(1);
var1=var(1); var2=var(2);
u1=u(1); u2=u(2);

x1=randn(1,m).*sqrt(var1)+u1; %Generate vector of RVs from first mixing model
x2=randn(1,m).*sqrt(var2)+u2;

%Create a vector that will 'select' a mixing model. p is an indicator variable
p=rand(1,m);
p(p<p1)=0;
p(p~=0)=1;

y=x1.*(1-p)+x2.*p;