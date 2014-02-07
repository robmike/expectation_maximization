% 1-d gaussian (normal with mean mu, variance var
% mu and var must be row vectors and x a column vector
function f=gaussian(mu,var,x)
l=length(x);
% Augment x, mu and var into  Nx2 matrices
k=length(var);
x=repmat(x,1,k);
mu=repmat(mu,l,1); var=repmat(var,l,1);

f=1./sqrt(2*pi*var).*exp(-0.5.*(x-mu).^2./var);
