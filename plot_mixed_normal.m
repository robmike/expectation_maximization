p1=0.25;
p2=0.75;
u1=0;
u2=3;
var1=1;
var2=0.5;

x=-10:0.01:10;

g=@(x,mu,var) 1./sqrt(2*pi*var).*exp(-0.5.*(x-mu).^2./var);

f=p1.*g(x,u1,var1)+p2.*g(x,u2,var2);
plot(x,f);

print -dpng mixed_gaussian_pdf.png
