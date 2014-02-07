% Return the value
function y=gmm_prob(u1,u2,var1,var2,p1,p2,x)
y=p1.*gaussian(u1,var1,x)+p2.*gaussian(u2,var2,x);
