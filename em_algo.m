% Find the expectation maximization of a bi-modal Gaussian mixing model.
% With few modificatons the program can be extended to arbritary number of
% modes. (Change the repmats from 2->J) where J is number of modes.
% Warning: the EM algorithm depends on reasonable guesses for the initial
% parameters. The algorithm will fail if initial guesses are not good.
% Take note that x is required to be a column vector while P, u and var
% must be row vectors.
% Plots are produced.

actual.P=[0.25 0.75];
actual.u=[0 3];
actual.var=[1 0.5];

eps=1e-5; %Allowable error for convergence
eps=eps.*ones(1,6); %Error vector (for each component of theta)
max_iterations=10000;
N=1000; %Number of samples

L=zeros(1,100); %Pre-allocate memory

%Generate samples
x=gmm(actual.P,actual.u,actual.var,N)'; %Transpose to get column vector

% Random Initial guesses (If the guess is not good it will not converge!)
%--------------------------
u0_max=5;
u01=rand*u0_max-u0_max/2; u02=rand*u0_max-u0_max/2;
var0_max=5; var01=rand*var0_max; var02=rand*var0_max;
p=rand;
this.P=[p 1-p];
this.u=[u01 u02];  %Generate random numbers from -a/2..a/2
this.var=[var01 var02];
%--------------------------

%Main loop
%--------------------------
tic;
for i=1:max_iterations

last=this; % Backup last guesses so we can compare against them later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expectation step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% w is an Nx2 matrix, where N is the number of samples, corresponding to 
% the probability that sample x_n belongs to the first or second mixture.
% p_x is the probability matrix (IxJ) of each sample, x_i, 
% occuring given that the sample comes from mixture model j and
% parameterized by current estimates of the mean and variance for mixture 
% model j. J=2

p_x=gaussian(this.u,this.var,x);
%Probability of RV X coming from mixture model j and having value x
p_x_and_mm=repmat(last.P,N,1).*p_x;
w=p_x_and_mm./repmat(sum(p_x_and_mm,2),1,2); %Sum is row-wise (index j)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximization step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find new estimates
%-------------------
this.P=real(1/N*sum(w,1)); %Sum is column wise (index i)
this.u=real(1./(N.*this.P).*sum(w.*repmat(x,1,2),1));
x_c=real(repmat(x,1,2)-repmat(this.u,N,1)); %x_c=x-estimateof{mean(x)}
this.var=real(1./(N.*this.P).*sum(w.*x_c.^2,1)); %Sum against index i (columnwise);
%-------------------


this.var(this.var==0)=0.001; % Ensure that variance is not zero

phi(i)=this; % Store the guesses so we can plot mixture pdf at each step
%Exit out of loop if error<epsilon 

% Store the incomplete likelihood for this iteration
L(i)=sum(log(sum(repmat(this.P,length(x),1).*gaussian(this.u,this.var,x),2)));

error=abs([this.P this.u this.var]-[last.P last.u last.var]);
	if (i~=1 && max(error<eps))
		break; 
    end
    
pause(0); % Hack. Makes Matlab more responsive to Ctrl-C in Windows.
end %End 'for i=1:max_iterations'

toc
total_iterations=i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display EM-predictions and error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e1=real([this.P-actual.P this.u-actual.u this.var-actual.var]);
e2=real([fliplr(this.P)-actual.P fliplr(this.u)-actual.u fliplr(this.var)-actual.var]);
Predicted_values=real([this.P this.u this.var])
% Error depends on ordering of predicted values (i.e. which one is first, 
% mixing model 1 or 2). Choose the error value which is smallest
% (corresponding to the 'correct' ordering.
if(sum(e1.^2)<sum(e2.^2))
    em_error=e1
else em_error=e2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the mixture model at selected iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_plots=6; %Make sure this is even
r=-6:0.01:8; r=r'; % x coordinates of plot. 'gaussian' takes col vector
step_size=floor(total_iterations/num_plots);
%Generate range of iterations to plot. Logarithmically space because
%convergence does not happen linearly -> alot of change at beginning,
%little at end.
indexes=round(logspace(0,log10(total_iterations),6))
for i=indexes
        subplot(3,2,find(indexes==i));
        f=repmat(phi(i).P,length(r),1).*gaussian(phi(i).u,phi(i).var,r);
        f=sum(f,2);
        plot(r,f);
        title(['After ' int2str(i) ' iterations']);
        axis([-6 8 0 0.45]);
end

% Compare plots of predicted and actual PDFs
f1=sum(repmat(this.P,length(r),1).*gaussian(this.u,this.var,r),2);
f2=sum(repmat(actual.P,length(r),1).*gaussian(actual.u,actual.var,r),2);
figure;
plot(r,f1,'r:', r,f2,'b--');
title('EM-predicted vs. Actual GMM PDFs');
legend('Em-predicted','Actual',2);

%figure;
%plot(r,abs(f1-f2)); title('Error between predicted and actual PDFs');

figure; plot(L(1:total_iterations));
title('Log-likelihood of incomplete data')
xlabel('Iteration number');