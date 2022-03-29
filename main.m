clear all
clc




%%initialization 
alpha = 0.36;
beta = 0.99;
delta = 0.025;
gamma = 2;
rho = 0.91;
sigma = 0.1;
n_z = 7;
n_a = 500;
a_min = 0;
a_max = 200;
n_dis = 700;




%% grid
% asset
a_prime = linspace(a_min,a_max,n_a)';
% labor skill shock
[z, p_z] = mytauchen(0,rho,sigma,n_z);
% constrcut G_matrix
z = exp(z);
%Amat = repmat(a_prime,n_z,1);
%Zmat = repmat(z,1,n_a)';
%Zmat = reshape(Zmat,[3500,1]);
%Gmat = [Amat, Zmat];
Zmat = repmat(z,1,n_a)';
Amat = repmat(a_prime,1,n_z);

% initial guess r
r = 1/beta -1;
% recover capital rental rate and wage
r_k = r + delta;
w = (r_k / alpha)^(alpha/(alpha-1)) * (1-alpha);

%endogenous grid method
%define marginal utility function 
up = @(x) (x).^(-gamma);
inv_up = @(y) y.^(-1/gamma);
%define Endogenous function
C0 = @(cp_0,r,pz) inv_up(beta * (1+r) * up(cp_0)*pz');
A0 = @(A_prime,z,c0,r,w) 1/(1+r)*(c0+A_prime-z.*w);


% note Amat is a_prime
% initial guess for c_prime
cp_0 = (1+r)*Amat + w*Zmat ;
dif = Inf;
tol = 1e-5;
iter = 0;
maxiter = 100000;

while dif > tol && iter < maxiter
    % using the guess, get the current consumption
    %c0_test = C0(cp_0,r,p_z);
    c0 = ((beta*(1+r)*p_z*cp_0'.^(-gamma)).^(-1/gamma))';
    % get the current asset
    a0 = Amat - w*Zmat + c0;
    % how to deal with negative a? When borrowing constraint binds
    ind_bind = a0 < 0;
    c_bind = ind_bind.*a0*(1+r)+ind_bind.*(w*Zmat);
    %a0 = A0(Amat', Zmat', c0, r, w);
    % If any a0 is negative replace it with a_min
    %a0 = max(a0,a_min);
    cp_next = zeros(n_a,n_z)+c_bind;
    for p = 1:n_z
    cp_next(:,p) = interp1(a0(:,p),c0(:,p),Amat(:,p),'spline');
    end
    
    dif = norm(cp_next - cp_0);
    if rem(iter,1000) == 0
      fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dif]);
    end

    cp_0 = cp_next;
    iter = iter+1;


end
    






