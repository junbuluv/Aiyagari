clear all
clc




%%initialization 
alpha = 0.36;
beta = 0.985;
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
Zmat = repmat(z,1,n_a)';
Amat = repmat(a_prime,1,n_z);

Amat_dis = repmat(a_prime,1,n_z);



upper_r = 1/beta - 1;
lower_r = 0;
% initial guess r
r = (upper_r+lower_r)/2;
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
cp_0 = (r)*Amat + w*Zmat ;
dif = Inf;
tol = 1e-7;
iter = 0;
maxiter = 100000;


tic;
while dif > tol && iter < maxiter
    % using the guess, get the current consumption
    c0 = ((beta*(1+r)*(cp_0.^(-gamma))*p_z).^(-1/gamma));
    % get the current asset
    a0 = (Amat - w*Zmat + c0)/(1+r);
    % how to deal with negative a? When borrowing constraint binds 
    % get index for negative asset
    ind_bind = a0 <= a_min;
    c_bind = ind_bind.*a0*(1+r)+ind_bind.*(w*Zmat) - a_min;
    %Interpolation to get updated policy function
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
    toc


a0 = max(a0,0);

% get the policy function
for e = 1:n_z
    aux_func(:,e) = interp1(a0(:,e),Amat(:,e),"spline","extrap");
end
    policy_func(:,e) = aux_func
end




Amat_new = (1+r)*Amat_dis + w*Zmat - cp_0;
Zmat = repmat(z,1,n_a)';
Zmat = reshape(Zmat,[3500,1]);
Gmat = [Amat, Zmat];









