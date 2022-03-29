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
Zmat = repmat(z,1,n_a);
Amat = repmat(a_prime,1,n_z)';

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
cp_0 = r*Amat + w*Zmat ;
cp_0 = cp_0';

while dif > crit && iter < maxiter
    % using the guess, get the current consumption
    c0 = C0(cp_0,r,p_z);
    % get the current asset
    a0 = A0(Amat', Zmat', c0, r, w);
    % If any a0 is negative replace it with a_min
    %a0 = max(a0,a_min);
    cp_next = zeros(n_a,n_z);
    for nz = 1:n_z
    cp_next(:,nz) = interp1(a0(:,nz),c0(:,nz),Amat(nz,:)','spline');
    end

    dif = norm(cp_next - c0);
    cp_0 = cp_next; 



end
    


euler_rhs = beta*R*(c_0).^(-gamma);
c_tilde = euler_rhs.^(-1/gamma);

a_star = (c_tilde + Amat - w*Zmat)./R;

