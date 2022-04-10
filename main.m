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




upper_r = 1/beta - 1;
lower_r = 0;
bisection_dif = Inf;
bisection_tol = 10;

%while bisection_dif > bisection_tol
%% outer loop
%% grid
% asset
a_prime = linspace(a_min,a_max,n_a)';
% labor skill shock
[z, p_z] = mytauchen(0,rho,(sigma^2/(1-rho^2)),n_z);
% constrcut G_matrix
z = exp(z);
Zmat = repmat(z,1,n_a)'; % z matrix 500 x 7 (z1,z2,...z7 columns)
A_pZmat = meshgrid(a_prime,z)'; % meshgrid of a' and z
%% inner loop simulation grid
% initial guess r
r = (upper_r+lower_r)/2;

%while bisection_dif > bisection_tol 
%% inner loop

% recover capital rental rate and wage
r_k = r + delta;
k_d =((r_k)/alpha)^(1/(alpha-1));


w = (r_k / alpha)^(alpha/(alpha-1)) * (1-alpha);
M_f = (1+r)*A_pZmat + w * Zmat; % total available resource
ap_endo = zeros(n_a,n_z); % initialize ap_endo grid
c_0 = M_f - A_pZmat; % initial guess of consumption 
env_old = (1+r)*c_0.^(-gamma); % initial guess of envelope condition

% convergence criterion
criter = 1e-9;
dif = Inf;
iter = 0;

while dif > criter
c_g = (beta* env_old*p_z').^(-1/gamma); %c_g choice
M_g = c_g + A_pZmat; % market resource choice
a_endo = (M_g - w* Zmat)./(1+r); % a_endogenous grid
% ap_endogenous using interpolation sample point(M_g), Query point(M_f) and function (A_prime)
for i = 1:n_z
ap_endo(:,i) = max(0,interp1(M_g(:,i),A_pZmat(:,i),M_f(:,i),"spline"));
end
% get c_endogenous (updated c_current)
c_endo = M_f - ap_endo;
% get new envelope condition
env_new = (1+r)*c_endo.^(-gamma);

dif = norm(env_old - env_new);
    if rem(iter,250) == 0
      fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dif]);
    end
    iter = iter +1;
env_old = env_new;
end


%%% get Asset future using policy function from EGM
A_dis = linspace(a_min,a_max,n_dis)';

Ag_dis = meshgrid(A_dis,z)';

AP_dis = meshgrid(zeros(n_dis,1),z)';

for i = 1:n_z
AP_dis(:,i) = interp1(A_pZmat(:,i),ap_endo(:,i),A_dis');
end


% construct the tensor grid and vector of future asset 

Amat_dis = meshgrid(A_dis,z)';
Amat_dis = reshape(Amat_dis,[n_dis*n_z,1]);
Zmat_dis = repmat(z,1,n_dis)';
Zmat_dis = reshape(Zmat_dis,[n_dis*n_z,1]);
G_grid = [Amat_dis, Zmat_dis];
ap_grid = reshape(AP_dis,[n_dis*n_z,1]);





%% get the index for below
%% construct gamma matrix working in progress
P_a = zeros(n_z*n_dis, n_z*n_dis);
Interval_G = G_grid(2,1) - G_grid(1,1);
for m = 1: n_z*n_dis
for n = 1: n_z*n_dis
if mod(n,n_dis) ==1
    P = 1 - (ap_grid(m,1)-G_grid(n,1))/Interval_G;
    P_a(m,n) = min(1.0, max(0.0, P));
elseif ( ap_grid(m,1) >= G_grid(n,1) && ap_grid(m,1) <= G_grid(n+1,1) ) 
    P_a(m,n) =  1 - (ap_grid(m,1) - G_grid(n,1))/Interval_G;
elseif ( ap_grid(m,1) >= G_grid(n-1,1) && ap_grid(m,1) <= G_grid(n,1))
    P_a(m,n) = (ap_grid(m,1) - G_grid(n-1,1))/Interval_G;
elseif mod(n,n_dis) == 0
    P1 = (ap_grid(m,1) - G_grid(n-1))/Interval_G;
    P_a(m,n) =  min(1.0, max(0.0, P1));
else
    P_a(m,n) = 0;
end 
end
end



mu_init = ones(n_z*n_dis,1)./n_dis;




    







