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
[z, p_z] = mytauchen(0,rho,(sigma^2/(1-rho^2)),n_z);
% constrcut G_matrix
z = exp(z);
Zmat = repmat(z,1,n_a)'; % z matrix 500 x 7 (z1,z2,...z7 columns)
A_pZmat = meshgrid(a_prime,z)'; % meshgrid of a' and z
%% inner loop simulation grid



%% outer loop
upper_r = 1/beta - 1;
lower_r = 0;
% initial guess r
r = (upper_r+lower_r)/2;


%while bisection_dif > bisection_tol 
%% inner loop

% recover capital rental rate and wage
r_k = r + delta;
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
    if rem(iter,500) == 0
      fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dif]);
    end
    iter = iter +1;
env_old = env_new;
end



%% 
A_dis = linspace(a_min,a_max,n_a);
Amat_dis = repmat(A_dis,1,n_z)';
Zmat_dis = repmat(z,1,n_a)';
Zmat_dis = reshape(Zmat_dis,[n_a*n_z,1]);
G_grid = [Amat_dis, Zmat_dis];
ap_grid = reshape(ap_endo,[n_a*n_z,1]);



%% construct gamma matrix working in progress


%% coordinate
%((mod(m-1,n_a)+1) , floor((m/(n_a+1))+1))
m = 1:1:n_a*n_z;
COO = [(mod(m-1,n_a)+1)' , floor((m/(n_a+1))+1)'];
temp = zeros(n_a*n_z,2);

for m = 1:n_a*n_z
if m == 1
temp(m,1) = 1 - (ap_grid(m,1) - G_grid(m,1))/(G_grid(m+1,1) - G_grid(m,1));
elseif m == n_z*n_a
    if (ap_grid(m,1) >= G_grid(m-1,1) && ...
        ap_grid(m,1) <= G_grid(m,1))
temp(m,:) = (ap_grid(m,1) - G_grid(m-1,1))/(G_grid(m,1) - G_grid(m-1,1)) ;
    else
        temp(m,:) = 0;
    end
elseif (ap_grid(m ,1) >= G_grid(m,1) && ...
        ap_grid(m ,1) <= G_grid(m+1,1))
temp(m,1) = 1 - (ap_grid(m,1) - G_grid(m,1))/(G_grid(m+1,1) - G_grid(m,1)) ;
elseif (ap_grid(m,1) >= G_grid(m-1,1) && ...
        ap_grid(m,1) <= G_grid(m,1))
temp(m,2) = (ap_grid(m,1) - G_grid(m-1,1))/(G_grid(m,1) - G_grid(m-1,1)) ;
else
temp(m,:) = 0;
end
end




m_test = 1:1:3500;

COO_row = [mod(m_test -1, n_a)+1]';
COO_column= [floor((m_test/(n_a+1))+1)]';
test = sparse([COO_row, COO_row ] ,[COO_column , COO_column] ...
   ,[temp(:,1), temp(:,2)]);


init_mu = ones(n_a,1);
interp1(G_grid(1:500,1),G_grid(1:500,1), ap_grid(1:500,1),'previous')
interp1(G_grid(1:500,1),G_grid(1:500,1), ap_grid(1:500,1),'next')



