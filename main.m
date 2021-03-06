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
lower_r = 1e-6;
bisection_dif = Inf;
bisection_tol = 1e-8;

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
r = (upper_r+lower_r)/2;
bisection_iter = 0;
lss=exp((0.5*sigma^2)/(1.0-rho^2));

%% Loop start (Bisection)
while bisection_dif > bisection_tol 

% recover capital rental rate and wage
r_k = r + delta;
%k_d =((r_k)/alpha)^(1/(alpha-1));
k_d =((r_k*lss^(alpha-1.0))/alpha)^(1/(alpha-1));
w=(1-alpha)*k_d^alpha*lss^(-alpha);

%w = (r_k / alpha)^(alpha/(alpha-1)) * (1-alpha);

% EGM initialization
M_f = (1+r)*A_pZmat + w * Zmat; % total available resource
ap_endo = zeros(n_a,n_z); % initialize ap_endo grid
c_0 = M_f - A_pZmat; % initial guess of consumption 
env_old = (1+r)*c_0.^(-gamma); % initial guess of envelope condition

% convergence criterion
criter = 1e-5;
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


idx_x = zeros(n_dis, n_z);
 for i = 1:n_dis
     for j = 1:n_z 
          idx_x(i, j) = bisection_search(1,Ag_dis,AP_dis(i,j),n_dis);
     end
 end
 for j = 1:n_z-1
        idx_x(:,j+1) = idx_x(:,j+1) + j*n_dis;
 end

 idx_x = reshape(idx_x,[n_dis*n_z,1]);

    P_a = zeros(n_dis*n_z, n_dis*n_z);
    for i = 1:n_dis*n_z
        if rem(idx_x(i,:), n_dis) == 0
            P_a(idx_x(i,1),i) = 1; 
        else
            P_a(idx_x(i,1),i) = 1 - (ap_grid(i,1) - Amat_dis(idx_x(i,1)))/...
                (Amat_dis(idx_x(i,1)+1) - Amat_dis(idx_x(i,1)));
            P_a(idx_x(i,1)+1,i) = (ap_grid(i,1) - Amat_dis(idx_x(i,1)))/...
                (Amat_dis(idx_x(i,1)+1) - Amat_dis(idx_x(i,1)));
        end
    end
    
P_a = sparse(P_a);

tol_mu = 10e-6;
dif_mu =Inf;
mu_init = ones(n_dis*n_z,1)/(n_dis);
mu_iter = 0;

%% Need to fix
while dif_mu > tol_mu

  mu_ap = P_a*mu_init;

  mu_ap = reshape(mu_ap,[n_dis, n_z]);

  mu_ap_zp = mu_ap*p_z;
  
  mu_update = reshape(mu_ap_zp, [n_dis*n_z, 1]);
  
  dif_mu=norm(mu_init-mu_update,Inf);
  if rem(mu_iter,500) ==0
      fprintf('mu convergence, iteration: %3i, Norm: %2.6f \n',[mu_iter,dif_mu]);
  end
        mu_init(:,:)=mu_update;
        mu_iter=mu_iter+1;
end


k_s = (Amat_dis'*mu_init)./n_z;

cap_dif = norm(k_s-k_d);

if (k_s > k_d)
    upper_r = r;
elseif (k_s < k_d)
    lower_r = r;
elseif abs(cap_dif) < bisection_tol
    break;
end
 
r = (upper_r + lower_r) / 2;
bisection_dif = upper_r - lower_r;
fprintf('Bisection convergence, iteration: %3i, Norm: %2.6f \n',[bisection_iter,bisection_dif]);
bisection_iter = bisection_iter +1;
end


%Plotting the asset policy & consumption policy functions
figure
plot(A_dis, AP_dis(:,1), 'LineWidth',2)
hold on
plot(A_dis, AP_dis(:,2),'LineWidth',2)    
hold on
plot(A_dis, AP_dis(:,3),'LineWidth',2)
hold on
plot(A_dis, AP_dis(:,4),'LineWidth',2)
hold on
plot(A_dis, AP_dis(:,5),'LineWidth',2)
hold on
plot(A_dis, AP_dis(:,6),'LineWidth',2)
hold on
plot(A_dis, AP_dis(:,7),'LineWidth',2)
hold on
plot(A_dis, A_dis,'LineWidth',2,'Color','black')
hold off
legend ('Ap_1', 'Ap_2', 'Ap_3', 'Ap_4', 'Ap_5', 'Ap_6', 'Ap_7', '45 degree line','Location','northwest')


%Get consumption policy function

c_dis = zeros(n_dis,n_z);
for i = 1:n_z
c_dis(:,i) = interp1(A_pZmat(:,i),c_endo(:,i),A_dis');
end


figure
plot(A_dis, c_dis(:,1),'LineWidth',2)
hold on
plot(A_dis, c_dis(:,2),'LineWidth',2)    
hold on
plot(A_dis, c_dis(:,3),'LineWidth',2)
hold on
plot(A_dis, c_dis(:,4),'LineWidth',2)
hold on
plot(A_dis, c_dis(:,5),'LineWidth',2)
hold on
plot(A_dis, c_dis(:,6),'LineWidth',2)
hold on
plot(A_dis, c_dis(:,7),'LineWidth',2)
hold off
legend ('C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7','Location','northwest')

%Plotting the stationary distribution of assets
asset_dist = Amat_dis.*mu_init;
histogram(asset_dist, 100)


histrange=40;
histgraph=zeros(histrange,1);
for histg=1:n_dis*n_z
    for histgr=1:histrange
        if a_max/histrange*(histgr-1)<=ap_grid(histg)&& ap_grid(histg)<a_max/histrange*(histgr)
            histgraph(histgr,1)=histgraph(histgr,1)+mu_update(histg,1);
        end
    end
end
histgraph=histgraph/7;

xaxis=linspace(a_min,a_max,histrange)
bar(xaxis,histgraph,1),xlabel('asset'),ylabel('probability')
title('asset distribution')