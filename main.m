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
bisection_tol = 1e-7;

while bisection_dif > bisection_tol
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
A_dis = linspace(a_min,a_max,n_dis);
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


%guess initial mu
mu_old=ones(n_dis,n_z)/(n_dis*n_z)*p_z;

criter_mu=1e-5;
dif_mu=Inf;
iter_mu=0;

%check where belongs
mu_idxlow=zeros(n_dis,n_z);
for z=1:n_z
for a=1:n_dis
    for b=1:n_a
    if AP_dis(a,z)>=Ag_dis(b,z)
        mu_idxlow(a,z)=mu_idxlow(a,z)+1;
    end
    end
end
end
mu_idxhigh=mu_idxlow+1;

mu_new=zeros(n_dis,n_z);

while dif_mu>criter_mu && iter_mu < 50
%Pz    
for mu_newz=1:n_z
    mu_newa=1;
    while mu_newa>min(mu_idxlow(:,mu_newz))
        mu_new(mu_newa,n_z)=0;
        mu_newa=mu_newa+1;
    end
    find_alow=find(mu_idxlow(:,mu_newz)==mu_newa);
    mu_new(mu_newa,mu_newz)=mu_old(min(find_alow):1:max(find_alow),mu_newz)'*(1-(AP_dis(min(find_alow):1:max(find_alow),mu_newz)-A_dis(mu_newa))/(A_dis(mu_newa+1)-A_dis(mu_newa)));
    mu_newa=mu_newa+1;
    while mu_newa<=max(mu_idxlow(:,mu_newz))
        find_ahigh=find(mu_idxlow(:,mu_newz)==mu_newa-1);
        find_alow=find(mu_idxlow(:,mu_newz)==mu_newa);
        mu_new(mu_newa,mu_newz)=mu_old(min(find_alow):1:max(find_alow),mu_newz)'*(1-(AP_dis(min(find_alow):1:max(find_alow),mu_newz)-A_dis(mu_newa))/(A_dis(mu_newa+1)-A_dis(mu_newa)))...
            +mu_old(min(find_ahigh)+1:1:max(find_ahigh)+1,mu_newz)'*((AP_dis(min(find_ahigh):1:max(find_ahigh),mu_newz)-A_dis(mu_newa-1))/(A_dis(mu_newa)-A_dis(mu_newa-1)));
    mu_newa=mu_newa+1;
    end
    find_ahigh=find(mu_idxlow(:,mu_newz)==mu_newa-1);
    mu_new(mu_newa,mu_newz)=mu_old(min(find_ahigh):1:max(find_ahigh),mu_newz)'*((AP_dis(min(find_ahigh):1:max(find_ahigh),mu_newz)-A_dis(mu_newa-1))/(A_dis(mu_newa)-A_dis(mu_newa-1)));
    mu_newa=mu_newa+1;
    while mu_newa<=n_dis
        mu_new(mu_newa,mu_newz)=0;
        mu_newa=mu_newa+1;
    end
end

%apply Z transition
mu_new=(mu_new)*p_z;

dif_mu=norm(mu_old-mu_new);
if dif_mu>criter_mu
mu_old=(mu_old+mu_new)/2;
end
iter_mu=iter_mu+1;
end

test_mu = reshape(mu_new,[n_dis*n_z,1]);

k_s = test_mu'*G_grid(:,1);
k_d = (r_k/(alpha))^(1/(alpha-1));
bisection_dif = norm(k_s - k_d)


if bisection_dif > 0   
    upper_r = r;
elseif bisection_dif < 0
    lower_r = r;
else
    return
end



end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%DONT LOOK AT THIS BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% construct gamma matrix working in progress

P_a = zeros(n_z*n_dis, n_z*n_dis);
Interval_G = G_grid(2,1) - G_grid(1,1);

for n = 1: n_z*n_dis
for m = 1: n_z*n_dis
if n == 1
    P_a(m,n) = 0;
elseif ( ap_grid(m,1) >= G_grid(n,1) && ap_grid(m,1) <= G_grid(n+1,1) ) 
    P_a(m,n) =  (1 - (ap_grid(m,1) - G_grid(n,1))/Interval_G);
elseif ( ap_grid(m,1) >= G_grid(n-1,1) && ap_grid(m,1) <= G_grid(n,1) ) 
    P_a(m,n) = ((ap_grid(m,1) - G_grid(n-1,1))/Interval_G);
elseif n == n_z*n_dis
    P_a(m,n) = 0;
else
    P_a(m,n) = 0;
end 
end
end



%% p_Z
tol_mu = 10e-6;
dif_mu = Inf;
mu_init = ones(n_dis,n_z)*(1/n_dis);
mu_init = reshape(mu_init,[1,n_dis*n_z]);

while dif_mu > tol_mu
mu_ap = mu_init * P_a;
mu_ap = reshape(mu_ap, [n_dis,n_z]);
mu_ap_z = mu_ap*p_z;
mu_update = reshape(mu_ap_z, [n_dis*n_z, 1]);
dif_mu=norm(mu_init-mu_update,Inf);
    mu_init(:,:)=mu_update;
end



A_vec = G_grid(1:end,1)';
A_vec = repmat(A_vec,length(ap_grid),1);
AP_vec = repmat(ap_grid,1,length(ap_grid));
test = 1 - ((A_vec - AP_vec) / Interval_G);
test2 = (A_vec - AP_vec) / Interval_G;



idx1_1 = G_grid(1:n_dis) ==  interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(1:n_dis),'previous');
idx1_2 = G_grid(1:n_dis) ==  interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(n_dis+1:2*n_dis),'previous');
idx1_3 = G_grid(1:n_dis) ==  interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(2*n_dis+1: 3*n_dis),'previous');
idx1_4 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(3*n_dis+1: 4*n_dis),'previous');
idx1_5 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(4*n_dis+1: 5*n_dis),'previous');
idx1_6 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(5*n_dis+1: 6*n_dis),'previous');
idx1_7 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(6*n_dis+1: 7*n_dis),'previous');
p_a = vertcat(idx1_1, idx1_2, idx1_3, idx1_4, idx1_5, idx1_6, idx1_7);
p_a = repmat(p_a,1,n_z);

test1 = test.*p_a;
test2 = test2 .*p_a;
mu = test1 + test2 ;


idx1_1 = G_grid(1:n_dis) ==  interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(1:n_dis),'nearest');
idx1_2 = G_grid(1:n_dis) ==  interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(n_dis+1:2*n_dis),'nearest');
idx1_3 = G_grid(1:n_dis) ==  interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(2*n_dis+1: 3*n_dis),'nearest');
idx1_4 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(3*n_dis+1: 4*n_dis),'nearest');
idx1_5 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(4*n_dis+1: 5*n_dis),'nearest');
idx1_6 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(5*n_dis+1: 6*n_dis),'nearest');
idx1_7 = G_grid(1:n_dis) == interp1(G_grid(1:n_dis,1), G_grid(1:n_dis,1), ap_grid(6*n_dis+1: 7*n_dis),'nearest');
p_b = vertcat(idx1_1, idx1_2, idx1_3, idx1_4, idx1_5, idx1_6, idx1_7);




