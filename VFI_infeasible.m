clear all
clc




%%initialization 
alpha = 0.36;
beta = 0.985;
delta = 0.025;
gamma = 2;
rho = 0.91;
sigma = 0.1;
N_z = 7;
N_a = 50;
a_min = 0;
a_max = 200;
N_dis = 700;
%% grid
% asset
a_prime = linspace(a_min,a_max,N_a)';
% labor skill shock
[z,p_z] = mytauchen(0,rho,sigma,N_z);
% constrcut G_matrix
z = exp(z);
Zmat = z;
Amat = a_prime;
Tmat = p_z;
model_params=struct("alpha",alpha,"beta",beta,"delta",delta,"gamma",gamma);



upper_r = 1/beta - 1;
lower_r = 0;
% initial guess r
r0 = (upper_r+lower_r)/2;
% recover capital rental rate and wage
r_k = r0 + delta;
w = (r_k / alpha)^(alpha/(alpha-1)) * (1-alpha);

%Iteration parameters
tol=10e-6;
maxits=1000;
dif=10; 
its=0; 

%Grid characteristics
%Grid
%Initialize the value funcion:
v0=zeros(N_a*N_z,1);
vp=zeros(N_a*N_z,1);
ap=zeros(N_a*N_z,1);

%I'm going to constuct the interpolant using functions from the CompEcon
%toolbox. This toolbox is frequently used and is available in Matlab and
%Julia (possibly Python as well, but I am not certain).
fun=fundef({'spli',[Amat(1);Amat(end)],N_a,1},{'spli',[Zmat(1);Zmat(end)],N_z,1});
basis_mat=funbas(fun);
nodes=gridmake(Amat,Zmat);
a_nodes=nodes(:,1);
z_nodes=nodes(:,2);
N_nodes=length(a_nodes);
p_stack = reshape(Tmat',[1,49]);
p_stack_test = repmat(p_stack,N_a,1);
p_stack_final = reshape(p_stack_test,[N_a*N_z*N_z,1]);
z_nodes_idx=zeros(N_nodes,1);
for i = 1:N_z
z_nodes_idx(z_nodes == Zmat(i),1) = i;
end

oh=ones(N_z,1);

v0_newton=zeros(N_a*N_z,1);
coeff_newton=basis_mat\v0_newton;
vp_newton=zeros(N_a*N_z,1);
ap_newton=zeros(N_a*N_z,1);

dif=10.0;
its=0;
its_newton=0;
tic;
while dif > tol && its < maxits 
    for i=1:N_nodes
        a0=a_nodes(i,1);
        z0=z_nodes(i,1);
        z_idx=z_nodes_idx(i,1);
        Tvec=Tmat(z_idx,:);
        a1=fminbnd(@(x) valfun(x,fun,coeff_newton,a0,z0,Zmat,Tvec,oh,model_params,r0,w),a_min,a_max);
        vp_newton(i,1)=-valfun(a1,fun,coeff_newton,a0,z0,Zmat,Tvec,oh,model_params,r0,w);
        ap_newton(i,1)=a1; 
    end

   ap_stack = repmat(ap_newton,N_z,1);
   z_stack = repmat(z_nodes,N_z,1);
   basis_mat_ap = funbas(fun,[ap_stack z_stack]);
  

    dif=norm(vp_newton-v0_newton,Inf)
    if its <its_newton
            coeff_newton(:)=basis_mat\vp_newton;
    else
        update_mat=model_params.beta*[p_stack_final' * basis_mat_ap];  
        coeff_newton(:)=coeff_newton-(basis_mat-update_mat)\(basis_mat*coeff_newton-vp_newton);
    end
    v0_newton(:)=vp_newton;
    its=its+1;
end
toc;

