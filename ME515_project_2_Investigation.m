clc; 
clear;
close all;
format short e

b = input('Enter beam height b (m): ');
t = input('Enter thickness t (m): ');

E = 200e9;
nu = 0.3;

q = input('Enter distributed load q (N/m): ');

aspect_ratios = [2 3 5 8 10 15 20 25 30 40 50];

nx = 300;
ny = 6;

uy_FEM = zeros(size(aspect_ratios));
uy_EB  = zeros(size(aspect_ratios));
uy_Timo= zeros(size(aspect_ratios));

for i = 1:length(aspect_ratios)
    
    ratio = aspect_ratios(i);
    a = ratio * b;
    
    % FEM solution
    [uy_A, U_fem, ~, ~, ~] = solve_fem(nx, ny, a, b, t, E, nu, q);
    
    uy_FEM(i) = uy_A;
    
  
    A = t*b;
    I = t*b^3/12;
    G = E/(2*(1+nu));
    kappa = 5/6;
    
    x = a/2;
    
    % Euler-Bernoulli
    uy_EB(i) = (q*x^2/(24*E*I))*(6*a^2 - 4*a*x + x^2);
    
    % Timoshenko
    uy_Timo(i) = uy_EB(i) + (q/(kappa*G*A))*(a*x - x^2/2);
end


% Displacement comparison
figure;
plot(aspect_ratios, abs(uy_FEM),'o-','LineWidth',2); hold on;
plot(aspect_ratios, abs(uy_EB),'--g','LineWidth',2);
plot(aspect_ratios, abs(uy_Timo),'--r','LineWidth',2);
xlabel('Aspect Ratio (L/b)');
ylabel('|uy(A)|');
title('Displacement Comparison');
legend('FEM','Euler-Bernoulli','Timoshenko');
grid on;

% Error plot
err_EB = abs((uy_FEM - uy_EB)./uy_EB)*100;
err_T  = abs((uy_FEM - uy_Timo)./uy_Timo)*100;

figure;
plot(aspect_ratios, err_EB,'o-','LineWidth',2); hold on;
plot(aspect_ratios, err_T,'s-','LineWidth',2);
xlabel('Aspect Ratio');
ylabel('Error (%)');
title('Error vs Aspect Ratio');
legend('Error vs EB','Error vs Timoshenko');
grid on;

U_FEM  = zeros(size(aspect_ratios));
U_EB   = zeros(size(aspect_ratios));
U_Timo = zeros(size(aspect_ratios));

for i = 1:length(aspect_ratios)
    
    ratio = aspect_ratios(i);
    L = ratio * b;
    
    % Fem solution
    [uy_A, U_fem, ~, ~, ~] = solve_fem(nx, ny, L, b, t, E, nu, q);
    
    uy_FEM(i) = uy_A;
    U_FEM(i)  = U_fem; 
    
    A = t*b;
    I = t*b^3/12;
    G = E/(2*(1+nu));
    kappa = 5/6;
    
    x = L/2;
    
   
    % Euler - Bernoulli
    U_EB(i) = q^2 * L^5 / (40 * E * I);

    % Timoshenko
    U_shear = q^2 * L^3 / (6 * kappa * G * A);
    U_Timo(i) = U_EB(i) + U_shear;
end

% Strain Energy comparison
figure;
plot(aspect_ratios, U_FEM,'o-','LineWidth',2); hold on;
plot(aspect_ratios, U_EB,'--g','LineWidth',2);
plot(aspect_ratios, U_Timo,'--r','LineWidth',2);

xlabel('Aspect Ratio (L/b)');
ylabel('Strain Energy');
title('Strain Energy Comparison');
legend('FEM','Euler-Bernoulli','Timoshenko');
grid on;

% Error plot
errU_EB = abs((U_FEM - U_EB)./U_EB)*100;
errU_T  = abs((U_FEM - U_Timo)./U_Timo)*100;

figure;
plot(aspect_ratios, errU_EB,'o-','LineWidth',2); hold on;
plot(aspect_ratios, errU_T,'s-','LineWidth',2);

xlabel('Aspect Ratio');
ylabel('Error (%)');
title('Energy Error vs Aspect Ratio');
legend('Error vs EB','Error vs Timoshenko');
grid on;


a = 10*b;

[~,~,node_coords,elements,d] = solve_fem(nx, ny, a, b, t, E, nu, q);

D = (E/(1 - nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

[Ux, Uy, strain_gp, stress_gp, gp_coords] = ...
    compute_fields(node_coords, elements, d, D);

xq = linspace(min(node_coords(:,1)), max(node_coords(:,1)), 80);
yq = linspace(min(node_coords(:,2)), max(node_coords(:,2)), 40);
[Xq,Yq] = meshgrid(xq,yq);

Uy_i = griddata(node_coords(:,1), node_coords(:,2), Uy, Xq, Yq);
sigx = griddata(gp_coords(:,1), gp_coords(:,2), stress_gp(:,1), Xq, Yq);
tauxy = griddata(gp_coords(:,1), gp_coords(:,2), stress_gp(:,3), Xq, Yq);

epsx = griddata(gp_coords(:,1), gp_coords(:,2), strain_gp(:,1), Xq, Yq);
gamxy = griddata(gp_coords(:,1), gp_coords(:,2), strain_gp(:,3), Xq, Yq);

figure;
contourf(Xq,Yq,Uy_i,20,'LineColor','none'); colorbar
title('Total Deformation (Uy)'); axis equal

figure;
contourf(Xq,Yq,sigx,20,'LineColor','none'); colorbar
title('\sigma_x'); axis equal

figure;
contourf(Xq,Yq,tauxy,20,'LineColor','none'); colorbar
title('\tau_{xy}'); axis equal

figure;
contourf(Xq,Yq,epsx,20,'LineColor','none'); colorbar
title('\epsilon_x'); axis equal

figure;
contourf(Xq,Yq,gamxy,20,'LineColor','none'); colorbar
title('\gamma_{xy}'); axis equal

A = t*b;
I = t*b^3/12;
G = E/(2*(1+nu));

x = a/10;
y_vals = linspace(-b/2, b/2, 100);

for i = 1:length(y_vals)
    
    y = y_vals(i);
    
    M = q*(a-x)^2/2;
    V = q*(a-x);
    
    sigma_x = M*y/I;
    tau_xy = (3/2)*(V/A)*(1 - 4*y^2/b^2);
    
    sigma_vm(i) = sqrt(sigma_x^2 + 3*tau_xy^2);
    
    eps_x = sigma_x/E;
    gamma_xy = tau_xy/G;
    
    eps_eq(i) = sqrt(2/3*(eps_x^2 + gamma_xy^2/2));
end

figure;
plot(sigma_vm, y_vals,'LineWidth',2);
xlabel('\sigma_{vm}'); ylabel('y');
title('Von Mises Stress (Analytical)'); grid on;

figure;
plot(eps_eq, y_vals,'LineWidth',2);
xlabel('\epsilon_{eq}'); ylabel('y');
title('Equivalent Strain (Analytical)'); grid on;


a = 10*b;          % or any aspect ratio
x = a/10;

A = t * b;

y_vals = linspace(-b/2, b/2, 100);
tau_analytical = zeros(size(y_vals));

for i = 1:length(y_vals)
    
    y = y_vals(i);
    
    V = q * (a - x);
    
    tau_analytical(i) = (3/2) * (V/A) * (1 - (4*y^2)/(b^2));
end

figure;
plot(tau_analytical, y_vals,'LineWidth',2);
xlabel('\tau_{xy}');
ylabel('y');
title('Shear Stress Distribution (Analytical)');
grid on;


tau_fem = zeros(size(y_vals));

for i = 1:length(y_vals)
    tau_fem(i) = compute_shear(x, y_vals(i), ...
                              node_coords, elements, d, D);
end

hold on;
plot(tau_fem, y_vals,'--r','LineWidth',2);
legend('Analytical','FEM');

function [uy_A, U, node_coords, elements, d] = solve_fem(nx, ny, a, b, t, E, nu, q)

D = (E/(1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

nnode = (nx+1)*(ny+1);
nelem = nx*ny;
ndof = 2*nnode;

dx = a/nx;
dy = b/ny;

node_coords = zeros(nnode,2);
c = 1;
for j = 0:ny
    for i = 0:nx
        node_coords(c,:) = [i*dx, j*dy];
        c = c + 1;
    end
end

elements = zeros(nelem,4);
e = 1;
for j = 1:ny
    for i = 1:nx
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);
        elements(e,:) = [n1 n2 n3 n4];
        e = e + 1;
    end
end

K = sparse(ndof, ndof);
F = zeros(ndof,1);

gp = [-1 1]/sqrt(3);

for el = 1:nelem
    
    nodes = elements(el,:);
    coords = node_coords(nodes,:);
    
    Ke = zeros(8,8);
    Fe = zeros(8,1);
    
    for i = 1:2
        for j = 1:2
            
            xi  = gp(i);
            eta = gp(j);
            
            dN_dxi  = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];
            
            J = [dN_dxi*coords(:,1), dN_deta*coords(:,1);
                 dN_dxi*coords(:,2), dN_deta*coords(:,2)];
             
            detJ = det(J);
            
            dN = J \ [dN_dxi; dN_deta];
            dN_dx = dN(1,:);
            dN_dy = dN(2,:);
            
            B = zeros(3,8);
            B(1,1:2:end) = dN_dx;
            B(2,2:2:end) = dN_dy;
            B(3,1:2:end) = dN_dy;
            B(3,2:2:end) = dN_dx;
            
            Ke = Ke + B' * D * B * detJ * t;
        end
    end
    
    % Load on top edge
    if abs(coords(3,2)-b)<1e-6 && abs(coords(4,2)-b)<1e-6
        for xi = gp
            eta = 1;
            N = 0.25 * [(1-xi)*(1-eta)
                        (1+xi)*(1-eta)
                        (1+xi)*(1+eta)
                        (1-xi)*(1+eta)];
            
            dx_dxi = [-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)] * coords(:,1);
            dy_dxi = [-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)] * coords(:,2);
            
            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);
            
            for k = 1:4
                Fe(2*k) = Fe(2*k) - q * N(k) * Jedge;
            end
        end
    end
    
    dof = reshape([2*nodes-1; 2*nodes],1,[]);
    K(dof,dof) = K(dof,dof) + Ke;
    F(dof) = F(dof) + Fe;
end

fixed = find(abs(node_coords(:,1)) < 1e-8);
fixed_dofs = reshape([2*fixed-1 2*fixed]',1,[]);
free = setdiff(1:ndof, fixed_dofs);

d = zeros(ndof,1);
d(free) = K(free,free)\F(free);

U = 0.5 * d' * K * d;

[~, idx] = min(abs(node_coords(:,1)-a/2));
uy_A = d(2*idx);

end

function tau = compute_shear(x, y, node_coords, elements, d, D)

tau = NaN;

for el = 1:size(elements,1)
    
    nodes = elements(el,:);
    coords = node_coords(nodes,:);
    
    % Check if point lies in element (bounding box)
    if x >= min(coords(:,1)) && x <= max(coords(:,1)) && ...
       y >= min(coords(:,2)) && y <= max(coords(:,2))
   
        % -------- LOCAL COORDINATES (APPROX MAPPING) --------
        x1 = coords(1,1); x2 = coords(2,1);
        y1 = coords(1,2); y4 = coords(4,2);
        
        xi  = 2*(x - (x1+x2)/2)/(x2-x1);
        eta = 2*(y - (y1+y4)/2)/(y4-y1);
        
        % Clamp values (avoid numerical overflow)
        xi = max(min(xi,1),-1);
        eta = max(min(eta,1),-1);
        
        % -------- SHAPE FUNCTION DERIVATIVES --------
        dNdxi  = 0.25 * [-(1-eta),(1-eta),(1+eta),-(1+eta)];
        dNdeta = 0.25 * [-(1-xi),-(1+xi),(1+xi),(1-xi)];
        
        J = [dNdxi*coords(:,1), dNdeta*coords(:,1);
             dNdxi*coords(:,2), dNdeta*coords(:,2)];
        
        dN = J \ [dNdxi; dNdeta];
        
        dN_dx = dN(1,:);
        dN_dy = dN(2,:);
        
        % -------- B MATRIX --------
        B = zeros(3,8);
        B(1,1:2:end) = dN_dx;
        B(2,2:2:end) = dN_dy;
        B(3,1:2:end) = dN_dy;
        B(3,2:2:end) = dN_dx;
        
        % -------- ELEMENT DISPLACEMENT --------
        dof = reshape([2*nodes-1;2*nodes],1,[]);
        d_elem = d(dof);
        
        stress = D * B * d_elem;
        
        tau = stress(3);   % τ_xy
        
        return
    end
end

% fallback
tau = 0;

end
function [Ux, Uy, strain_gp, stress_gp, gp_coords] = compute_fields(node_coords, elements, d, D)

nnode = size(node_coords,1);
nelem = size(elements,1);

Ux = d(1:2:end);
Uy = d(2:2:end);

gp = [-1 1]/sqrt(3);

strain_gp = [];
stress_gp = [];
gp_coords = [];

for el = 1:nelem
    
    nodes = elements(el,:);
    coords = node_coords(nodes,:);
    
    dof = reshape([2*nodes-1; 2*nodes],1,[]);
    d_elem = d(dof);
    
    for i = 1:2
        for j = 1:2
            
            xi = gp(i);
            eta = gp(j);
            
            dN_dxi  = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];
            
            J = [dN_dxi*coords(:,1), dN_deta*coords(:,1);
                 dN_dxi*coords(:,2), dN_deta*coords(:,2)];
             
            dN = J \ [dN_dxi; dN_deta];
            dN_dx = dN(1,:);
            dN_dy = dN(2,:);
            
            B = zeros(3,8);
            B(1,1:2:end) = dN_dx;
            B(2,2:2:end) = dN_dy;
            B(3,1:2:end) = dN_dy;
            B(3,2:2:end) = dN_dx;
            
            strain = B * d_elem;
            stress = D * strain;
            
            % Gauss point location (approx)
            N = 0.25 * [(1-xi)*(1-eta)
                        (1+xi)*(1-eta)
                        (1+xi)*(1+eta)
                        (1-xi)*(1+eta)];
            
            x_gp = N' * coords(:,1);
            y_gp = N' * coords(:,2);
            
            strain_gp = [strain_gp; strain'];
            stress_gp = [stress_gp; stress'];
            gp_coords = [gp_coords; x_gp y_gp];
        end
    end
end

end