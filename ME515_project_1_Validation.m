clc; clear; close all;

format short e

b = input('Enter beam height b (m): ');
t = input('Enter beam thickness t (m): ');
nx = input('No. of elements in x- direction nx: ');
ny = input('No. of elements in y- direction ny: ');

ratio = 50;
a = ratio * b;

E = 200e9;
nu = 0.3;

q = input('Enter distributed load q (N/m): ');

x_A = a/2;
y_A = 0;

D = (E/(1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

tol = 0.001;
max_iter = 20;

A = t * b;
I = t * b^3 / 12;
G = E / (2*(1+nu));
kappa = 5/6;

x = x_A;

v_EB = (q * x^2 / (24 * E * I)) * (6*a^2 - 4*a*x + x^2);
v_shear = (q / (kappa * G * A)) * (a*x - x^2/2);
v_Timo = v_EB + v_shear;

U_EB = q^2 * a^5 / (40 * E * I);
U_shear = q^2 * a^3 / (6 * kappa * G * A);
U_Timo = U_EB + U_shear;

results = [];
uy_vals = [];
U_vals = [];

err_EB_list = [];
err_Timo_list = [];
errU_EB_list = [];
errU_Timo_list = [];

prev_uy = NaN;
prev_U = NaN;

fprintf(' iter   ny    nx    elems      |uy(A)|            U_fem           rel_uy         rel_U\n');
fprintf('-------------------------------------------------------------------------------------------\n');

for it = 1:max_iter
    
    [uy_A, U_fem, t_elem, node_coords, elements, d] = ...
        solve_fem(nx, ny, a, b, t, D, q, x_A, y_A);
    
    uy_vals(end+1) = abs(uy_A);
    U_vals(end+1) = U_fem;
    results(end+1) = t_elem;
    
    % Relative convergence
    if it == 1
        rel_uy = NaN;
        rel_U = NaN;
    else
        rel_uy = abs((uy_vals(end)-prev_uy)/uy_vals(end));
        rel_U = abs((U_vals(end)-prev_U)/U_vals(end));
    end
    
    fprintf('%4d  %4d  %4d  %7d  %.6e  %.6e  %.6e  %.6e\n', ...
        it, ny, nx, t_elem, uy_vals(end), U_vals(end), rel_uy, rel_U);
    
    % -------- ERROR STORAGE --------
    err_EB_list(end+1) = abs((uy_vals(end) - v_EB)/v_EB)*100;
    err_Timo_list(end+1) = abs((uy_vals(end) - v_Timo)/v_Timo)*100;
    
    errU_EB_list(end+1) = abs((U_vals(end) - U_EB)/U_EB)*100;
    errU_Timo_list(end+1) = abs((U_vals(end) - U_Timo)/U_Timo)*100;
    
    if it > 1 && rel_uy < tol && rel_U < tol
        break;
    end
    
    prev_uy = uy_vals(end);
    prev_U = U_vals(end);
    
    nx = nx + 50;
    ny = ny + 1;
end

fprintf('\nAnalytical Results:\n');
fprintf('Euler-Bernoulli uy(A) = %.6e\n', v_EB);
fprintf('Timoshenko uy(A)     = %.6e\n', v_Timo);
fprintf('FEM uy(A)            = %.6e\n', uy_vals(end));

fprintf('\nStrain Energy Comparison:\n');
fprintf('Euler-Bernoulli U = %.6e\n', U_EB);
fprintf('Timoshenko U      = %.6e\n', U_Timo);
fprintf('FEM U             = %.6e\n', U_vals(end));

figure;

subplot(1,2,1)
plot(results, uy_vals, 'o-','LineWidth',2); hold on;
yline(abs(v_EB),'--g','EB');
yline(abs(v_Timo),'--r','Timoshenko');
xlabel('Elements'); ylabel('|uy(A)|')
title('Displacement Convergence');
legend('FEM','Euler-Bernoulli','Timoshenko')
grid on

subplot(1,2,2)
plot(results, U_vals, 'o-','LineWidth',2); hold on;
yline(U_EB,'--g','EB');
yline(U_Timo,'--r','Timoshenko');
xlabel('Elements'); ylabel('Strain Energy')
title('Energy Convergence');
legend('FEM','Euler-Bernoulli','Timoshenko')
grid on

figure;

subplot(1,2,1)
plot(results, err_EB_list, 'o-','LineWidth',2); hold on;
plot(results, err_Timo_list, 's-','LineWidth',2);
xlabel('Elements');
ylabel('Error in Displacement (%)');
title('Displacement Error Convergence');
legend('Error vs EB','Error vs Timoshenko');
grid on;

subplot(1,2,2)
plot(results, errU_EB_list, 'o-','LineWidth',2); hold on;
plot(results, errU_Timo_list, 's-','LineWidth',2);
xlabel('Elements');
ylabel('Error in Energy (%)');
title('Energy Error Convergence');
legend('Error vs EB','Error vs Timoshenko');
grid on;

function [uy_A, U, t_elem, node_coords, elements, d] = ...
         solve_fem(nx, ny, a, b, t, D, q, xA, yA)

nnode = (nx+1)*(ny+1);
nelem = nx*ny;
ndof = 2*nnode;

dx = a/nx;
dy = b/ny;

node_coords = zeros(nnode,2);
id = 1;
for j = 0:ny
    for i = 0:nx
        node_coords(id,:) = [i*dx, j*dy];
        id = id + 1;
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
    
    nodes_elem = elements(el,:);
    coords = node_coords(nodes_elem,:);
    
    Ke = zeros(8,8);
    Fe = zeros(8,1);
    
    % ---- stiffness ----
    for i = 1:2
        for j = 1:2
            xi  = gp(i);
            eta = gp(j);

            dN_dxi  = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];

            J = [dN_dxi*coords(:,1), dN_deta*coords(:,1);
                 dN_dxi*coords(:,2), dN_deta*coords(:,2)];

            detJ = det(J);
            invJ = inv(J);

            dN_dx = invJ(1,1)*dN_dxi + invJ(1,2)*dN_deta;
            dN_dy = invJ(2,1)*dN_dxi + invJ(2,2)*dN_deta;

            B = zeros(3,8);
            for k = 1:4
                B(:,2*k-1:2*k) = [dN_dx(k) 0;
                                  0 dN_dy(k);
                                  dN_dy(k) dN_dx(k)];
            end

            Ke = Ke + B' * D * B * detJ * t;
        end
    end

    % ---- LOAD (CORRECT POSITION) ----
    if abs(coords(3,2)-b)<1e-6 && abs(coords(4,2)-b)<1e-6
        
        for xi = gp
            eta = 1;

            N = 0.25 * [(1-xi)*(1-eta)
                        (1+xi)*(1-eta)
                        (1+xi)*(1+eta)
                        (1-xi)*(1+eta)];

            dNdxi = 0.25 * [-(1-eta) (1-eta) (1+eta) -(1+eta)];

            dx_dxi = dNdxi * coords(:,1);
            dy_dxi = dNdxi * coords(:,2);

            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);

            for k = 1:4
                Fe(2*k) = Fe(2*k) - q * N(k) * Jedge;
            end
        end
    end

    % ---- ASSEMBLY (CRITICAL) ----
    dof = reshape([2*nodes_elem-1; 2*nodes_elem],1,[]);
    
    K(dof,dof) = K(dof,dof) + Ke;
    F(dof) = F(dof) + Fe;  
end


fixed = find(abs(node_coords(:,1)) < 1e-8);
fixed_dofs = reshape([2*fixed-1 2*fixed]',1,[]);

free = setdiff(1:ndof, fixed_dofs);

d = zeros(ndof,1);
d(free) = K(free,free) \ F(free);
U = 0.5 * d' * K * d;

[~, idx] = min(abs(node_coords(:,1)-xA));
uy_A = d(2*idx);

t_elem = nelem;

end