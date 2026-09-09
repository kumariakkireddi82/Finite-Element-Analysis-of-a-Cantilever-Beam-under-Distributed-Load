clc; clear; close all;

format short e

b = input('Enter beam height b (m): ');
t = input('Enter thickness t (m): ');
rho = input('Enter density rho (kg/m^3): ');
nx = input('Initial elements in x-direction nx: ');
ny = input('Elements in y-direction ny (from static): ');

ratio = 50;
a = ratio * b;

E = 200e9;
nu = 0.3;

D = (E/(1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

freq5_list = [];
elements_list = [];

max_iter = 10;

for it = 1:max_iter
    
    [K, M, node_coords, ~] = assemble_system(nx, ny, a, b, t, D, rho);
    
    % Apply BC (cantilever)
    fixed = find(abs(node_coords(:,1)) < 1e-8);
    fixed_dofs = reshape([2*fixed-1 2*fixed]',1,[]);
    free = setdiff(1:size(K,1), fixed_dofs);
    
    Kf = K(free,free);
    Mf = M(free,free);
    
    % Solve eigenvalue problem
    num_modes = 5;
    
    [V, D_eig] = eigs(Kf, Mf, num_modes, 'smallestabs');
    omega = sqrt(diag(D_eig));
    freq = omega / (2*pi);
    [freq, idx] = sort(freq);
    V = V(:,idx);

    % Store 5th mode
    freq5_list(end+1) = freq(5);
    elements_list(end+1) = nx*ny;
    
    fprintf('Iter %d → Elements: %d → 5th freq = %.4f Hz\n', ...
        it, nx*ny, freq(5));
    
    nx = nx + 20;   % refine mesh
end

figure;
plot(elements_list, freq5_list,'o-','LineWidth',2)
xlabel('Number of Elements')
ylabel('5th Natural Frequency (Hz)')
title('Convergence of 5th Mode')
grid on

nx_final = nx;
[K, M, node_coords, elements] = assemble_system(nx_final, ny, a, b, t, D, rho);

fixed = find(abs(node_coords(:,1)) < 1e-8);
fixed_dofs = reshape([2*fixed-1 2*fixed]',1,[]);
free = setdiff(1:size(K,1), fixed_dofs);

Kf = K(free,free);
Mf = M(free,free);

num_modes = 5;

[V, D_eig] = eigs(Kf, Mf, num_modes, 'smallestabs');

omega = sqrt(diag(D_eig));
freq = omega/(2*pi);

[freq, idx] = sort(freq);
V = V(:,idx);

fprintf('\nFirst 5 Natural Frequencies (Hz):\n');
disp(freq)

for mode = 1:5
    
    d = zeros(size(K,1),1);
    d(free) = V(:,mode);
    
    Ux = d(1:2:end);
    Uy = d(2:2:end);
    
    scale = 0.1*max(a,b)/max(abs(Uy));
    
    deformed = node_coords + scale*[Ux Uy];
    
    figure;
    patch('Faces',elements,'Vertices',node_coords,...
          'EdgeColor','k','FaceColor','none');
    hold on;
    
    patch('Faces',elements,'Vertices',deformed,...
          'EdgeColor','r','FaceColor','none');
    
    title(['Mode Shape ', num2str(mode), ...
          '  (f = ', num2str(freq(mode)), ' Hz)'])
    axis equal; grid on;
end


for mode = 1:5
    
    d = zeros(size(K,1),1);
    d(free) = V(:,mode);
    
    Ux = d(1:2:end);
    Uy = d(2:2:end);
    
    scale = 0.3 * max(a,b) / max(abs(Uy));
    
    figure;
    
    for t = linspace(0, 2*pi, 40)
        
        factor = sin(t);
        
        Z = scale * factor * Uy(:);   % FIXED
        
        clf
        
        patch('Faces',elements,...
              'Vertices',[node_coords(:,1), node_coords(:,2), Z],...
              'FaceVertexCData',Z,...   % FIXED
              'FaceColor','interp',...
              'EdgeColor','k');
        
        title(['3D Mode ', num2str(mode), ...
              '  f = ', num2str(freq(mode),'%.2f'), ' Hz'])
        
        xlabel('X'); ylabel('Y'); zlabel('Displacement')
        
        % FIXED axis
        zmax = max(abs(Z));
        if zmax < 1e-12
            zmax = 1e-6;
        end
        
        axis([0 a -b b -zmax zmax])
        
        view(30,30)
        grid on
        
        camlight
        lighting gouraud
        
        drawnow
        pause(0.05)
    end
end



function [K, M, node_coords, elements] = assemble_system(nx, ny, a, b, t, D, rho)

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
M = sparse(ndof, ndof);

gp = [-1 1]/sqrt(3);

for el = 1:nelem
    
    nodes = elements(el,:);
    coords = node_coords(nodes,:);
    
    Ke = zeros(8);
    Me = zeros(8);
    
    for i = 1:2
        for j = 1:2
            
            xi = gp(i);
            eta = gp(j);
            
            N = shape_func(xi, eta);
            [dNdxi, dNdeta] = shape_deriv(xi, eta);
            
            J = [dNdxi*coords(:,1), dNdeta*coords(:,1);
                 dNdxi*coords(:,2), dNdeta*coords(:,2)];
             
            detJ = det(J);
            
            dN = J \ [dNdxi; dNdeta];
            dN_dx = dN(1,:);
            dN_dy = dN(2,:);
            
            B = zeros(3,8);
            B(1,1:2:end) = dN_dx;
            B(2,2:2:end) = dN_dy;
            B(3,1:2:end) = dN_dy;
            B(3,2:2:end) = dN_dx;
            
            Ke = Ke + B'*D*B * detJ * t;
            
            % Mass matrix
            Nmat = zeros(2,8);
            for k = 1:4
                Nmat(:,2*k-1:2*k) = [N(k) 0; 0 N(k)];
            end
            
            Me = Me + rho * t * (Nmat'*Nmat) * detJ;
        end
    end
    
    dof = reshape([2*nodes-1; 2*nodes],1,[]);
    
    K(dof,dof) = K(dof,dof) + Ke;
    M(dof,dof) = M(dof,dof) + Me;
end

end

function N = shape_func(xi, eta)
N = 0.25 * [(1-xi)*(1-eta)
            (1+xi)*(1-eta)
            (1+xi)*(1+eta)
            (1-xi)*(1+eta)];
end

function [dNdxi, dNdeta] = shape_deriv(xi, eta)
dNdxi  = 0.25 * [-(1-eta) (1-eta) (1+eta) -(1+eta)];
dNdeta = 0.25 * [-(1-xi) -(1+xi) (1+xi) (1-xi)];
end