clear; clc; close all

% Parameters
nx = 64;
ny = nx;
Lx = 2*pi;
Ly = 1.0;

dx  = Lx / nx;
dy  = Ly / (ny - 1);
dyi = 1/dy;

x = (0:nx-1)*dx;
y = linspace(0,Ly,ny);

n = 1;
m = 2;

% Exact solution and RHS
A_coef = -1 / (n^2 + (2*pi*m/Ly)^2);

pext = zeros(nx,ny);
rhsp = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        pext(i,j) = A_coef * sin(n*x(i)) * cos(2*pi*m*y(j)/Ly);
        rhsp(i,j) =        sin(n*x(i)) * cos(2*pi*m*y(j)/Ly);
    end
end

% FFT in x
rhspc = fft(rhsp,[],1);
rhspc = rhspc(1:nx/2+1,:);

kx  = (0:nx/2) * (2*pi/Lx);
kx2 = kx.^2;

% Preallocate
pc = zeros(nx/2+1, ny);

% Loop over Fourier modes
for i = 1:nx/2+1

    % --- Allocate tridiagonal system including ghost cells ---
    N = ny + 2;      % j = 0 ... ny+1
    a = zeros(N,1);
    b = zeros(N,1);
    c = zeros(N,1);
    d = zeros(N,1);

    % --- Interior rows: j = 1 ... ny ---
    for j = 1:ny
        a(j+1) =  dyi^2;
        b(j+1) = -2*dyi^2 - kx2(i);
        c(j+1) =  dyi^2;
        d(j+1) = rhspc(i,j);
    end

    % --- Neumann BC at bottom: p(1) - p(0) = 0 ---
    % row j=0 in ghost system
    a(1) = 0;
    b(1) = -1;
    c(1) = +1;
    d(1) = 0;

    % --- Neumann BC at top: p(ny+1) - p(ny) = 0 ---
    % row j = ny+1
    a(ny+2) = +1;
    b(ny+2) = -1;
    c(ny+2) =  0;
    d(ny+2) =  0;

    % --- Pressure fix for kx=0 mode ---
    if i == 1
        % Enforce p(1)=0  → row j=1 in extended system = interior j=1
        b(2) = 1;
        c(2) = 0;
        d(2) = 0;
    end

    % === TDMA ===
    % Forward sweep
    for j = 2:N
        factor = a(j)/b(j-1);
        b(j) = b(j) - factor*c(j-1);
        d(j) = d(j) - factor*d(j-1);
    end

    % Back substitution
    sol = zeros(N,1);
    sol(N) = d(N)/b(N);
    for j = N-1:-1:1
        sol(j) = (d(j) - c(j)*sol(j+1)) / b(j);
    end

    % Extract interior j = 1:ny from sol(j+1)
    pc(i,:) = sol(2:ny+1).';
end

% Reconstruct full spectrum
p_hat_full = zeros(nx,ny);
p_hat_full(1:nx/2+1,:) = pc;
for i = 2:nx/2
    p_hat_full(nx-i+2,:) = conj(pc(i,:));
end

% Inverse FFT in x
p = real(ifft(p_hat_full,[],1));

% L2 error
err = p - pext;
l2norm = sqrt(mean(err(:).^2));
disp(['L2 error = ', num2str(l2norm)]);

% Plot results from Matlab
figure(3)
subplot(2,4,1)
contourf(x,y,rhsp', 30, 'EdgeColor', 'none');
hold on
title('RHS - Matlab'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

subplot(2,4,2)
contourf(x, y, p', 30, 'EdgeColor', 'none');
hold on
title('Numerical - Matlab'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

subplot(2,4,3)
contourf(x, y, pext', 30, 'EdgeColor', 'none');
title('Analytical'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

subplot(2,4,4)
contourf(x, y, (p-pext)', 30, 'EdgeColor', 'none');
title('Error - Matlab'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off