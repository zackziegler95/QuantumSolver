% Constants for the simulation
hbar = 1;
m = 1;
w = 5;

L = 10; % length of the space
d = 0.25; % spatial step
dt = 0.3; % time step
N = L/d-1;

% Crank-Nicolson methods
% A*psi(n+1) = B*psi(n)
% psi(n+1) = A^-1*B*psi(n)
% Because it's in 2D, each matrix is now N^2 x N^2... they get pretty big
r = -hbar/(2*m*sqrt(-1))*dt/(d^2);
s = dt/(sqrt(-1)*hbar);

VV_x = zeros(N^2, 1);
VV_y = zeros(N^2, 1);

% Potential matrix
VV = sparse(N^2, N^2);

for x_i = 1:N
  x = (x_i-1)*d;
  for y_i = 1:N
    y = (y_i-1)*d;
    ind_x = y_i+(x_i-1)*N;
    ind_y = x_i+(y_i-1)*N;
    if x < 0.4*L || x > 0.6*L
      #VV(ind, ind) = 0;
    else
      #VV(ind, ind) = c;
    end
  end
end

A = diag(ones(N, 1)*2*r, 0) + diag(ones(N-1, 1)*(-r/2), -1)+ diag(ones(N-1, 1)*(-r/2), 1);

% Left hand side
AA = sparse(N^2, N^2);
for n = 1:N
  AA((n-1)*N+1:n*N, (n-1)*N+1:n*N) = sparse(A);
  if n != 1
    AA((n-1)*N+1:n*N, (n-2)*N+1:(n-1)*N) = -r/2*speye(N);
  end
  if n != N
    AA((n-1)*N+1:n*N, n*N+1:(n+1)*N) = -r/2*speye(N);
  end
end
MM = inv(speye(N^2)+AA-s/2*VV)*(speye(N^2)-AA+s/2*VV);


next_psi = zeros(N*N,1);
for x_i = 1:N
  x = (x_i)*d;
  for y_i = 1:N
    y = (y_i)*d;
    ind = y_i + (x_i-1)*N;
    c0 = 0.1;
    width = (c0*L);
    kx0 = 2.1;
    #ky0 = 1;
    eex = sqrt(2/width)*(x-0.3*L);
    #eey = sqrt(2/width)*(y-L/2);
    %next_psi(ind) = 0.5*e^(-(eex^2)/2+sqrt(-1)*kx0*eex)*sin(pi/L*y);
    next_psi(ind) = (sin(pi/L*x)+sin(2*pi/L*x))*sin(pi/L*y);
  end
end

figure;
xyval = d:d:L-d;
zplot = surf(xyval, xyval, reshape(abs(next_psi).^2, [N, N]), 'EdgeColor', 'None', 'facecolor', 'interp');
axis([d, L-d, d, L-d, 0, max(abs(next_psi).^2)])
xlabel('x')
ylabel('y')
view(10, 30)

psi_star = zeros(N^2, 1);
for t_i = 1:200
  t = t_i*dt;
  next_psi = MM*next_psi;
  
  set(zplot, 'ZData', reshape(abs(next_psi).^2, [N, N]));
  set(zplot, 'CData', reshape(abs(next_psi).^2, [N, N]));
  drawnow
end
