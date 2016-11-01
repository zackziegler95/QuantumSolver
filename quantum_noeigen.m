% Constants for the simulation
hbar = 5;
m = 1;
w = 0.8;

L = 10; % length of the space
d = 0.01; % spatial step
dt = 0.04; % time step
N = L/d-1;

% Crank-Nicolson methods
% A*psi(n+1) = B*psi(n)
% psi(n+1) = A^-1*B*psi(n)
r = -hbar/(2*m*sqrt(-1))*dt/(d^2);
s = dt/(sqrt(-1)*hbar);


V = zeros(N, 1);
c = 235000; % Pretty arbitrary

% Set the potential here, some examples are included
% The walls implicity form an infiite square well
for x = 1:N
  %V(x) = 1/2*m*w^2*(x*d-L/2)^2; % harmonic oscillator
  %V(x) = w*x; % gravitational potential
  %if x > 0.5*N && x < 0.6*N % Finite wall
  %  V(x) = c;
  %end
end

% The Schrödinger equation
A = diag(ones(N, 1)*2*(r+1), 0) + diag(ones(N-1, 1)*(-r), -1)+ diag(ones(N-1, 1)*(-r), 1) - s*diag(V);
B = diag(ones(N, 1)*2*(1-r), 0) + diag(ones(N-1, 1)*(r), -1)+ diag(ones(N-1, 1)*(r), 1) + s*diag(V);
M = inv(A)*B; % This will be applied at each time step

a = ones(N, 1)*2*(r+1) - s*V;
b = ones(N, 1)*(-r);

% Sets the initial wave function
next_psi = zeros(N, 1);
for x_i = 1:N
  x = (x_i)*d;
  next_psi(x_i) = sin(pi/L*x)+sin(2*pi/L*x); % First two eigenfunctions of the infinite square well
  
  % Moving wave packet
  %{
  kx = 25;
  width = (0.005*L);
  eex = sqrt(2/width)*(x-0.2*L);
  next_psi(x_i) = e^(-(eex^2)/2-kx*sqrt(-1)*eex);
  %}
end

figure;
xval = d:d:L-d;
hLine = plot(xval, abs(next_psi).^2);
axis([d, L-d, 0, 1.5*max(abs(next_psi).^2)])


for t_i = 1:300
  t = t_i*dt;
  next_psi = M*next_psi;
  set(hLine, 'YData', abs(next_psi).^2);
  drawnow
end
