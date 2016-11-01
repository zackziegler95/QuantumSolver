% Constants for the simulation
hbar = 5;
m = 1;
w = 0.8;
c = 50;

L = 10; % Spatial length
d = 0.2; % Space step
dt = 0.03; % Time step
N = L/d-1;

% Set the potential
VV = sparse(N^2, N^2);
for x_i = 1:N
  x = (x_i-1)*d;
  for y_i = 1:N
    y = (y_i-1)*d;
    ind = y_i+(x_i-1)*N;
    if x < 0.4*L || x > 0.6*L
      #VV(ind, ind) = 0;
    else
      #VV(ind, ind) = c;
    end
  end
end

% Evaluate the eigenvectors and eigenvalues
tic
[evec_sorted, eval_sorted, start] = get_stat_states_2d(hbar, m, N, d, VV);
toc
clear VV

eval_sorted(1:10) % print the first 10 eigenvalues
figure;
surf(reshape(abs(evec_sorted(:, 1)+evec_sorted(:, 2)).^2, [N, N]));
xlabel('x')
ylabel('y')

%{
% Testing the analytical solution for the infinite square well
real_eval = zeros(N, N);
for xn = 1:N
  for yn = 1:N
    real_eval(yn, xn) = pi^2*hbar^2/(2*m*L^2)*(xn^2+yn^2);
  end
end
real_eval = sort(reshape(real_eval, [N^2, 1]));
real_eval(1:10)
%}

initial_state = zeros(N*N,1);
for x_i = 1:N
  x = (x_i)*d;
  for y_i = 1:N
    y = (y_i)*d;
    ind = y_i + (x_i-1)*N;
    %initial_state(ind) = (sin(pi/(L)*x))*sin(pi/(L)*y);
    c0 = 0.2;
    width = (c0*L);
    kx0 = -0.3/(c0*d);
    %ky0 = 1;
    eex = sqrt(2/width)*(x-0.2*L);
    %eey = sqrt(2/width)*(y-L/2);
    initial_state(ind) = 0.5*e^(-(eex^2)/2+sqrt(-1)*kx0*eex)*sin(pi/L*y);
    %initial_state(ind) = e^(-(eex^2)/2)*(1)*e^(-(eey^2)/2);
    %initial_state(ind) = (1 + 1/sqrt(2)*2*eex + 1/sqrt(8)*(4*eex^2 -2))*e^(-(eex^2)/2)*(1 + 1/sqrt(2)*2*eey)*e^(-(eey^2)/2);
    %initial_state(ind) = e^(-m*w/hbar*((x-L/2)^2+(y-L/2)^2)/2); + 1/sqrt(2)*2*sqrt(m*w/hbar)*(x-L/2)*e^(-(sqrt(m*w/hbar)*(x-L/2))^2/2);
    %initial_state(ind) = (32*eee^5 - 160*eee^3 + 120*eee)*e^(-(eee^2)/2);
    %initial_state(ind) = sin(pi/L*x)+0.0001*sin(2*pi/L*x);
  end
end
initial_state /= norm(initial_state);

% C_n for each eigenvector
c_n = zeros(N^2-start, 1);
for x_i = 1:N^2-start
  c_n(x_i) = dot(initial_state, evec_sorted(:, x_i));
end
norm(c_n)
c_n(1:20);

max = 0.04*d;
figure;
xyval = d:d:L-d;
zplot = surf(xyval, xyval, reshape(abs(initial_state).^2, [N, N]), 'EdgeColor', 'None', 'facecolor', 'interp');
#axis([d, L-d, d, L-d, 0, max])
axis([d, L-d, d, L-d])
xlabel('x')
ylabel('y')
view(2)

% Then go!
for t_i = 1:50
  t = t_i*dt;
  
  next_psi = zeros(N^2, 1);
  for x_i = 1:N^2-start
    next_psi += c_n(x_i)*evec_sorted(:, x_i)*e^(-sqrt(-1)*eval_sorted(x_i)*t/hbar);
  end
  
  set(zplot, 'ZData', reshape(abs(next_psi).^2, [N, N]));
  set(zplot, 'CData', reshape(abs(next_psi).^2, [N, N]));
  drawnow
end
