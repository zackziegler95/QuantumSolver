% Constants for the simulation
hbar = 5;
m = 1;
w = 3;

L = 10; % Length of space
d = 0.01; % spatial step
dt = 4*pi/w/100; % time step
N = L/d-1;

% Second derivative
d2_op = sparse(diag(ones(N, 1)*(-2), 0) + diag(ones(N-1, 1), -1)+ diag(ones(N-1, 1), 1));

V = sparse(N, N);

c = 1;
% A few different potentials
for x_i = 1:N
  x = (x_i)*d;
  %V(x_i,x_i) = e^((x-L/2)^2/(2*(L/4)^2));
  V(x_i,x_i) = 1/2*m*w^2*abs(x-L/2)^2; % Harmonic oscillator. Classic.
  %V(x,x) = w*x;
  %{
  V(x_i,x_i) = c;
  if x > L/4 && x < 3*L/4
    V(x_i,x_i) = c/2;
  end
  if x > 3*L/8 && x > 5*L/8
    V(x_i,x_i) = 0;
  end
  %}
end

% Hamiltonian
H = -hbar^2/(2*m*d^2)*d2_op + V;
clear d2_op


[eigenvectors, eigenvalues] = eigs(H, N, 'sa');
clear H
[eval_sorted, I] = sort(diag(eigenvalues));

% Need to get the right eigenvalues with the eigenvectors
evec_sorted = zeros(N, N);
for i = 1:N
  evec_sorted(:, i) = eigenvectors(:, I(i));
end
%eval_sorted

clear I
clear eigenvalues
clear eigenvectors

%{
% Test for scaling of energies
enums = 1:10;
plot(enums, eval_sorted(enums), 'LineWidth', 2, enums, (enums).^2*3.1415^2*hbar^2/(2*m*L^2), 'LineWidth', 2, enums, (enums-0.5)*hbar*w, 'LineWidth', 2);
legend('Actual', 'Inf. Sq. Well', 'Harmonic osc.')
%}

figure;
hold on
colors = ['blue'; 'red'; 'green'; 'cyan'; 'red'; 'red'; 'red'; 'brown'; 'brown'; 'brown'; 'brown';];
for i = 1:4
  plot(eval_sorted(i)+5000*evec_sorted(:, i).^2, 'color', colors(i), 'LineWidth', 2)
end
plot(diag(V))

clear V

% Set the initial conditions
initial_state = zeros(N, 1);
for x_i = 1:N
  x = (x_i)*d;
  initial_state(x_i) = e^(-(x-L/2)^2);
  %initial_state(x_i) = sin(pi/L*x)+sin(2*pi/L*x);
  %initial_state(x_i) = e^(-(x-0.2*L)^2+sqrt(-1)*1*(x-0.2*L));
  %initial_state(x_i) = sin(pi/L*x);
  %c0 = 0.2;
  %kx = 25;
  %width = (0.005*L);
  %eex = sqrt(2/width)*(x-0.2*L);
  %initial_state(x_i) = e^(-(eex^2)/2-kx*sqrt(-1)*eex);
end

%initial_state = evec_sorted(:, 3) + 2i*evec_sorted(:, 4)+ 8*evec_sorted(:, 5)+ (2.3+3i)*evec_sorted(:, 6)+ evec_sorted(:, 7);
initial_state /= norm(initial_state);
%plot(abs(initial_state).^2)

% Set c_n for each eigenvalue/eigenvector pair
c_n = zeros(N, 1);
for x_i = 1:N
  c_n(x_i) = dot(initial_state, evec_sorted(:, x_i));
end



figure('Position',[100 100 1200 600]);

#start_button = uicontrol('Style','pushbutton', 'String','Start','Position',[1000,500,70,25]);
#stop_button = uicontrol('Style','pushbutton', 'String','Stop','Position',[1000,440,70,25]);
#ha = axes('Units','pixels','Position',[50,60,900,500]);
#xlim([d, L-d])
#ylim([0, 8*max(abs(initial_state).^2)]);
#set(f, 'CurrentAxes', ha2);
#set(f, 'Visible', 'on');

xval = d:d:L-d;
ydata = plot(xval, zeros(N, 1));

axis([d, L-d, 0, 1.5*max(abs(initial_state).^2)]);

% Propagate the wavefunction. Also calculate the expectation value of x, for example
t_tot = 100;
exp_overt = zeros(t_tot, 1); % Expectation value of x as a funciton of t
for t_i = 1:t_tot
  t = t_i * dt;
  result = zeros(N, 1);
  
  % Step by propagating the eigenvectors
  for x_i = 1:N
    result += c_n(x_i)*evec_sorted(:, x_i)*e^(-sqrt(-1)*eval_sorted(x_i)*t/hbar);
  end
  
  exp_overt(t_i) = dot(result, ([d:d:L-d]').*result);
  
  set(ydata, 'YData', abs(result).^2);
  drawnow
end

figure;
plot([1:t_tot], exp_overt)