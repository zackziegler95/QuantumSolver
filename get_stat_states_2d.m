function [evec_sorted, eval_sorted, start] = get_stat_states_2d(hbar, m, N, dx, VV)
  % Get the sorted eigenvectors and eigenvalues of the Hamiltonian

  r = -hbar^2/(2*m*dx^2);

  % Second derivative operator
  A = diag(ones(N, 1)*(-4*r), 0) + diag(ones(N-1, 1)*(r), -1)+ diag(ones(N-1, 1)*(r), 1);

  % Form the Hamiltonian. Need to be careful converting the NxN space into a N^2 by N^2 matrix
  H = sparse(N^2, N^2);
  for n = 1:N
    H((n-1)*N+1:n*N, (n-1)*N+1:n*N) = sparse(A);
    if n != 1
      H((n-1)*N+1:n*N, (n-2)*N+1:(n-1)*N) = r*speye(N);
    end
    if n != N
      H((n-1)*N+1:n*N, n*N+1:(n+1)*N) = r*speye(N);
    end
  end

  H = H + VV;
  clear VV
  clear A
  
  start = 0;
  
  opts.issym = 1; % Symmetric matrix, so it's faster
  num = N^2;
  [eigenvectors, eigenvalues] = eigs(H, num, 'sa', opts);
  clear H
  
  [eval_sorted, I] = sort(diag(eigenvalues));

  % Get the correct eigenvalues for the eigenvectors
  evec_sorted = zeros(N^2, N^2-start);
  for i = 1:N^2-start
    evec_sorted(:, i) = eigenvectors(:, I(i));
  end

  clear I
  clear eigenvalues
  clear eigenvectors
  
end