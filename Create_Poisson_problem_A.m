function [ A ] = Create_Poisson_problem_A( N )
  A = zeros(N^2,N^2);                                 %make empty matrix
  % Create the archtypical matrix A for an N x N Poisson problem (5-point
  % stencil.
  if(N == 1)                                          %N = 1 is just 4 
      A = 4;
  end
  % Set the diagonal
  if N >= 2                                           %if N is 3 or greater 
    temp =  4 * ones(N^2);                            %prepare diagonal of 4
    new_diag = temp(1,:);
    A = A + diag(new_diag, 0);
  % Set the entries of the first sub and super diagonals
    temp =  (-1) * ones((N^2) - N);                    
    new_diag = temp(1,:);
    A = A + diag(new_diag, -N);
    A = A + diag(new_diag,  N);
  % Set the other off-diagonal entriesagonl
    temp = (-1) * ones((N^2) - 1);                    %prepare -1 diagonals
    new_diag = temp(1,:);                             %according to 7.1
    A = A + diag(new_diag, -1);                       %add lower diagonal 
    A = A + diag(new_diag,  1);                       %add higher di
  end