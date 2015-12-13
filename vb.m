function [F,mu_N,la_N,a_N,b_N] =vb(x,be,mu_0,la_0,a_0,b_0)

max_iter = 100;
N = length(x);

% Initial values for posterior
a_N = a_0; 
b_N = b_0;
mu_N = 0;

% Initialize FE
F = -inf;

for i = 1:max_iter
    
    la_N = (la_0 + be*N)*a_N/b_N;
    %la_N = (la_0 + N)*a_N/b_N;
    
    mu_N = ((la_0*mu_0 + be*N*mean(x))/(la_0 + be*N));
    %mu_N = ((la_0*mu_0 + N*mean(x))/(la_0 + N));
    
    a_N = (a_0 + be*N/2);
    %a_N = (a_0 + N/2);
    
    
    %b_N = b_0 + (1/2)*(be*sum((x-mu_N).^2) + la_0*((mu_N-mu_0)^2));
    b_N = b_0 + (1/2)*(sum((x-mu_N).^2) + la_0*((mu_N-mu_0)^2));
    

%     la_N = (la_0 + N)*a_N/b_N;
%     mu_N = ((la_0*mu_0 + N*mean(x))/(la_0 + N));
%   
%     a_N = (a_0 + *N/2);
%     b_N = b_0 + (1/2)*(sum((x-mu_N).^2) + la_0*((mu_N-mu_0)^2));

    
    F_old = F;
    F = free_energy(x,N,be,a_N, b_N, la_N, mu_N,a_0,b_0,mu_0,la_0);
    
    % Convergence criterion
    if (F - F_old < 10e-4), break; end
    if (i == max_iter), warning('Reached %d iterations',max_iter); end

end
[F] = free_energy(x,N,be,a_N, b_N, la_N, mu_N,a_0,b_0,mu_0,la_0);
