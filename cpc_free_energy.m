%% Explore dynamics of beta
clear
% Define priors
a_0 = 1;
b_0 = 1;
mu_0 = 0;
la_0 = 1;

N = 1e5;
x = normrnd(5,0.5,N,1);

betas = [0.5 1 10];

% Put data and sufficient statistics into F function
for b = 1:length(betas)
    be = betas(b);
    [F(b),mu_N(b),la_N(b),a_N(b),b_N(b)]  = vb(x,be,mu_0,la_0,a_0,b_0);
    [X, J(b),H(b)] = free_energy(x,N,be,a_N(b), b_N(b), la_N(b), mu_N(b),a_0,b_0,mu_0,la_0);
end
J
H
mu_N

all_t = sqrt(1./(a_N./b_N))

 prior = normrnd(mu_0,1/la_0,D,1);

% Plot: Free energy versus beta
figure(101)
plot(betas,F,'.--','LineWidth',2)


% Plot approx posterior and data distribution:
    sig_hat = sqrt(1/(a_N(g)/b_N(g)));
    [n1,x1] = hist(normrnd(mu_N(g),sig_hat,D,1),M);
    [n2,x2] = hist(x,M);
    [n3,x3] = hist(prior,M);
    
    bar(x2,n2/N/diff(x2(1:2)))
    bar(x3,n3/D/diff(x3(1:2)))
    bar(x1,n1/D/diff(x1(1:2)))
    aa = get(gca,'child');
    
    set(aa(1),'FaceColor','b','EdgeColor','none');
    set(aa(2),'FaceColor','k','EdgeColor','none');
    set(aa(3),'FaceColor','r','EdgeColor','none');


