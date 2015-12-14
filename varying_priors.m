%% Explore dynamics of beta
clear

addpath('~/Dropbox/Tools'));

% Define priors
a_0_all = [1 1 5];
b_0_all = [5 5 1];
mu_0_all = [-50 0 -30];
la_0_all = [1 1 1];

N = 30;
x = normrnd(20,0.5,N,1);

be = 1;

% Put data and sufficient statistics into F function
for b = 1:length(mu_0_all)
    mu_0 = mu_0_all(b);
    la_0 = la_0_all(b);
    a_0 = a_0_all(b);
    b_0 = b_0_all(b);
    [F(b),mu_N(b),la_N(b),a_N(b),b_N(b)]  = vb(x,be,mu_0,la_0,a_0,b_0);
    [X, J(b),H(b)] = free_energy(x,N,be,a_N(b), b_N(b), la_N(b), mu_N(b),a_0,b_0,mu_0,la_0);
end

D = 100000;
x = normrnd(20,0.5,D,1);
M = 50;
tau_all = sqrt(1./(a_0_all./b_0_all));
sig_hat = sqrt(1./(a_N./b_N));
%sig_0_all = sqrt(1./la_0_all.*tau_all);

prior1 = normrnd(mu_0_all(1),tau_all(1),D,1);
prior2 = normrnd(mu_0_all(2),tau_all(2),D,1);
prior3 = normrnd(mu_0_all(3),tau_all(3),D,1);
% Plot: Free energy versus beta
% figure(101)
% plot(betas,F,'.--','LineWidth',2)

% Plot approx posterior and data distribution:

% figure(101);
% clf
% hold on;
% [n1,x1] = hist(normrnd(mu_N(1),sig_hat(1),D,1),M);
% [n2,x2] = hist(x,M);
% [n3,x3] = hist(prior1,M);
% 
% bar(x2,n2/D/diff(x2(1:2)))
% bar(x3,n3/D/diff(x3(1:2)))
% bar(x1,n1/D/diff(x1(1:2)))
% aa = get(gca,'child');
% 
% set(aa(1),'FaceColor','none','EdgeColor','b');
% set(aa(2),'FaceColor','none','EdgeColor','k');
% set(aa(3),'FaceColor','none','EdgeColor','r');
% xlim([-20 40]);
% ylim([0 1]);
% 
figure(101);
clf
hold on;
[n1,x1] = hist(normrnd(mu_N(2),sig_hat(2),D,1),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior2,M);

bar(x2,n2/D/diff(x2(1:2)))
bar(x3,n3/D/diff(x3(1:2)))
bar(x1,n1/D/diff(x1(1:2)))
aa = get(gca,'child');

set(aa(1),'FaceColor','none','EdgeColor','b');
set(aa(2),'FaceColor','none','EdgeColor','k');
set(aa(3),'FaceColor','none','EdgeColor','r');
xlim([-40 40]);
ylim([0 1]);
purty_plot(101,'good_prior');

figure(102);
clf
hold on;
[n1,x1] = hist(normrnd(mu_N(3),sig_hat(3),D,1),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior3,M);

bar(x2,n2/D/diff(x2(1:2)))
bar(x3,n3/D/diff(x3(1:2)))
bar(x1,n1/D/diff(x1(1:2)))
aa = get(gca,'child');

set(aa(1),'FaceColor','none','EdgeColor','b');
set(aa(2),'FaceColor','none','EdgeColor','k');
set(aa(3),'FaceColor','none','EdgeColor','r');
xlim([-40 40]);
ylim([0 1]);
purty_plot(102,'bad_prior');

keyboard

% x = -100:100;
% 
% prior = normpdf(x,0,1);
% likelihood = normpdf(x,5,0.5);
% posterior = prior.*likelihood;
% posterior = posterior./sum(posterior);



