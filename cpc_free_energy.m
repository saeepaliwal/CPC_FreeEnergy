% Explore dynamics of beta
clear

% Define priors
a_0 = 1;
b_0 = 1;
mu_0 = 0;
la_0 = 1;

% Observed data
N = 10;
x = normrnd(5,0.5,10,1);

alphas = [0.5 15 100];

% Put data and sufficient statistics into F function
for b = 1:length(alphas)
    be = alphas(b);
    [F(b),mu_N(b),la_N(b),a_N(b),b_N(b)]  = vb(x,be,mu_0,la_0,a_0,b_0);
    [X, J(b),H(b)] = free_energy(x,N,be,a_N(b), b_N(b), la_N(b), mu_N(b),a_0,b_0,mu_0,la_0);
end


% Plot it all up!
D = 10000;
N = D;

x = normrnd(5,0.5,N,1);
M = 10;
sig_hat = sqrt(1./(a_N./b_N));
prior = normrnd(mu_0,1/la_0,D,1);

% Plot approx posterior and data distribution:
figure(101);
clf
hold on;
final_distro(:,1) = normrnd(mu_N(1),sig_hat(1),D,1);
[n1,x1] = hist(final_distro(:,1),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior,M);

bar(x2,n2/N/diff(x2(1:2)))
bar(x3,n3/D/diff(x3(1:2)))
bar(x1,n1/D/diff(x1(1:2)))
aa = get(gca,'child');

set(aa(1),'FaceColor','b','EdgeColor','none');
set(aa(2),'FaceColor','k','EdgeColor','none');
set(aa(3),'FaceColor','r','EdgeColor','none');
xlim([-6 10]);
ylim([0 1]);
purty_plot(101,'low_beta','png');

figure(102);
clf
hold on;
final_distro(:,2) = normrnd(mu_N(2),sig_hat(2),D,1);
[n1,x1] = hist(final_distro(:,2),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior,M);

bar(x2,n2/N/diff(x2(1:2)))
bar(x3,n3/D/diff(x3(1:2)))
bar(x1,n1/D/diff(x1(1:2)))
aa = get(gca,'child');

set(aa(1),'FaceColor','b','EdgeColor','none');
set(aa(2),'FaceColor','k','EdgeColor','none');
set(aa(3),'FaceColor','r','EdgeColor','none');
xlim([-6 10]);
ylim([0 1]);
purty_plot(102,'normal_beta','png');
figure(103);
clf
hold on;
final_distro(:,3) = normrnd(mu_N(3),sig_hat(3),D,1);
[n1,x1] = hist(final_distro(:,3),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior,M);

bar(x2,n2/N/diff(x2(1:2)))
bar(x3,n3/D/diff(x3(1:2)))
bar(x1,n1/D/diff(x1(1:2)))
aa = get(gca,'child');

set(aa(1),'FaceColor','b','EdgeColor','none');
set(aa(2),'FaceColor','k','EdgeColor','none');
set(aa(3),'FaceColor','r','EdgeColor','none');
xlim([-6 10]);
ylim([0 1]);
purty_plot(103,'high_beta','png');

figure(104)
bar(F)
ylabel('-\Delta F');
aa = get(gca,'child');
set(aa(1),'FaceColor',[0 0.5 0],'EdgeColor','none');
set(gca,'XTickLabel',{'Low alpha','Normal alpha','Large alpha'})
purty_plot(104,'fe_diff','png');

%Entropy
y = -20:0.1:20;
likelihood_info = 0.5*log2(2*pi*0.25) + 0.5;
baseline = 0.5*log2(2*pi) + 0.5;
for k = 1:3
    information(:,k) = 0.5*log2(2*pi*(sig_hat(k)^2)) + 0.5;
end
