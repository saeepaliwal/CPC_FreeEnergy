function [F,J,H] = free_energy(x,N,be,a_N, b_N, la_N, mu_N,a_0,b_0,mu_0,la_0)
% F = < ln p(D|mu,tau) + ln p(mu|tau) + ln(p(tau)) - ln q(mu) - ln q(tau) >_q

mult=10^(round(log10(1/a_N)));
a_N = a_N*mult;
b_N = b_N*mult;

% < ln p(D|mu,tau)>_q
T1 = (N/2)*(digamma(a_N) - log(b_N) - log(2*pi)) - 0.5*(a_N/b_N)*sum((x.^2) - 2*(mu_N*x) + (mu_N^2) + (1/la_N)); % likelihood

% <ln q(mu)>_q_mu
T2 = 0.5 + 0.5*log(2*pi/la_N);

% <ln q(tau) >_q_tau
T3 = a_N - log(b_N) + gammaln(a_N) + (1-a_N)*digamma(a_N);

% <ln p(mu|tau)>_q
T4 = (N/2)*(log(la_0) + digamma(a_N) - log(b_N) - log(2*pi)) - 0.5*la_0*(a_N/b_N)*((mu_N^2) + (1/la_N) - 2*mu_N*mu_0 + (mu_0^2));

% <ln p(tau)>_q
T5 = a_0*log(b_0)-gammaln(a_0) + (a_0 - 1)*(digamma(a_N) - log(b_N)) - b_0*(a_N/b_N);

% F = < ln p(D|mu,tau)  - ln q(mu) - ln q(tau) + ln p(mu|tau) + ln(p(tau)) >_q
% Likelihood
J = T1;

% KL (constraint)
H = - T2 - T3 + T4 + T5;

F = J + (1/be)*H;


 