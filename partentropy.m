function [hKS,hP] = partentropy(varargin)
% PARTENTROPY calculate partitioned entropy
%   [hP, hKS] = PARTENTROPY(t,Y,k,taumax)
%   [hP, hKS] = PARTENTROPY(t,Y,k,taumax,ntimes)
%   input:
%   t - time signal, N x 1 vector
%   Y - system trajectory, N x dim vector
%   k - number of clusters
%   taumax - greatest value of tau for the statistics h_P(tau)
%   ntimes - how many time to repeat max time of leaving the initial
%   cluster for estimating hKS (default is 6)
%   output

t = varargin{1,1};
Y = varargin{1,2};
k = varargin{1,3};
taumax = varargin{1,4};
ntimes = 6;
if nargin > 4
    ntimes = varargin{1,5};
end


h = mean(diff(t)); %time step
[N, dim] = size(Y);

[idx,C] = kmeans(Y,k);

if dim >= 2
    figure;
    hold on
    plot(Y(:,1),Y(:,2),'.-k');
    voronoi(C(:,1),C(:,2));
    xlabel('y_1');
    ylabel('y_2');
end

% example illustration
Nbins = k;
p = histcounts(idx,Nbins);
rhos = p/sum(p); %determine rho

tau_span = h:h:taumax;
nv = length(tau_span);

Nleave = N - ceil(tau_span(end)/h);
tleave = zeros(Nleave,1); %mean time of leaving clusters

entrs = zeros(1,nv);
for ti = 1:nv
    tau = tau_span(ti);
    ns = ceil(tau/h);
    idxold = idx(1:N-ns);
    idxnew = idx(ns + 1:N);
    cngs = idxold ~= idxnew; %points where idx changed

    idleave = (tleave == 0) & cngs(1:Nleave); %poins which left their cluster at tau
    tleave(idleave) = tau;

    sumet = 0;
    for i = 1:k
        %get all events of i cluster
        ni = (idxold == i) &  cngs;
        edges = 0.5:k+0.5;
        n = histcounts(idxnew(ni),edges);
        idsz = n==0;
        n(idsz) = []; %remove zero entries
        n = n/sum(n);
        hi = -rhos(i)*sum(n.*log(n)); %get i-th entropy calculation
        sumet = sumet + hi;
    end
    entrs(ti) = sumet;
end



%times of leaving cluster
taulmax = max(tleave);


%then, polynomial fit to line
nmax = ceil(taulmax/h);
nt = length(tau_span);
ids = nmax:min((ntimes+1)*nmax,nt);
ts = tau_span(ids);
es = entrs(ids);
E = [ts', ones(length(ids),1)];
V = es';
h = E \ V; % LSM

hKS = h(1); %KS entropy estimation
hP = entrs;

figure;
hold on
plot(tau_span,entrs);
plot([taulmax taulmax],[min(entrs), max(entrs)],'--r');
plot(ts, ts*h(1) + h(2),'-','LineWidth',1);
legend('$h_p$','$\tau_{0}$',['fit $h_{KS} = ',num2str(hKS),'$'],'Interpreter','latex');
xlabel('$\tau$','Interpreter','latex');
ylabel('$h_p$','Interpreter','latex');

end






