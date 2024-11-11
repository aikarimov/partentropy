close all
dat = readmatrix("Person51.xlsx"); %rec1, 26, male

t = dat(:,1);
N = length(t);
y = dat(:,3);
y = sgolayfilt(y,2,7); %smooth signal

figure(2); hold on
plot(t,y);
xlabel('t');
ylabel('x');
grid

%append dimension
dim = 3;
tau = 0.01;
h = mean(diff(t)); %time step
Ntau = ceil(tau/h);

Y = zeros(N - (dim-1)*Ntau,dim);
for d = 1:dim
    Y(:,d) = y((d-1)*Ntau + 1: N - (dim - d)*Ntau);
end

N = N - (dim-1)*Ntau; %reduce length

figure(3); hold on
plot(Y(:,1),Y(:,2),'.-k');
grid

k = 64; %clusters
taumax = 3; %max time for the experiment
ntimes = 20; %number of max times of leaving cluster for calculating the slope hKS
hKS = partentropy(t,Y,k,taumax,ntimes)