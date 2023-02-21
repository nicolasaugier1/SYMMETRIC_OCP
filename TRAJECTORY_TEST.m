% optimal control 
% minimal time double integrator with state constraint
clear all
close all
mset clear

d=6; % degree of relaxation
x0 = [0; 0;1]; u0 = 0; % initial condition at the north pole

%parameters of the dynamics
M=3;
E=1;
alpha=atan(M/E); 
k=2*(M^2+E^2)^(0.5);
Tf=2*pi/k;
eps=0.5;

mpol t
mpol x 3
mpol v



m = meas([t;x;v]); %occupation measure

mpol T
mpol xT 3
mT = meas([T;xT]); % terminal measure

g = mmon([t;x],d-1);
gT = mmon([T;xT],d-1);
assign([t;x;v],[0;x0;u0]);
g0 = double(g);


f = [- k*cos(alpha)*x(2);k*cos(alpha)*x(1)-v*k*sin(alpha)*x(3);v*k*sin(alpha)*x(2)]; %vectorfield

 
%% optimal traj
TTT=linspace(0,2*pi/k,1000)';
U1=[];
Uoptbis=[];
Uopttime=[];
Uoptbistime=[];

for j=1:size(TTT)
    if TTT(j)-(pi/k-acos(1/(tan(alpha))^2)/k)<0
      U1=[U1,-1];
    else U1=[U1,+1];
    end
end 

dt = TTT(2)-TTT(1);

t_s1 = pi/k-acos(1/(tan(alpha))^2)/k;
[~,idx] = min(abs(TTT-t_s1));
clear u 
u = ones(size(TTT));
u(idx:end) = -1;

%Definition of a good trajectory

ybon = [0;0;1];

for i = 1:length(TTT)-1
    ybon(:,i+1) = ybon(:,i)+dt*[-k*cos(alpha)*ybon(2,i);
                k*cos(alpha)*ybon(1,i)-u(i)*k*sin(alpha)*ybon(3,i);
                u(i)*k*sin(alpha)*ybon(2,i)];
end

%Definition of a bad trajectory
ybad2=0;
for i = 1:length(TTT)/2
    ybad2(i+1) = ybon(2,i+1)
end

for i = length(TTT)/2:length(TTT)-1
    ybad2(i+1) = -ybon(2,i+1)
end

% Good first order moment
MOMx=0;

for i=1:length(TTT)-1 
    MOMx = MOMx+ybon(2,i)*dt;                
end

% Bad first order moment of bad solution ybad2

MOMxbad=0;

for i=1:length(TTT)-1 
    MOMxbad = MOMxbad+ybad2(i)*dt;                
end


%% SOLVE


tic
P = msdp(min(mass(m)),x(1)^2+x(2)^2+x(3)^2==1,-(mom(x(2))-MOMx)<=eps,(mom(x(2))-MOMx)<=eps,xT-[0;0;-1]==0,T*(T-Tf)<=0,t*(t-Tf)<=0,v^2<=1,0 == mom(diff(g,t)+diff(g,x)*f) - mom(gT) + mom(g0)); 
[status,obj] = msol(P);
toc

%the problem is feasible 

tic
Q = msdp(min(mass(m)),x(1)^2+x(2)^2+x(3)^2==1,-(mom(x(2))-MOMxbad)<=eps,(mom(x(2))-MOMxbad)<=eps,xT-[0;0;-1]==0,T*(T-Tf)<=0,t*(t-Tf)<=0,v^2<=1,0 == mom(diff(g,t)+diff(g,x)*f) - mom(gT) + mom(g0)); 
[statusq,objq] = msol(Q);
toc

%the problem is not feasible


