close all
clear all

% Yalmip modeling for symmetric OCP for the qubit: dual version on polynomials.


% Parameters
n = 3;        % Dimension of the state space
d = 8;        % Degree of the test-function polynomials 

M=3; % effective bound on the control
E=1; % size of the drift term
alpha=atan(M/E); %angle: useful for computations
k=2*(M^2+E^2)^(0.5); %idem
Tf = 2*pi/k;
xin = [0;0;1];
xf= [0;0;-1];

%  Qubit Dynamics on the Bloch sphere: real dynamics with 3 variables
f_qubit = @(x,u)([-k*cos(alpha)*x(2);k*cos(alpha)*x(1)-u*k*sin(alpha)*x(3);u*k*sin(alpha)*x(2)] );

% Yalmip sdpvar 
x = sdpvar(n,1); % the state
u = sdpvar(1,1); % the control
t = sdpvar(1,1); % time
T = sdpvar(1,1); % final time 

% Definition of symmetry-adapted monomial lists
v=monolist([t;x;u],d/2);
v1=[];
v2=[];
for k=1:length(v)
    if mod(degree(v(k),x(1))+degree(v(k),x(2))+degree(v(k),u),2)==0
  v2=[v2;v(k)];;
  v1=[v1];
    else
  v2=[v2]
  v1=[v1;v(k)];
    end
end


wx=monolist(x,(d-2)/2);
w1x=[];
w2x=[];
for k=1:length(wx)
    if mod(degree(wx(k),x(1))+degree(wx(k),x(2)),2)==0
  w2x=[w2x;wx(k)];
  w1x=w1x;
    else
  w2x=w2x;
  w1x=[w1x;wx(k)];
    end
end


wtx=monolist([t;x],d-1);
wTx=monolist([T;x],d-1);

wu=monolist(u,(d-2)/2);
w1u=[];
w2u=[];
for k=1:length(wu)
    if mod(degree(wu(k),u),2)==0
  w2u=[w2u;wu(k)];
  w1u=w1u;
    else
  w2u=w2u;
  w1u=[w1u;wu(k)];
    end
end

Wt=monolist([t],(d-2)/2);
WT=monolist([T],(d-2)/2-1);

% Definition of the blocks of symmetry-adapted matrices associated with SOS multipliers

Qg = sdpvar(1,length(wtx));

Q1 = sdpvar(length(v1)); 
Q2 = sdpvar(length(v2));

Q1u = sdpvar(length(w1u));
Q2u = sdpvar(length(w2u));

Q1x = sdpvar(length(w1x));
Q2x = sdpvar(length(w2x));

Qt = sdpvar(length(Wt));
QT = sdpvar(length(WT));

% Vector field
f = f_qubit(x,u);

% State and control constraints
cx = 1 - (x(1)^2+x(2)^2+x(3)^2); %equality constraint
cu = 1 - u^2;
ct = t*(Tf - t);
cT = T*(Tf - T);

% Polynomials
g= Qg*wtx;
g0 = replace(g,[t;x],[0;xin]);
gg=Qg*wTx;
gT = replace(gg,[x],[xf]);




% Sum of squares multipliers

s1 = v1'*Q1*v1+v2'*Q2*v2;
su=w1u'*Q1u*w1u+w2u'*Q2u*w2u;
sx=w1x'*Q1x*w1x+w2x'*Q2x*w2x;
st=Wt'*Qt*Wt;
sT =WT'*QT*WT;

% Positivity condition: A\phi+1 \geq 0 on X
% and \phi(T,x)\leq 0 on the target set K

constraint = [coefficients(jacobian(g,t)+jacobian(g,x)*f+1 -ct*st - cx*sx - cu*su -s1,[t;x;u])==0; Q1>=0; Q2>=0;Q1u>=0; Q2u>=0;Qt>=0]; % the constraint of cx is here a sphere constraint ||x||=1
constraint= [constraint;coefficients(-gT-cT*sT,[T])==0;QT>=0]; %constraint at the final time

% Objective 
obj = g0;

SDPsolver = lower('sedumi');
[Fd,objd,X,free] = dualize(constraint,-obj)

tic 
optimize(Fd,-objd);
toc
objf=-double(objd);

