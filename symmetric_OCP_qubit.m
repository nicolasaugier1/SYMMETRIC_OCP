% optimal control with symmetry
% minimal time qubit
% Gloptipoly modeling with moments

clear all
close all
mset clear


%% Parameters

d = 10; % degree of relaxation
x0 = [0; 0;1]; % initial condition at the north pole
u0 = 0; 
M=3;
E=1;
alpha=atan(M/E); 
k=2*(M^2+E^2)^(0.5);
Tf=2*pi/k;


% First optimization variables
mpol t
mpol x 3
mpol u

mpol T
mpol xT 3

% Second optimization variables
mpol tb
mpol z 3
mpol v

mpol Tb
mpol zTb 3


% Third optimization variables
mpol tbb
mpol zb 3
mpol vb


mpol Tbb
mpol zTbb 3

% Vector fields

f = [- k*cos(alpha)*x(2);k*cos(alpha)*x(1)-u*k*sin(alpha)*x(3);u*k*sin(alpha)*x(2)]; % vector field 1
fbis = [- k*cos(alpha)*z(2);k*cos(alpha)*z(1)-v*k*sin(alpha)*z(3);v*k*sin(alpha)*z(2)]; % vector field 2
fbisb = [- k*cos(alpha)*zb(2);k*cos(alpha)*zb(1)-vb*k*sin(alpha)*zb(3);vb*k*sin(alpha)*zb(2)]; %vectorfield 3

%% Definition of measures and test polynomial functions

% First optimization step (A_1)
m = meas([t;x;u]); %first occupation measure
mT = meas([T;xT]); % first terminal measure

% Test functions
g = mmon([t;x],d-1);
gT = mmon([T;xT],d-1);
assign([t;x;u],[0;x0;u0]);
g0 = double(g);

% Second optimization step (A_2)
mbis = meas([tb;z;v]); % second occupation measure
mTb = meas([Tb;zTb]); % second terminal measure

gbis = mmon([tb;z],d-1);
gTb = mmon([Tb;zTb],d-1);
assign([tb;z;v],[0;x0;u0]);
g0bis = double(gbis);

% Third optimization step
mbisb = meas([tbb;zb;vb]); %third occupation measure
mTb = meas([Tbb;zTbb]); % third terminal measure

gbisb = mmon([tbb;zb],d-1);
gTbb = mmon([Tbb;zTbb],d-1);
assign([tbb;zb;vb],[0;x0;u0]);
g0bisb = double(gbisb);

%% List of moments for step (A_1)

Ot=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if mod(j+q+l-3,2)==1 && j+q+l+p+r-5<=d;
                    Ot=[Ot,mom(t^(p-1)*x(1)^(j-1)*x(2)^(q-1)*x(3)^(r-1)*u^(l-1))];
                  end
                end
            end  
        end
    end
end

OT=[];

for j=1:d+1
    for q=1:d+1
            for r=1:d+1
                for p=1:d+1
                   if mod(j+q-2,2)==1 && p+j+q+r-4<=d;
                   OT=[OT,mom(T^(p-1)*xT(1)^(j-1)*xT(2)^(q-1)*xT(3)^(r-1))];
                   end
                end
            end
    end
end

%List of all moments

Otot=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if j+q+l+p+r-5<=d;
                  Otot=[Otot,mom(t^(p-1)*x(1)^(j-1)*x(2)^(q-1)*x(3)^(r-1)*u^(l-1))];
                  end
                end
            end  
        end
    end
end





OTtot=[];

for j=1:d+1
    for q=1:d+1
            for r=1:d+1
                for p=1:d+1
                   if p+j+q+r-4<=d;
                   OTtot=[OTtot,mom(T^(p-1)*xT(1)^(j-1)*xT(2)^(q-1)*xT(3)^(r-1))];
                   end
                end
            end
    end
end


%% List of moments for step 2

List=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if mod(j+q+l-3,2)==0 & j+q+l+p+r-5<=d
                  List=[List,(1-2*rand(1))*mom(tb^(p-1)*z(1)^(j-1)*z(2)^(q-1)*z(3)^(r-1)*v^(l-1))];
                  end
                end
            end  
        end
    end
end

Listede1=ones(size(List));

Otb=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if mod(j+q+l-3,2)==1 && j+q+l+p+r-5<=d;
                  Otb=[Otb,mom(tb^(p-1)*z(1)^(j-1)*z(2)^(q-1)*z(3)^(r-1)*v^(l-1))];
                  end
                end
            end  
        end
    end
end

 OTb=[];
 for j=1:d+1
     for q=1:d+1
            for r=1:d+1
                 for p=1:d+1
                    if mod(j+q-2,2)==1 && p+j+q+r-4<=d;
                    OTb=[OTb,mom(Tb^(p-1)*zTb(1)^(j-1)*zTb(2)^(q-1)*zTb(3)^(r-1))];
                    end
                 end
           end
    end
 end
 
%% Costcomputation without symmetry with SDP (Q_k)
%tic
%P = msdp(min(mass(m)),x(1)^2+x(2)^2+x(3)^2==1,T*(T-Tf)<=0,t*(t-Tf)<=0,u^2<=1,xT-[0;0;-1]==0,0 == mom(diff(g,t)+diff(g,x)*f) - mom(gT) + mom(g0)); 
%[status,obj] = msol(P);
%toc
%% Costcomputation with symmetry with SDP (I_k^G)

tic
P = msdp(min(mass(m)),x(1)^2+x(2)^2+x(3)^2==1,Ot==0,OT==0,T*(T-Tf)<=0,t*(t-Tf)<=0,u^2<=1,(xT'-[0;0;-1]')*(xT-[0;0;-1])<=0,0 == mom(diff(g,t)+diff(g,x)*f) - mom(gT) + mom(g0)); 
[status,obj] = msol(P);
toc


%% Computation with linear functional and symmetry (A_2), SDP (Z_k^G)

tic
Q = msdp(min(List*Listede1'),mass(mbis)<=obj+0.04,Otb==0,OTb==0,z(1)^2+z(2)^2+z(3)^2==1,(Tf-tb)*tb>=0,(Tf-Tb)*Tb>=0, v^2<=1, (zTb'-[0;0;-1]')*(zTb-[0;0;-1])<=0,0 == mom(diff(gbis,tb)+diff(gbis,z)*fbis) - mom(gTb) + mom(g0bis)); 
[statusquo,objbis] = msol(Q);
toc 

%% Definition of symmetry invariant moments

B=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if 2*j+2*q+2*l+p+r-8<=d;
                  B=[B,mom(tb^(p-1)*z(1)^(2*(j-1))*z(2)^(2*(q-1))*z(3)^(r-1)*v^(2*(l-1)))];
                  end
                end
            end  
        end
    end
end

BT=[];
for j=1:d+1
    for q=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if 2*j+2*q+p+r-6<=d;
                  BT=[BT,mom(Tb^(p-1)*zTb(1)^(2*(j-1))*zTb(2)^(2*(q-1))*zTb(3)^(r-1))];
                  end
                end
            end  
    end
end

%% TRAJECTORY RECONSTRUCTION WITH CD KERNELS AND THEORETICAL OPTIMAL TRAJECTORIES

b = mmon([tb;z],d/2);
b1 = mmon([tb;z(1)],d/2);
b2 = mmon([tb;z(2)],d/2);
b3 = mmon([tb;z(3)],d/2);
%bv = mmon([tb;v],d/2);

G = double(mom(b*b'));
G1 = double(mom(b1*b1'));
G2 = double(mom(b2*b2'));
G3 = double(mom(b3*b3'));
%Gv = double(mom(bv*bv'));

TT = linspace(0,obj,1e3)';

[XXL1,p1] = momgraph(G1,TT);
[XXL2,p2] = momgraph(G2,TT);
[XXL3,p3] = momgraph(G3,TT);
%[XXLv,pv] = momgraph(Gv,TT);


 
%% optimal controls and trajectories

TTT=linspace(0,2*pi/k,1000)';
U1=[];
U2=[];
U3=[];
U4=[];

for j=1:size(TTT)
    if TTT(j)-(pi/k-acos(1/(tan(alpha))^2)/k)<0
      U1=[U1,-1];
    else U1=[U1,+1];
    end
end 
for j=1:size(TTT)
    if TTT(j)-(pi/k-acos(1/(tan(alpha))^2)/k)<0
      U2=[U2,1];
    else U2=[U2,-1];
    end
end 
for j=1:size(TTT)
    if TTT(j)-(pi/k+acos(1/(tan(alpha))^2)/k)<0
      U3=[U3,-1];
    else U3=[U3,+1];
    end
end 

for j=1:size(TTT)
    if TTT(j)-(pi/k+acos(1/(tan(alpha))^2)/k)<0
      U4=[U4,1];
    else U4=[U4,-1];
    end
end 

% Optimal trajectories 

dt = TTT(2)-TTT(1);

t_s1 = pi/k-acos(1/(tan(alpha))^2)/k;
[~,idx] = min(abs(TTT-t_s1));
clear u 
u = ones(size(TTT));
u(idx:end) = -1;


y1 = [0;0;1];

for i = 1:length(TTT)-1
    y1(:,i+1) = y1(:,i)+dt*[-k*cos(alpha)*y1(2,i);
                k*cos(alpha)*y1(1,i)-u(i)*k*sin(alpha)*y1(3,i);
                u(i)*k*sin(alpha)*y1(2,i)];
end

clear u
u = ones(size(TTT));
u(1:idx) = -1;

y2 = [0;0;1];

for i = 1:length(TTT)-1
    y2(:,i+1) = y2(:,i)+dt*[-k*cos(alpha)*y2(2,i);
                k*cos(alpha)*y2(1,i)-u(i)*k*sin(alpha)*y2(3,i);
                u(i)*k*sin(alpha)*y2(2,i)];
end



t_s2 = pi/k+acos(1/(tan(alpha))^2)/k;
[~,idx] = min(abs(TTT-t_s2));

clear u
u = ones(size(TTT));
u(idx:end) = -1;

y3 = [0;0;1];

for i = 1:length(TTT)-1
    y3(:,i+1) = y3(:,i)+dt*[-k*cos(alpha)*y3(2,i);
                k*cos(alpha)*y3(1,i)-u(i)*k*sin(alpha)*y3(3,i);
                u(i)*k*sin(alpha)*y3(2,i)];
end

clear u
u = ones(size(TTT));
u(1:idx) = -1;

y4 = [0;0;1];

for i = 1:length(TTT)-1
    y4(:,i+1) = y4(:,i)+dt*[-k*cos(alpha)*y4(2,i);
                k*cos(alpha)*y4(1,i)-u(i)*k*sin(alpha)*y4(3,i);
                u(i)*k*sin(alpha)*y4(2,i)];
end 

%% List of moments for control reconstruction with SDP \tilde{Z_k}

Otbb=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if 2*j+2*q+2*l+p+r-8<=d;
                  Otbb=[Otbb,mom(tbb^(p-1)*zb(1)^(2*(j-1))*zb(2)^(2*(q-1))*zb(3)^(r-1)*vb^(2*(l-1)))];
                  end
                end
            end  
        end
    end
end


OTbb=[];
for j=1:d+1
    for q=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if 2*j+2*q+p+r-6<=d;
                    OTbb=[OTbb,mom(Tbb^(p-1)*zTbb(1)^(2*(j-1))*zTbb(2)^(2*(q-1))*zTbb(3)^(r-1))];
                  end
                end
            end  
    end
end

Listb=[];
for j=1:d+1
    for q=1:d+1
        for l=1:d+1
            for r=1:d+1
                for p=1:d+1
                  if j+q+l+p+r-5<=d
                    Listb=[Listb,(1-2*rand(1))*mom(tbb^(p-1)*zb(1)^(j-1)*zb(2)^(q-1)*zb(3)^(r-1)*vb^(l-1))];
                  end
                end
            end  
        end
    end
end

Listede1b=ones(size(Listb));
%% Control reconstruction with SDP \tilde{Z_k}
tic
Qb = msdp(min(Listb*Listede1b'),mass(mbisb)<=obj+0.04,Otbb==double(B),OTbb==double(BT),(Tf-tbb)*tbb>=0,(Tf-Tbb)*Tbb>=0, vb^2<=1, (zTbb'-[0;0;-1]')*(zTbb-[0;0;-1])<=0,0 == mom(diff(gbisb,tbb)+diff(gbisb,zb)*fbisb) - mom(gTbb) + mom(g0bisb)); 
[statusquob,objbisb] = msol(Qb);
toc 
%% Trajectories obtained with SDP \tilde{Z_k}

bvb = mmon([tbb;vb],d/2);
Gvb = double(mom(bvb*bvb'));

TT = linspace(0,obj,1e3)';
[XXLvb,pvb] = momgraph(Gvb,TT);

Yinp = [0;0;1];

%solution with symmetry
for i = 1:length(XXLvb)-1
    Yinp(:,i+1) = Yinp(:,i)+dt*[-k*cos(alpha)*Yinp(2,i);
                k*cos(alpha)*Yinp(1,i)-XXLvb(i)*k*sin(alpha)*Yinp(3,i);
                XXLvb(i)*k*sin(alpha)*Yinp(2,i)];
end


%% PLOTS



figure(1)
 plot(TT,XXL1.^2,'-k','linewidth',3);
 hold on
  plot(TTT,y1(1,:).^2,'--','linewidth',2);
  hold on
plot(TTT,y3(1,:).^2,'--','linewidth',2);
xlabel('Time t','Fontsize',12);
ylabel('Trajectories','Fontsize',12);
legend('1st component reconstruction','Squared 1st component of optimal trajectory','Squared 1st component of optimal trajectory','Fontsize',12)


 figure(2)
 plot(TT,XXL2.^2,'-k','linewidth',3);
 hold on
  plot(TTT,y1(2,:).^2,'--','linewidth',2);
  hold on
plot(TTT,y3(2,:).^2,'--','linewidth',2);
xlabel('Time t','Fontsize',12);
ylabel('Trajectories','Fontsize',12);
legend('2nd component reconstruction','Squared 2nd component of optimal trajectory','Squared 2nd component of optimal trajectory','Fontsize',12)


 figure(3)
 plot(TT,XXL3,'-k','linewidth',3);
 hold on
 plot(TTT,y1(3,:),'--','linewidth',2);
 legend('Optimal trajectories')
 hold on
plot(TTT,y3(3,:),'--','linewidth',2);
xlabel('Time t','Fontsize',12);
ylabel('Trajectories','Fontsize',12);
legend('3rd component reconstruction','Optimal trajectory','Optimal trajectory','Fontsize',12)
hold on
 

figure(4)
plot(TT,XXLvb,'-k','linewidth',3);
hold on
plot(TTT,U1,'--','linewidth',3.5);
 hold on
plot(TTT,U3,'--','linewidth',3.5);
 hold on
legend('Rebuilt optimal control','u1(t)','u2(t)','Fontsize',12)
xlabel('Time t','Fontsize',12);
ylabel('Optimal controls','Fontsize',12);
hold on

figure(5)
plot3(Yinp(1,:),Yinp(2,:),Yinp(3,:),'-k','linewidth',3);
hold on
p1=plot3(y1(1,:),y1(2,:),y1(3,:),'--','linewidth',3);
hold on
p2=plot3(y2(1,:),y2(2,:),y2(3,:),'--','linewidth',3);
hold on
p3=plot3(y3(1,:),y3(2,:),y3(3,:),'--','linewidth',3);
hold on
p4=plot3(y4(1,:),y4(2,:),y4(3,:),'--','linewidth',3);
hold on
[X,Y,Z] = sphere;
plot3(0, 0 ,-1,'r*','linewidth',7)
plot3(0,0,+1,'r*','linewidth',7)
surf(X,Y,Z)
axis equal
shading interp 
xlabel('x1','Fontsize',14);
ylabel('x2','Fontsize',14);
zlabel('x3','Fontsize',14);
legend('Rebuilt trajectory x(t)','x^1(t)','x^2(t)','x^3(t)','x^4(t)','Fontsize',14)
hold on
 

