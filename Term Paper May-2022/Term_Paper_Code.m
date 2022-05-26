clc;
clear;

%% Import data
nodes = readmatrix('nodesreduced.xlsx');
triangles = readmatrix('trianglesreduced.xlsx');
edge = readmatrix('edgesreduced.xlsx');
%% Parameter assumptions
C = 1;
Cg = 1;
T = sqrt(L/1.56);
w = 2*pi/T;
k = 2*pi/L;
L = 1;
H = 0.04;
ang = 120;
K = zeros(length(nodes),length(nodes));
F = zeros(length(nodes),1);
%% K and F matri formulation
tic
syms x y
for i = 1:length(triangles(:,1))
    x1 = nodes(triangles(i,2),2);
    x2 = nodes(triangles(i,3),2);
    x3 = nodes(triangles(i,4),2);
    y1 = nodes(triangles(i,2),3);
    y2 = nodes(triangles(i,3),3);
    y3 = nodes(triangles(i,4),3);
    coord = [x1 x2 x3;y1 y2 y3];  
    [temp, order] = sort(coord(2,:));
    coord = coord(:,order);
    x1 = coord(1,1);
    x2 = coord(1,2);
    x3 = coord(1,3);
    y1 = coord(2,1);
    y2 = coord(2,2);
    y3 = coord(2,3);
    ai = x2*y3-x3*y2;
    bi = y2-y3;
    ci = x3-x2;
    aj = x3*y1-y3*x1;
    bj = y3-y1;
    cj = x1-x3;
    ak = x1*y2-x2*y1;
    bk = y1-y2;
    ck = x2-x1;
%     if (x3-x1)==0
%         xc = x1;
%         yc = y2;
%     elseif (y3-y1)==0
%         xc = x2;
%         yc = y1;
%     else
%         xc1 = (y1 - y2 - ((y3-y1/x3-x1))*x1 + ((x1-x3)/(y3-y1))*x2);
%         xc2 = ((x1-x3)/(y3-y1))-((y3-y1)/(x3-x1));
%         xc = xc1/xc2;
%         yc = y1 + ((y3-y1)/(x3-x1))*(xc-x1) ;
%     end
    dy11 = y1;
    dy12 = y2;
    dx11 = x1+((x2-x1)/(y2-y1))*(y-y1);
    dx12 = x1+((x3-x1)/(y3-y1))*(y-y1);
    dy21 = y2;
    dy22 = y3;
    dx21 = x2+((x2-x3)/(y2-y3))*(y-y2);
    dx22 = x1+((x3-x1)/(y3-y1))*(y-y1);
    A = 0.5*(det([x1 y1 1;x2 y2 1;x3 y3 1]));

    N11 = ai^2+(x^2*bi^2)+(y^2*ci^2)+(2*ai*x*bi)+(2*ai*y*ci)+(2*x*y*bi*ci);
    N12=(ai*aj)+(ai*x*bj)+(ai*y*cj)+(aj*x*bi)+(x^2*bi*bj)+(x*y*bi*cj)+(aj*y*ci)+(x*y*bj*ci)+(y^2*cj*ci);
    N13=(ai*ak)+(ai*x*bk)+(ai*y*ck)+(ak*x*bi)+(x^2*bi*bk)+(x*y*bi*ck)+(ak*y*ci)+(x*y*bk*ci)+(y^2*ck*ci);
    N22= aj^2+(x^2*bj^2)+(y^2*cj^2)+(2*aj*x*bj)+(2*aj*y*cj)+(2*x*y*bj*cj);
    N23=(aj*ak)+(aj*x*bk)+(aj*y*ck)+(ak*x*bj)+(x^2*bj*bk)+(x*y*bj*ck)+(ak*y*cj)+(x*y*bk*cj)+(y^2*ck*cj);
    N33= ak^2+(x^2*bk^2)+(y^2*ck^2)+(2*ak*x*bk)+(2*ak*y*ck)+(2*x*y*bk*ck);
    N_new = [N11 N12 N13 N22 N23 N33];
    F1 = int(N_new,x,dx11,dx12);
    F1 = int(F1,y,dy11,dy12);
    F2 = int(N_new,x,dx21,dx22);
    F2 = int(F2,y,dy21,dy22);
    F1 = double(F1);
    F2 = double(F2);
    F_total = double(F1)+double(F2);
    F_total_new = F_total./(4*A^2);
    K12 = C*Cg*k^2*(F_total_new);
    K12 = [K12(1,1) K12(1,2) K12(1,3);K12(1,2) K12(1,4) K12(1,5);K12(1,3) K12(1,5) K12(1,6)];
    delN = (1/(2*A))*[bi+ci bj+cj bk+ck];
    K11 = -(C*Cg)* (delN'*delN);
    K1 = K11 + K12;
    if(K1(1,1)<0)
        K1 = K1*(-1);
    end
    for j = 1:3
        K(triangles(i,j+1),triangles(i,2)) = K(triangles(i,j+1),triangles(i,2)) + K1(j,1);
        K(triangles(i,j+1),triangles(i,3)) = K(triangles(i,j+1),triangles(i,3)) + K1(j,2);
        K(triangles(i,j+1),triangles(i,4)) = K(triangles(i,j+1),triangles(i,4)) + K1(j,3);        
    end
    
end
toc
%% fluid domain boundary 
tic
syms theta
for  i = 1:length(edge(:,1))
    x1 = nodes(edge(i,1),2);
    x2 = nodes(edge(i,2),2);
    x3 = nodes(edge(i,3),2);
    y1 = nodes(edge(i,1),3);
    y2 = nodes(edge(i,2),3);
    y3 = nodes(edge(i,3),3);
    % C = sqrt((9.81/k)*tanh(k*((x1+x2)/2)));
    % Cg = 0.5*(1+((2*k*((x1+x2)/2))/(sinh(2*k*((x1+x2)/2)))));
    ai = x2*y3-x3*y2;
    bi = y2-y3;
    ci = x3-x2;
    aj = x3*y1-y3*x1;
    bj = y3-y1;
    cj = x1-x3;
    A = 0.5*abs(det([x1 y1 1;x2 y2 1;x3 y3 1]));
    r1 = sqrt(x1^2+y1^2);
    r2 = sqrt(x2^2+y2^2);
    angle1 = atan(y1/x1);
    angle2 = atan(y2/x2);
    Ni = ai + bi*r1*cos(theta) + ci*r1*sin(theta);
    Nj = aj + bj*r2*cos(theta) + cj*r2*sin(theta);
    N = [Ni Nj];
    K4 = sqrt(r1*r2)*1i*C*Cg*k*int(N'*N,theta,angle1,angle2);
    K4 = double(K4);
    K4 = K4./(4*A^2);
    F_without_phi = 2*K4;
    phi_incident = [H*0.5*exp(1i*k*r1*cosd(ang)); H*0.5*exp(1i*k*r2*cosd(ang))];
    F_final = imag(F_without_phi*phi_incident);
    K(edge(i,1),edge(i,1)) = K(edge(i,1),edge(i,1)) + K4(1,1);
    K(edge(i,1),edge(i,2)) = K(edge(i,1),edge(i,2)) + K4(1,2);
    K(edge(i,2),edge(i,1)) = K(edge(i,2),edge(i,1)) + K4(2,1);
    K(edge(i,2),edge(i,2)) = K(edge(i,2),edge(i,2)) + K4(2,2);
    F(edge(i,1),1) = F(edge(i,1),1) + F_final(1,1);
    F(edge(i,2),1) = F(edge(i,2),1) + F_final(2,1);
end
toc
%% Shore Data import
edge_shore = readmatrix('edgesshore.xlsx');
K(edge_shore(:,1)',:) = [];
K(:,edge_shore(:,1)') = [];
F(edge_shore(:,1)',:) = [];
%% Solve phi and eta
phi = K\F;
phi_final = phi';
for i=1:length(edge_shore(:,1))
    phi_final = [phi_final(1:(edge_shore(i,1)-1)) 0 phi_final(edge_shore(i,1):end)];
end
phi_final = phi_final';
eta = 1i*(w/9.81)*phi_final;
%% eta conversions
etamagnitude = abs(eta);
eta_phase_radian = angle(eta);
%% eta plots
X = nodes(:,2);
Y = nodes(:,3);
% [x,y]=meshgrid(X,Y);
% Z1 = zeros(length(nodes),length(nodes));
% Z2 = zeros(length(nodes),length(nodes));
% for j=1:length(Z1(:,1))
%     Z1(j,j) = etamagnitude(j,1);
% end
% for j=1:length(Z2(:,1))
%     Z2(j,j) = eta_phase_radian(j,1);
% end
% contour(x,y,cos(Z2).*Z1)

%% surface plot for wave
Z = etamagnitude.*cos(eta_phase_radian.*(180/pi));
% Z = reshape(Z,[length(nodes),length(nodes)]);
j=1;
x = zeros(length(triangles(:,1)),3);
y = zeros(length(triangles(:,1)),3);
z = zeros(length(triangles(:,1)),3);
for i = 1:length(triangles(:,1))
    for j=1:3
       x(triangles(i,j+1),j) = nodes(triangles(i,j+1),2);
       y(triangles(i,j+1),j) = nodes(triangles(i,j+1),3);
       z(triangles(i,j+1),j) = Z(triangles(i,j+1),1);
    end
end

surf(x,y,z)
