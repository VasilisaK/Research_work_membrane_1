function membrane_solution
clear, clc
PI = 3.14159265;

% experimental data
M = [0 0 1;
0 -10 0.8;
0 -20 0.3;
0 -30 0.1;
0 -40 0.0;
sqrt(50) -sqrt(50) 1;
sqrt(200) -sqrt(200) 0.6;
sqrt(450) -sqrt(450) 0.1;
sqrt(800) -sqrt(800) -0.1;
20 0 0.7;
30 0 0.3;
40 0 0.1;
sqrt(50) sqrt(50) 0.9;
sqrt(200) sqrt(200) 0.5;
sqrt(450) sqrt(450) 0.0;
sqrt(800) sqrt(800) -0.2;
0 10 0.6;
0 20 0.3;
0 30 0.1;
0 40 -0.1;
-10 0 0.6;
-20 0 0.3;
-30 0 0.1;
-40 0 0.0;
10 0 1.3]; 

x = M(1:16,1);
y = M(1:16,2);
hexp = M(1:16,3);

r = sqrt(x.^2 + y.^2);
Fi = atan2(y,x)*180/PI;
Fi_rad = Fi.*PI./180;

% polar coordinates
[rho, Teta] = calc_cyl(r,Fi_rad);

% SLAE, model coefficients
[coeffs] = solve_slau_m1(rho, Teta, hexp);

A1 = coeffs(1)
B1 = coeffs(2)
D0 = coeffs(3)

% solution plot
draw_solution_field_m1_2(A1,B1,D0);

% SLAE, model coefficients
[coeffs] = solve_slau_m2(rho, Teta, hexp);

A1 = coeffs(1)
A2 = coeffs(2)
B1 = coeffs(3)
B2 = coeffs(4)
D0 = coeffs(5)

% solution plot
draw_solution_field_m2_2(A1,A2,B1,B2,D0);


end

% conformal mapping, new polar coordinates
function [rho, Teta] = calc_cyl(r,Fi_rad)

r0 = 3.95; % domain under cargo radius
Rm = 50; % membrane radius
a = 10; % distance between membrane center and cargo center

Qplus = Rm^2 + a^2 - r0^2 + sqrt( (Rm^2 + a^2 - r0^2)^2 - 4*a^2*Rm^2 );
Qminus = Rm^2 + a^2 - r0^2 - sqrt( (Rm^2 + a^2 - r0^2)^2 - 4*a^2*Rm^2 );
k = Qplus/(2*a*Rm);
rho_up = 4*a^2*r.^2 + 4*a*Qminus*r.*cos(Fi_rad) + Qminus^2;
rho_down = 4*a^2*r.^2 + 4*a*Qplus*r.*cos(Fi_rad) + Qplus^2;
rho = k*sqrt( rho_up./rho_down );
Teta_up = r*sqrt( ( Rm^2+a^2-r0^2 )^2 - 4*Rm^2*a^2 ).*sin(Fi_rad);
Teta_down = a*(Rm^2+r.^2) + r*(Rm^2+a^2-r0^2).*cos(Fi_rad);
Teta = atan2( Teta_up,Teta_down );

end

function [rho, Teta] = calc_cyl_2d(r,Fi_rad)

r0 = 3.95; % domain under cargo radius
Rm = 50; % membrane radius
a = 10; % distance between membrane center and cargo center

Qplus = Rm^2 + a^2 - r0^2 + sqrt( (Rm^2 + a^2 - r0^2)^2 - 4*a^2*Rm^2 );
Qminus = Rm^2 + a^2 - r0^2 - sqrt( (Rm^2 + a^2 - r0^2)^2 - 4*a^2*Rm^2 );
k = Qplus/(2*a*Rm);

rho = zeros(size(r));
Teta = zeros(size(Fi_rad));

for i = 1:size(rho,1)
    for j = 1:size(Teta,2)
        rho_up = 4*a^2*r(i,j)*r(i,j) + 4*a*Qminus*r(i,j)*cos(Fi_rad(i,j)) + Qminus^2;
        rho_down = 4*a^2*r(i,j)*r(i,j) + 4*a*Qplus*r(i,j)*cos(Fi_rad(i,j)) + Qplus^2;
        rho(i,j) = k*sqrt( rho_up/rho_down );
        Teta_up = r(i,j)*sqrt( ( Rm^2+a^2-r0^2 )^2 - 4*Rm^2*a^2 )*sin(Fi_rad(i,j));
        Teta_down = a*(Rm^2+r(i,j)*r(i,j)) + r(i,j)*(Rm^2+a^2-r0^2)*cos(Fi_rad(i,j));
        Teta(i,j) = atan2( Teta_up,Teta_down );
    end
end


end


function [coeffs] = solve_slau_m1(rho, Teta, hexp)

% create system matrix and right part 
M = [0 0 0; 0 0 0; 0 0 0];
B = [0; 0; 0;];

for i=1:size(rho)

        M(1, 1) = M(1, 1) + (rho(i) - 1/rho(i))^2*cos(Teta(i))^2;
        M(1, 2) = M(1, 2) + (rho(i) - 1/rho(i))^2*cos(Teta(i))*sin(Teta(i));
        M(1, 3) = M(1, 3) - log(rho(i))*(rho(i) - 1/rho(i))*cos(Teta(i));
        M(2, 1) = M(1, 2);
        M(2, 2) = M(2, 2) + (rho(i) - 1/rho(i))^2*sin(Teta(i))^2;
        M(2, 3) = M(2, 3) - log(rho(i))*(rho(i) - 1/rho(i))*sin(Teta(i));
        M(3, 1) = M(1, 3);
        M(3, 2) = M(2, 3);
        M(3, 3) = M(3, 3) + log(rho(i))^2;
            
        B(1) = B(1) + (rho(i) - 1/rho(i))*hexp(i)*cos(Teta(i));
        B(2) = B(2) + (rho(i) - 1/rho(i))*hexp(i)*sin(Teta(i));
        B(3) = B(3) - hexp(i)*log(rho(i));

end
%disp(M);
%disp(B);

% matrix system conditional number
c = cond(M)

% SLAE solution, model coefficients
coeffs = linsolve(M,B);

end

function [coeffs] = solve_slau_m2(rho, Teta, hexp)

% create system matrix and right part 
M = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
B = [0; 0; 0; 0; 0];

for i=1:size(rho)

        M(1, 1) = M(1, 1) + (rho(i) - 1/rho(i))^2*cos(Teta(i))^2;
        M(1, 2) = M(1, 2) + (rho(i) - 1/rho(i))*(rho(i)^2 - 1/rho(i)^2)*cos(Teta(i))*cos(2*Teta(i));
        M(1, 3) = M(1, 3) + (rho(i) - 1/rho(i))^2*cos(Teta(i))*sin(Teta(i));
        M(1, 4) = M(1, 4) + (rho(i) - 1/rho(i))*(rho(i)^2 - 1/rho(i)^2)*cos(Teta(i))*sin(2*Teta(i));       
        M(1, 5) = M(1, 5) - log(rho(i))*(rho(i) - 1/rho(i))*cos(Teta(i));
        M(2, 1) = M(1, 2);
        M(2, 2) = M(2, 2) + (rho(i)^2 - 1/rho(i)^2)^2*cos(2*Teta(i))^2;
        M(2, 3) = M(2, 3) + (rho(i) - 1/rho(i))*(rho(i)^2 - 1/rho(i)^2)*cos(2*Teta(i))*sin(Teta(i));
        M(2, 4) = M(2, 4) + (rho(i)^2 - 1/rho(i)^2)^2*cos(2*Teta(i))*sin(2*Teta(i));
        M(2, 5) = M(2, 5) - log(rho(i))*(rho(i)^2 - 1/rho(i)^2)*cos(2*Teta(i));
        M(3, 1) = M(1, 3);
        M(3, 2) = M(2, 3);
        M(3, 3) = M(3, 3) + (rho(i) - 1/rho(i))^2*sin(Teta(i))^2;
        M(3, 4) = M(3, 4) + (rho(i) - 1/rho(i))*(rho(i)^2 - 1/rho(i)^2)*sin(Teta(i))*sin(2*Teta(i));
        M(3, 5) = M(3, 5) + - log(rho(i))*(rho(i) - 1/rho(i))*sin(Teta(i));
        M(4, 1) = M(1, 4);
        M(4, 2) = M(2, 4);
        M(4, 3) = M(3, 4);
        M(4, 4) = M(4, 4) + (rho(i)^2 - 1/rho(i)^2)^2*sin(2*Teta(i))^2;
        M(4, 5) = M(4, 5) - log(rho(i))*(rho(i)^2 - 1/rho(i)^2)*sin(2*Teta(i));
        M(5, 1) = M(1, 5);
        M(5, 2) = M(2, 5);
        M(5, 3) = M(3, 5);
        M(5, 4) = M(4, 5);
        M(5, 5) = M(5, 5) + log(rho(i))^2;
            
        B(1) = B(1) + (rho(i) - 1/rho(i))*hexp(i)*cos(Teta(i));
        B(2) = B(2) + (rho(i)^2 - 1/rho(i)^2)*hexp(i)*cos(2*Teta(i));
        B(3) = B(3) + (rho(i) - 1/rho(i))*hexp(i)*sin(Teta(i));
        B(4) = B(4) + (rho(i)^2 - 1/rho(i)^2)*hexp(i)*sin(2*Teta(i));
        B(5) = B(5) - hexp(i)*log(rho(i));

end

%disp(M);
%disp(B);

% matrix system conditional number
c = cond(M)

% SLAE solution, model coefficients
format short
coeffs = linsolve(M,B);

end


function draw_solution_field_m1_2(A1,B1,D0)

PI = 3.14159265;
Nx = 101;
Ny = 101;
xn = linspace(-50,50,Nx);
yn_t = linspace(-50,50,Ny);
yn = yn_t';
coords2d_r = zeros(size(xn,2),size(yn,1));
coords2d_fi = zeros(size(xn,2),size(yn,1));

for i = 1:size(xn,2)
    for j = 1:size(yn,1)
        coords2d_fi(i,j) = atan2(yn(j,1),xn(1,i))*PI/180;
    end
end
    
for i = 1:size(xn,2)
    for j = 1:size(yn,1)
        coords2d_r(i,j) = sqrt(xn(1,i)^2 + yn(j,1)^2);
    end
end

[rho2d, Teta2d] = calc_cyl_2d(coords2d_r,coords2d_fi);
h2d = zeros(size(coords2d_r));
for i = 1:size(rho2d,1)
    for j = 1:size(Teta2d,2)
        h2d(i,j) = (rho2d(i,j) - 1/rho2d(i,j))*(A1*cos(Teta2d(i,j)) + B1*sin(Teta2d(i,j))) + D0*log(1/rho2d(i,j));
    end
end


format long

Mexp = [0 0 1;
0 -10 0.8;
0 -20 0.3;
0 -30 0.1;
0 -40 0.0;
sqrt(50) -sqrt(50) 1;
sqrt(200) -sqrt(200) 0.6;
sqrt(450) -sqrt(450) 0.1;
sqrt(800) -sqrt(800) -0.1;
20 0 0.7;
30 0 0.3;
40 0 0.1;
sqrt(50) sqrt(50) 0.9;
sqrt(200) sqrt(200) 0.5;
sqrt(450) sqrt(450) 0.0;
sqrt(800) sqrt(800) -0.2;
0 10 0.6;
0 20 0.3;
0 30 0.1;
0 40 -0.1;
-10 0 0.6;
-20 0 0.3;
-30 0 0.1;
-40 0 0.0;
10 0 1.3]; 

xexp = Mexp(:,1);
yexp = Mexp(:,2);
hexp = Mexp(:,3);

%solution and experimental data plot
figure('Color','w')
set(gca,'FontSize',12)
a = gradient(0:0.005:0.1);
h = surf(xn, yn, h2d, 'AlphaData',a, 'FaceAlpha',.3);
set(h,'edgecolor','r','facecolor',[1 1 1])
xlim([-50 50])
ylim([-50 50])
zlim([0 1.5])

hold on
plot3(xexp(1:24),yexp(1:24),hexp(1:24),'.k','MarkerSize',20)

xlabel('x(sm)')
ylabel('y(sm)')
zlabel('h(sm)')
legend('h(x,y)', 'experiment',1)

%calculating errors
Arr_calc = zeros(1,25);
[X, Y] = meshgrid(xn, yn);

for i = 1:25
    Arr_calc(i) = griddata(X,Y,h2d,xexp(i),yexp(i),'natural');
end

Arr_calc_t = Arr_calc';

err = Arr_calc_t - hexp;
disp(err);
err_sq = sqrt((1/size(Arr_calc_t(1:16,1),1))*sum((Arr_calc_t(1:16,1) - hexp(1:16,1)).*(Arr_calc_t(1:16,1) - hexp(1:16,1))))
err_sq2 = sqrt((1/size(Arr_calc_t(17:24,1),1))*sum((Arr_calc_t(17:24,1) - hexp(17:24,1)).*(Arr_calc_t(17:24,1) - hexp(17:24,1))))

%plot discrepancy learning and testing
figure('Color','w')
set(gca,'FontSize',12)
plot3(xexp(1:16),yexp(1:16),err(1:16),'.','MarkerSize',20)
grid on
set(gca,'XTick',-50:25:50)
set(gca,'YTick',-50:25:50)
xlim([-50 50])
ylim([-50 50])
zlim([-0.3 0.3])
set(gca,'ZTick',-0.3:0.1:0.3)
xlabel('x(sm)')
ylabel('y(sm)')
zlabel('d(sm)')

hold on
plot3(xexp(17:24),yexp(17:24),err(17:24),'r.','MarkerSize',20)

legend('discrepancy learning','discrepancy testing',1)

end

function draw_solution_field_m2_2(A1,A2,B1,B2,D0)

PI = 3.14159265;
Nx = 101;
Ny = 101;
xn = linspace(-50,50,Nx);
yn_t = linspace(-50,50,Ny);
yn = yn_t';
coords2d_r = zeros(size(xn,2),size(yn,1));
coords2d_fi = zeros(size(xn,2),size(yn,1));

for i = 1:size(xn,2)
    for j = 1:size(yn,1)
        coords2d_fi(i,j) = atan2(yn(j,1),xn(1,i))*PI/180;
    end
end
    
for i = 1:size(xn,2)
    for j = 1:size(yn,1)
        coords2d_r(i,j) = sqrt(xn(1,i)^2 + yn(j,1)^2);
    end
end

[rho2d, Teta2d] = calc_cyl_2d(coords2d_r,coords2d_fi);
h2d = zeros(size(coords2d_r));
for i = 1:size(rho2d,1)
    for j = 1:size(Teta2d,2)
        h2d(i,j) = (rho2d(i,j) - 1/rho2d(i,j))*(A1*cos(Teta2d(i,j)) + B1*sin(Teta2d(i,j))) + (rho2d(i,j)*rho2d(i,j) - 1/(rho2d(i,j)*rho2d(i,j)))*(A2*cos(2*Teta2d(i,j)) + B2*sin(2*Teta2d(i,j))) + D0*log(1/rho2d(i,j));
    end
end


format long

Mexp = [0 0 1;
0 -10 0.8;
0 -20 0.3;
0 -30 0.1;
0 -40 0.0;
sqrt(50) -sqrt(50) 1;
sqrt(200) -sqrt(200) 0.6;
sqrt(450) -sqrt(450) 0.1;
sqrt(800) -sqrt(800) -0.1;
20 0 0.7;
30 0 0.3;
40 0 0.1;
sqrt(50) sqrt(50) 0.9;
sqrt(200) sqrt(200) 0.5;
sqrt(450) sqrt(450) 0.0;
sqrt(800) sqrt(800) -0.2;
0 10 0.6;
0 20 0.3;
0 30 0.1;
0 40 -0.1;
-10 0 0.6;
-20 0 0.3;
-30 0 0.1;
-40 0 0.0;
10 0 1.3]; 

xexp = Mexp(:,1);
yexp = Mexp(:,2);
hexp = Mexp(:,3);

%solution and experimental data plot
figure('Color','w')
set(gca,'FontSize',12)
a = gradient(0:0.005:0.1);
h = surf(xn, yn, h2d, 'AlphaData',a, 'FaceAlpha',.3);
set(h,'edgecolor','r','facecolor',[1 1 1])
xlim([-50 50])
ylim([-50 50])
zlim([0 1.5])

hold on
plot3(xexp(1:24),yexp(1:24),hexp(1:24),'.k','MarkerSize',20)

xlabel('x(sm)')
ylabel('y(sm)')
zlabel('h(sm)')
legend('h(x,y)', 'experiment',1)

%calculating errors
Arr_calc = zeros(1,25);
[X, Y] = meshgrid(xn, yn);

for i = 1:25
    Arr_calc(i) = griddata(X,Y,h2d,xexp(i),yexp(i),'natural');
end

Arr_calc_t = Arr_calc';

err = Arr_calc_t - hexp;
disp(err);
err_sq = sqrt((1/size(Arr_calc_t(1:16,1),1))*sum((Arr_calc_t(1:16,1) - hexp(1:16,1)).*(Arr_calc_t(1:16,1) - hexp(1:16,1))))
err_sq2 = sqrt((1/size(Arr_calc_t(17:24,1),1))*sum((Arr_calc_t(17:24,1) - hexp(17:24,1)).*(Arr_calc_t(17:24,1) - hexp(17:24,1))))

%plot discrepancy learning and testing 
figure('Color','w')
set(gca,'FontSize',12)
plot3(xexp(1:16),yexp(1:16),err(1:16),'.','MarkerSize',20)
grid on
set(gca,'XTick',-50:25:50)
set(gca,'YTick',-50:25:50)
xlim([-50 50])
ylim([-50 50])
zlim([-0.3 0.3])
set(gca,'ZTick',-0.3:0.1:0.3)
xlabel('x(sm)')
ylabel('y(sm)')
zlabel('d(sm)')

hold on
plot3(xexp(17:24),yexp(17:24),err(17:24),'r.','MarkerSize',20)

legend('discrepancy learning','discrepancy testing',1)

end
