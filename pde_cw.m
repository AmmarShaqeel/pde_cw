clearvars
close all


D = 1.2e-7;
G = 2.4e-2*1.2;
T_inital = 20;
T_ref = 4;

r_max = 0.01;
t_max = 1000;
n_max=1000; 
Nr=101;
Nt=101;

r=linspace(0,r_max,Nr); % Define range of r values
t=linspace(0,t_max,Nt); % Define range of t values
[R,T]=meshgrid(r,t); % Defines grids for r and t
Uk=zeros(Nt,Nr,n_max); % Initialise Tk

% Calculate Ny by Nx grid of Tk for each value of k
for n=1:n_max
%  Uk(:,:,n)= 1./R .* 2 .* (-1)^n .* r_max./(pi*n) * ((T_ref - T_inital) + r_max^2 .* G./D.*(n.*pi)^2).*sin(pi.*n.*R./r_max)*exp(-1*(pi.*n./R)^2 .* D .* T) + T_ref + r_max^2.*G./(6.*D) - R^2*G/(6*D);
Uk(:,:,n)= 2.*r_max./(pi.*R) .* (((-1)^n)/n .* sin(n*pi.*R/r_max) .* exp(-(n^2.*pi^2)*D.*T/r_max^2))*((T_ref-T_inital)+ r_max^2 * G/(D*n^2*pi^2));
end

% Sum k along the 3rd dimension to obtain k
U=sum(Uk,3) + T_ref + r_max^2.*G /(6.*D) - R.^2 .* G /(6.*D);
figure(1);
mesh(R,T,U);
xlabel('radius of apple (m)');
ylabel('Time (s)');
zlabel('temperature (c)');
title('Graph when G = 2.4e-3');
colorbar;

G = 0;

for n=1:n_max
%  Uk(:,:,n)= 1./R .* 2 .* (-1)^n .* r_max./(pi*n) * ((T_ref - T_inital) + r_max^2 .* G./D.*(n.*pi)^2).*sin(pi.*n.*R./r_max)*exp(-1*(pi.*n./R)^2 .* D .* T) + T_ref + r_max^2.*G./(6.*D) - R^2*G/(6*D);
Uk(:,:,n)= 2.*r_max./(pi.*R) .* (((-1)^n)/n .* sin(n*pi.*R/r_max) .* exp(-(n^2.*pi^2)*D.*T/r_max^2))*((T_ref-T_inital)+ r_max^2 * G/(D*n^2*pi^2));
end

% Sum k along the 3rd dimension to obtain k
U=sum(Uk,3) + T_ref + r_max^2.*G/(6.*D) - R.^2 .* G/(6.*D);
figure(2);
mesh(R,T,U);
xlabel('Radius of apple (m)');
ylabel('Time (s)');
zlabel('Temperature (c)');
title('Graph when G =  0');

colorbar;



