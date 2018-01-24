clearvars
close all


D = 1.2e-7;
G = 2.4e-3;
T_inital = 20;
T_ref = 4;

r_max = 0.01;
t_max = 100;
n_max=100; % Define maximum value of k
Nr=101; % Define number of x values between 0 and L
Nt=101; % Define number of y values between 0 and L

r=linspace(0,r_max,Nr); % Define range of x values
t=linspace(0,t_max,Nt); % Define range of y values
[R,T]=meshgrid(r,t); % Defines grids for x and y
Uk=zeros(Nt,Nr,n_max); % Initialise Tk

% Calculate Ny by Nx grid of Tk for each value of k
for n=1:n_max
 Uk(:,:,n)=2.*r_max./(pi.*R) .* (((-1)^n)/n .* sin(n*pi.*R/r_max) .* exp((n^2.*pi^2)*D.*T/r_max^2)) + T_ref + r_max^2.*G/(6.*D) - R.^2 .* G/(6.*D);
end

% Sum Tk along the 3rd dimension to obtain k
U=sum(Uk,3);
figure
mesh(R,U,T)
xlabel('radius of apple (m)')
ylabel('Time (s)')
zlabel('temperature (c)')
% axis([0 L 0 L 0 1.5*A])
% colorbar
