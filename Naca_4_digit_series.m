% compute NACA 4 Digit Airfoil
% Written by: Ayaz Shariff
% Started: 28/10/2021
% End Date: 28/10/2021

clear all;
clc;

%% User Inputs


% Type Of Airfoil
typeNACA = '2412';

% Extract Values From Type
Minit = str2double(typeNACA(1));
Pinit = str2double(typeNACA(2));
Tinit = str2double (typeNACA(3:4));

% Number of Grid Points
gridpts = 500;

% Constants
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = - 0.1015;   % open trailing edge
% a4 = - 0.1036; % closed trailing edge

%% Calculations

% Actual Percentage Values Of Airfoil Properties
M = Minit/100;
P = Pinit/10;
T = Tinit/100;

% Airfoil Grid

x = linspace(0,1,gridpts)'; 

%% Camber and Gradient
yc = ones (gridpts,1);
dyc_dx = ones (gridpts,1);
theta = ones (gridpts,1);
for i = 1:1:gridpts
    if (x(i) >=0 && x(i) < P)
        yc(i)     = (M/P^2)*((2*P*x(i))-(x(i)^2));
        dyc_dx(i) =  ((2*M)/(P^2))*(P-x(i));
    elseif (x(i) >= P && x(i)<=1)
        yc(i)     = (M/(1-P)^2)*(1-(2*P) + (2*P*x(i))-(x(i)^2));
        dyc_dx(i) = ((2*M)/((1-P)^2))*(P-x(i));
    end
    theta(i) = atan(dyc_dx(i));
end

%% Thickness Distributions
yt = ones (gridpts,1);

for i = 1:1:gridpts
    term0 = a0*sqrt(x(i));
    term1 = a1*x(i);
    term2 = a2*x(i)^2;
    term3 = a3*x(i)^3;
    term4 = a4*x(i)^4;
    
    yt(i) = 5*T*(term0+term1+term2+term3+term4);
end

%% Upper Surface Points
xu = ones(gridpts,1);
yu = ones(gridpts,1);

for i = 1:1:gridpts
    xu (i) = x(i) - yt(i)*sin(theta(i));
    yu (i) = yc(i) + yt(i)*cos(theta(i));
end

%% Lower Surface Points
xl = ones(gridpts,1);
yl = ones(gridpts,1);

for i = 1:1:gridpts
    xl(i) = x(i) + yt(i)*sin(theta(i));
    yl(i) = yc(i) - yt(i)*cos(theta(i));
end

%% Plot The Airfoil (With Lines)
f1 = figure(1);
hold on; grid on;
axis equal;
plot(xu,yu, 'r-');
plot(xl,yl, 'k-');
    


