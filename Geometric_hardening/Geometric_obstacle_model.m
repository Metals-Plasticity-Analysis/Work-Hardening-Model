%This script is for Geomtric obstacle model . The script is
%based on Kocks-Mecking equation dp/de=M(K-K2*p). Here, 'p' is
% forest dislocation density. 'M' is Taylor factor. K2 is related to dynamic recovery. K is
% related to dislocation storage due to dispersoids/precipitates/grain boundaries. 

clc;
clear;

%Provide the path of the experimental data. Experimental data should be
%true plastic strain versus true stress in a 'txt' file.

data = load('example_1.txt');
exp_strain = data(:,1);
exp_stress = data(:,2);

%Input Material parameters
%Al
b=2.86E-10; %Burgers vector
G=27000;    %Shear modulus
orientationfactor=2.3; %Taylor factor
alpha=0.3; %alpha parameter

%Specify the total plastic strain till which simulations should be
%performed

total_plastic_strain=0.0745;

%define the strain increment for performing the simulations
h=0.0005; %a small value should always be used for accurate numerical integartion to be performed later

%Calculate the total number of steps in the calculation
total_steps=total_plastic_strain/h;

%Intialize variables for storing the output variables such as dislocation
%density and flow stress

Rho=zeros(1,(total_steps+1)); %initialize dislocation density
stress=zeros(1,(total_steps+1));%initialize flow stress
strain=zeros(1,(total_steps+1)); %intialize strain

%set initial counter
dynamic=1;

%Give an intial guess of K1 , K2, intial dislocation density and yield
%strength

                    
 k_2=2.0;             %Initial guess of K2 %1.0;
 k=8.69E16;            %Initial guess of K %8.69E14; 
 %Keep changing K1 and K2 to get good match between experimental and
 %simulated stress-strain curves
 
 Rho(dynamic)=1E11; %intial dislocation density
 yield=333;          %yield stength of the material
 stress(dynamic)=yield+orientationfactor*alpha*G*b*sqrt(Rho(dynamic));
 
while dynamic <=total_steps
    
              y_1= (k-k_2*Rho(dynamic))*orientationfactor*h;
              
              y_2= (k-k_2*(Rho(dynamic)+0.5*y_1))*orientationfactor*h;
              
              y_3= (k-k_2*(Rho(dynamic)+0.5*y_2))*orientationfactor*h;
              
              y_4= (k-k_2*(Rho(dynamic)+y_3))*orientationfactor*h;
              
              Rho(dynamic+1)= Rho(dynamic)+((y_1+2*y_2+2*y_3+y_4)/6);
              
              strain(dynamic+1)=strain(dynamic)+h;
              stress(dynamic+1)=yield+orientationfactor*alpha*G*b*sqrt(Rho(dynamic+1));
              dynamic=dynamic+1;
end

%plot forest dislocation density with plastic strain
plot(strain,Rho,'linewidth',2)
title('Forest dislocation density')
xlabel('True Strain')
ylabel('Dislocation density (m-2)')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20,'fontweight','bold')
set(gcf,'color','w');
set(gca,'linewidth',2);

figure ()

%compare simulated flow stress with plastic strain with experimental flow
%stress data
plot(strain,stress,'linewidth',2) %simulated
hold on
plot(exp_strain,exp_stress,'linewidth',2,'LineStyle','--') %experimental
title('Comparison of experimental and simulated flow curves')
xlabel('True Strain')
ylabel('True Stress (MPa)')
set(gca,'FontSize',20,'fontweight','bold')
set(gcf,'color','w');
set(gca,'linewidth',2);
legend('Simulated','Experimental','Location','Southeast','Orientation','vertical','fontsize',20);



