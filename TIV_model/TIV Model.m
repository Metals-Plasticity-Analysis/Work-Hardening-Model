%Tensile flow and work hardening behaviour of type 316L(N) austenitic stainless steel 
%in the framework of one-internal-variableand two-internal-variable approaches

clc;
clear;

data = load('example_1.txt');
exp_strain = data(:,1);
exp_stress = data(:,2);

%Input Material parameters
%Al
b=2.86E-10; %Burgers vector
G=27000;    %Shear modulus
orientationfactor=3.1;      %Taylor factor
alpha=0.3;

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
Rhom=zeros(1,(total_steps+1));   %intialize mobile dislocation density
meanfreepath=zeros(1,(total_steps+1)); %intialize mean free path

%set initial counter
dynamic=1;

%Give an intial guess of K1 , K2, intial dislocation density and yield
%strength

 kl=30;          %Initial guess of K1
 k2=15;           %Initial guess of K2 
 L_S=1E-6;        %saturation value of mean free path;
 
 Rho(dynamic)=1E11;                   %intial forest dislocation density
 Rhom(dynamic)=1E11;                  %intial mobile dislocation density
 meanfreepath(dynamic)=10E-6;         %Intial mean free path
 yield=333;                           %yield stength of the material
 stress(dynamic)=yield+orientationfactor*alpha*G*b*sqrt(Rho(dynamic));
 
 
while dynamic <=total_steps
    
              y_1= -kl*(meanfreepath(dynamic)-L_S)*h;
              y_2= -kl*((meanfreepath(dynamic)+0.5*y_1)-L_S)*h;
              y_3= -kl*((meanfreepath(dynamic)+0.5*y_2)-L_S)*h;
              y_4= -kl*((meanfreepath(dynamic)+y_3)-L_S)*h;
    
              meanfreepath(dynamic+1)= meanfreepath(dynamic)+((y_1+2*y_2+2*y_3+y_4)/6);
              
              z_1=(orientationfactor/b)*((1/L_S)-(1/(meanfreepath(dynamic))))*h;
              z_2=(orientationfactor/b)*((1/L_S)-(1/(meanfreepath(dynamic))))*h;
              z_3=(orientationfactor/b)*((1/L_S)-(1/(meanfreepath(dynamic))))*h;
              z_4=(orientationfactor/b)*((1/L_S)-(1/(meanfreepath(dynamic))))*h;
              
             Rhom(dynamic+1)= Rhom(dynamic)+((z_1+2*z_2+2*z_3+z_4)/6);
              
             m_1= orientationfactor*((1/(b*meanfreepath(dynamic)))-k2*Rho(dynamic))*h;
             m_2= orientationfactor*((1/(b*meanfreepath(dynamic)))-k2*(Rho(dynamic)+0.5*m_1))*h;
             m_3= orientationfactor*((1/(b*meanfreepath(dynamic)))-k2*(Rho(dynamic)+0.5*m_2))*h;
             m_4= orientationfactor*((1/(b*meanfreepath(dynamic)))-k2*(Rho(dynamic)+m_3))*h;
             
              Rho(dynamic+1)= Rho(dynamic)+((m_1+2*m_2+2*m_3+m_4)/6);
              strain(dynamic+1)=strain(dynamic)+h;
              stress(dynamic+1)=yield+orientationfactor*alpha*G*b*sqrt(Rho(dynamic+1)); 
             dynamic=dynamic+1;
end

 %plot forest dislocation density with plastic strain
plot(strain,Rho,'linewidth',2)
hold on
plot(strain,Rhom,'linewidth',2)
title('Forest and mobile dislocation density')
xlabel('True Strain')
ylabel('Dislocation density (m^-^2)')
set(gca, 'YScale', 'log')
set(gca,'FontSize',20,'fontweight','bold')
set(gcf,'color','w');
set(gca,'linewidth',2);
legend('Forest','Mobile','Location','Southeast','Orientation','vertical','fontsize',20);
figure ()

plot(strain,meanfreepath,'linewidth',2)
title('Mean free path')
xlabel('True Strain')
ylabel('Mean free path (micro-metre)')
set(gca,'FontSize',20,'fontweight','bold')
set(gcf,'color','w');
set(gca,'linewidth',2);
legend('Mean free path','Location','Southeast','Orientation','vertical','fontsize',20);
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
 
