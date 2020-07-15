%LOCAL STRAIN HARDENING AND NONUNIFORMITY
%OF PLASTIC DEFORMATION 

clc;
clear;

data = load('example_1.txt');
exp_strain = data(:,1);
exp_stress = data(:,2);

%Input Material parameters
%Al
b=2.86E-10; %Burgers vector
b2=b^2;
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

%set initial counter
dynamic=1;

%Give an intial guess of C1 , C2,C3, C4 intial dislocation density and yield
%strength
C1=5.5E-5;
C2=5;
C3=0.04;
C4=15;

 Rho(dynamic)=1E11;                   %intial forest dislocation density
 Rhom(dynamic)=1E11;                  %intial mobile dislocation density
 yield=333;                           %yield stength of the material
 stress(dynamic)=yield+orientationfactor*alpha*G*b*sqrt(Rho(dynamic));
 
 while dynamic <=total_steps
     
     y1=orientationfactor*(((C1/b2)*( Rho(dynamic)/Rhom(dynamic)))-(C2*Rhom(dynamic))-((C3/b)*sqrt(Rho(dynamic))))*h;
     y2=orientationfactor*(((C1/b2)*( Rho(dynamic)/(Rhom(dynamic)+0.5*y1)))-(C2*(Rhom(dynamic)+0.5*y1))-((C3/b)*sqrt(Rho(dynamic))))*h;
     y3=orientationfactor*(((C1/b2)*( Rho(dynamic)/(Rhom(dynamic)+0.5*y2)))-(C2*(Rhom(dynamic)+0.5*y2))-((C3/b)*sqrt(Rho(dynamic))))*h;
     y4=orientationfactor*(((C1/b2)*( Rho(dynamic)/(Rhom(dynamic)+y3)))-(C2*(Rhom(dynamic)+y3))-((C3/b)*sqrt(Rho(dynamic))))*h;
     
      Rhom(dynamic+1)= Rhom(dynamic)+((y1+2*y2+2*y3+y4)/6);
     
      z1=orientationfactor*((C3/b)*(sqrt(Rho(dynamic)))-(C4*Rho(dynamic))+(C2*Rhom(dynamic+1)))*h;
      z2=orientationfactor*((C3/b)*(sqrt(Rho(dynamic)+0.5*z1))-(C4*(Rho(dynamic)+0.5*z1))+(C2*Rhom(dynamic+1)))*h;
      z3=orientationfactor*((C3/b)*(sqrt(Rho(dynamic)+0.5*z2))-(C4*(Rho(dynamic)+0.5*z2))+(C2*Rhom(dynamic+1)))*h;
      z4=orientationfactor*((C3/b)*(sqrt(Rho(dynamic)+z3))-(C4*(Rho(dynamic)+z3))+(C2*Rhom(dynamic+1)))*h;
      
       Rho(dynamic+1)= Rho(dynamic)+((z1+2*z2+2*z3+z4)/6);
       strain(dynamic+1)=strain(dynamic)+h;
       stress(dynamic+1)=yield+orientationfactor*alpha*G*b*sqrt(Rho(dynamic+1)); 
       dynamic=dynamic+1;
 end
 
 %plot forest and mobile dislocation density with plastic strain
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


