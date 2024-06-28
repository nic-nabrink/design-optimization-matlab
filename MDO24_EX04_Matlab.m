%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Truss optimization      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�//////�//////� 
% \     |     /
%  \    |    /
%   \   |   / 
%    \  |  /
%     \ | / 
%      \|/
%    F1/ \F2

% Clean up the workspace
clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%
%%  Problem Setup  %%
%%%%%%%%%%%%%%%%%%%%%

syms A1 A2;

%System parameters
F = 150000;    %load (F1=F2=F)
sigma_yield = 216;      %maximum allowed stress for steel
sigma_yield_2 = 140;    %maximum allowed stress for alloy
Rho = 7.9e-6;
l = 1000;

% objective function
m = Rho*(2*sqrt(2)*A1+A2)*l;

% constraints
s1 = F*(sqrt(2)*A1+A2)/(sqrt(2)*A1^2+2*A1*A2);
s2 = F*(sqrt(2)*A1)/(sqrt(2)*A1^2+2*A1*A2);
s3 = F*(A2)/(sqrt(2)*A1^2+2*A1*A2);

[G1,p,c]=solve(s1==sigma_yield,A2,'ReturnConditions',true); %%g1=solve(s1==sigma_yield,A2)
[G2,p,c]=solve(s2==sigma_yield_2,A2,'ReturnConditions',true); %%g2=solve(s2==sigma_yield,A2)
[G3,p,c]=solve(s3==sigma_yield,A2,'ReturnConditions',true); %%g3=solve(s3==sigma_yield,A2)

g=[G1,G2,G3]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Numerical Gradients   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem Setup
m=@(A1,A2) Rho*(2*sqrt(2)*A1+A2)*l;
g1=@(A1,A2) F*(sqrt(2)*A1+A2)/(sqrt(2)*A1^2+2*A1*A2)-sigma_yield;
g2=@(A1,A2) F*(sqrt(2)*A1)/(sqrt(2)*A1^2+2*A1*A2)-sigma_yield_2;

format long

curves(1).name = 'm';
curves(1).value = m;
curves(2).name = 'g1';
curves(2).value = g1;
curves(3).name = 'g2';
curves(3).value = g2;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Input Arguments   %%
%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO: vary the function input arguments to calculate all required 
%gradients.
%   4.3.4: Use Forward finite differences (ffd), a stepsize of 0.5 and no noise
%          on the data. Implement finite difference function from line 153. 
%   4.3.5: This time only evaluate the curve g1 and enable noise on the
%          data. How would you rate the quality of the result and what can
%          be done to improve prediction accuracy?

%Edit starts here:---------------------------------------------------------
method = 'ffd';       %Finite Difference method ('cfd','ffd','bfd')
curve = curves(2).value;       %Select curves to analyze {m,g1,g2} 
stepsize = 0.5;       %Determine stepsize
noise = 0;            %Enable (1) or disable (0) noise on the data
%Edit ends here:-----------------------------------------------------------

% Call Finite Differences Function  
if isequal(curve,curves)
    for i=1:3
        [grad_A1 grad_A2 data]=finiteDifferences(method,curve(i).value,stepsize);
        text1=['d',curve(i).name,'/','d','A1','  =  ',num2str(grad_A1)];
        text2=['d',curve(i).name,'/','d','A2','  =  ',num2str(grad_A2)];
        disp(text1)
        disp(text2)
        disp(' ')
    end
end

%%%%%%%%%%%%%%%%%%%%
%%  Contour Plot  %%
%%%%%%%%%%%%%%%%%%%%
if noise == 1
    g(1)= @(A1) -(9*2^(1/2)*A1^2 - 6250*2^(1/2)*A1)/(18*A1 - 6250)...
        + 7.5*sin(A1*0.5)+5*sin(A1);
end

figure;

fcontour(m,'LineWidth',1,'Levelstep',3)
set(gca,'fontsize',15)
colormap([0 0 0])
hold on 

for i=1:3                           %%for i=1:3
    fplot(g(i),'k','LineWidth',3);  %%fplot(g(i),'k','LineWidth',3);
end

axis([0,1000,0,1000])
xlabel('A_1,A_3 (mm^2)','Fontsize',20)
ylabel('A_2 (mm^2)','Fontsize',20)
hold on

%% Plotting of Result (for g1 only) %%
if isequal(curve,g1)
    [grad_A1 grad_A2 data]=finiteDifferences(method,curve,stepsize);
    f1=@(A1) -(9*2^(1/2)*A1^2 - 6250*2^(1/2)*A1)/(18*A1 - 6250); % 2-dimensional constraint
    f2=@(A1) -(9*2^(1/2)*A1^2 - 6250*2^(1/2)*A1)/(18*A1 - 6250)...
    + 7.5*sin(A1*0.5)+5*sin(A1);                                 % 2-dimensional constraint with noise
    
    if noise == 0                                                % Select datacurve depending on user input
        f=f1;
    else 
        f=f2;
    end
    
    X=data;                                                      % X-Coordinates of Sample points
    Y=[f(data(1)),f(data(2))];                                   % Y-Coordinates of Sample points
    plot(X,Y,'*','MarkerSize',10,...
        'MarkerEdgeColor','red')
    hold on;
    p=polyfit(X,Y,1);
    p1=@(A1) p(2)+p(1)*A1;
    fplot(p1,'r','LineWidth',1.5);
end

%%Function Declaration
function [gradx1, gradx2, sample] = finiteDifferences(method,output_func,delta_x)
    x1=513.7;
    x2=394.37;

    switch method
        %%%%%%%%%%%%%%%%%   Forward finite differences:   %%%%%%%%%%%%%%%%%
        case 'ffd'
            %TODO: Complete the function routine to calculate the Forward
            %finite differences for both variables x1,x2 and return the 
            %sample points.

            %Edit starts here:---------------------------------------------
            gradx1= (output_func(x1+delta_x,x2)-output_func(x1,x2))/delta_x; 
            gradx2= (output_func(x1,x2+delta_x)-output_func(x1,x2))/delta_x; 
            sample= [x1+delta_x,x1];                                         
            %Edit ends here:-----------------------------------------------
                       
        %%%%%%%%%%%%%%%%%   Central finite differences:   %%%%%%%%%%%%%%%%%    
        case 'cfd'
            gradx1=(output_func(x1+delta_x,x2)-output_func(x1-delta_x,x2))/(2*delta_x);
            gradx2=(output_func(x1,x2+delta_x)-output_func(x1,x2-delta_x))/(2*delta_x);
            sample=[x1-delta_x,x1+delta_x];
            
        %%%%%%%%%%%%%%%%%   Backward finite differences:  %%%%%%%%%%%%%%%%%
        case 'bfd'
            gradx1=(output_func(x1,x2)-output_func(x1-delta_x,x2))/(delta_x);
            gradx2=(output_func(x1,x2)-output_func(x1,x2-delta_x))/(delta_x);
            sample=[x1-delta_x,x1];
            
    end
end