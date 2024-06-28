close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%
%%  Problem Setup  %%
%%%%%%%%%%%%%%%%%%%%%


%System parameters
F = 150000;         %load (F1=F2=F)
sigma_yield = 216;  %maximum allowed stress for steel
Rho = 7.9e-3;       %density of material (steel)

% objective function
fun = @(x) Rho*(2*sqrt(2)*x(1)+x(2));

% TODO: Initialize the algorithm with all required input arguments. Look
% into the official documentation for support: 
% https://de.mathworks.com/help/optim/ug/fmincon.html

%Edit starts here:---------------------------------------------------------
x0=[900,900];
A = [];
b = [];
Aeq = [];
beq = [];
lb=[1 1];
ub=[1000 1000];
%Edit ends here:-----------------------------------------------------------
nonlcon=@(x)constraints(x,F,sigma_yield);
options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','OutputFcn',@TrussOutFn,'MaxFunEvals',100)

TrussPlot
% Output result   
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

% TODO: complete the constraints function, that yields all nonlinear
% constraints. Use x,F and sigma_yield as input arguments

%Edit starts here:---------------------------------------------------------
function [c,ceq] = constraints(x,F,sigma_yield)
    c(1) = F*(sqrt(2)*x(1)+x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2))-sigma_yield;
    c(2) = F*(sqrt(2)*x(1))/(sqrt(2)*x(1)^2+2*x(1)*x(2))-sigma_yield;
    c(3) = F*(x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2))-sigma_yield;
    ceq =  [];
end
%Edit ends here:-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%
%%  Contour Plot   %%
%%%%%%%%%%%%%%%%%%%%%
function [] = TrussPlot
    syms A1 A2;

    %System parameters
    F = 150000;   
    sigma_yield = 216;   
    Rho=7.9e-6;

    % objective function
    m=Rho*(2*sqrt(2)*A1+A2);


    g1=-(9*2^(1/2)*A1^2 - 6250*2^(1/2)*A1)/(18*A1 - 6250);
    g2=(3125*2^(1/2))/9 - (2^(1/2)*A1)/2;
    g3=-(9*2^(1/2)*A1^2)/(18*A1 - 6250);
    g=[g1,g2,g3]; %%g=[g1,g2,g3]

    f1=figure;
    movegui(f1,'west');
    fcontour(m,'LineWidth',1,'Levelstep',0.3e-2)
    set(gca,'fontsize',15)
    colormap([0 0 0])
    hold on 
    for i=1:3                         
        fplot(g(i),'k','LineWidth',3);  
    end
    axis([0,1000,0,1000])
    xlabel('A_1,A_3 (mm^2)','Fontsize',20)
    ylabel('A_2 (mm^2)','Fontsize',20)
end


%%%%%%%%%%%%%%%%%%%%%%
%% Convergence Plot %%
%%%%%%%%%%%%%%%%%%%%%%
function stop = TrussOutFn(x,optimValues,state)
    stop = false;
    
    persistent history
    
    history.fval = [];
    history.x = [];
    history.iter =[];
    
    switch state
        case 'init'
            f2=figure(2);
            movegui(f2,'east');
            xlabel('Iteration nr.');
            title('Optimization');
            hold on
            
            figure(1)
            title('Solution Sequence computed by fmincon');
            hold on
        case 'iter'
             figure(2);  
             xlabel('Iteration nr.');
             title('Mass after Optimization');
             hold on
             plot(optimValues.iteration,optimValues.fval);
             set(gca,'fontsize',15)
             H1 = line(optimValues.iteration,optimValues.fval);
             ylabel('m (kg)','Fontsize',20);
             set(H1,'LineStyle','-','Color','r','Marker','*','MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor','r');
             hold off;
            
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
              history.fval = [history.fval; optimValues.fval];
              history.x = [history.x; x];
              history.iter =[history.iter; optimValues.iteration];

            pause(.3)
            figure(1);
            plot(x(1),x(2),'r*','MarkerSize',10);
       
         case 'done'
             hold off
            
        otherwise
     end
end



