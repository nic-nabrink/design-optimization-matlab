close all;
clear;
clc;
warning('off')
%%%%%%%%%%%%%%%%%%%%%
%%  Problem Setup  %%
%%%%%%%%%%%%%%%%%%%%%

syms x_1 x_2

%TODO: Setup the Optimization Problem in matrix notation

%Edit starts here:---------------------------------------------------------
x=[x_1;x_2];
Q=eye(2);
e=[2; 1];
f=1/2*x'*Q*x + e'*x;
Cons=[1,1; -1 1; -1 -1; 1 0; 0 1; 0 -1];
d=[-5; -2; 0; -5; -2; -1];
%Edit ends here:-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%
%%     Plotting    %%
%%%%%%%%%%%%%%%%%%%%%

fcontour(f,'LineWidth',1,'Levelstep',3)
set(gca,'fontsize',15);
axis([-8 8 -8 8]);
colormap([0 0 0]);
hold on;

for k=1:size(Cons,1) %Plotting of the constraints 
    g_i=Cons(k,:)*x+d(k);
    g=solve(g_i,x_2);
    if isempty(g) %For upper and lower bounds wrt x1
        line([-d(k)/Cons(k,1) -d(k)/Cons(k,1)], [-8 8],'Color','black','Linewidth',1.5)
    else
        fplot(g,'k','Linewidth',1.5);
    end
    axis([-8 8 -8 8])
end

%%%%%%%%%%%%%%%%%%%%%%
%%     Algorithm    %%
%%%%%%%%%%%%%%%%%%%%%%

alpha = 1;  %Initialization of step size
mu=[];      %Initialization of Lagrange Multiplier
s=[];       %Initialization of search direction
x=[5;0];    %Start Vector

%TODO: Define a vector 'active', that contains the indices of all currently
%active constraints [i_1;i_2]

%Edit starts here:---------------------------------------------------------
constr_dist=Cons*x+d;   %calculate the distance of x from all constraints
active=find(~constr_dist);  %find all elements with constr_dist=0;
%Edit ends here:-----------------------------------------------------------

plot(x(1),x(2),'g*','MarkerSize',10);

%Abort if start vector is infeasible
if any(constr_dist > 0)
    disp('Start vector is infeasible!')
    return
end 

for k=1:20 %Iteration limit
    
    pause(1)
    grad=Q*x+e;         %Gradient calculation of current iteration
    C=Cons(active,:);   %Set up C-matrix with active constraints

    %Lagrange-Newton Equation:

    rhs=vertcat(grad,zeros(size(active,1),1));      %Right hand side of Lagrange Newton Equations
    sol=-inv([Q C';C zeros(size(active,1))])*rhs;   %Solve for search direction and Lagrange multiplier mu
    idc = abs(sol)< 1e-14; sol(idc) = 0;            %Numerical errors
    s=horzcat(s,sol(1:2));                          %Stepsize s
    addmu=[sol(3:end);zeros(size(Cons,1)-size(sol(3:end),1),1)];
    mu=horzcat(mu,addmu);                           %Lagrange multiplier mu
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: Complete the conditions for the following cases:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %"The new design is feasible" 
    %Edit starts here:---------------------------------------------------------
    if(any(s(:,k))~=0 && max(Cons*(x+s(:,k))+d)<=0)
        %A new iteration step can be made.
        alpha=1;
    %Edit ends here:-----------------------------------------------------------
    
    %"The new design is infeasible"
    %Edit starts here:---------------------------------------------------------
    elseif (any(s(:,k))~=0 && max(Cons*(x+s(:,k))+d)>0)
        %Plot infeasible design
        x_viol = x+s(:,k);
        plot(x_viol(1),x_viol(2),'r*','MarkerSize',10);
        %Add most violated constraint to set of active constraints
        [~, m] = max(Cons*(x+s(:,k))+d);
        active=[active;m];
        %Calculate stepsize
        alpha = -(Cons(m,:)*x+d(m))/(Cons(m,:)*s(:,k));
    %Edit ends here:-----------------------------------------------------------
    
    %"The design remains the same but does not satisfy the KKT conditions"
    elseif (all(s(:,k)==0) && min(mu(:,k))<0)
        %Remove inequality constraint with smallest mu from the set of
        %active costraints 
        [~, m] = min(mu(:,k));
        active(m)=[];
    
    %"The design remains the same and does satisfy the KKT conditions"
    %Edit starts here:---------------------------------------------------------
    elseif (all(s(:,k)==0) && min(mu(:,k))>0)
    %Edit ends here:-----------------------------------------------------------
        result1=num2str(x(1));
        result2=num2str(x(2));
        it=num2str(k);
        text=['Covergence at ','[',result1,' ',result2,']',' after ',it, ' iterations.'];
        disp(text)
        break
    end   
    
    %set new x according to stepsize
    x = x + alpha*s(:,k);
    plot(x(1),x(2),'g*','MarkerSize',10);
    
end