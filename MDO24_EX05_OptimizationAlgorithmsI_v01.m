%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Two Bar Truss        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/|\
%/| \
%/|  \
%/|   \
%/|    \
%/|     \   | F
%/|      \  |
%/|       \ |
%/|________\|
%/|
%/|

clear;
close all;
clc;
method='Descent'; % 'Descent'for Steepest Desceent Algorithm, 'Newton' for Modified Newton's Method 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Plotting        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1=linspace(0,400);
A2=linspace(0,400);
[X,Y]=meshgrid(A1,A2);
p = psi(X,Y);
p(p>2)=2; %Remove too high values for plotting reasons
contour(X,Y,p,'LineWidth',1,'Levelstep',0.05)
set(gca,'fontsize',15)
hold on 

g1 = sqrt(2)*-X+362.03;
plot(X,g1,'k','LineWidth',3.0);
axis([0,400,0,400])
xlabel('A_1 (mm^2)','Fontsize',20)
ylabel('A_2 (mm^2)','Fontsize',20)
colormap([0 0 0])
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Steepest Descent          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Select start vector A^(0) %%

%TODO: Set the termination criteria and the start vector according to the
%exercise

%Edit starts here:---------------------------------------------------------
epsilon=1E-9 ;     %Convergence Test
k_max=25 ;       %Convergence Test
lb= [1 1];       %Lower Bound of Design Space
A_0=[100,50] ;      %Start vector
%Edit ends here:-----------------------------------------------------------

%Catching error because of missing/infeasible initialization
if isempty(epsilon) == 1 || isempty(k_max) == 1 || isempty(lb) == 1 || g(A_0(1),A_0(2))> -1e-3 || 1 == any(A_0 < lb)
    disp('Missing/infeasible initialization!')
    return
end 

k=0;              %Initialize counter variable
plot(A_0(1),A_0(2),'r*');
hold on
delta_psi=psi(A_0(1),A_0(2))-0; %Initialize termination criterium
A_prev=A_0; %initialize solution vector

text= ['Iteration =',num2str(k),'   A =',num2str(A_prev),'   displacement =',num2str(psi(A_prev(1),A_prev(2))),...
    '   delta_psi =',num2str(abs(delta_psi))];
    disp(text)

%TODO: Implement a while loop to repeat algorithm steps as long as
%one of the two termination criteria is not fulfilled.

%Edit starts here:---------------------------------------------------------
while (abs(delta_psi)>epsilon && k<k_max)
%Edit ends here:-----------------------------------------------------------
    
    pause(1) % Pause to track algorithm progress
    
    %% Set search direction s^(k)=?grad(psi)*(A^(k)) %%
    %TODO: Calculate the numerical gradients gradA1,gradA2
    % via central finite differences and a delta_x of 0.1.

    %Edit starts here:---------------------------------------------------------
    delta_x=0.1;
    gradA1=(psi(A_prev(1)+delta_x,A_prev(2)) - psi(A_prev(1)-delta_x,A_prev(2)))/2*delta_x;
    gradA2=(psi(A_prev(1),A_prev(2)+delta_x) - psi(A_prev(1),A_prev(2)-delta_x))/2*delta_x;
    %Edit ends here:-----------------------------------------------------------
        
    if isequal(method, 'Descent')
                
        %TODO: Calculate the search direction from the normalized
        %Gradients (Euclidian norm)
        
        %Edit starts here:-------------------------------------------------
        s=[gradA1;gradA2];
        s_norm = -s/sqrt(s' * s);  %Normalize Gradients
        %Edit ends here:---------------------------------------------------
        
        %TODO: Only Bonus Question! Calculate a search direction s_norm in
        %case you are using a Modified Newton's Method 
    elseif isequal(method,'Newton')
        %Using an Approximation of the Hessian
        H= zeros(2,2);
        H(1,1)=0;
        H(2,2)=0;
        H(1,2)=0;
        H(2,1)=H(1,2);
        s = -1*inv(H)*[0;0]; 
        s_norm = -s / sqrt(s' * s);  %Normalize Gradients
    end
    
    %%  Compute step size a^(k) %%
    
   % Bracketing (Nothing to do here! -> Jump right into Approximation)
    bracket = [0;100/3]; %bracket initialization
    bracketLimit= 200; %set limit for bracketing iterations
    s_norm=s_norm';
    A=zeros(bracketLimit,2);
    A(1,:)=A_prev;
    A(2,:)=A_prev+bracket(2)*s_norm;
    
    %First Bracketing Step
    i=0;loopLimit= 1e2;
    while i<loopLimit
        if psi(A(2,1),A(2,2)) < psi(A(1,1),A(1,2)) && ...
                g(A(2,1),A(2,2))< -1e-3 && 1==all(A(2,:) >= lb)
            break
        end
        %Reduce step size if initial step increases function value
        bracket = [0;bracket(2)*0.5];
        A(2,:)=A_prev+bracket(2)*s_norm;
        i=i+1;
    end
    
    %Prevent infinite while-loop because of zero/incorrect gradient
    if i==loopLimit;disp('Could not find a solution!');return;end
    
    bStep = 2;
    i=3;
    while i <= bracketLimit
        b=bracket(i-1)+bracket(i-1)*bStep;
        bracket=vertcat(bracket,b);
        A(i,:)=A_prev+bracket(i)*s_norm;
        %Checks if the bracketing was succesfull and all designs
        %are within the feasible domain
        if (psi(A(i-2,1),A(i-2,2)) > psi(A(i-1,1),A(i-1,2))... 
         && psi(A(i-1,1),A(i-1,2)) < psi(A(i,1),A(i,2)))...
         && g(A(i,1),A(i,2))< -1e-3 && 1==all(A(i,:) >= lb)
            break
        end
        i=i+1;
        %If no feasible third bracketing step is found, reduce bracketing
        %step size
        if i>bracketLimit || g(A(i-1,1),A(i-1,2)) > -1e-3 || 1 == any(A(i-1,:) < lb)
            bStep = bStep*0.5;
            bracket = [bracket(1);bracket(2)];
            i=3;
        end
    end

    % Approximation
    A=zeros(3,2); %initialize A
    alpha = [bracket(i-2) bracket(i-1) bracket(i)]'; %set sample points according to bracketing result
    
    %TODO: Approximate the objective function along the search direction
    %with a second order polynomial. Use the first order derivative to
    %calculate the stepsize a_star
    
    %Edit starts here:---------------------------------------------------------
    A= A_prev + alpha*s_norm ;%calculate sample points for different alphas
    u= psi(A(:,1),A(:,2)) ; %calculate meta objective funcion value at sample points
    B=[1 alpha(1) alpha(1)^2; ...
        1 alpha(2) alpha(2)^2; ...
        1 alpha(3) alpha(3)^2]; %set up matrix B
    a=B\u; %solve polynomial interpolation for coefficients a_i
    a_star = - a(2)/(2*a(3)); %calculate approximated stepsize a_star
    %Edit ends here:-----------------------------------------------------------
    
    A_curr=A_prev+a_star*s_norm; %calculate solution vector for the current time step
    plot (A_curr(1),A_curr(2),'r*');
    x=[A_prev(1),A_curr(1)];
    y=[A_prev(2),A_curr(2)];
    line (x,y,'Color','red');
    hold on
    
    delta_psi=psi(A_curr(1),A_curr(2))-psi(A_prev(1),A_prev(2)); %test for convergence
    A_prev=A_curr;  %set A(k-1)=A(k) for next iteration
    k=k+1;
    
    text= ['Iteration =',num2str(k),'   A =',num2str(A_prev),'   displacement =',num2str(psi(A_prev(1),A_prev(2))),...
    '   delta_psi =',num2str(abs(delta_psi))];
    disp(text)

end


%% Function Declaration %%
function [u] = psi(A_1,A_2)
    C_1=13.47;
    C_2=4.76;
    C_3=362.03;
    r=0.5;
    u= C_1./(A_1) + C_2./(A_2)- r./(sqrt(2).*(A_1)+A_2-C_3);
end

function [m] = g(A_1,A_2)
    C_3=362.03;
    m = sqrt(2).*(A_1)+A_2-C_3;
end