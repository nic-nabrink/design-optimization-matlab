%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Truss optimization      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%?//////?//////? 
% \     |     /
%  \    |    /
%   \   |   / 
%    \  |  /
%     \ | / 
%      \|/
%    F1/ \F2

% Clean up the workspace
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%
%  Problem Setup  %%%
%%%%%%%%%%%%%%%%%%%%%

syms A1 A2;

%System parameters
F = 150000;               %load (F1=F2=F)
sigma_yield = 216;        %maximum allowed stress for steel
Rho = 7.9e-6;
l = 1000;
lb = ones(2,1)*1;         %lower bounds on cross sectional areas
ub = ones(2,1)*1000;      %upper bounds on cross sectional areas

% objective function
m = Rho*l*(2*A1*sqrt(2)+A2);

% constraints
s1 = F*(sqrt(2)*A1+A2)/(sqrt(2)*A1^2+2*A1*A2);
s2 = F*(sqrt(2)*A1)/(sqrt(2)*A1^2+2*A1*A2);
s3 = F*(A2)/(sqrt(2)*A1^2+2*A1*A2);

%TODO: solve all constraints for A2. g(A1,A2)=< 0 --> A2 = f(A1)
%(Hint: you can use the solve(eqn, x,'ReturnConditions',true) function from
%the symbolic math toolbox)

%Edit starts here:---------------------------------------------------------
[g1,p,c]=solve(s1==sigma_yield,A2,'ReturnConditions',true); 
[g2,p,c]=solve(s2==sigma_yield,A2,'ReturnConditions',true);
[g3,p,c]=solve(s3==sigma_yield,A2,'ReturnConditions',true);
%Edit ends here:-----------------------------------------------------------

%TODO: save all constraints in one constraint-vector g to make them easily
%accessible for plotting later.

%Edit starts here:---------------------------------------------------------
g=[g1;g2;g3]; 
%Edit ends here:-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%
%  Contour Plot  %%%
%%%%%%%%%%%%%%%%%%%%
figure(1);
fcontour(m,'LineWidth',1,'Levelstep',0.3e+1)
set(gca,'fontsize',15)
colormap([0 0 0])
ratio=0;

%TODO: calculate the gradient of the objective function for the direction of decreasing
%mass in the contour plot. (remember: m = Rho*l*(2*A1*sqrt(2)+A2))

%Edit starts here:---------------------------------------------------------
dmA1 = Rho*l*(2*sqrt(2)); 
dmA2 = Rho*l;
dmA = -[dmA1; dmA2];
dmA = dmA/norm(dmA);
%Edit ends here:-----------------------------------------------------------

aspect_y=1.5;

if exist('dmA') ~=0
    x = [0.7 0.7 + 0.1*dmA(1,1)];
    y = [0.7 0.7 + 0.1*dmA(2,1)*aspect_y];
    str='$$ -\nabla f$$';
    a = annotation('textarrow',x,y,'LineWidth',2,'String',str,'Interpreter','Latex');
    a.FontSize=16;
end    
hold on 

axis([lb(1),ub(1),lb(2),ub(2)])
xlabel('A_1,A_3 (mm^2)','Fontsize',20)
ylabel('A_2 (mm^2)','Fontsize',20)

%TODO: complete the following for-loop to plot all constraints inside the
%contour plot. Hint: You can use the fplot() function

%Edit starts here:---------------------------------------------------------
for i=1:3                          
    fplot(g(i),'k','LineWidth', 3)
end
%Edit ends here:-----------------------------------------------------------

pbaspect([aspect_y 1 1])