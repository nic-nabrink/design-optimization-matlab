%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Cantilever beam      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/|            
%/|            |F
%/|-------------
%/| 
%/| 


% Clean up the workspace
clear;
close all;
clc;

% Seed random number generator with current time
rng ('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Problem Setup      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System parameters
E        = 210000;      % [MPa]
F        = 10;          % [N]

L_lb     = 900;         % [mm]
L_ub     = 1100;

I_lb     = 4218.75;     % [mm^4]
I_ub     = 33750;


% Plotting of System Behaviour
syms L I;
u_surf=@(L,I) (F*L.^3)./(3*E*I);

figure(1)
fsurf(u_surf,[L_lb L_ub I_lb I_ub],'FaceColor','w','FaceAlpha',0.)
view(-120,30)
ylabel('I (mm^4)');
xlabel('L (mm)');
zlabel('u (mm)');
hold on

%%%%%%%%%%%%%%%%%%%%%%%
%%     Sampling      %%
%%%%%%%%%%%%%%%%%%%%%%%

nParameters=2;
samplingMethod='Full Factorial' %'Monte Carlo','Full Factorial'

switch samplingMethod
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     3.4.1 Full Factorial Sampling      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Full Factorial'
    
    %TODO: set the number of Levels according to the task and find an
    %expression for the number of samples as a function of nLevels and 
    %nParameters.
        
    %Edit starts here:-----------------------------------------------------    
    nLevels= 4; %(Assuming same number of levels for each parameter)  
    nSamples= nLevels^nParameters;% Formula from Lecture
    %--> Go to Line 116 to distribute the sample points according to R.
    %Edit ends here:-------------------------------------------------------
    
    % Initialize FullFact Matrix with the correct dimensions
    FullFact=ones(nSamples,nParameters);

    % Actual implementation of the Full Factorial Table
    for j=1:nParameters

        % counter variable determines level at which sample is evaluated
        counter=1;
        for i = 1:nSamples
            if counter>nLevels %Check if counter exceeds nLevels and reset to 1
                counter=1;
            end
            FullFact(i,j)=counter;  %Assign current cell with counter value
            dimRemaining=nParameters-j; 

            if mod(i,nLevels^(dimRemaining))==0 
                counter=counter+1;
            end
        end
    end
    
    % Condense samples to [0,1]
    R=1/(nLevels-1).*(FullFact-ones(size(FullFact)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    3.4.3 Monte Carlo Sampling     %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Monte Carlo'   
        
    %TODO: Use the rand() function to create an a matrix R containing 30
    %sample points in 2 coordinates 

    %Edit starts here:-----------------------------------------------------
    nSamples = 30; 
    R= rand(nSamples,nParameters);
    %Edit ends here:-------------------------------------------------------
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     3.4.1 Full Factorial Sampling      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: Use the matrix R (already defined for 'Full Factorial' case) to 
%distribute the sample points within [L_lb,L_ub];[I_lb,I_ub]..

%Edit starts here:---------------------------------------------------------
L= L_lb + R(:,1)*(L_ub - L_lb); % enter vector with distributed points for L
I= I_lb + R(:,2)*(I_ub - I_lb); % enter vector with distributed points for I
%Edit ends here:-----------------------------------------------------------

% System response
u=(F *L.^3)./(3 *E *I);

% Complete DOE Table
DOE_table=[u L I]; 

%Matlab itself provides a Design of Experiment Toolbox, which allows to
%generate and evealuate different Experiment Designs 
%(Including Fractional Factorial and LHS) with very little effort.
%Look into:
%https://de.mathworks.com/help/stats/design-of-experiments-1.html


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     3.4.2 Correlation Coefficients     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(samplingMethod,'Full Factorial')
    
    % For I=4218.75 mm^4 only a fraction of the Full Factorial Table is
    % required. However for I=[I_lb,I_ub] all sampling points are relevant.

    %TODO: Complete the following lines of code to extract the relevant 
    %sample points.

    %Edit starts here:---------------------------------------------------------
    x= 1:nLevels:nSamples; %Hint: Use the ":" operator to define the increment between the sample points
    relevantData= DOE_table(x,:); %Vary this line for 3.4.2.b)
    %Edit ends here:-----------------------------------------------------------

    N=size(relevantData,1);

    %TODO: Calculate Mean, standard Deviation and Covariance in order to 
    %compute the Pearson Correlacion Coefficient.

    %Edit starts here:---------------------------------------------------------
    Mean= 1/N*sum(relevantData,1);
    stdDeviation= 1/N*(sum(relevantData.^2,1)-N*Mean);
    Covariance= 1/N*(sum(relevantData(:,1).*relevantData(:,2))-N*Mean(1)*Mean(2));
    %Edit ends here:-----------------------------------------------------------
    
    if Covariance~=0
        PearsonCoefficient=Covariance/(stdDeviation(1)*stdDeviation(2))
    end
end

%Plotting of Sample Points
scatter3(L,I,u,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r');
hold on

Y=DOE_table(:,1); 
X_L=[ones(nSamples,1) DOE_table(:,2)]; 
X_I=[ones(nSamples,1) DOE_table(:,3)]; 

if exist('Y') ~=0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    3.4.4 Regression Model      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %TODO: Calculate a 2-Dimensional Regression Model based on
    %L and I

    %Edit starts here:---------------------------------------------------------
    % Define X and Y:
    X = [ones(nSamples,1) DOE_table(:,2) DOE_table(:,3)];
    Y = DOE_table(:,1);
    
    % Calculate Regression coefficient:
    beta= (X'*X)^(-1)*X'*Y; %Use formula from the lecture

    syms L I;

    % Calculate Displacement based on L,I and beta:
    u_reg=@(L,I) beta(1) + beta(2)*L + beta(3)*I; %Complete the function handle to calculate u;
    %Edit ends here:-----------------------------------------------------------
end

if exist('u_reg') ~=0
    figure(1)
    fsurf(u_reg,[L_lb L_ub I_lb I_ub],'FaceAlpha',0.5,'EdgeColor','none')
    view(-120,30)
    ylabel('I (mm^4)');
    xlabel('L (mm)');
    zlabel('u (mm)');
end



