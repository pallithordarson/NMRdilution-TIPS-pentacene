
% Before running this file, the data needs to be loaded on to the command prompt
%In this case use load nmrdilutiondata;

% Sets the initial start value for the parameter(s)
start1 = [1.2];

ttb=size(DA); %measures size of raw data matrix DA (raw UV-Vis, NMR... to be fitted).
tth=ttb(1); %number of line in raw data = number of titration point in experiment. should ideally be at least > 10.

%Note that number of columns in DA = number of spectra data points (in UV and fluorescence)
%or number of proton resonances in NMR used. 

hostA=repmat(DA(1,:),tth,1); %Creates a new matrix containing all the initial raw data values (in raw 1, before 
%any guest is added). The size of this maxtrix should be the same as DA.

DA=DA(:); %This operation changes the matrix DA into a single column vector. This is 
%essential for the fminsearch optimisation function to work.
DA3=DA;
DA(end+1)=ttb(2); %adds the ttb(2) = number of columns in original DA to the DA column
%vector. This is essential so the DA matrix can be regenerated on the
%"other side" of the fminsearch function.

%%
%This segment execute the fminsearch optimisation process. It is set up as loop 
%that will run up to 10 times or until there is no significant changes in
%the sum of residual squares (ss) between the raw and fitted data which
%fminsearch is trying to make = 0;

options = optimset('TolFun',1e-18,'TolX',1e-18); %Sets the tolerance options for fminsearch

ssold=1e60; %initial (artificially high) value for sum of squares. Stops
%the loop from terminating before the first run.
itn=0; %initial counter for loop function
ssdiff=1e-8;%criteria for loop to stop where ssdiff=ss(this cycle) - ss(last cycle)
tic;%Starts a clock to time the process
while itn<10 %start of loop - itn<10 means it will stop after 10 iterations.

    %The next line executes the fminsearch optimisation.
    %INPUTS for fminsearch are:
    
    %'xx#to#fitbind$$$$' = fitting function that contains the program that calculates
    %naming convention for fitting function:
    %xx = uv or nmr or flu...
    %#to# = stochiometry, e.g. 1to1, 1to2, 2to1,...]\
    %$$$ = optional suffix for variants such as "stat" = statistical 1:2
    %model
    
    %the fitted data (isotherms) and then compares the outcome with the raw
    %data to calculate the sum of squares (ss) which fminsearch tries to
    %minimize.
    
    %start1=input of the initial guess(ss) for the parameters to be fitted
    %(K1, K2....)
    %options = the options for fminsearch as set above.
    %initial = a two column vector with initial host and guest
    %concentrations
    %DA = raw data to be fitted
%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
[results, fval, exitflag, output] = fminsearch('nmraggfitbind',start1,options,initial,DA);
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

%OUTPUT from fminsearch
%results = optimised parameters = K1, K2....
%fval = value of sum of squares (ss, target = 0)
%exitflag = messages from fminsearch about success or otherwise
%output = string containing information about the performance of fminsearch
%such as no of iterations

%loop now checks the outcome from fminsearch:

criter=(ssold-fval)/ssold; %measures difference between this cycle and the last
iter(itn+1)=output.iterations; %saves the no of iterations from fminsearch as = iter
if abs(criter) <=ssdiff; %if the ss-difference between this and the last cycle is 
    %less than (<) the preset criteria above
    break % the loop terminates.
end
ssold=fval; %if loop didn't terminate, ssold (previous cycle) is now reassigned 
%as ss from this cycle (fval) 
start1=results; %The initial guess for the next cycle is now changed to the 
%best fit of parameters for this cycle.
itn=itn+1; %increments the number of cycles completed (until 10)
%start1=results;
end %end of loop.

fittime=toc; %records the time the optimisation took as fittime


%%
%Post-processing - part 1 - calculating the final spectra 

tic; %starts a clock to measure the time for all the post-processing processes.

para=results; %takes the best fit from loop agove (results) and assigns it as = para

%Fitting function 'xx#to#fitbind$$$$' is now executed once more with the
%best fit from optimisation loop above.
%Input = para = K1, K2...., initial and DA = initial concentration and raw
%data as before

%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
[ss, EA, H, HG, HG2, RR]=nmraggfitbind(para,initial,DA);
ss=ss./10000; %corrects for a factor of 10000 that was used to make ss more sensitive in the fitting process.
HGG2=[H+HG/2 HG/2+HG2];
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)
EAspare=EA;
%outputs:

%ss = sum of squares from fitting function
%EA = vector or matrix of Y(DeltaHG), Y(DeltaHG2),... values. 
%E.g. for 1:1 binding and: 
%UV; EA = [DeltaEpsilonWavelength1,DeltaEpsilonWavelength2,....]
%NMR; EA = [DeltadeltainppmorHzProton#1, DeltadeltainppmorHzProton#2,...]
%see e.g. the first term on the right of the "=" in eq. 14
%and 15 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
%HG = Calculated concentration of the HG species from the fitting process
%RR = The calculated matrix of residuals (Raw data - Fitted data) from the
%fitting process.
DAmat=DA;
DA=DA(1:end-1); %removes the "dummy" ttb(2) (see above) again from 
%raw data column DA
DA=reshape(DA,tth,ttb(2)); %Regenerates the original matrix of 
%raw data (DA) as a m x n matrix with tth = number of rows and
%ttb(2) = number of columns

%% 
%Post-processing - part 2 - calculates key statistics 

ndata=length(DA); %measures the number of experimental datapoints in raw data
%ndata = number of titration point x number of spectra in global fit
npara=length(results); %measures the binding of binding constants fitted (usually 1 or 2)
nEA=length(EA(:)); %measures the number of datapoints in global fit fitted 
%(= no of columns in original DA x number of species).
nEA2=length(EA);%measures the number of spectra in global fit fitted (= no of columns in original DA).
deFree=ndata-npara-nEA; %deFree = degree's of freedom 
SE=(sqrt(ss/(deFree))); %Standard devitation of the y-estimate (Eq 52 in 
%P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
chiSQ=(sqrt(ss/(deFree-1))); %Chi-square = Eq 53 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.

%The next section is to calculate the uncertainty for the linear least square process
%that was used in the fitting function program to obtain EA (spectral shifts) by matrix
%division (linear regression).

%Next we calculate the covariance of fit (covfit). This parameter is less
%sensitive to no of parameters than chiSQ and SE - see Eq 54 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
covRR=cov(RR(:)); %Takes the covariance of the residual matrix after converting it to a column vector
covDA=cov(DA(:));%Takes the covariance of the raw data (-initial values)after converting it to a column vector
covFit=covRR./covDA; %calcules the covfit according to Eq 54 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
%

%%

%Data preparation for Monte Carlo process
%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
CA=HGG2*EA; %Calculates the fitted data as EA*HG (spectral shifts x concentration of host-guest)
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

UDA=DA; %Creates a copy of the raw data UDA = DA.
UCA=CA; %Creates a copy of the fitted data UCA=CA

initial2=initial(:); %Converts the concentration matrix into a column vector (required for nlinfit) called initial2
initial2(tth*1+1:tth*1+nEA)=EA(:); %Adds the EA parameters to the initial2 column to pass on to nlinfit  
initial2(tth*1+1+nEA:tth*1+nEA+nEA2)=DA(1,:); %Adds the first row of DA to initial2
initial2(end+1)=nEA; %Adds the number of EA datapoints (nEA) to initial2
initial2(end+1)=nEA2; %Adds the columns of EA datapoints (nEA) to initial2
initial2(end+1)=tth; %Add the number of rows in DA = number of titration point to initial2
DA2=DA(:); %Converts the raw data to a new column vector = DA2



%% The Monte Carlo approach

% First we set the % error on the inputs. These can be changed but I picked
% values that seem reasonable from my tests.
herror=[0.02]; % This is the % error on the host concentration (htot or H0)
%gerror=[0.01]; % This is the % error on the guest concentration (Ltot or G0)
% N.b. initial = [htot Ltot]
dataerror=[0.001]; %This is the % errror on the 
sizeDA=ttb(2); % If doing a global fit sizeDA > 1 

cycln=20000; %Number of iteration
tic %times the process

% The above % errors are now used to create the "random" multipliers
% randn = normally distributed random numbers with sigma = +/- 1.

Dataerror=randn(tth*sizeDA,cycln)*dataerror;    
%Herror=repmat(1+randn(1,cycln)*herror,tth,1);
Herror=1+(randn(tth,cycln)*herror/2)+(repmat(randn(1,cycln)*herror/2,tth,1));
%n.b. the error on G0 (Gerror) is split between a "systemic" and "randon" 
%error while on H0 (Herror) it is all "systemic" and on data (Dataerror)
% is all "random".

% This loop now fits the cycln x number of input data (here 200) to the
% model and collects all the key data.
oldinitial=initial;
for i=1:cycln
start1=para; % we start always with the "correct" value as a gues.
initial=[oldinitial.*Herror(:,i)]; 
DA=DA2+Dataerror(:,i);
DA22=DA;
DA(end+1)=ttb(2);

% The next part is exactly as the loop above
ssold=1e60;
itn=0; %initial counter for loop function
ssdiff=1e-8;%criteria for loop to stop where ssdiff=ss(this cycle) - ss(last cycle)

while itn<10
    
[results, fval, exitflag, output] = fminsearch('nmraggfitbind',start1,options,initial,DA);
criter=(ssold-fval)/ssold; %measures difference between this cycle and the last
iter(itn+1)=output.iterations; %saves the no of iterations from fminsearch as = iter
if abs(criter) <=ssdiff; %if the ss-difference between this and the last cycle is 
    %less than (<) the preset criteria above
    break % the loop terminates.
end
ssold=fval; %if loop didn't terminate, ssold (previous cycle) is now reassigned 
%as ss from this cycle (fval) 
start1=results; %The initial guess for the next cycle is now changed to the 
%best fit of parameters for this cycle.
itn=itn+1; %increments the number of cycles completed (until 10)
%start1=results;
end %end of loop.

% This part of the loop collects the initial values
initialall(:,i)=initial(:,1);
Datumall(:,i)=DA22(:,1);
resultsall(i,:)=results; % collects all the results from the Monte Carlo
[ss, EA, HG, HG2, RR]=nmraggfitbind(results,initial,DA); % Calc once more the outputs according
%to the results from the Monte Carlo in each cycle
HGG3=[H+HG2/2, HG2/2+HG];
EAall(i,:)=EA(:)';
CA=HGG3*EA;
CA=CA(:);
yall(:,i)=CA; % gathers all the calc. outputs from the Monte Carlo
end

% post-processing after the cycln (200) Monte Carlo calculations are
% finished

% Rank all "results" and find the 2.5% lowest and 97.5% highest values
% (the 2.5% and 97.5% percentile). The difference = 95% confidence Interval
perspace2=prctile(resultsall,[2.5 97.5]); 
perspace1=prctile(resultsall,[15.8 84.1]); 
perspace2=perspace2(:)'; % creates a row vector (in case > 1 parameter)
perspace1=perspace1(:)';
EAspace2=prctile(EAall,[2.5 97.5]); 
EAspace1=prctile(EAall,[15.8 84.1]); 
EAper2 = 100*(1-bsxfun(@rdivide, EAspace2, EAspare(:)'));
EAper1 = 100*(1-bsxfun(@rdivide, EAspace1, EAspare(:)'));
% stores the results in a structure for the future
%Fitstuff={'initialstuff',initialall,'yall',yall,'perspace',perspace2,perspace1,resultsall,herror,'DAerror',dataerror,EAper2,EAper1,EAall};
% times how long the Monte Carlo took.

fitMC=toc;



    