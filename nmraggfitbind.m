function [ss, EA, H, HG, HG2, RR, CA3]=nmraggfitbind(para,initial,DA) %ss is the outcome of this function, it
%%
%               FIT BINDING
%
%
%(C) Dr. Pall Thordarson
%School of Chemistry
%UNSW
%AUSTRALIA
%p.thordarson@unsw.edu.au
%
%Please cite: P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323 
%when using this program.
%
%A program for determining binding constants from titration experiments in
%supramolecular chemistry
%
%This sub-program calculates the binding isotherms based on the 
%parameters (K1, K2,...) provided, either directly or more imporantly by the 
%fminsearch optimisation program. 

%The inputs for this function are
%"para", "initial" and "DA" are the parameters to be fitted
%(K1, K2...), total host and guest concentrations (initial) and 
%raw experimentla data (DA) as explained in more details in 
%the function that calls this program (note para = start1 in that program).

%The output of this function are 
%ss = sum of squares, EA = vector or matrix of Y(DeltaHG), Y(DeltaHG2),...
%values, HG (and HG2 etc..) are the calculated HostGuest complex
%concentration and RR = residual matrix. See the function that calls this
%program for more details. 

%%
%This section extracts the relevant values from the initial input arguments
%above.

%Extracts total host and guest concentrations:
htot = initial(:,1);%extracts the Ltot = total conc. of guest 
%Ltot = initial(:,2);%extracts the htot = total conc. of host 

%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
KE=para(1); %Extracts the first binding constant = K1 (K11)from the input "para" (start1 in program calling this function).
%rho=1; %Extracts the first binding constant = K1 (K11)from the input "para" (start1 in program calling this function).
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

%Converts the raw data (DA) column vector back to the original DA matrix
ttb=DA(end); %extracts the terminal value ttb = number of columns in original matrix
DA=DA(1:end-1); %removes ttb from the DA column vector. 
tth=length(DA)/ttb; %calculates the no of rows in original raw data matrix
%DA as tth = length of DA column vector / ttb (no of columns).

DA=reshape(DA,tth,ttb);%Regenerates the original matrix of 
%raw data (DA) as a m x n matrix with tth = number of rows and
%ttb = number of columns

%%
%This section calculates the required concentration(s) of host-guest
%complex(es).


%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)  
%The function uv1to2bb below solves the cubic equation Eq 17 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323 
%The solution to this cubic equation (b) is the free guest [G] concentration


%Now that the free guest concentration [G]=b has been calcualted we can
%calcuate the concentrations of the 1:1 host-guest complex (HG) and the 1:2
%host-guest*2(HG2) complex as derived from Eq 17 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323 

%Because this is UV-Vis data we use a factor of htot = total concentration
%see Eq 35 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323 
LLL=htot.*KE;

sqr = (((4.*LLL+1)));
nemi = ((2.*LLL.*LLL));
toppi = (((1+2.*LLL)-(sqrt(sqr))));
harx = (toppi./nemi);
hh=harx;   


alfa=hh;
xxxx=alfa.*LLL;
svigi=(1-xxxx);
gamm=((2.*alfa.*xxxx)./(svigi));
epsi=((alfa.*(xxxx.^2))./(((svigi).^2)));


H=alfa;
HG2=epsi;
HG=gamm;

H=alfa;
HG2=epsi;
HG=gamm;

%HGG=[H HG HG2];
HGG=[H+HG/2 HG/2+HG2]; %Gathers the 1:1 host-guest and 1:2 host-guest*2 values to a
%matrix HGG for concentrations.
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)


%%
%This section performs the linear regression least-square process 

%to find the EA = vector or matrix of Y(DeltaHG), Y(DeltaHG2),...values
%by performing a matrix divison.
% See for instance Eq 4.51, page 140 in "Practical Data Analysis in Chemistry" by M. Maeder and Y.-M. Neuhold,
%Elsevier, 2007 (Vol 26 in the series Data handling in Science and
%Technology) for more info.

UA=DA; %Assigns UA = DA = Raw data matrix.
EA=HGG\UA;
%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)  
%EA=nonneg(HGG',UA'); %Performs the matrix division = linear regression by least-squares
RR=UA-HGG*EA; %Calculates the residual matrix RR as the difference between
%UA = raw data matrix and HGG*EA - calculated data which is the matrix 
%product (HGG x EA) of hostguest complex matrix (HGG) and the EA = vector of
%Y(DeltaHG)... values.
CA=HGG*EA; %Calculates the data again as the matrix product HGG x EA so it can 
%be passed back to the original function as one of the outputs of this
%function.
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)
CA3=CA;
%%
%This final segment calculates the sum of squares (ss) which is the target
%function for the optimisation function fminsearch to make as small as
%possible (=0)
rr=RR(:); %Converts the residual matrix RR from above to a column vector = rr

ss=sum(rr.^2)*10000; %calculates the sum of squares of rr (residual matrix).
%A factor of 10000 was added to increase the sensitivity of the fitting
%process.
%end of this function.