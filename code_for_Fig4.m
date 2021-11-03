%% loading the matlab file containing parameter values
clear all
load('xxxxx.mat');
shortlist = xxxxx; % just the variable name being used in this code

%% chosing the parameter set  and defing the values of parameters 

i = 64; % parameter set row number

ga = shortlist.Prod_of_A(i);  gb = shortlist.Prod_of_B(i);  gc = shortlist.Prod_of_C(i);   %production rates 
ka = shortlist.Deg_of_A(i);   kb = shortlist.Deg_of_B(i);   kc = shortlist.Deg_of_C(i);    %degradationtion rates 

xab=shortlist.Inh_of_AToB(i);   xba=shortlist.Inh_of_BToA(i);       % fold-change of Inhibition of B(A) by A(B)
nab=shortlist.Num_of_AToB(i);   nba=shortlist.Num_of_BToA(i);       % # of binding sites of Inhibition of B(A) by A(B)
a0b=shortlist.Trd_of_AToB(i);   b0a=shortlist.Trd_of_BToA(i);       % threshold of Inhibition of B(A) by A(B)

xbc=shortlist.Inh_of_BToC(i);   xcb=shortlist.Inh_of_CToB(i);       % fold-change of Inhibition of C(B) by B(C)
nbc=shortlist.Num_of_BToC(i);   ncb=shortlist.Num_of_CToB(i);       % # of binding sites of Inhibition of C(B) by B(C)
b0c=shortlist.Trd_of_BToC(i);   c0b=shortlist.Trd_of_CToB(i);       % threshold of Inhibition of C(B) by B(C)

xca=shortlist.Inh_of_CToA(i);   xac=shortlist.Inh_of_AToC(i);       % fold-change of Inhibition of A(C) by C(A)
nca=shortlist.Num_of_CToA(i);   nac=shortlist.Num_of_AToC(i);       % # of binding sites of Inhibition of A(C) by C(A)
c0a=shortlist.Trd_of_CToA(i);   a0c=shortlist.Trd_of_AToC(i);       % threshold of Inhibition of A(C) by C(A)

%% Noise, alpha and other epigenetic feedback parameters ------------------

N = 0;                                 % kept 0 for fig4 simulations
beta=10;                               % scaling factor by which the threshold term in ODE is divided by

tau=0.001;                             % time step of sampling
days = 50;                             % number of days simulation is run (a.u.)
sampling_factor = 100;                 % how many times in an hour of the above arbitrary day, sampling is done
maxstep=days*24*sampling_factor;       % days*hrs_in_a_day*time_step_per_hr gives total time steps taken (a.u.)
maxtimes=100;                          % total number of cells in population

%% defining alpha values for phase plots

% defining the ranges for which the specific alpha values are to be varied
alpha_ba_ar = 0;
alpha_ab_ar = 0:0.05:0.1;
alpha_bc_ar = 0;
alpha_cb_ar = 00:0.05:0.1;
alpha_ca_ar = 0;
alpha_ac_ar = 0;

iterations = 30; % number of repetetions to take an average at end to account for various initial conditions

% sizes below are to be changed based on which alpha value is varied, iterations and how long simulation is being run for
PA_final = zeros(size(alpha_cb_ar,2),iterations,days+1);  
PB_final = zeros(size(alpha_cb_ar,2),iterations,days+1);
PC_final = zeros(size(alpha_cb_ar,2),iterations,days+1);
    
PA_state = zeros(iterations,days+1);
PB_state = zeros(iterations,days+1);
PC_state = zeros(iterations,days+1);

% Initializing alpha values - will be chosen later inside loop accordingly
alpha_ba = 0;
alpha_bc = 0;
alpha_ab = 0;
alpha_ac = 0;
alpha_cb = 0;
alpha_ca = 0;

% variables used to define fraction of time for which specific/selective
% alpha values will be given epigenetic feedback for
y = 0.2;

%% Outer loops to first provide values of alpha_1, alpha_2 and iteration number

for o = 1:size(alpha_ac_ar,2)
    
    alpha_ab = alpha_ab_ar(o);
    alpha_cb = alpha_cb_ar(o);
    
for p=1:iterations
    
%% initialization the matrices

% node expression values
a = zeros(maxtimes,maxstep+1);  
b = zeros(maxtimes,maxstep+1);
c = zeros(maxtimes,maxstep+1);

% updating threshold values with epigenetic feedback
newb0a=zeros(maxtimes,maxstep+1);
newa0b=zeros(maxtimes,maxstep+1);
newb0c=zeros(maxtimes,maxstep+1);
newc0b=zeros(maxtimes,maxstep+1);
newc0a=zeros(maxtimes,maxstep+1);
newa0c=zeros(maxtimes,maxstep+1);

% initial conditions of threhold values for starting simulation
newb0a(:,1)=b0a;
newa0b(:,1)=a0b;
newb0c(:,1)=b0c;
newc0b(:,1)=c0b;
newc0a(:,1)=c0a;
newa0c(:,1)=a0c;

for j=1:maxtimes
   
% initial conditions of node expression for starting simulation
a(:,1) = (ga/ka)*rand;
b(:,1) = (gb/kb)*rand;
c(:,1) = (gc/kc)*rand;
    
for i=1:maxstep-1
    
%% Defining the ODE

Hills_ab = (1+xab*(a(j,i)/newa0b(j,i))^nab)/(1+(a(j,i)/newa0b(j,i))^nab);
Hills_ba = (1+xba*(b(j,i)/newb0a(j,i))^nba)/(1+(b(j,i)/newb0a(j,i))^nba);

Hills_bc = (1+xbc*(b(j,i)/newb0c(j,i))^nbc)/(1+(b(j,i)/newb0c(j,i))^nbc);
Hills_cb = (1+xca*(c(j,i)/newc0b(j,i))^ncb)/(1+(c(j,i)/newc0b(j,i))^ncb);

Hills_ca = (1+xcb*(c(j,i)/newc0a(j,i))^nca)/(1+(c(j,i)/newc0a(j,i))^nca);
Hills_ac = (1+xac*(a(j,i)/newa0c(j,i))^nac)/(1+(a(j,i)/newa0c(j,i))^nac);

% add noise to the ODES every 2000 steps
   if mod(i,500)==0
    rnumber3 = normrnd(0,0.25,[1,3]);
    a(j,i+1) = a(j,i) + tau*(ga*Hills_ba*Hills_ca - ka*a(j,i)) + rnumber3(1,1)*N;
    b(j,i+1) = b(j,i) + tau*(gb*Hills_ab*Hills_cb - kb*b(j,i)) + rnumber3(1,2)*N;
    c(j,i+1) = c(j,i) + tau*(gc*Hills_ac*Hills_bc - kc*c(j,i)) + rnumber3(1,3)*N;
   else
    a(j,i+1) = a(j,i) + tau*(ga*Hills_ba*Hills_ca - ka*a(j,i));
    b(j,i+1) = b(j,i) + tau*(gb*Hills_ab*Hills_cb - kb*b(j,i));
    c(j,i+1) = c(j,i) + tau*(gc*Hills_ac*Hills_bc - kc*c(j,i));
   end

% controlling the time when epigenetic feedback is switched ON and OFF with
% if conditions and the predefined variable y. Arbitrary values of 
% y were taken here as an example

   if i > maxstep*0
       newb0c(j,i+1) = newb0c(j,i) + tau*(b0c-newb0c(j,i)-0*b(j,i))/beta;    % NO epigenetic feedback
   else
   newb0c(j,i+1) = newb0c(j,i) + tau*(b0c-newb0c(j,i)-alpha_bc*b(j,i))/beta; % epigenetic feedback
   end
   

   if i > maxstep*0
       newb0a(j,i+1) = newb0a(j,i) + tau*(b0a-newb0a(j,i)-0*b(j,i))/beta;    % NO epigenetic feedback
   else
   newb0a(j,i+1) = newb0a(j,i) + tau*(b0a-newb0a(j,i)-alpha_ba*b(j,i))/beta; % epigenetic feedback
   end
   
    if i > maxstep*y
       newc0b(j,i+1) = newc0b(j,i) + tau*(c0b-newc0b(j,i)-0*c(j,i))/beta;    % NO epigenetic feedback
   else
   newc0b(j,i+1) = newc0b(j,i) + tau*(c0b-newc0b(j,i)-alpha_cb*c(j,i))/beta; % epigenetic feedback
   end
   

   if i > maxstep*y
       newa0b(j,i+1) = newa0b(j,i) + tau*(a0b-newa0b(j,i)-0*a(j,i))/beta;    % NO epigenetic feedback
   else
   newa0b(j,i+1) = newa0b(j,i) + tau*(a0b-newa0b(j,i)-alpha_ab*a(j,i))/beta; % epigenetic feedback
   end
   
   if i > maxstep*y1
       newc0a(j,i+1) = newc0a(j,i) + tau*(c0a-newc0a(j,i)-0*c(j,i))/beta;    % NO epigenetic feedback
   else
   newc0a(j,i+1) = newc0a(j,i) + tau*(c0a-newc0a(j,i)-alpha_ca*c(j,i))/beta; % epigenetic feedback
   end
   

   if i > maxstep*0
       newa0c(j,i+1) = newa0c(j,i) + tau*(a0c-newa0c(j,i)-0*a(j,i))/beta;    % NO epigenetic feedback
   else
   newa0c(j,i+1) = newa0c(j,i) + tau*(a0c-newa0c(j,i)-alpha_ac*a(j,i))/beta; % epigenetic feedback
   end
   
% if expression value goes below 0 due to perturbations
if a(j,i)<0
   a(j,i)=0;
end
if b(j,i)<0
   b(j,i)=0;
end  
if c(j,i)<0
   c(j,i)=0;
end  

end
end

%% statistical analysis

% defining counters for population corresponding to three single-positive
% states
a_state=zeros(1,days+1);
b_state=zeros(1,days+1);
c_state=zeros(1,days+1);

% checking state of every cell at end of every day and updating population
% of that state
for i=0:days
    for j=1:maxtimes
        if (a(j,i*24*sampling_factor+1)/b(j,i*24*sampling_factor+1))>1 && (a(j,i*24*sampling_factor+1)/c(j,i*24*sampling_factor+1))>1
           a_state(i+1)=a_state(i+1)+1; 
        elseif (b(j,i*24*sampling_factor+1)/a(j,i*24*sampling_factor+1))>1 && (b(j,i*24*sampling_factor+1)/c(j,i*24*sampling_factor+1))>1
           b_state(i+1)=b_state(i+1)+1;
        elseif (c(j,i*24*sampling_factor+1)/a(j,i*24*sampling_factor+1))>1 && (c(j,i*24*sampling_factor+1)/a(j,i*24*sampling_factor+1))>1
           c_state(i+1)=c_state(i+1)+1;
        end
    end
end

% calculate the population distribution
pa_state=zeros(1,days+1);
pb_state=zeros(1,days+1);
pc_state=zeros(1,days+1);
day=zeros(1,days+1);

% calculate population percentages
for i=1:days+1
    pa_state(i)=100*a_state(i)/(a_state(i)+b_state(i)+c_state(i));
    pb_state(i)=100*b_state(i)/(a_state(i)+b_state(i)+c_state(i));
    pc_state(i)=100*c_state(i)/(a_state(i)+b_state(i)+c_state(i));
    day(i)=i-1;
end

% for a fixed value of alpha_1, data corresponding to (all alpha_2,day)
PA_state(p,:) = pa_state;
PB_state(p,:) = pb_state;
PC_state(p,:) = pc_state;

end

% data corresponding to all (alpha_1,alpha_2,day)
PA_final(o,:,:) = PA_state;
PB_final(o,:,:) = PB_state;
PC_final(o,:,:) = PC_state;
    
end

%% Averaging over the number of iterations

A_mean = mean(PA_state(:,end-1));
B_mean = mean(PB_state(:,end-1));
C_mean = mean(PC_state(:,end-1));

A_std = std(PA_state(:,end-1));
B_std = std(PB_state(:,end-1));
C_std = std(PC_state(:,end-1));

%% Final steps
% WE RAN THE CODE FOR VARYING VALUE OF y (i.e. varying fraction of time for
% which epigenetic feedback is witched ON)AT COMBINATIONS OF DIFFERENT 
% alpha_1=alpha_2 VALUES FOR BISTABLE AS WELL AS TRISTABLE PARAMETER SETS.
% TOGETHER TEHY WERE THEN PLOTTED AS FIG 4D, 4F ANF 4G