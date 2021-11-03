% SAME CODE IS USED IN FIG2, FIG3, FIGS2, FIGS3, FIGS4 BUT WITH CHANGING 
% COMBINATIONS OF ALPHA_1 AND ALPHA_2 CORRESPONDING TO THE PARTICULAR SIMULATION

%% loading the matlab file containing parameter values
clear all
load('xxxxx.mat');
shortlist = xxxxx; % just the variable name being used in this code

%% chosing the parameter set  and defing the values of parameters 

i = 64; % parameter set row number

ga = shortlist.Prod_of_A(i);  gb = shortlist.Prod_of_B(i);  gc = shortlist.Prod_of_C(i);   %production rates 
ka = shortlist.Deg_of_A(i);  kb = shortlist.Deg_of_B(i);  kc = shortlist.Deg_of_C(i);      %degradationtion rates 

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

N = 5;                                 % noise arbitrary magnitude
beta=100;                              % scaling factor by which the threshold term in ODE is divided by

tau=0.01;                              % time step of sampling
days = 100;                            % number of days simulation is run (a.u.)
sampling_factor = 100;                 % how many times in an hour of the above arbitrary day, sampling is done
maxstep=days*24*sampling_factor;       % days*hrs_in_a_day*time_step_per_hr gives total time steps taken (a.u.)
maxtimes=50;                           % total number of cells in population

%% defining alpha loops for phase plots

% defining the ranges for which the specific alpha values are to be varied
alpha_ba_ar = 0;
alpha_ab_ar = 0:0.05:0.1;
alpha_bc_ar = 0;
alpha_cb_ar = 00:0.05:0.1;
alpha_ca_ar = 0;
alpha_ac_ar = 0;

% sizes below are to be changed based on which alpha value is varied and how long simulation is being run for
PA_final = zeros(size(alpha_ab_ar,2),size(alpha_cb_ar,2),days+1);  
PB_final = zeros(size(alpha_ab_ar,2),size(alpha_cb_ar,2),days+1);
PC_final = zeros(size(alpha_ab_ar,2),size(alpha_cb_ar,2),days+1);
    
PA_state = zeros(size(alpha_cb_ar,2),days+1);
PB_state = zeros(size(alpha_cb_ar,2),days+1);
PC_state = zeros(size(alpha_cb_ar,2),days+1);

% Initializing alpha values - will be chosen later inside loop accordingly
alpha_ba = 0;
alpha_bc = 0;
alpha_ab = 0;
alpha_ac = 0;
alpha_cb = 0;
alpha_ca = 0;

%% Outer loops to first provide values of alpha_1 and alpha_2 

for o = 1:size(alpha_ab_ar,2)
    
alpha_ab = alpha_ab_ar(o); % chosing value of bifurcation parameter 1 or alpha_1
    
parfor p = 1:size(alpha_cb_ar,2)
        
alpha_cb = alpha_cb_ar(p); % chosing value of bifurcation parameter 2 or alpha_2
       
%% initializing the matrices

% node expression values
a = zeros(maxtimes,2);  
b = zeros(maxtimes,2);
c = zeros(maxtimes,2);

% updating threshold values with epigenetic feedback
newb0a=zeros(maxtimes,2);
newa0b=zeros(maxtimes,2);
newb0c=zeros(maxtimes,2);
newc0b=zeros(maxtimes,2);
newc0a=zeros(maxtimes,2);
newa0c=zeros(maxtimes,2);

% variable to store state of cells at every day
a_final = zeros(maxtimes,days);
b_final = zeros(maxtimes,days);
c_final = zeros(maxtimes,days);


%% Loop for solving the Ordinary Differential Equations for one set of values od alpha_1 and alpha_2

for j=1:maxtimes
    
% initial conditions of node expression for starting simulation
a(j,1) = (ga/ka)*rand;
b(j,1) = (gb/kb)*rand;
c(j,1) = (gc/kc)*rand;

% initial conditions of threhold values for starting simulation
newb0a(j,1)=b0a;
newa0b(j,1)=a0b;
newb0c(j,1)=b0c;
newc0b(j,1)=c0b;
newc0a(j,1)=c0a;
newa0c(j,1)=a0c;

k = 1; % counter for number of days

for i=1:maxstep

%% Defining the ODE

Hills_ab = (1+xab*(a(j,1)/newa0b(j,1))^nab)/(1+(a(j,1)/newa0b(j,1))^nab);
Hills_ba = (1+xba*(b(j,1)/newb0a(j,1))^nba)/(1+(b(j,1)/newb0a(j,1))^nba);

Hills_bc = (1+xbc*(b(j,1)/newb0c(j,1))^nbc)/(1+(b(j,1)/newb0c(j,1))^nbc);
Hills_cb = (1+xcb*(c(j,1)/newc0b(j,1))^ncb)/(1+(c(j,1)/newc0b(j,1))^ncb);

Hills_ca = (1+xca*(c(j,1)/newc0a(j,1))^nca)/(1+(c(j,1)/newc0a(j,1))^nca);
Hills_ac = (1+xac*(a(j,1)/newa0c(j,1))^nac)/(1+(a(j,1)/newa0c(j,1))^nac);

%add noise to the ODES every 500 steps
if mod(i,500)==0
 rnumber3 = normrnd(0,1,[1,3]);
 a(j,2) = a(j,1) + tau*(ga*Hills_ba*Hills_ca - ka*a(j,1)) + rnumber3(1,1)*N;
 b(j,2) = b(j,1) + tau*(gb*Hills_ab*Hills_cb - kb*b(j,1)) + rnumber3(1,2)*N;
 c(j,2) = c(j,1) + tau*(gc*Hills_ac*Hills_bc - kc*c(j,1)) + rnumber3(1,3)*N;
else
 a(j,2) = a(j,1) + tau*(ga*Hills_ba*Hills_ca - ka*a(j,1));
 b(j,2) = b(j,1) + tau*(gb*Hills_ab*Hills_cb - kb*b(j,1));
 c(j,2) = c(j,1) + tau*(gc*Hills_ac*Hills_bc - kc*c(j,1));
end

% taking note of cell state at end of every day
if mod(i,2400)==0
 a_final(j,k) = a(j,2);
 b_final(j,k) = b(j,2);
 c_final(j,k) = c(j,2);
 k = k + 1;
end

% epigenetic feedback updating rules
newb0a(j,2) = newb0a(j,1) + tau*(b0a-newb0a(j,1)-alpha_ba*b(j,1))/beta;
newa0b(j,2) = newa0b(j,1) + tau*(a0b-newa0b(j,1)-alpha_ab*a(j,1))/beta;
newb0c(j,2) = newb0c(j,1) + tau*(b0c-newb0c(j,1)-alpha_bc*b(j,1))/beta;
newc0b(j,2) = newc0b(j,1) + tau*(c0b-newc0b(j,1)-alpha_cb*c(j,1))/beta;
newc0a(j,2) = newc0a(j,1) + tau*(c0a-newc0a(j,1)-alpha_ca*c(j,1))/beta;
newa0c(j,2) = newa0c(j,1) + tau*(a0c-newa0c(j,1)-alpha_ac*a(j,1))/beta;

% if expression value goes below 0 due to perturbations
if a(j,1)<0
   a(j,1)=0;
end
if b(j,1)<0
   b(j,1)=0;
end  
if c(j,1)<0
   c(j,1)=0;
end   

% update for next iteration
a(j,1) = a(j,2);
b(j,1) = b(j,2);
c(j,1) = c(j,2);
newb0a(j,1) = newb0a(j,2);
newb0c(j,1) = newb0c(j,2);
newa0b(j,1) = newa0b(j,2);
newa0c(j,1) = newa0c(j,2);
newc0b(j,1) = newc0b(j,2);
newc0a(j,1) = newc0a(j,2);

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
for i=1:days-2
for j=1:maxtimes
    if a_final(j,i)>b_final(j,i) && a_final(j,i)>c_final(j,i)
       a_state(i) = a_state(i)+1; 
    elseif b_final(j,i)>a_final(j,i) && b_final(j,i)>c_final(j,i)
       b_state(i) = b_state(i)+1;
    elseif c_final(j,i)>b_final(j,i) && c_final(j,i)>a_final(j,i)
       c_state(i) = c_state(i)+1;
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

% for a fixed value of alpha_1, data corresponding to all (alpha_2,day)
PA_state(p,:) = pa_state;
PB_state(p,:) = pb_state;
PC_state(p,:) = pc_state;

end

% data corresponding to all (alpha_1,alpha_2,day)
PA_final(o,:,:) = PA_state;
PB_final(o,:,:) = PB_state;
PC_final(o,:,:) = PC_state;
    
end

%% averaging node expression for values of alpha_1 and alpha_2 to accomodate for noise

n = size(alpha_ab_ar,2); % nxn matrix is size of alpha values of 1 and 2 (we generally choose both same size)
A = zeros(n); %A expression vlaue matrix
B = zeros(n); %B expression vlaue matrix
C = zeros(n); %C expression vlaue matrix

for i = 1:n
    for j = 1:n
        
        A(i,j) = mean(PA_final(i,j,20:end-10));
        B(i,j) = mean(PB_final(i,j,20:end-10));
        C(i,j) = mean(PC_final(i,j,20:end-10));
        
    end
end
%% Plotting the phase plots

figure()
graphA = pcolor(alpha_ab_ar,alpha_cb_ar,A');
% caxis([0 100]);
title('Population percentage of A');
xlabel('alpha value of epigenetic feedback from A-|B');
ylabel('alpha value of epigenetic feedback from C-|B');
pbaspect([1 1 1]);

figure()
graphB = pcolor(alpha_ab_ar,alpha_cb_ar,B');
% caxis([0 100]);
title('Population percentage of B');
xlabel('alpha value of epigenetic feedback from A-|B');
ylabel('alpha value of epigenetic feedback from C-|B');
pbaspect([1 1 1]);

figure()
graphC = pcolor(alpha_ab_ar,alpha_cb_ar,C');
% caxis([0 100]);
title('Population percentage of C');
xlabel('alpha value of epigenetic feedback from A-|B');
ylabel('alpha value of epigenetic feedback from C-|B');
pbaspect([1 1 1]);

%% Plotting the poulation percentage dynamics for a desired pair of values (alpha_1,alpha_2)

day=zeros(1,days+1);
for i=1:days+1
    day(i)=i-1;
end

f = 1; % alpha_1 value index
g = 1; % alpha_2 value index
pa(1,:) = PA_final(f,g,:);
pb(1,:) = PB_final(f,g,:);
pc(1,:) = PC_final(f,g,:);

figure()
hold on
grid(gca,'minor')
grid on
plot(day,pa,'LineWidth',3);
hold on
plot(day,pb,'LineWidth',3);
hold on
plot(day,pc,'LineWidth',3);
hold on
legend('A','B','C');
title('b0c b0a feedback alpha  0.1');
xlabel('Time'); 
ylabel('% population');
ylim([0 100]);