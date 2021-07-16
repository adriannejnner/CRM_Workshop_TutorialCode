%% Loading and plotting the data

load('Tumour_growth_data.mat') %load in experiment data

%set up matlab variable for the data
data.time = experiment_time;                                               %experiment time points
data.control = control_tumour_volume_mean;                                 %control experiment tumour volume mean  
data.std_control = control_tumour_volume_std;                              %control experiment tumour volume standard deviation
data.treatment = treatment_tumour_volume_mean;                             %immune experiment tumour size mean  
data.std_treatment = treatment_tumour_volume_mean;                         %immune experiment tumour volume standard deviation

%plot control fit
figure
hold on 
errorbar(data.time,data.control,data.std_control,'LineWidth',2)
xlabel('Time (days)')
ylabel('Tumour volume')
legend('Tumour volume measurement')
set(gca,'FontSize',16)
title('Control')

%plot immune experiment fit
figure
hold on 
errorbar(data.time,data.treatment,data.std_treatment,'Color',[0.51, 0.78, 0.95],'LineWidth',2)
ylabel('Tumour volume')
xlabel('Time (days)')
set(gca,'FontSize',16)
ylim([0 4500])
legend('Tumour volume measurement')
title('Immune treatment')

%% running one simulation of the model with parameter values and comparing to the data

%set parameter values
r = 0.2629;                                                       %tumour cell replication rate
K = 7.6252e+03;                                                   %tumour cell carrying capacity
kappa = 0.0026;                                                   %T cell tumour cell killing rate
a = 7.8451;                                                       %T cell activation rate
d = 189.4870;                                                     %T cell death rate

S0 = 100;                                                                   %initial tumour volume
T0 = 0;                                                                  %initial number of T cells

%simulate control experiment model
[controlmodel_S] = simulated_controlmodel(r,K,S0,data); 

%simulate full model    
[fullmodel_S] = simulated_fullmodel(r,K,kappa,a,d,S0,T0,data);
    
%plot control model vs data
figure
hold on 
errorbar(data.time,data.control,data.std_control,'LineWidth',2)
plot(data.time,controlmodel_S,'k:','LineWidth',2)
xlabel('Time (days)')
ylabel('Tumour volume')
legend('Tumour volume measurement','Tumour cells, S(t)')
set(gca,'FontSize',16)
title('Control')

%plot immune experiment model vs data
figure
hold on 
errorbar(data.time,data.treatment,data.std_treatment,'Color',[0.51, 0.78, 0.95],'LineWidth',2)
plot(data.time,fullmodel_S,'k:','LineWidth',2)
ylabel('Tumour volume')
xlabel('Time (days)')
set(gca,'FontSize',16)
ylim([0 4500])
legend('Tumour volume measurement','Tumour cells, S(t)')
title('Immune treatment')

%% creating virtual patients

% parameter distribution means
r_mu = r;
K_mu = K;
kappa_mu = kappa;
a_mu = a;
d_mu = d;

%parameter distribution variances
r_sig = 0.1*r;
K_sig = 0.1*K;
kappa_sig = 0.1*kappa;
a_sig = 0.1*a;
d_sig = 0.1*d;

%initial number of patients
N = 50;                                                                    

%initialising the counter for the virtual patient loop
M = 1;

%setting up figures to plot patient control and immune growths
fig1 = figure; 
fig2 = figure;
while M<N
    
    %sample from the normal distributions
    r = normrnd(r_mu,r_sig);
    K = normrnd(K_mu,K_sig);
    kappa = normrnd(kappa_mu,kappa_sig);
    a = normrnd(a_mu,a_sig);
    d = normrnd(d_mu,d_sig);
    
    %simulate control experiment model
    [controlmodel_S] = simulated_controlmodel(r,K,S0,data);
    g_control = (controlmodel_S - (data.control+3*data.std_control+data.control-3*data.std_control)./2).^2-(data.control+3*data.std_control-(data.control+3*data.std_control+data.control-3*data.std_control)./2).^2;
    
    %simulate full model
    [fullmodel_S] = simulated_fullmodel(r,K,kappa,a,d,S0,T0,data);
    g_treatment = (fullmodel_S - (data.treatment+3*data.std_treatment+data.treatment-3*data.std_treatment)./2).^2-(data.treatment+3*data.std_treatment-(data.treatment+3*data.std_treatment+data.treatment-3*data.std_treatment)./2).^2;
    g = [g_control; g_treatment]';
    
    %check that patient tumour growths lie within data standard deviation
    if isempty(find(g>0))==1
        
        patients(M,:) = [r,K,kappa,a,d];                                    % if patient growth satisfies requirement, add their parameters to the patient list
        M = M+1;
        
        %plot patient growth
        figure(fig1)
        hold on 
        plot(data.time,controlmodel_S,':','Color',[0.5 0.5 0.5])
        
        figure(fig2)
        hold on
        plot(data.time,fullmodel_S,':','Color',[0.5 0.5 0.5])
    end
        
        fprintf('Patient #%d\n', M);
end

% label x and y axis of fig1 and fig 2 and plot data over the top
figure(fig1)
hold on 
errorbar(data.time,data.control,data.std_control,'LineWidth',2)
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)
title('Control patient dynamics')

figure(fig2)
hold on
errorbar(data.time,data.treatment,data.std_treatment,'Color',[0.51, 0.78, 0.95],'LineWidth',2)
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)
title('Patient dynamics under treatment')

%% Run a genetic algorithm to determine the optimial initial drug concentration for the cohort

% initialise drug related parameters (PAC-1)
p.delta = 0.11;                                                             % killing rate of tumour cells by drug
p.eta = 33;                                                                 % half-effect of drug
p.rho = 0.02;                                                               % decay rate of drug

% run a genetic algorithm to determine optimal initial drug for the cohort
[optimalparam]  = ga(@(x)cohort_response_model(x,p,patients,data),1,[],[],[],[],[1],[100]);

for j = 1:length(patients)                              

    r = patients(j,1);
    K = patients(j,2);
    kappa = patients(j,3);
    a = patients(j,4);
    d = patients(j,5);

    initial_drug = optimalparam;
    [fullmodel_S] = simulated_fullmodel_with_drug(r,K,kappa,a,d,p,data,initial_drug);
    Tumour_burden_day20_patients(j) = fullmodel_S(end);                 % fullmodel_S is a vector of tumour volumes, we want the last tumour volume on day 20

end

figure
hist(Tumour_burden_day20_patients)
xlabel('Tumour volume day 20')
set(gca,'FontSize',16)









