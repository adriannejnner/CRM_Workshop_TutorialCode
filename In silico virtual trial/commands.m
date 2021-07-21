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
r_sig = 0.2*r;
K_sig = 0.2*K;
kappa_sig = 0.2*kappa;
a_sig = 0.2*a;
d_sig = 0.2*d;

%initial number of patients
N = 30;                                                                    

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
l1 = plot(data.time,controlmodel_S,':','Color',[0.5 0.5 0.5]);
l2 = errorbar(data.time,data.control,data.std_control,'Color',[0, 0.45, 0.75],'LineWidth',2);
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)
title('Control patient dynamics')
legend([l1 l2],'Inidividual Patient','Data')

figure(fig2)
hold on
l1 = plot(data.time,fullmodel_S,':','Color',[0.5 0.5 0.5]);
l2 = errorbar(data.time,data.treatment,data.std_treatment,'Color',[0.51, 0.78, 0.95],'LineWidth',2);
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)
title('Patient dynamics under treatment')
legend([l1 l2],'Inidividual Patient','Data')

%% Simulate chemotherapy treatment of patient 1

% initialise drug related parameters (PAC-1)
p.delta = 0.5;                                                             % killing rate of tumour cells by drug
p.eta = 50;                                                                 % half-effect of drug
p.rho = 0.2;                                                               % decay rate of drug

initial_drug = 100;                                                         % Initial concentration of drug

day_of_second_dose = 1;                                                     % Day the second dose is administered

%Simulate the tumour immune model with the chemotherapy drug
[time, chemomodel_S, chemomodel_T, chemomodel_drug] = simulated_fullmodel_with_drug(r,K,kappa,a,d,p,data,initial_drug,day_of_second_dose);

% Plotting the outputs of the model
figure
hold on
l1 = plot(time,chemomodel_S,'LineWidth',2);
l2 = plot(time,chemomodel_T,'LineWidth',2);
l3 = plot(time,chemomodel_drug,'LineWidth',2);
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)
title('Patient dynamics under treatment')
legend([l1 l2 l3],{'Tumour cells','Immune cells','Drug'})

%% Run a genetic algorithm to determine the optimial day for the second dose for the cohort

% run a genetic algorithm to determine optimal initial drug for the cohort
rng default % For reproducibility
[optimalparam]  = ga(@(x)cohort_response_model(x,p,patients,data),1,[],[],[],[],[1],[20]);

%set the day of the second dose as the optimal parameter determined by the
%genetic algorithm
day_of_second_dose = optimalparam;

%Simulate the cohort under this optimal protocol and plot the individual
%patient dynamics
fig3 = figure;
hold on
for j = 1:length(patients)                              

    r = patients(j,1);
    K = patients(j,2);
    kappa = patients(j,3);
    a = patients(j,4);
    d = patients(j,5);

    [time,fullmodel_S] = simulated_fullmodel_with_drug(r,K,kappa,a,d,p,data,initial_drug,day_of_second_dose);
    Tumour_burden_day20_patients(j) = fullmodel_S(end);                 % fullmodel_S is a vector of tumour volumes, we want the last tumour volume on day 20

    plot(time,fullmodel_S,':','LineWidth',1,'Color',[0.75 0.5 0.5])
    
end

figure(fig3)
ylabel('Tumour volume (mm^3')
xlabel('Time (days)')
set(gca,'FontSize',18)

% Plot a histogram of the tumour volume of the cohort on day 20
figure
hist(Tumour_burden_day20_patients)
xlabel('Tumour volume day 20')
ylabel('Frequency')
title('Histogram of cohort tumour size day 20')
set(gca,'FontSize',18)

%% Run a genetic algorithm to determine the optimial day for the second dose for each patient


%loop through the patients and determine each patients optimal day of
%second dosage
for j = 1:length(patients)                              

    patient_number = j;
    
    rng default % For reproducibility
    [optimalparam]  = ga(@(x)individual_response_model(x,p,patients,patient_number,data),1,[],[],[],[],[1],[20]);
    
    optimal_day_vec(j) = optimalparam;
    
    patient_number
end

[sorted_optimal_day_vec, sort_index] = sort(optimal_day_vec);

%setting up color map
coloring_matrix=cbrewer('seq', 'PuBu', 9);
rgrid = linspace(min(patients(:,1)),max(patients(:,1))+1e-8,10);

%Plotting each patient's optimal day for the second injection and colouring
%the bar by the patients tumour proliferation rate r
figure
b = bar(sorted_optimal_day_vec);
b.FaceColor = 'flat';
for j = 1:length(patients)
    
    %determine which colour to plot the bar based on patients parameter r
    patient_rval = patients(sort_index(j),1);
    color_index = find(rgrid>patient_rval,1)-1;
    
    %changing the colour of the bar to match the r value
    b.CData(j,:) = coloring_matrix(color_index,:);
    
end
ylabel('Optimal day of second injection')
xlabel('Patient number')
set(gca,'FontSize',18)
colormap(coloring_matrix)
c = colorbar;
set(c,'ticklabels',{'0.2026','0.2047','0.2069','0.2090','0.2111','0.2132'})
title('Individualised protocol')
ylim([8 10])








