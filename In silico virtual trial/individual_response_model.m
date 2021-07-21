function [Tumour_burden_day20] = individual_response_model(x,p,patients,patient_number,data)

% recall genetic algorithm evaluates model on x parents each time, so x is
% a vector of parent initial dosages to try

% we set the parameter values for which patient we are optimising
r = patients(patient_number,1);
K = patients(patient_number,2);
kappa = patients(patient_number,3);
a = patients(patient_number,4);
d = patients(patient_number,5);

initial_drug = 100;
    
for i = 1:length(x)
    
    day_of_second_dose = x(i);
    [time, model_S,model_T,model_drug] = simulated_fullmodel_with_drug(r,K,kappa,a,d,p,data,initial_drug,day_of_second_dose);
    Tumour_burden_day20_patients = model_S(end);                 % fullmodel_S is a vector of tumour volumes, we want the last tumour volume on day 20

    Tumour_burden_day20(i) = Tumour_burden_day20_patients;
end

end