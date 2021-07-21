function [Tumour_burden_day20] = cohort_response_model(x,p,patients,data)

% recall genetic algorithm evaluates model on x parents each time, so x is
% a vector of parent initial dosages to try

for i = 1:length(x)
    
    % we then evaluate each patients tumour size for each parent initial dose x
    for j = 1:length(patients)                              

        r = patients(j,1);
        K = patients(j,2);
        kappa = patients(j,3);
        a = patients(j,4);
        d = patients(j,5);
        
        initial_drug = 100;
        day_of_second_dose = x(i);
        [time, model_S,model_T,model_drug] = simulated_fullmodel_with_drug(r,K,kappa,a,d,p,data,initial_drug,day_of_second_dose);
        Tumour_burden_day20_patients(j) = model_S(end);                 % fullmodel_S is a vector of tumour volumes, we want the last tumour volume on day 20
        
    end
    Tumour_burden_day20(i) = sum(Tumour_burden_day20_patients)
end
end