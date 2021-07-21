function [time, model_S,model_T,model_drug] = simulated_fullmodel_with_drug(r,K,kappa,a,d,p,data,initial_drug,day_of_second_dose)

S0 = 100;                                                                   %initial tumour volume
T0 = 0;                                                                  %initial number of T cells

initialcondition = [S0,T0,initial_drug];
sol1 = ode45(@fullmodel,[0 day_of_second_dose],initialcondition);
sol2 = ode45(@fullmodel,[day_of_second_dose 20],sol1.y(:,end)+[0; 0; initial_drug]);

model_S = [sol1.y(1,:),sol2.y(1,:)];
model_T = [sol1.y(2,:),sol2.y(2,:)];
model_drug = [sol1.y(3,:),sol2.y(3,:)];

time = [sol1.x,sol2.x];

function dydt = fullmodel(t,y,Z)
   S = y(1);
   T = y(2);
   drug = y(3);
   
   dS = r*S*(1-S/K)-kappa*S*T-p.delta*drug/(drug+p.eta)*S;
   dT = a*S-d*T-p.delta*drug/(drug+p.eta)*T;
   ddrug=  -p.rho*drug;

   dydt = [dS;dT;ddrug];

end

end