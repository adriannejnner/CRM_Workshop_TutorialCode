function [controlmodel_S] = simulated_controlmodel(r,K,S0,data)

initialcondition = [S0];
sol = ode45(@controlmodel,[0 data.time(end)],initialcondition);
controlmodel_S = deval(sol,data.time,1);

function dydt = controlmodel(t,y,Z)
   S = y(1);

   dS = r*S*(1-S/K);

   dydt = [dS];

end

end