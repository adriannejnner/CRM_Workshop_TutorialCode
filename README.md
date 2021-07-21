To demonstrate the generation and use of a virtual clinical trial, we will consider a model for the interaction between a cancer drug, a tumour population and the immune system. Let S(t) be the number of cancer cells at day t and T(t) be the number of immune cells. Cancer cells are proliferating at a logistic rate with proliferation constant r and carrying capacity K. Cancer cells also undergo apoptosis through contact with immune cells at a rate κ. We assume immune cells are recruited at a rate proportional to the amount of tumour cells at rate a and decay at a rate d. Consider the concentration of drug, [drug], which is an anti-cancer drug (such as chemotherapy). The drug is introduced into the system both initially and at time I_D. This drug has an ability to kill cancer cells and has a half effect η and rate δ. If we assume the drug decays from the system at a rate ρ, then all these assumptions give the following system of ODEs:

dS/dt=rS(1-S/K)-κST-δ[drug]/([drug]+η) S, 
dT/dt=aS-dT, 
d[drug]/dt=-ρ[drug] 
A schematic summary of the system can be found in Fig. 1.  

To generate the virtual cohort, we first take parameter values that match to data for the control tumour growth (immune-absent system) and tumour growth with the presence of the immune system. This is achieve by setting the appropriate parts of the model to zero (i.e. the drug and the immune population in the immune-absent case) and fitting the parameters in the model simultaneously. We fit r,K,κ,a and d using least-squares non-linear fitting and the parameters can be found in Table 1. 
To create the virtual patients, we then sampled parameters for N virtual patients from a normal distribution centered at μ with variance σ, rejecting any negative parameters. We took μ to be the set of parameters obtained from fitting. To only include patients with realistic tumour growths, we confirm that these patient’s parameters result in tumour growth which is within three standard deviations of the data mean (this involved simulating the model in the absence of the drug). Accepted patients were then considered part of the cohort and any patients whose dynamics did not lie within the physiological ranges were rejected and parameters resampled. 
	With this patient distribution we then examined the impact of different dosage protocols. Consider an initial concentration of drug C_0 is administer initially and on day I_D, we used a genetic algorithm to determine first the optimal day for administration for the whole cohort and then the optimal day of administration for each individual patient. We found that on the whole the therapy is robust when administer on day 9.  
	
![image](https://user-images.githubusercontent.com/48768705/126488377-68c1ed7b-30bd-4a59-8f0e-b07c3b13f254.png)

