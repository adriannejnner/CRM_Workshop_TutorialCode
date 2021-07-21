# CRM_Workshop_TutorialCode

To demonstrate the generation and use of a virtual clinical trial, we will consider how the interaction between a cancer drug and the immune system. Let S(t) be the number of cancer cells at day t,  and T(t) be the number of immune cells. Cancer cells are proliferating at a logistic rate with proliferation constant r and carrying capacity K. Cancer cells also undero apoptosis through contact with immune cells at a rate κ. We assume immune cells are recruited at a rate proportional to the amount of tumour cells at rate a and decay at a rate d. Consider the concentration of drug [drug] which is an anti-cancer drug (such as chemotherapy) is introduced into the system. This drug has an ability to kill cancer cells at a mckelis menton rate with half effect η and rate δ. If we assume the drug decays from the system at arate ρ, this gives the following model

dS/dt=rS(1-S/K)-κST-δ[drug]/([drug]+η) S,
dT/dt=aS-dT,
d[drug]/dt=-ρ[drug],

A schematic summary of the system can be found in Fig. 1(a). 
To generate the virtual cohort, we first take parameter values that match to data. We then sample parameters for N virtual patients from a normal distribution centered at μ with variance σ, rejecting any negative parameters. To only include patients with realistic tumour growths, we confirm that these patient’s parameters result in tumour growth which is within three standard deviations of the data mean. Accepted patients are then considered part of the cohort and any patients whose dynamics do not lie within the physiological ranges are rejected and parameters resampled. 
	With this patient distribution we can then examine the impact of different dosage protocols. Consider an initial concentration of drug C_0 is administer in N injections given τ days apart. We can simulate the impact of varying the number of injections between 1 to 15 for each patient with possible days apart from 1 to 7. 
	
