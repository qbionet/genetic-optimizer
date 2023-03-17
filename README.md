# genetic-optimizer

This repository provides the MATLAB code for implementing and deploying the optimizer developed and featured in _A blueprint for a synthetic genetic feedback optimizer_. 
* Dynamical models featured in the manuscript are included in the folder **ODE models**. 
* Parameters can be set by running the files ==optimizer_parameters.m==, growth_parameters.m, and process_parameters.m, all located in the folder **Parameters**. 
* For custom application examples, the Simulink file ClosedLoop_custom can be modified by creating an Sfunction (custom_Sfun) containing the dynamics of the regulated process, alongside with the parameters defined in process_parameters.m (sample examples are provided in ClosedLoop_FFL.slx and ClosedLoop_reporter.slx). 
* Alternatively, one may include the closed loop dynamics in wrapper_ODE.m by including the ODE models of the regulated process directly together with the optimizer.
