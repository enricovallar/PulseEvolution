# PulseEvolution

This is a thesis project by Vallar Enrico. 
The thesis supervisors are prof. M. Santagiustina and prof. A. Galtarossa.
University of Padua, DEI, Bachelor's degree in electronic engineering.

It consists in a MATLAB app called PulseEvolution made with App Designer.

PulseEvolution simulates the propagation of pulses in optical fibers 
by solving the NLSE using the Split Step Fourier Method. 
A GUI allows you to easily configure the simulated fiber and pulse parameters.
Furthermore, the application allows you to study the evolution 
of the pulse at different distances from the origin, 
displaying the trend of the power and other useful graphs.
It is possible to generate optical solitons or to compare up to three pulses.

The main function used in this program is based on the Split Step Fourier Method presented in "Nonlinear Fiber Optics" by Govind P. Agrawal.

I have uploaded the MATLAB code of PulseEvolution in this repo.

-Fiber.m is a MATLAB class used to simulate optical fibers.

-Pulse.m is a MATLAB class used to simulate pulses which propagate in optical fibers and whose evolution is governed by the NLSE.
 
-holdAX.m is a MATLAB function used to speed up axis manipulation.

-PulseEvolution.mlapp is PulseEvolution app code. Part of this code has been automatically written by MATLAB AppDesigner.
 
-PulseEvolution_3.8.mlappinstall is the installer of PulseEvolution v.3.8. MATLAB is required.
