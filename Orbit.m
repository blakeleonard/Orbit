% Blake Leonard	2009

% University of Missouri
% Computational Physics


% Orbit - Program to compute the orbit of a comet using four separate methods.  Completes after one full period.

% print header

fprintf('Orbit - Program to compute the orbit of a comet, completes after 1 full period\n\n');


%* Set initial position and velocity of the comet.

r0 = input('Enter initial radial distance (AU): ');  

v0 = input('Enter initial tangential velocity (AU/yr): ');

r = [r0 0];  v = [0 v0];

state = [ r(1) r(2) v(1) v(2) ];   % Used by R-K routines


%* Set physical parameters (mass, G*M)

GM = 4*pi^2;      % Grav. const. * Mass of Sun (au^3/yr^2)

mass = 1.;        % Mass of comet 

adaptErr = 1.e-3; % Error parameter used by adaptive Runge-Kutta

time = 0;

signcount = 0; oldth = 0;

maxvmag = 0; minvmag = 1000000000;

TotKin = 0; TotPot = 0; AvgKin = 0; AvgPot = 0; period = 0;


%* Loop until 1 orbit complete using specified numerical method.

nStep = 1000000000;

tau = input('Enter time step (yr): '); 

NumericalMethod = menu('Choose a numerical method:', ...
       'Euler','Euler-Cromer','Runge-Kutta','Adaptive R-K');

for iStep=1:nStep  

  %* Record position and energy for plotting.

  rplot(iStep) = norm(r);           % Record position for polar plot

  thplot(iStep) = atan2(r(2),r(1));

  tplot(iStep) = time;

  kinetic(iStep) = .5*mass*norm(v)^2;   % Record energies

  potential(iStep) = - GM*mass/norm(r);

  TotKin = TotKin + kinetic(iStep);

  TotPot = TotPot + potential(iStep);	


  if (thplot(iStep)*oldth < 0)      % Break loop when sign has changed twice, one orbit has completed

	signcount = signcount + 1;

	if (signcount == 2)

		period = time;          % Get time for one orbit    

		AvgKin = TotKin/iStep;  % Compute Average Kinetic Energy

		AvgPot = TotPot/iStep;  % Compute Average Potential Energy

		break;
	end
  end

  oldth = thplot(iStep);
  

  %* Calculate new position and velocity using desired method.

  if( NumericalMethod == 1 )

    accel = -GM*r/norm(r)^3;  
 
    r = r + tau*v;             % Euler step

    v = v + tau*accel; 

    time = time + tau;  

 
  elseif( NumericalMethod == 2 )

    accel = -GM*r/norm(r)^3;  
 
    v = v + tau*accel; 

    r = r + tau*v;             % Euler-Cromer step

    time = time + tau;   

 
  elseif( NumericalMethod == 3 )

    state = rk4(state,time,tau,'gravrk',GM);

    r = [state(1) state(2)];   % 4th order Runge-Kutta

    v = [state(3) state(4)];

    time = time + tau;   


  else

    [state time tau] = rka(state,time,tau,adaptErr,'gravrk',GM);

    r = [state(1) state(2)];   % Adaptive Runge-Kutta

    v = [state(3) state(4)];

  end
  

  if ( norm(v) > maxvmag )        % If velocity is max, point is perihelion

	maxvmag = norm(v);

	q = norm(r);                % Get distance from sun at perihelion

  end


  if ( norm(v) < minvmag )        % If velocity is min, point is aphelion

	minvmag = norm(v);

	Q = norm(r);                % Get distance from sun at aphelion

  end


end

a = (q + Q)/2;                 % compute semi-major axis

e = - (q/a - 1);               % compute eccentricity

L = mass*v0*r0;                % Compute Angular Momentum

E = AvgKin + AvgPot;           % Compute Total Average Energy


theor_e = (1 + (2*E*L^2)/(GM^2*mass^3)) ^(.5);   % compute e from book equation

e_error = abs(100 * (theor_e - e)/theor_e);

kepler_a = (((GM)/(4*pi^2))*period^2) ^ (1/3);   % compute a from kepler’s law

a_error = abs(100 * (kepler_a - a)/kepler_a);

virial_error = abs(100 * ((AvgKin + AvgPot/2)/AvgKin)); % compute adherence to virial theorem


fprintf('Period (yr): %g \n', period);
fprintf('Semi-Major Axis (AU): %g \n', a);
fprintf('Semi-Major Axis derived from Kepler’s 3rd Law (AU): %g \n', kepler_a);
fprintf('Percent Error in Kepler’s 3rd Law: %g \n', a_error);
fprintf('Perihelion Distance (AU): %g \n', q);
fprintf('Eccentricity: %g \n', e);
fprintf('Theoretical Eccentricity: %g \n', theor_e);
fprintf('Percent Error in eccentricity: %g \n', e_error);
fprintf('Percent Adherence to Virial Theorem: %g \n', virial_error);


%* Graph the trajectory of the comet.

figure(1); clf;  % Clear figure 1 window and bring forward
polar(thplot,rplot,'+');  % Use polar plot for graphing orbit
xlabel('Distance (AU)');  grid;
pause(1)   % Pause for 1 second before drawing next plot

z = input('Hit Enter for Next: ');


%* Graph the energy of the comet versus time.

figure(2); clf;   % Clear figure 2 window and bring forward
totalE = kinetic + potential;   % Total energy
plot(tplot,kinetic,'-.',tplot,potential,'--',tplot,totalE,'-')
legend('Kinetic','Potential','Total');
xlabel('Time (yr)'); ylabel('Energy (M AU^2/yr^2)');

z = input('Hit Enter to End: ');
