%PLEASE SUBMIT MODIFICATIONS AND IMPROVEMENTS!

%File: particle_diff_8.m (MATLAB)
%Version: 0.0
%Author: Stewart Nash
%Date: March 1, 2015
%Description: 2-Dimensional simulation of non-interacting particles originating in region with constant velocity. A reflecting boundary is present. A deflecting ring is present in the center.

%Note: Particle index means a consecutive list of integers starting with 1.

%>>>>Conversion Factors<<<<%
%cf_degree - Degree conversion factor
%cf_radian - Radian conversion factor

cf_degree = 360;
cf_radian = 2 * 3.14159265359; 
cf_second = 1;
cf_nanosecond = 10e-9;
cf_meter = 1;
cf_micrometer = 10e-6;
cf_micron = cf_micrometer;

%>>>>Configuration<<<<%

s_FRAMERATE = 25;
s_FLAGNUMBER = 1;

%>>>>Flags<<<<%

f_withinring = 1;

%>>>>Constants<<<<%
%DEGREES - Degree conversion factor
%RADIANS - Radian conversion factor

%MAX_PARTICLES - Maximum number of particles
%GRID_SIZE - Grid size for square grid
%TERMINAL_SPEED - Maximum particle speed
%TIME_STEPS - Number of time steps to run
%TIME_INCREMENT  - Number of time steps per frame (to advance particle position)
%DIMENSION - Dimension of space
%RING_RADIUS - Radius of central ring

C_PI =  3.14159265358979323846264338327950288419716939937510582;
DEGREES = 360;
RADIANS = 2 * C_PI;

MAX_PARTICLES = 50;
GRID_SIZE = 500 * cf_micrometer;
TERMINAL_SPEED = 100 * cf_micrometer / cf_nanosecond;
TIME_STEPS = 50000;
TIME_INCREMENT  = 1 * cf_nanosecond;
DIMENSION = 2;
RING_RADIUS = 100 * cf_micrometer;
WAVELENGTH = 0.75 * cf_micrometer;
cf_megahertz = 1;

temperature = 273.15;
molarmass = 0.08547;
rydberg = 8.315;

%Create video file and main figure window.
aviobj = avifile('example_08.avi', 'compression', 'None');
main_figure = figure;
set(main_figure, 'Position', [50, 50, 720, 880])

%Create a random number of particles and set size of grid (a 2-by-2 square matrix).
particlenumber = MAX_PARTICLES;
length = TERMINAL_SPEED;
domainsize = [ceil(GRID_SIZE / cf_micrometer), ceil(GRID_SIZE / cf_micrometer)];
domaincenter = domainsize / 2.0;
%Convert to nanoseconds
timeincrement = TIME_INCREMENT / cf_nanosecond;

%>>>>Variables<<<<%
%domain - Grid
%domain_2 - Grid which will be displayed, a copy of 'domain'
%particleflags - flags for particles indexed by particle number; flag description given by flag name
%particleposition - position of particles indexed by particle number
%particleangle - angle of particles velocity vector which is indexed by particle number
%particlespeed - speed of particles (magnitude of velocity vector) which is indexed by particle number
%particlevelocity - velocity of particles indexed by particle number

%ringradius - Holds value of the constant RING_RADIUS
%ringangle - Holds offset from zero radians that the field starts to vary sinusoidally in the theta direction
%ringfields - Holds maximum magnitude of coupling strength or field
%ringwavelength - Holds the wavelength of the field radiation
%ringwavevector - Holds the wavevector which describes the field radiation
%ringvalues - Holds the preceding ring values in an array

%x - Value for x-axis of a plot of coupling values in integer time steps
%y_current - Immediate value of coupling strength at individual time instances
%y_average - Average value of coupling strength at individual time instances
%counter - counter of current time step

%Allocate variables.
particleflags = zeros(particlenumber, s_FLAGNUMBER);
particleposition = zeros(particlenumber, DIMENSION);
particleangle = zeros(particlenumber, 1);
%particlespeed = zeros(particlenumber, 1);
particlevelocity = zeros(particlenumber, DIMENSION);

ringradius = RING_RADIUS / cf_micrometer;
ringangle = 0;
ringfields = 100 * cf_megahertz;
ringwavelength = WAVELENGTH / cf_micrometer;
ringwavevector = 1 / ringwavelength;
ringvalues = [ringradius, ringangle, ringfields, ringwavelength, ringwavevector];

x = 1 : TIME_STEPS;
y_current = zeros(TIME_STEPS, 1);
y_average = y_current;
counter = 0;

%Initialize variables to random values.
particleposition = rand(particlenumber, DIMENSION) * GRID_SIZE;
%Convert to micrometers
particleposition = particleposition / cf_micrometer;
%Set flags
%The f_withingring flag is set if the particle is within the radius of the ring from the center point.
for j = 1 : particlenumber
	if (sqrt((particleposition(j, 1) - domainsize(1) / 2)^2 + (particleposition(j, 2) - domainsize(2) / 2)^2) < ringradius)
		particleflags(j, f_withinring) = 1;
	end
end
particlespeed = (rand(particlenumber, 1) * TERMINAL_SPEED);
particlespeed = sqrt(-2 * rydberg * temperature * log(1 - particlespeed / length) / molarmass);
%Convert units to micrometers per nanosecond
particlespeed = particlespeed * cf_nanosecond / cf_micrometer;
particleangle = (rand(particlenumber, 1) * 2 * C_PI);
%Compute additive coupling strength
couplestrength = 0;
for j = 1 : particlenumber
    if particleflags(j, f_withinring) == 0
        k = [particleposition(j, 1), particleposition(j, 2)];
        temporary = ringfield(ringradius, ringangle, ringfields, ringwavelength, ringwavevector, k, domainsize);
        temporary = real(temporary);
        temporary = abs(temporary);
        couplestrength = couplestrength + temporary;
    end
end
counter = counter + 1;
y_current(counter) = couplestrength;
y_average(counter) = sum(y_current) / counter;

%Update velocity vector and note location of particles in domain.
particlevelocity = updatevelocityring(particlespeed, particleangle, particlenumber, particleflags);
%Display image and save to video file.
domain_image = ringimage(particlenumber, particleposition, domainsize, ringvalues, particleflags);
%mkdir('./results_08')
subplot(3, 1, [1 2])
imshow(domain_image)
subplot(3, 1, 3)
plot(x, y_current, 'y', x, y_average, 'b')
mystring = strcat('Current Coupling Strength: ', num2str(y_current(counter), 5),', Average Coupling Strength: ', num2str(y_average(counter), 5));
xlabel(mystring)
%filename = strcat('./results_08/', num2str(counter));
%print(filename, '-dpng')
aviobj = addframe(aviobj, main_figure);

%Iterate the previous process through a number of time steps.
for i = 1 : TIME_STEPS - 1
	[particleposition, particleangle] = updatepositioncirc(particleflags, particlespeed, particleangle, particleposition, timeincrement, particlenumber, ringradius, domainsize);
	%Update velocity as the particle angle may have been changed.
	particlevelocity = updatevelocityring(particlespeed, particleangle, particlenumber, particleflags);
	%Compute additive coupling strength
	couplestrength = 0;
    for j = 1 : particlenumber
        if particleflags(j, f_withinring) == 0
            k = [particleposition(j, 1), particleposition(j, 2)];
            temporary = ringfield(ringradius, ringangle, ringfields, ringwavelength, ringwavevector, k, domainsize);
            temporary = real(temporary);
            temporary = abs(temporary);
            couplestrength = couplestrength + temporary;
        end
    end
    counter = counter + 1;
    y_current(counter) = couplestrength;
    y_average(counter) = sum(y_current) / counter;
    if (mod(i, s_FRAMERATE) == 0)
        domain_image = ringimage(particlenumber, particleposition, domainsize, ringvalues, particleflags);
		subplot(3, 1, [1 2])
		imshow(domain_image)
		subplot(3, 1, 3)
		plot(x, y_current, 'y', x, y_average, 'b')
		mystring = strcat('Current Coupling Strength: ', num2str(y_current(counter), 5),', Average Coupling Strength: ', num2str(y_average(counter), 5));
		xlabel(mystring)
%		filename = strcat('./results_08/', num2str(counter));
%		print(filename, '-dpng')
    	aviobj = addframe(aviobj, main_figure);
    end
end

%Close video file.
aviobj = close(aviobj);
