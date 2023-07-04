# Introduction #
Our code uses Radiative transfer simulation to fit the SED and ALMA continuum images to obtain the best-fit parameter of the Protoplanetary disk.

# Flow chart #
## Radmc3d_Simulation.py
1. Set the model parameters values(amin,amax,grainP,grainNum,grainLam,ifkappa) for preparing dust opacity files
2. Set the grid of model space
3. Set your model parameters values (flaring, hGas100, hGasLocation, mDisk, dustToGas, turbulence ) for preparing the input file of the dust density
4. Set auxiliary files for radmc3d(radmc.inp,dust_opac.inp,wavelength_micron.inp,stars.inp)
5. Run the thermal simulation to get the dust temperature distribution (radmc3d mctherm)
6. Simulate the continuum image, and compare the model brightness profile with the observed brightness profile
7. Update the dust surface density (see section 4.2 in Li+ 2023MNRAS.518.6092L)
8. Obtain the converaged dust surface density
9. Simulate the finall SED and continuum image, and calculate the chi2_SED and chi2_image

## Radmc3d_Simulation_mult.py
Multi-process programs can be utilized to concurrently determine the optimal fitting results for multiple models within the parameter space.

# Requirements #

## python env ##
python3.X

## package ##
Our code contains a large amount of code from the DHSHARP and RADMC3D packages.

DSHARP can be downloaded through:
        
        pip3 install dsharp_opac

Radmc3d can be downloaded through:
        https://github.com/dullemond/radmc3d-2.0
        
manual of RADMC-3D can be found:
        https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html
        
# How to use our code #
## 0.Example ##
You can execute the Radmc3d_Simulation.py file.

The DS_Tau_b6avgf.dat file represents the radial distribution of observed flux, while the DS_Tau.txt file represents the observed Spectral Energy Distribution (SED).

## 1.Configuring the environment: This determines the program's working directory ##
Related parameters: pathIn

Location of the folder for dust_kappa_x.inp files

Tip: pathIn must end in '/'

Related parameters: pathOut
Location of the folder for output files

The observation files for sed and image should be placed in the pathOut location.

Tip: pathOut must end in '/'

The SED observation file consists of three columns. The first column represents the wavelength of each observation [um], the second column contains the corresponding measured values at each wavelength [Jy], and the third column represents the measurement errors [Jy].

The radial profile of flux observation file consists of three columns. The first column represents the distance to the young star [AU], the second column contains the corresponding measured values at each radii [mJy/beam], and the third column represents the measurement errors [mJy/beam].
                  
## 2.Set parameters ##
There are four primary parameter categories: dust properties, grid settings, disk characteristics, and parameters for calculating the Spectral Energy Distribution (SED) and continuous emission images of the dust.
    
### Grain
amin: Minimum size of dust [cm]

amax: Maximum size of dust [cm]

grainP: -3.5

grainNum: The quantity of dust in logarithmic space 

grainLam: The wavelength range of the absorption and scattering coefficients of dust
    
### Grid
rin,rout: The maximum and minimum values of grid points on the radial direction of the disk. This value can exceed the inner and outer radii of the disk, but our program sets the range of radial grid points to be the same as the inner and outer radii of the disk.

x(y,z)bound: Range in the x(y,z) direction. 

nx(y,z): The number of grid points in each dimension. Corresponds to x (y, z) bound, which can be a list. Using a list allows for individual specification of the number of grid points within a given range.
    
### Disk ###
hGas100,hGasLocation: The scale height of the gas at a specific location is determined. The default value is set as 0.1 and 100, representing a scale height of 10 au at a distance of 100 au.

flaring: flaring index of the disk. 

**Please note that when the value is set to 0.1, it corresponds to a flaring index of 1.1.**

surfaceDensityP: The power-law exponent of surface density is less significant during the fitting of continuous dust images, as the density can be adjusted accordingly.

mDisk: This parameter represents the combined mass of the protoplanetary disk, encompassing both the dust and gas components.

dustToGas: The parameter signifies the ratio between the mass of dust and the mass of gas within the protoplanetary disk, with a default value set at 0.01.

turbulence: turbulence 
    
### SED and Image ###

SEDLam: This parameter specifies the wavelength(s) used for calculating the Spectral Energy Distribution (SED). The comparison between the observed SED and the model-calculated SED during iterations determines the updated disk mass. 

seds,flux: These two parameters refer to the observed SED and the observed flux extracted radially along the disk.   

incl, PA, dpc:  The inclination, position angle, and distance of the protoplanetary disk.

pixelNum, sizeau, pixelSize: During the fitting process in the code, a set of fitted images will be generated. These three parameters correspond to the number of pixels in the image, the total size of the image in astronomical units (AU), and the size of each pixel in AU.

imageLam,beamMajor,beamMinor: The band of the generated image and the size of the beam.

**Please note that:**

**1. Multiple band images can be generated, but it is important to ensure that the number of elements in the imageLam, beamMajor, and beamMinor lists are consistent and matching.**

**2. The fitted wavelength(s) should be placed as the last element(s) in the list.**
                    
### Star and wavelength ###
see Radmc3d manual 

Teff: Effective temperature of young stars.

mstar: The Mass of Young Stars

rstar: The radius of a young star.

wavelength:  The wavelength of light emitted by young stars 
    
        
**Please take note that these parameters will generate two files: star.inp and wavelength_micron.inp. If you plan to modify the star.inp file, it is crucial to ensure that the wavelength specified in that file matches the wavelength in the wavelength_micron.inp file. It is important to maintain consistency between the two files; otherwise, radmc3d may encounter errors during execution.**

    
### radmc.inp control file parameter
The following is a list of the parameter settings we used in radmc3d 

These parameters are from the radmc3d package, and I will copy the meaning of the parameters here.

nphot = 10000000: The number of photon packages used for the thermal Monte Carlo simulation. 

nphot_scat = 10000000: The number of photon packages for the scattering Monte Carlo simulations, done before image-rendering.

nphot_spec = 200000: The number of photon packages for the scattering Monte Carlo simulations, done during spectrum-calculation. This is actually the same functionality as for nphot_scat, but it is used (and only used) for the spectrum and SED calculations. The reason to have a separate value for this is that for spectra you may not need as many photon packages as for imaging, because you anyway integrate over the images. Many of the annoying ‘stripe noise’ in images when using insufficiently large nphot_scat will cancel each other out in the flux calculation. So nphot_spec is usually taken smaller than nphot_scat.

scattering_mode_max = 1: When radmc3d reads the dust opacity files it checks if one or more of the opacity files has scattering opacity included. If yes, the scattering_mode will automatically be set to 1. It will also check if one or more includes anisotropic scattering. If yes, the scattering_mode will automatically be set to 2. But the user may nevertheless want to exclude anisotropic scattering or exclude scattering altogether (for instance for testing purposes, or if the user knows from experience that the scattering or anisotropic nature of scattering is not important for the problem at hand). Rather than editing the opacity files to remove the scattering and/or Henyey-Greenstein g-factors, you can limit the value that radmc3d is allowed to make scattering_mode by setting the variable scattering_mode_max. If you set scattering_mode_max=0 then no matter what opacity files you have, scattering will not be treated. If you set scattering_mode_max=1, then no matter what opacity files you have, scattering will be treated in an isotropic way.

camera_min_dr = 0.1 [Fine-tuning only]: Fine-tuning parameter for recursive subpixeling, for spherical coordinates, assuring that not too fine subpixeling would slow down the rendering of images or spectra too much.

istar_sphere = 1: If 0 (=default), then all stars are treated as point-sources. If 1, then all stars are treated as finite-size spheres. This mode is more accurate and more realistic, but the applications are a bit more restricted. Such finite-size stars are (for technical reasons) not always allowed anywhere in the model. But for problems of circumstellar disks and envelopes in spherical coordinates, it is recommended to set this to 1. Typically, if a star is outside the grid (in spherical coordinates this can also be at the origin of the coordinate system, as long as the inner radius of the coordinate system is larger than the stellar radius!) the use of the finite-size star mode is always possible. But if the star is on the grid, there are technical limitations.

setthreads:  Number of CPU cores used simultaneously.
    
## 3.Function  ##
createDSHARPGrain: Generate a customized chemical composition (DHSHARP), a specific number, and a designated size of dust particles. Consuming a lot of time. 

**If you want to skip this step, please set ifkappa to 0 and ensure that there are kappa files in the pathIn folder.**

               
calculateDustMassWeight: Allocate the disk quality among individuals based on the size of the dust particles.

writeSpatialGrid:  Create grid and generate amr_grid.inp file.

calculatePho: Generate a density distribution of the protoplanetary disk with a vertical layered structure, organized according to the size of dust particles as per the specified settings. 

writeDustDensity: Generate density.inp file.

createSedLam: Specify the radmc3d code to calculate the wavelength of the SED and generate the camera_wavelength_micron.inp file. 

makeRadSet,makeOpacity,makeStarandWavelength: Generate radmc.inp, dustOpac.inp, stars.inp and wavelength_micron.inp.

iteration: Main Iteration Code

calculateImage: Generate the image.out and .fits files

extractedFlux: Compute the average distribution of flux along the radial direction based on the angle of the principal axis direction.

calculateDustMass: Calculate the mass of dust in the disk. The principle is to calculate the mass of dust from the data in the density. inp file.

calculateSED: Calculate SED and generate spectrum.out and sed.txt files.

## 4.Output ##
    
There are many generated files, so if you don't want to generate a certain file, you can comment it out.
Saved to the current folder or 'iteration_warehouse' folder. 

Files for each step(iteration_warehouse):

1.dust_density_x.inp

2.flux_R_x_x.txt

3.fx_x:  fitting coefficient 

4.time_section_x

The last step of the flux and sed files:

1.flux_R_x.txt

2.sed_x.txt


Chi2 with observed values:

1.ChiFlux.txt

2.ChiSed.txt
    
# multiprocess #

If you want to run a multiprocess program to search for parameter space, you can use Generate_Sigma_mult.py program. Need a multiprocess package, and the download method is:
            
            pip3 install multprocess
            


**Please note that multi process programs should be set based on the number of CPU cores. The required number of cores is poolNum * setthreads. For example, if the parameters are pool (5) and setthread = 8, the required cores are 5 * 8=40.**


    








 





        
