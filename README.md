# Introduction #
Our code is developed to fit the SED and ALMA continuum images of protoplanetary disks.

# Flow chart #
## Radmc3d_Simulation.py
1. Set the model parameter values (for amin, amax, grainP, grainNum, grainLam, ifkappa) for calculating dust opacities.
2. Set the grid of model space (number of radial and vertical grid points).
3. Set the model parameter values (flaring, hGas100, hGasLocation, mDisk, dustToGas, turbulence), and assume an initial dust surface density, and prepare the input file of the dust density (i.e. dust_density.inp). The equations are taken from the equation 10 and 11 in Liu+ 2022, A&A, 668, A175.
   
   The density of dust is calculated by:
   $$\rho(R, z, a)=\frac{\Sigma(R, a)}{\sqrt{2 \pi} h(R, a)} \exp \left[-\frac{1}{2}\left(\frac{z}{h(R, a)}\right)^2\right]$$
   $`\Sigma(R, a)`$ and $`h(R,a)`$ is:
   $$\Sigma(R, a)=\Sigma_0(a)\left(\frac{R}{100 \mathrm{AU}}\right)^{-\gamma}$$
   $$h(R, a)=H_{\mathrm{gas}}\left(1+\frac{\mathrm{St}}{\alpha_{\mathrm{turb}}} \frac{1+2 \mathrm{St}}{1+\mathrm{St}}\right)^{-1 / 2}$$
   where $`St`$ is Stokes number $`= \frac{\pi}{2} \frac{\rho_{\text {grain }} a}{\sum_{\text {gas }}(R)}`$. $`\Sigma_0(a)`$ is related to the size of the dust, and the mass ratio obtained for each dust of different sizes is:
   $$f\left(a_j\right)=\frac{\int_{a_0}^{a_1} \frac{4 \pi}{3} \rho_{\text {grain }} n(a) a^3 \mathrm{~d} a}{\int_{a_{\text {lower }}}^{a_{\text {uppe }}} \frac{4 \pi}{3} \rho_{\text {grain }} n(a) a^3 \mathrm{~d} a}$$

5. Prepare auxiliary files for radmc3d(radmc.inp,dust_opac.inp,wavelength_micron.inp,stars.inp).
6. Run the thermal simulation to get the dust temperature distribution (radmc3d mctherm).
7. Simulate the continuum image, and compare the model radial intensity profile with the observed intensity profile.
8. Update the dust surface density (see section 4.2 in Li+ 2023, MNRAS, 518, 6092L).
9. interate the steps from 3 to 7, until the dust surface density is converaged.
10. Once the dust density is converaged (typically with 12 iterations), simulate the final SED and continuum image, and calculate the chi2_SED and chi2_image.

## Radmc3d_Simulation_mult.py
Multi-process programs can determine the fitting results for multiple models within the parameter space.

# Requirements #

## python env ##
python3.X

## package ##
Our code contains a large amount of code from the DHSHARP and RADMC3D packages.

dsharp_opac can be downloaded through:
        
        pip3 install dsharp_opac

Radmc3d can be downloaded through:
        https://github.com/dullemond/radmc3d-2.0
        
manual of RADMC-3D can be found:
        https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html
        
# How to use our code #
## 0.Example ##
You can execute the Radmc3d_Simulation.py file.

The DS_Tau_b6avgf.dat file gives the radial intensity profile of DS Tau that is extracted from the ALMA band6 continuum image, while the DS_Tau.txt file is the observed Spectral Energy Distribution (SED) of DS Tau.

## 1.Configuring the environment: This determines the program's working directory ##
Related parameters: pathIn

Location of the folder for dust_kappa_x.inp files

Tip: pathIn must end in '/'

Related parameters: pathOut
Location of the folder for output files

The observation files for sed and image should be placed in the pathOut location.

Tip: pathOut must end in '/'

The SED file consists of three columns. The first column represents the wavelength [um], the second column contains the observed flux densities [Jy], and the third column represents the errors [Jy].

The radial intensity profile file consists of three columns. The first column represents the distance to the central star [AU], the second column contains the corresponding intensity [mJy/beam], and the third column represents the errors [mJy/beam].
                  
## 2.Set parameters ##
There are four primary parameter categories: dust properties, grid settings, disk characteristics, and parameters for calculating the Spectral Energy Distribution (SED) and continuous emission images of the dust.
    
### Grain
amin: Minimum size of dust [cm]

amax: Maximum size of dust [cm]

grainP: -3.5

grainNum: The number of dust grain sizes. We typically take 32.

grainLam: The wavelength range of the absorption and scattering coefficients of dust. It is typically from 0.1 micron to 10000 micron.
    
### Grid
rin,rout: The maximum and minimum values of grid points on the radial direction of the disk. This value can exceed the inner and outer radii of the disk, but our program sets the range of radial grid points to be the same as the inner and outer radii of the disk.

x(y,z)bound: Range in the x(y,z) direction. 

nx(y,z): The number of grid points in each dimension. Corresponds to x (y, z) bound, which can be a list. Using a list allows for individual specification of the number of grid points within a given range.
    
### Disk ###
hGas100,hGasLocation: The scale height of the gas at a specific location is determined. The default value is set as 0.1 and 100, representing a scale height of 10 au at a distance of 100 au.

flaring: flaring index of the disk. 

**Please note that when the value is set to 0.1, it corresponds to a flaring index of 1.1.**

surfaceDensityP: The **initial** power-law exponent of the dust surface density. Since the dust surface density is iteraterd during the fitting process, it is not important. We typically take -0.5.

mDisk: This parameter is the total mass of the disk, i.e.  dust + gas mass, assuming the dust-to-gas ratio of 0.01.

dustToGas: The parameter is the ratio between the mass of dust and the mass of gas within the protoplanetary disk, with a default value set at 0.01.

turbulence: turbulence level. default:1e-4
    
### SED and Image ###

SEDLam: This parameter specifies the wavelength(s) used for calculating the Spectral Energy Distribution (SED). The comparison between the observed SED and the model-calculated SED during iterations determines the updated disk mass. 

seds,flux: These two parameters refer to the observed SED and the observed flux extracted radially along the disk.   

incl, PA, dpc:  The inclination, position angle, and distance of the protoplanetary disk.

pixelNum, sizeau: During the fitting process in the code, a set of fitted images will be generated. These three parameters correspond to the number of pixels in the image, the total size of the image in astronomical units (AU).

imageLam,beamMajor,beamMinor,beamPad: The band of the simulated image, and the size of the beam, position angle of the beam.

**Please note that:**

**1. You have the flexibility to assign distinct beam sizes for images in different bands. Simply add the desired elements to the respective lists. .**

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

# Acknowledgements
If you want to use this code, please:

1. Contact the author Yao Liu or DafaLi via email yliu@pmo.ac.cn or dfli@pmo.ac.cn.
2. Cite Liu+ 2022, A&A, 668, A175 and Li+ 2023, MNRAS, 518, 6092L.

    








 





        
