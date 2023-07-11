import numpy as np
import os
import radmc3dPy
from radmc3dPy import *
from astropy import constants as c
from scipy.interpolate import interp1d
from itertools import islice
from astropy.io import fits
import time
import multiprocess as mul


def createSEDLam(pathOut, SEDLam):
    file_path = '%scamera_wavelength_micron.inp' % pathOut
    with open(file_path, "w") as file:
        file.write("%s\n" % (len(SEDLam)))
        for lam in SEDLam:
            file.write("%s\n" % lam)


def makeRadSet(pathOut,
               nphot=10000000,
               nphot_scat=10000000,
               nphot_spec=200000,
               scattering_mode_max=1,
               camera_min_drr=0.1,
               modified_random_walk=1,
               istar_sphere=1,
               setthreads=8
               ):
    # define path
    file_path = "%sradmc.inp" % pathOut

    # create radmc.inp
    with open(file_path, "w") as file:
        file.write("nphot = %s\n" % nphot)
        file.write("nphot_scat = %s\n" % nphot_scat)
        file.write("nphot_spec =  %s\n" % nphot_spec)
        file.write("scattering_mode_max =  %s\n" % scattering_mode_max)
        file.write("camera_min_drr =  %s\n" % camera_min_drr)
        file.write("modified_random_walk =  %s\n" % modified_random_walk)
        file.write("istar_sphere =  %s\n" % istar_sphere)
        file.write("setthreads = %s\n" % setthreads)

    print("set_radmc.inp：{%s}" % pathOut)


def makeOpacity(pathOut, grainNum):
    file_path = "%sdustopac.inp" % pathOut
    with open('dustopac.inp', 'w') as wfile:
        # File format
        wfile.write('%-15s %s\n' % ('2', 'Format number of this file'))
        # Number of dust species
        wfile.write('%-15s %s\n' % (str(grainNum), 'Nr of dust species'))
        # Separator
        wfile.write('%s\n' % '============================================================================')
        for idust in range(grainNum):
            # Dust opacity will be read from a file
            wfile.write('%-15s %s\n' % ('1', 'Way in which this dust species is read'))

            wfile.write('%-15s %s\n' % ('0', '0=Thermal grain, 1=Quantum heated'))

            # Dustkappa filename extension
            wfile.write('%s %s %s\n' % ((idust + 1), '    ', 'Extension of name of dustkappa_***.inp file'))
            # Separator
            wfile.write('%s\n' % '----------------------------------------------------------------------------')


def makeStarandWavelength(pathOut, rstar, mstar, wavelength, Teff):
    fname = "%sstars.inp" % pathOut
    with open(fname, 'w') as wfile:
        wfile.write('%d\n' % 2)
        wfile.write('%d %d\n' % (1, len(wavelength)))
        wfile.write('%.9e %.9e %.9e %.9e %.9e\n' % (rstar, mstar,
                                                    0, 0, 0))

        wfile.write('%s\n' % ' ')
        for ilam in range(len(wavelength)):
            wfile.write('%.9e\n' % wavelength[ilam])
        wfile.write('%s\n' % ' ')
        wfile.write('%.9e\n' % (-Teff))

    fname = "%swavelength_micron.inp" % pathOut
    with open(fname, 'w') as wfile:
        print('Writing ' + fname)
        with open(fname, 'w') as wfile:
            wfile.write('%d\n' % len(wavelength))
            for ilam in range(len(wavelength)):
                wfile.write('%.9e\n' % wavelength[ilam])


def createDSHARPGrain(pathOut, lam, amid, a1, a2, p, num, na=500):
    # os.chdir('/Users/lidafa/Desktop/dsharp/dsharp_opac-master 2/notebooks')

    # from disklab import opacity
    import dsharp_opac as opacity

    a = np.logspace(a1, a2, 500)
    # lam = np.logspace(-5, 0, 160)
    density_water = 0.92
    density_silicates = 3.30
    density_troilite = 4.83
    density_organics = 1.50

    # default values

    constants_default = [
        opacity.diel_warrenbrandt08(),
        opacity.diel_draine2003('astrosilicates'),
        opacity.diel_henning('troilite'),
        opacity.diel_henning('organics'),
    ]

    densities_default = np.array([
        density_water,
        density_silicates,
        density_troilite,
        density_organics,
    ])

    fv_default = np.array([
        0.3642,
        0.1670,
        0.0258,
        0.4430])
    rhos_default = (fv_default * densities_default).sum()
    rhos_default = (fv_default * densities_default).sum()
    d_def = opacity.diel_mixed(constants_default, fv_default, rule='Bruggeman')
    res_def = opacity.get_opacities(a, lam, rho_s=rhos_default, diel_const=d_def, extrapol=True,
                                    extrapolate_large_grains=True)

    ###############Output dustkappa_x.inp#############
    os.chdir('%s' % pathOut)

    lam * 1e-4

    s = a ** 0.5
    # s        = (a/a[0])**0.5
    ia = np.abs(a - a2).argmin()  # amax: cm
    s[ia + 1:] = 0.
    s = s / s.sum()

    print('using amax = {} cm'.format(amid))

    # average the opacities with that size distribution

    k_a_avg = (res_def['k_abs'].T * s).sum(1)
    k_s_avg = (res_def['k_sca'].T * s).sum(1)
    g_avg = (res_def['g'].T * s).sum(1)
    write_radmc3d_dustkappa_from_array('%s' % (num), lam, k_a_avg, k_s_avg, g_avg)


def write_radmc3d_dustkappa_from_array(name, lam, k_abs, k_sca, g=None, path='.'):
    """
    Write out the opacity such that it can be read from
    RADMC-3D. This function writes the file dustkappa_[name].inp which
    writes only the absorption opacity, scattering opacity, and the
    assymetry parameter g as function of wavelength into the file.

    Arguments:
    ----------

    name : str
        name to be put in the filename: dustkappa_[name].inp

    lam : array-like
        wavelength grid in cm (will be writtin in micron)

    k_abs, k_sca, g : array-like
        absorption and scattering opacity and asymmetry factor on array `lam`

    Keywords:
    ---------

    path : str
        path where to write the file, default: current dir

    Output:
    -------
    writes out the file dustkappa_[name].inp
    """
    filename = 'dustkappa_' + name + '.inp'

    if g is None:
        data = np.array([
            lam * 1e4,
            k_abs,
            k_sca
        ]).T
    else:
        data = np.array([
            lam * 1e4,
            k_abs,
            k_sca,
            g
        ]).T

    header = '3\n{:d}\n'.format(len(lam))
    np.savetxt(filename, data, header=header, comments='')


def calculateDustMassWeight(a, grainP, mDust):
    in_pwgt = (grainP - 2) / 3 + 2
    mass_per_grain = (1.675 * 4.0 * np.pi / 3.) * (a) ** 3
    in_dum = (mass_per_grain / mass_per_grain[0]) ** in_pwgt
    in_weights = in_dum / sum(in_dum)
    # layermass = mDust * in_weights

    return in_weights


def writeSpatialGrid(fname, crd_sys, act_dim, xi, yi, zi, old=False):
    """Writes the wavelength grid to a file (e.g. amr_grid.inp).

    Parameters
    ----------

    fname : str, optional
            File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used.

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used
    """

    #
    # Write the spatial grid for radmc3d
    #
    nx = len(xi) - 1
    ny = len(yi) - 1
    nz = len(zi) - 1
    nxi = nx + 1
    nyi = ny + 1
    nzi = nz + 1
    if not old:
        # if fname == '':
        # fname = 'amr_grid.inp'
        # fname = '%s/amr_grid.inp' %path

        print('Writing ' + fname)
        with open(fname + 'amr_grid.inp', 'w') as wfile:
            # Format number
            wfile.write('%d\n' % 1)
            # AMR style (0=regular  NO AMR)
            wfile.write('%d\n' % 0)
            # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if crd_sys == 'car':
                wfile.write('%d\n' % 0)
            # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if crd_sys == 'sph':
                wfile.write('%d\n' % 100)
            # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if crd_sys == 'cyl':
                wfile.write('%d\n' % 200)
            # Gridinfo
            wfile.write('%d\n' % 0)

            # Active dimensions
            wfile.write('%d %d %d \n' % (act_dim[0], act_dim[1], act_dim[2]))
            # Grid size (x,y,z or r,phi,theta, or r,phi,z)
            wfile.write('%d %d %d \n' % (nx, ny, nz))
            for i in range(nxi):
                wfile.write('%.9e\n' % xi[i])
            for i in range(nyi):
                wfile.write('%.9e\n' % yi[i])
            for i in range(nzi):
                wfile.write('%.9e\n' % zi[i])
        wfile.close()
    # #
    # # Write the spatial grid for radmc
    # #
    # else:
    #
    #     fname = 'radius.inp'
    #     with open(fname, 'w') as wfile:
    #
    #         print('Writing ' + fname)
    #         x = np.sqrt(xi[1:] * xi[:-1])
    #         wfile.write("%d\n" % nx)
    #         wfile.write(" \n")
    #         for i in range(nx):
    #             wfile.write("%.7e\n" % x[i])
    #
    #     fname = 'theta.inp'
    #     with open(fname, 'w') as wfile:
    #         print('Writing ' + fname)
    #         wfile.write("%d 1\n" % (ny / 2))
    #         wfile.write(" \n")
    #         for i in range(int(ny / 2)):
    #             wfile.write("%.7e\n" % y[i])


def gridx(xbound, nx, xres_nlev=0, xres_nspan=0, xres_nstep=0):
    act_dim = [1, 1, 0]
    nxi1 = [i + 1 for i in nx]
    if len(nxi1) > 1:
        nxi = sum(nxi1)
        nx = nxi - 1
        xi = xbound[0] * (xbound[1] / xbound[0]) ** (
                np.arange(nxi1[0], dtype=np.float64) / float(nxi1[0]))
        for ipart in range(1, len(nxi1) - 1):
            dum = xbound[ipart] * (xbound[ipart + 1] / xbound[ipart]) ** (
                    np.arange(nxi1[ipart], dtype=np.float64) / float(nxi1[ipart]))
            xi = np.append(xi, dum)

        ipart = len(nxi1) - 1
        dum = xbound[ipart] * (xbound[ipart + 1] / xbound[ipart]) ** (
                np.arange(nxi1[ipart], dtype=np.float64) / float(nxi1[ipart] - 1))
        xi = np.append(xi, dum)
        x = np.sqrt(xi[0:nx] * xi[1:nx + 1])
    else:
        if act_dim[0] == 1:
            nxi = nxi1[0]
            xi = xbound[0] * (xbound[1] / xbound[0]) ** (
                    np.arange(nxi, dtype=np.float64) / float(nxi - 1.))
            nx = nxi - 1
            x = np.sqrt(xi[0:nx] * xi[1:nx + 1])
        else:
            x = [0.]
            xi = [0., 0., ]
            nx = 1
            nxi = 2

    # Refinement of the inner edge of the grid
    # This has to be done properly

    if xres_nlev > 0:
        ri_ext = np.array([xi[0], xi[xres_nspan]])
        for i in range(xres_nlev):
            dum_ri = ri_ext[0] + (ri_ext[1] - ri_ext[0]) * np.arange(xres_nstep + 1,
                                                                     dtype=np.float64) / float(
                xres_nstep)
            # print ri_ext[0:2]/au
            # print dum_ri/au
            ri_ext_old = np.array(ri_ext)
            ri_ext = np.array(dum_ri)
            ri_ext = np.append(ri_ext, ri_ext_old[2:])

        r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

        xi = np.append(ri_ext, xi[xres_nspan + 1:])
        x = np.append(r_ext, x[xres_nspan:])
        nx = x.shape[0]
        nxi = xi.shape[0]
    return x, xi, nx, nxi


def gridy(ybound, ny):
    act_dim = [1, 1, 0]
    nyi1 = [i + 1 for i in ny]
    if len(nyi1) > 1:

        # Check if we go to the full [0,pi] interval or only use the upper half-plane [0, pi/2]

        if ybound[len(ybound) - 1] != np.pi / 2.:
            nyi = sum(nyi1) + 1
            ny = nyi - 1
            yi = ybound[0] + (ybound[1] - ybound[0]) * (
                    np.arange(nyi1[0], dtype=np.float64) / float(nyi1[0]))

            for ipart in range(1, len(nyi1) - 1):
                # Now make sure that pi/2 will be a cell interface
                #
                # BUGFIX! 16-05-2012
                # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi1)>1)
                # This is now fixed
                if ybound[ipart] < np.pi / 2.:
                    dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            np.arange(nyi1[ipart], dtype=np.float64) / float(nyi1[ipart]))
                else:
                    if ybound[ipart] == np.pi / 2.:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                (np.arange(nyi1[ipart] + 1, dtype=np.float64)) / (float(nyi1[ipart])))
                    else:
                        dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                                (np.arange(nyi1[ipart], dtype=np.float64) + 1.) / float(nyi1[ipart]))

                yi = np.append(yi, dum)

            ipart = len(nyi1) - 1
            if len(nyi1) == 2:
                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                        (np.arange(nyi1[ipart] + 1, dtype=np.float64)) / (float(nyi1[ipart])))
            else:
                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                        (np.arange(nyi1[ipart], dtype=np.float64) + 1.) / float(nyi1[ipart]))

        else:
            nyi = sum(nyi1) + 1
            ny = nyi - 1
            yi = ybound[0] + (ybound[1] - ybound[0]) * (
                    np.arange(nyi1[0], dtype=np.float64) / float(nyi1[0]))
            for ipart in range(1, len(nyi1) - 1):
                # Now make sure that pi/2 will be a cell interface
                #
                # BUGFIX! 16-05-2012
                # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi1)>1)
                # This is now fixed
                if ybound[ipart] < np.pi / 2.:
                    dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            np.arange(nyi1[ipart], dtype=np.float64) / float(nyi1[ipart]))
                else:
                    dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                            (np.arange(nyi1[ipart] + 1, dtype=np.float64)) / (float(nyi1[ipart])))

                yi = np.append(yi, dum)

            ipart = len(nyi1) - 1

            if len(nyi1) == 2:
                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                        (np.arange(nyi1[ipart] + 1, dtype=np.float64)) / (float(nyi1[ipart])))
            else:
                dum = ybound[ipart] + (ybound[ipart + 1] - ybound[ipart]) * (
                        (np.arange(nyi1[ipart], dtype=np.float64) + 1.) / float(nyi1[ipart]))

        yi = np.append(yi, dum)
        y = 0.5 * (yi[0:ny] + yi[1:ny + 1])

    else:
        if act_dim[1] == 1:
            nyi = nyi1[0]
            yi = ybound[0] + (ybound[1] - ybound[0]) * (
                    np.arange(nyi, dtype=np.float64) / float(nyi - 1.))
            ny = nyi - 1
            y = 0.5 * (yi[0:ny] + yi[1:ny + 1])
        else:
            y = [0.]
            yi = [0., 0., ]
            ny = 1
            nyi = 2
    return y, yi, ny, nyi


def gridz(zbound, nz):
    act_dim = [1, 1, 0]
    nzi1 = [i + 1 for i in nz]
    if len(nzi1) > 1:
        nzi = sum(nzi1)
        nz = nzi - 1

        zi = zbound[0] + (zbound[1] - zbound[0]) * (
                np.arange(nzi1[0], dtype=np.float64) / float(nzi1[0]))
        for ipart in range(1, len(nzi1) - 1):
            dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                    np.arange(nzi1[ipart], dtype=np.float64) / float(nzi1[ipart]))
            zi = np.append(zi, dum)
        ipart = len(nzi1) - 1
        dum = zbound[ipart] + (zbound[ipart + 1] - zbound[ipart]) * (
                np.arange(nzi1[ipart], dtype=np.float64) / float(nzi1[ipart] - 1))
        zi = np.append(zi, dum)
        z = 0.5 * (zi[0:nz] + zi[1:nz + 1])
    else:
        if act_dim[2] == 1:
            nzi = nzi1[0]
            zi = zbound[0] + (zbound[1] - zbound[0]) * (
                    np.arange(nzi, dtype=np.float64) / float(nzi - 1))
            nz = nzi - 1
            z = 0.5 * (zi[0:nz] + zi[1:nz + 1])
        else:
            z = np.array([0.])
            zi = np.array([0., np.pi * 2.])
            nz = 1
            nzi = 2
    return z, zi, nz, nzi


def writeSpatialGrid(fname, crd_sys, act_dim, xi, yi, zi, old=False):
    """Writes the wavelength grid to a file (e.g. amr_grid.inp).

    Parameters
    ----------

    fname : str, optional
            File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used.

    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used
    """

    #
    # Write the spatial grid for radmc3d
    #
    nx = len(xi) - 1
    ny = len(yi) - 1
    nz = len(zi) - 1
    nxi = nx + 1
    nyi = ny + 1
    nzi = nz + 1
    if not old:
        # if fname == '':
        # fname = 'amr_grid.inp'
        # fname = '%s/amr_grid.inp' %path

        print('Writing ' + fname)
        with open(fname + 'amr_grid.inp', 'w') as wfile:
            # Format number
            wfile.write('%d\n' % 1)
            # AMR style (0=regular  NO AMR)
            wfile.write('%d\n' % 0)
            # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if crd_sys == 'car':
                wfile.write('%d\n' % 0)
            # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if crd_sys == 'sph':
                wfile.write('%d\n' % 100)
            # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if crd_sys == 'cyl':
                wfile.write('%d\n' % 200)
            # Gridinfo
            wfile.write('%d\n' % 0)

            # Active dimensions
            wfile.write('%d %d %d \n' % (act_dim[0], act_dim[1], act_dim[2]))
            # Grid size (x,y,z or r,phi,theta, or r,phi,z)
            wfile.write('%d %d %d \n' % (nx, ny, nz))
            for i in range(nxi):
                wfile.write('%.9e\n' % xi[i])
            for i in range(nyi):
                wfile.write('%.9e\n' % yi[i])
            for i in range(nzi):
                wfile.write('%.9e\n' % zi[i])
        wfile.close()

def calculateSigma(x,xi,y,yi,hrpivot,plsig1,mdisk):
    sig = np.power((x/hrpivot), plsig1)
    sigSqua = np.pi*(xi[1:]**2 - xi[:-1]**2)
    sigMtot = (sig * sigSqua).sum()
    sig = sig * mdisk / sigMtot

    return sig

def calculateDustPho(pathOut,x, y,hrdisk, hrpivot, plh, plsig1, grainSize, alpha, ilayer, mdisk, dustTogas,sig0,):
    nx = len(x)
    ny = len(y)
    rhodust = np.zeros([nx, ny])
    hrstore = np.zeros(nx)
    hfactor = np.zeros(nx)

    sigmaGas = calculateSigma(x, xi, y, yi, hrpivot, plsig1, mdisk)
    np.savetxt('sigmaGas.txt', np.transpose([x/natconst.au,sigmaGas]))
    for ix in range(nx):
        hrstore[ix] = hrdisk * (x[ix] / hrpivot) ** plh

    stnumber = 1.675 * grainSize * np.pi / sigmaGas / 2

    for ix in range(nx):
        hfactor[ix] = (1.00 + (stnumber[ix] / alpha) * (1.00 + 2.00 * stnumber[ix]) / (1.00 + stnumber[ix])) ** -0.5
        #hfactor[i] = np.power(ftmp,-0.5)

    np.savetxt('layer%s.txt' %(ilayer),np.transpose([x/natconst.au,stnumber,hrstore,hfactor]))

    sigmaDust = calculateSigma(x, xi, y, yi, hrpivot, plsig1, mdisk*dustTogas*sig0)
    np.savetxt('sigmaDust%s.txt' %ilayer, np.transpose([x / natconst.au, sigmaDust]))

    for ir in range(nx):
        hr = hrstore[ir] * hfactor[ir]
        rhodust[ir, :] = 1 /np.sqrt(2*np.pi) * sigmaDust[ir] * np.exp(-0.5 * ((np.pi/2-y[:]) / hr) ** 2) / (hr * x[ir])

    return rhodust

def calculateIterationDustPho(pathOut,x, y,hrdisk, hrpivot, plh, grainSize, alpha, ilayer, dustToGas,sig0,sigmaGas):
    nx = len(x)
    ny = len(y)
    rhodust = np.zeros([nx, ny])
    hrstore = np.zeros(nx)
    hfactor = np.zeros(nx)

    sigmaDust = sigmaGas * dustToGas * sig0
    for ix in range(nx):
        hrstore[ix] = hrdisk * (x[ix] / hrpivot) ** plh

    stnumber = 1.675 * grainSize * np.pi / sigmaGas / 2

    for ix in range(nx):
        hfactor[ix] = (1.00 + (stnumber[ix] / alpha) * (1.00 + 2.00 * stnumber[ix]) / (1.00 + stnumber[ix])) ** -0.5

    np.savetxt('layer%s.txt' %(ilayer),np.transpose([x/natconst.au,stnumber,hrstore,hfactor]))

    np.savetxt('sigmaDust%s.txt' % ilayer, np.transpose([x / natconst.au, sigmaDust]))

    for ir in range(nx):
        hr = hrstore[ir] * hfactor[ir]
        rhodust[ir, :] = 1 / np.sqrt(2 * np.pi) * sigmaDust[ir] * np.exp(-0.5 * ((np.pi / 2 - y[:]) / hr) ** 2) / (
                    hr * x[ir])
    return rhodust

def getCellVolume(x, y, z, xi, yi, zi):
    """Calculates the volume of grid cells."""
    nx = len(x)
    ny = len(y)
    nz = len(z)
    vol = np.zeros([nx, ny, nz], dtype=np.float64)
    diff_r3 = xi[1:] ** 3 - xi[:-1] ** 3
    diff_cost = np.cos(yi[:-1]) - np.cos(yi[1:])
    diff_phi = 2. * np.pi
    for ix in range(nx):
        for iy in range(ny):
            vol[ix, iy, :] = 1. / 3. * diff_r3[ix] * diff_cost[iy] * diff_phi

    return vol


def writeDustDensity(pathOut, sigmad, n):
    f = open('%sdust_density.inp' % pathOut, 'w')
    f.write('1\n')
    f.write('%s\n' % (len(sigmad[:, 0, 0, 0]) * len(sigmad[0, :, 0, 0]) * len(sigmad[0, 0, :, 0])))
    f.write('%s\n' % (len(sigmad[0, 0, 0, :])))
    for z in range(len(sigmad[0, 0, 0, :])):
        for j in range(len(sigmad[0, :, 0, 0])):
            for i in range(len(sigmad[:, 0, 0, 0])):
                f.write('%.9e\n' % sigmad[i, j, 0, z])
    f.close()
    os.system('cp dust_density.inp ./iteration_warehouse/dust_density_%s.inp' % n)


def extractedFlux(PA, incl, lad, beamMo, beamMi, rdisk, pixelNum, pixelsize, n):
    for la, beam1, beam2 in zip(lad, beamMo, beamMi):
        hdulist = fits.open('image_con_%s.fits' % (la))
        tdata = hdulist[0].data
        maxx = tdata.max()
        position = np.argmax(tdata)
        positiony = position // pixelNum
        positionx = position - positiony * pixelNum
        pa = PA
        R = np.arange((rdisk / pixelsize))
        flux1 = []
        flux3 = []
        y2 = []
        x2 = []
        for r in R:
            theta = np.linspace(0, 2 * np.pi, num=181, endpoint=True)
            # if theta[-1] != 2*np.pi:
            #     theta = np.append(theta,2*np.pi
            for i in theta:
                x = r * np.cos(i)
                y = np.cos(incl / 180. * np.pi) * r * np.sin(i)

                sin_pa = np.sin(pa / 180. * np.pi)
                cos_pa = np.cos(pa / 180. * np.pi)
                cos_pa_x = cos_pa * x
                cos_pa_y = cos_pa * y
                sin_pa_x = sin_pa * x
                sin_pa_y = sin_pa * y
                xx = cos_pa_x - sin_pa_y
                yy = sin_pa_x + cos_pa_y

                xx = xx + positionx
                yy = yy + positiony
                x2.append(xx)
                y2.append(yy)

                xx0 = int(xx)
                xx1 = xx0 + 1
                yy0 = int(yy)
                yy1 = yy0 + 1
                fluxlu = tdata[0][yy0][xx0]
                fluxld = tdata[0][yy0][xx1]
                fluxru = tdata[0][yy1][xx0]
                fluxrd = tdata[0][yy1][xx1]
                flux = fluxlu * (xx1 - xx) * (yy1 - yy) / (xx1 - xx0) / (yy1 - yy0) + fluxld * (xx - xx0) \
                       * (yy1 - yy) / (xx1 - xx0) / (yy1 - yy0) + fluxru * (yy - yy0) * (xx1 - xx) / (xx1 - xx0) / (
                               yy1 - yy0) \
                       + fluxrd * (xx - xx0) * (yy - yy0) / (xx1 - xx0) / (yy1 - yy0)
                flux1.append(flux)

            y2 = []
            x2 = []
            flux2 = sum(flux1) / len(flux1)
            flux1 = []
            flux3.append(flux2)

        R2 = (R) * pixelsize
        flux3 = [flux3[i] * 1000 * (np.pi / 4 * beam1 * beam2) / (3600 / np.pi * 180) ** 2 / np.log(2) for i in
                 range(len(flux3))]  # 2mJy
        # flux3 = [flux3[i] * 1000 for i in range(len(flux3))]
        cc = np.transpose([R2, flux3])
        np.savetxt('flux_R_%s_%s.txt' % (la, n), cc)
        # os.system('cp flux_R_%s.txt ./iteration_warehouse/flux_R_%s_%s.txt' %(la,la,n))


def calculateSED(incl, Pa, dpc, mdisk, n):
    # SED
    os.system('radmc3d spectrum incl %s posang %s loadlambda' % (incl, (Pa)))
    sed = analyze.readSpectrum()
    sed[:, 1] = sed[:, 1] / (dpc) ** 2
    sed1 = sed.copy()
    sed[:, 1] = sed[:, 1] * 1.0 * 1.e4 * 2.99792458e10 / sed[:, 0]
    sed1[:, 1] = sed1[:, 1] * 10 ** 23
    np.savetxt('sed_%s.txt' % (n), sed1)
    # np.savetxt('sed_%s.txt,erg' % (n), sed)
    # mdiskFactor = seds / sed1[-1, 1]
    # mdisk = mdisk * mdiskFactor

    # return mdisk


def calculateImage(incl, Pa, dpc, lamd, beamMa, beamMi, pixelNum, sizeau, beamPad):
    for lam, beam1, beam2, beamPa in zip(lamd, beamMa, beamMi, beamPad):
        ##计算image
        os.system('rm image_%s.fits' % lam)
        os.system('rm image_con_%s.fits' % lam)
        os.system('radmc3d image lambda %s incl %s posang %s npixx %s npixy %s sizeau %s secondorder' % (
            lam, incl, (Pa), pixelNum, pixelNum, sizeau))
        im = image.readImage()
        im.writeFits('image_%s.fits' % lam)
        cim = im.imConv(fwhm=[beam1, beam2], pa=beamPa, dpc=dpc)
        cim.writeFits('image_con_%s.fits' % (lam))


def iteration(pathOut, pathkappa, grainNum, incl, PA, beamPad, dpc, mdisk, dustToGas, lam, beamMa, beamMi, sed,
              flux, pixelNum, sizeau, pixelSize, rdisk,massWeight):
    n = 0
    os.chdir('%s' % pathOut)
    time_pho = []
    time_start = time.time()
    mdiskiter = []
    for isize in range(1, grainNum + 1):
        os.system('cp %sdustkappa_%s.inp %s' % (pathkappa, isize, pathOut))
    for n in range(13):
        os.system('cp sigmaGas.txt ./iteration_warehouse/sigmaGas_%s.txt' % n)
        for ig in range(grainNum):
            os.system('mv sigmaDust%s.txt ./iteration_warehouse/sigmaDust_%s_iter_%s.txt' % (ig, ig, n))
            os.system('mv layer%s.txt ./iteration_warehouse/layer_%s_iter_%s.txt' % (ig, ig, n))
        os.system('cp dust_density.inp ./iteration_warehouse/dust_density_%s.ing' % (n))
        # temperature
        os.system('radmc3d mctherm')
        calculateImage(incl, PA, dpc, lam, beamMa, beamMi, pixelNum, sizeau, beamPad)
        extractedFlux(PA, incl, lam, beamMa, beamMi, rdisk / natconst.au, pixelNum, pixelSize, n)
        flux3 = np.loadtxt('flux_R_%s_%s.txt' % (lam[-1], n))
        os.system('mv flux_R_%s_%s.txt ./iteration_warehouse' % (lam[-1], n))
        # flux3 = np.loadtxt('flux_R_%s_erg.txt')
        # flux3[:,1] = flux3[:,1]*1000
        fluxs = flux
        ## if need interp1d
        li1 = interp1d(flux3[:, 0] * natconst.au, flux3[:, 1], kind='cubic')
        li2 = interp1d(fluxs[:, 0] * natconst.au, fluxs[:, 1], kind='cubic')
        grid = analyze.readGrid()
        flux3n = li1(grid.x)
        fluxsn = li2(grid.x)

        fx = []
        for i in range(grid.nx):
            fx1 = fluxsn[i] / flux3n[i]
            if 0.7 <= fx1 <= 1.3:
                fx.append(fx1)
            elif fx1 >= 1.3:
                fx1 = 1.3
                fx.append(fx1)
            elif fx1 <= 0.7:
                fx1 = 0.7
                fx.append(fx1)
        np.savetxt('fx%s' % (n), fx)
        os.system('mv fx%s ./iteration_warehouse' % (n))
        sigGas = np.loadtxt('sigmaGas.txt')
        sigmaGas = sigGas[:, 1] * fx
        np.savetxt('sigmaGas.txt', np.transpose([x / natconst.au, sigmaGas]))
        for i in range(grainNum):
            phoList[:, :, 0, i] = calculateIterationDustPho(pathOut, x, y, hGas100, hGasLocation, flaring, aList[i], turbulence,i,dustToGas, massWeight[i], sigmaGas)
        writeDustDensity(pathOut, phoList, n)
        mdust = calculateDustMass(pathOut, grainNum)
        mdiskiter.append(mdust/natconst.ms)

        time_end = time.time()
        time_pho.append(time_end - time_start)
        np.savetxt('time_section%s' % (n), time_pho)
        os.system('mv time_section%s ./iteration_warehouse' % (n))

        if max(abs(fx - np.ones(len(fx)))) < 0.05:
            break
    os.system('cp sigmaGas.txt ./iteration_warehouse/sigmaGas_%s.txt' %n)
    for ig in range(grainNum):
        os.system('mv sigmaDust%s.txt ./iteration_warehouse/sigmaDust_%s_iter_%s.txt' %(ig,ig,n))
        os.system('mv layer%s.txt ./iteration_warehouse/layer_%s_iter_%s.txt' %(ig,ig,n))
    os.system('cp dust_density.inp ./iteration_warehouse/dust_density_%s.inp' % (n))
    os.system('cp massWeight ./iteration_warehouse/massWeight')
    os.system('cp mdustiter ./iteration_warehouse/mdustiter')

    os.system('cp ./iteration_warehouse/flux_R_%s_%s.txt ./flux_R_%s.txt' % (lam[-1], n, lam[-1]))
    fluxBestfit1 = np.loadtxt('flux_R_%s.txt' % (lam[-1]))
    fluxBestfit = np.interp(fluxs[:, 0], fluxBestfit1[:, 0], fluxBestfit1[:, 1])
    chiFlux = (fluxs[:, 1] - fluxBestfit) ** 2 / (fluxs[:, 1] * 0.1) ** 2
    np.savetxt('ChiFlux.txt', chiFlux)
    calculateSED(incl, PA, dpc, mdisk, n + 1, )
    sedBestfit = np.loadtxt('sed_%s.txt' % (n + 1))
    # sedBestfit = np.interp(fluxs[:,0],fluxBestfit1[:,0],fluxBestfit1[:,1])
    chiSed = (sed[:, 1] - sedBestfit[:, 1]) ** 2 / sed[:, 2] ** 2
    np.savetxt('ChiSed.txt', chiSed)

    # np.savetxt('mdustiter', mdiskiter)
    # os.system('mv mdust ./iteration_warehouse')


def calculateDustMass(pathInput, grainNum):
    os.chdir('%s' % pathInput)
    dustmass = 0
    grid = analyze.readGrid()
    vol = grid.getCellVolume()
    data = analyze.readData(ddens=True, binary=False)
    for i in range(grainNum):
        dustmass = vol * data.rhodust[:, :, :, i] + dustmass
    surf = np.zeros([grid.nx, grid.nz], dtype=np.float64)
    diff_r2 = (grid.xi[1:] ** 2 - grid.xi[:-1] ** 2) * 0.5
    diff_phi = grid.zi[1:] - grid.zi[:-1]
    for ix in range(grid.nx):
        surf[ix, :] = diff_r2[ix] * diff_phi

    mass = dustmass.sum(1).sum(0)
    # sigma = mass / surf

    # sigmaInput = np. interp(xTarget,grid.x.ravel()/natconst.au,sigma.ravel())
    return mass


##### set parameter #####

# path
pathIn = './'
pathOut = './'  # pathOut must end in '/'
os.chdir('%s' % pathOut)
os.system('mkdir ./iteration_warehouse')

# grain(cm)
amin = 1e-6  # 0.01um
amax = 0.3
grainP = -3.5
grainNum = 32
grainLam = np.logspace(-5, 0, 160)
ifkappa = 1

# grid (see Radmc3dPy)
rin = 0.1 * natconst.au
rout = 80 * natconst.au
xbound = [rin, 1 * natconst.au, rout]
nx = [100, 100]
ybound = [0., np.pi / 3., np.pi / 2., 2 * np.pi / 3., np.pi]  # Boundaries for the y grid
ny = [25, 50, 50, 25]  # Number of grid points in the second dimension (to switch off this dimension set it to 0)
zbound = [0, 0]
nz = [0]

# disk parameters
hGas100 = 0.07
hGasLocation = 100 * natconst.au
flaring = 0.05  # flaring = 1.1
surfaceDensityP = -0.5
mDisk = 0.14 * natconst.ms  # total mass = gas + dust
dustToGas = 0.01
turbulence = 1e-3

# SED and ALMA
SEDLam = [1329, 2855]  # band3 and band6 um
seds = np.loadtxt('%sDS_Tau.txt' % pathIn)  ##check target wavelength index
flux = np.loadtxt('%sDS_Tau_b6avgf.dat' % pathIn)
incl = 65.2
PA = 159.62 - 90
dpc = 158
pixelNum = 1200
sizeau = 240
pixelSize = sizeau / pixelNum
imageLam = [1329, ]  # um
beamMajor = [0.13, ]
beamMinor = [0.09, ]
beamPad = [20, ]
# star and wavelength
Teff = 3792
# L = 10**(-0.61)/natconst.ls
Lstar = 0
# blackbody
if Lstar != 0:
    rstar = (Lstar / (4 * np.pi * natconst.ss * T ** 4)) ** 0.5
rstar = 2 * natconst.rs
mstar = 0.58 * natconst.ms
wavelength = np.logspace(-1, 4, num=100)

##### create grain ######
aList = np.logspace(np.log10(amin), np.log10(amax), num=grainNum)
aLog10 = np.log10(aList)
aBoundary = np.zeros(grainNum + 1)
aBoundary[0:grainNum] = aLog10 - (aLog10[1] - aLog10[0]) / 2
aBoundary[-1] = aLog10[-1] + (aLog10[1] - aLog10[0]) / 2
if ifkappa == 1:
    for aIndex in range(grainNum):
        createDSHARPGrain(pathIn, grainLam, aList[aIndex], aBoundary[aIndex], aBoundary[aIndex + 1], grainP, aIndex + 1,
                          na=500)
massWeight = calculateDustMassWeight(aList, grainP, mDisk * dustToGas)
np.savetxt('massWeight',massWeight*0.14)

##### create grid ######
x, xi, nx, nxi = gridx(xbound, nx)
y, yi, ny, nyi = gridy(ybound, ny)
z, zi, nz, nzi = gridz(zbound, nz)
writeSpatialGrid(fname=pathOut, crd_sys='sph', act_dim=[1, 1, 0], xi=xi, yi=yi, zi=zi, old=False)

##### calculate pho ######
phoList = np.zeros([len(x), len(y), len(z), grainNum], dtype=np.float64)
for i in range(grainNum):
    phoList[:, :, 0, i] = calculateDustPho(pathOut,x, y, hGas100, hGasLocation, flaring, surfaceDensityP,aList[i],turbulence,i,mDisk,dustToGas, massWeight[i])
writeDustDensity(pathOut, phoList, 'ini')

##### radiative transfer ######
createSEDLam(pathOut, SEDLam)
makeRadSet(pathOut,
           nphot=10000000,
           nphot_scat=10000000,
           nphot_spec=200000,
           scattering_mode_max=1,
           camera_min_drr=0.1,
           modified_random_walk=1,
           istar_sphere=1,
           setthreads=8
           )

makeOpacity(pathOut, grainNum)

makeStarandWavelength(pathOut, rstar, mstar, wavelength, Teff)
iteration(pathOut, pathIn, grainNum, incl, PA, beamPad, dpc, mDisk, dustToGas, imageLam, beamMajor, beamMinor,
          seds, flux, pixelNum, sizeau, pixelSize, rout,massWeight)

