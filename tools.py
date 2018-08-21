"""
list of def:
add_revise_fits_keyword
blur_array
cal_R_gc
cal_angular_resolution
cal_chi2
cal_chi2_no_weight
cal_chi2_pearson
cal_dt_from_df_glon
cal_hii_electron_density
cal_hii_emission_measure
cal_hii_mass
cal_int_flux_density_jy
cal_ionising_Lyman_continuum_photon_count
cal_overlapping_uv_range
cal_radius
cal_si
cal_wavelength
change_beam_convol_miriad
change_coord_swarp
change_fits_with_a_constant
change_names
change_pix_size
compare_speed_template
convert_csv_to_reg_template
convert_jyperbeam_k
convert_spec_to_img
cut_fits
cut_fits_by_hand
cut_fits_downloaded
download_webpage
filter_fits_table_template
find_num_line_in_file
gauss_kern
gcd
get_def_list
get_fits_all_pix_statistics
is_nan
join_fits_tables
join_npy_files
join_votable_vstack
k_per_pc_to_physical
move_fits
operate_two_fits
pix_val_in_reg
pix_val_in_reg_old
replace_all
rm_blank_txt
si_map
smooth_fits_img
statistical_te_from_Balsar17
template_generate_ds9_reg
use_multiple_cores_tamplate
GLong_GLat_to_RA_DEC
RA_DEC_to_GLong_GLat
"""

from astropy.io import fits
from astropy import wcs
from astropy.table import vstack, hstack
from astropy.table import join
from astropy.table import Table as astropy_table
from os import system as s
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, Galactic, SkyCoord
from urllib.request import urlopen
import pyregion
import astropy.io.fits as pyfits
from astropy import constants as const
from scipy import mgrid
from scipy import signal
import time
import timeit
import random
from astroquery import nasa_ads
import os
import copy
from astropy import constants as const


def find_num_line_in_file(fname):
    # find out the number of lines in a file
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def operate_two_fits(fits1, fits2, fits_ou, operator):
    # read fits data
    fits1_open = fits.open(fits1)
    fits2_open = fits.open(fits2)
    fits1_data = fits1_open[0].data
    fits2_data = fits2_open[0].data
    fits1_head = fits1_open[0].header

    # operation
    if operator == '+':
        fits_ou_data = fits1_data + fits2_data
    elif operator == '-':
        fits_ou_data = fits1_data - fits2_data
    elif operator == '*':
        fits_ou_data = fits1_data * fits2_data
    elif operator == '/':
        fits_ou_data = fits1_data / fits2_data
    else:
        print('What are you going to do? + - * or /')

    # write
    out = fits.PrimaryHDU()
    out.data = fits_ou_data
    out.header = fits1_head
    out.writeto(fits_ou, clobber=True)
    return 0


def change_coord_swarp(fits_in, fits_ou):
    """
    change the coordinate of fits image using swarp
    :param list_img:
    :param fits_ou:
    :return:
    """
    s('swarp ' + fits_in + \
      ' -COMBINE N -IMAGEOUT_NAME ' + fits_ou + \
      ' -CELESTIAL_TYPE GALACTIC -PROJECTION_TYPE CAR -SUBTRACT_BACK N' + \
      ' -COPY_KEYWORDS TELESCOP,MWAVER,MWADATE,BTYPE,BUNIT,BMAJ,BMIN,BPA,FREQ' + \
      ' -MEM_MAX 1024 -COMBINE_BUFSIZE 1024 -XML_NAME swarp.xml')


def get_fits_all_pix_statistics(fits_in, extension_num, action):
    fits_open = fits.open(fits_in)
    fits_data = fits_open[extension_num].data

    # remove nan
    fits_data_nan_removed = fits_data[~np.isnan(fits_data)]
    if action == 'avg':
        return np.average(fits_data_nan_removed)
    elif action == 'sum':
        return np.sum(fits_data_nan_removed)
    elif action == 'std':
        return np.std(fits_data_nan_removed)
    else:
        return 0


def rm_blank_txt(txt_in, txt_ou):
    """
    remove the duplicated spaces in text
    :param txt_in: the input text file
    :param txt_ou: the output text file
    :return: 0
    """
    txt_in_open = open(txt_in, 'r').readlines()
    txt_ou_open = open(txt_ou, 'w')

    for line in txt_in_open:
        line_joined = " ".join(line.split())
        txt_ou_open.write(line_joined+'\n')
    return 0


def cal_si(data1, data2, freq1, freq2):
    # cal the spectral index
    # data1 and data2 can be integrated flux density or brightness temperature
    # freq1 < freq2
    si = np.log10(data2 / data1) / np.log10(freq2 / freq1)
    return si


def si_map(fits_in1, fits_in2, fits_ou, freq1, freq2):
    """
    get the spectral index map from two fits file
    :param fits_in1: the first fits image
    :param fits_in2: the second fits image
    :param fits_ou: the output si map
    :param freq1: the frequency of the first image
    :param freq2: the frequency of the second image
    :return: 0
    """
    fits_in1_open = fits.open(fits_in1, dtype='float64')
    fits_in2_open = fits.open(fits_in2, dtype='float64')

    fits_in1_head = fits_in1_open[0].header

    fits_in1_data = fits_in1_open[0].data
    fits_in2_data = fits_in2_open[0].data

    # get the factor from Jy/beam to K
    fac1 = convert_jyperbeam_k(v=freq1, bmaj=0., bmin=0., fits_img=fits_in1)
    fac2 = convert_jyperbeam_k(v=freq2, bmaj=0., bmin=0., fits_img=fits_in2)

    # convert the Jy/beam to K
    #data1_k = fits_in1_data * fac1
    #data2_k = fits_in2_data * fac2

	# do not convert
    data1_k = fits_in1_data
    data2_k = fits_in2_data

    ratio_k = np.abs(data2_k / data1_k)

    fits_ou_data = np.log10(ratio_k) / np.log10(freq2/freq1)

    # write the spectral index to a fits file
    out = fits.PrimaryHDU()
    out.data = fits_ou_data
    out.header = fits_in1_head
    out.writeto(fits_ou, clobber=1)
    return 0


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()


def blur_array(im, n, ny=None):
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im, g, mode='same')
    return(improc)


def smooth_fits_img(fits_in, fits_ou, n):
    # blur image
    # another function to smooth fits image is change_beam_convol_miriad
    fits_open = fits.open(fits_in, dtype='float64')
    fits_data = fits_open[0].data
    fits_data_smooth = blur_array(fits_data, n)
    out = fits.PrimaryHDU()
    out.data = fits_data_smooth
    out.header = fits_open[0].header
    out.writeto(fits_ou, clobber=1)


def change_beam_convol_miriad(fits_in, fits_ou, fits_tp=None, bmaj=None, bmin=None, bpa=None):
# def change_beam_convol_miriad(fits_in, fits_ou, fits_tp=None):
    """
    change the beam size of a fits image using convol from miriad
    this also changed the pixel size, the pixel size in the output fits image is the same with that in the input fits file
    WARNING: Convol does not do a very good job with masked images
    :param fits_in: the input fits image file
    :param fits_ou: the output fits image file with beam changed
    :param fits_tp: the template fits image file for providing the beam information
    :param bmaj: the BMAJ for the output image if fits_tp=None; arcsec
    :param bmin: the BMIN for the output image if fits_tp=None; arcsec
    :param bpa: the BPA for the output image if fits_tp=None; deg
    :return: 0
    """
    # read the beam from the template fits file
    if fits_tp != None:
        fits_tem_open=fits.open(fits_tp)
        fits_tem_head = fits_tem_open[0].header
        bmaj = str(fits_tem_head['BMAJ']*3600.)  # arcsec
        bmin = str(fits_tem_head['BMIN']*3600.)  # arcsec
        bpa = str(fits_tem_head['BPA'])  # deg

    # do the convert
    s('fits op=xyin in=' + fits_in + ' out=$HOME/fits_in.im')
    s('convol map=$HOME/fits_in.im fwhm='+str(bmaj)+','+str(bmin)+' pa='+str(bpa)+' out=$HOME/fits_ou.im options="final"')
    s('fits op=xyout in=$HOME/fits_ou.im out='+fits_ou)
    s('rm -r $HOME/fits_in.im $HOME/fits_ou.im')

    return 0


def add_revise_fits_keyword(fits_in, fits_ou, keyword, value, comment):
    """
    add or revise one keyword in fits header
    :param fits_in: the input fits file
    :param fits_ou: the output fits file with one keyword added or revised
    :param keyword: the name of the keyword; string
    :param value: the value of the keyword; string or float
    :param comment: the comments; string
    :return: 0
    """
    hdulist = fits.open(fits_in)
    prihdr = hdulist[0].header
    prihdr[keyword] = (value, comment)
    hdulist.writeto(fits_ou, clobber=1)
    hdulist.close()
    return 0


def cal_overlapping_uv_range(freq_s, freq_e, baseline_min=7.7, baseline_max=2864.):
    """
    calculate the overlapping uv-range from the starting frequency to the ending frequency
    :param freq_s: the starting freq, MHz
    :param freq_e: the ending freq, MHz
    :param baseline_min: m, the minimum baseline of the array, default is 7.7m for the MWA
    :param baseline_max: m, the maximum baseline of the array, default is 2864m for the MWA
    :return: the starting and ending uv-range, in umit of wave number (lamda)
    """
    # add units
    freq_s = freq_s * u.MHz
    freq_e = freq_e * u.MHz
    baseline_min = baseline_min * u.m
    baseline_max = baseline_max * u.m

    # calculate the uv range
    uv_s = (baseline_min / const.c * freq_e).decompose()
    uv_e = (baseline_max / const.c * freq_s).decompose()

    return uv_s, uv_e


def cal_int_flux_density_jy(pix_sum, bmaj, bmin, pix_size):
    """
    calculate the integrated flux density in unit of Jy from the pixel value in unit of Jy/beam
    :param pix_sum: jy/beam; the sum of all the pixel values within the source region
    :param bmaj: deg, the beam major axis
    :param bmin: deg, the beam minior axis
    :param pix_size: arcsec, the pixel size
    :return: the integrated flux density in unit of Jy
    """
    bmaj_deg = bmaj * u.deg
    bmin_deg = bmin * u.deg
    pix_size_arcsec = pix_size * u.arcsec  # arcsec

    # convert the unit
    pix_size_rad = pix_size_arcsec.to(u.rad)
    bmaj_rad = bmaj_deg.to(u.rad)
    bmin_rad = bmin_deg.to(u.rad)

    # calculate the source area and the beam area
    # beam和rad都是不量纲，是无量纲的。从另一个角度想，Jy/Beam的含义是一个beam大小的面积内有多少Jy，当你把各个pixel的值相加后，就变成这些pixel内共有多少Jy，就是你的源有多少Jy。也就是Jy总是对应着一个区域。我们说这个源多少Jy，暗含着说Jy/源面积。
    area_pix = pix_size_rad ** 2
    # fwhm_to_sigma = 1. / (8. * np.log(2.)) ** 0.5
    # area_beam = 2. * np.pi * bmaj_rad * bmin_rad * fwhm_to_sigma ** 2.
    area_beam = 1.1330900354567983 * bmaj_rad * bmin_rad
    ### area_beam = np.pi * bmaj_rad * bmin_rad   # this is wrong! do not use this!

    return float(pix_sum * area_pix / area_beam)


def cut_fits_by_hand(fits_in, fits_ou, start_x, start_y, end_x, end_y):
    """
    Make a cutout of a fits image by hand and save it to the given file
    :param fits_in: input filename of fits image
    :param fits_ou: output filename of fits image
    :param start_x: the ra or galactic longitude in degrees of the left-bottom 
    :param start_y: the dec or galactic latitude in degrees of the left-bottom 
    :param end_x: the ra or galactic longitude in degrees of the top right
    :param end_y: the dec or galactic latitude in degrees of the left-bottom
    """
    data = fits.open(fits_in)
    w = wcs.WCS(data[0].header)

    # define the rectangle u want to cut
    # the 1 is the ra_dec_order : bool, optional, When True will ensure that world coordinates are always given
    # and returned in as (ra, dec) pairs, regardless of the order of the axes specified by the in the CTYPE keywords.
    # Default is False.
    bottom_left = w.wcs_world2pix(start_x, start_y, 1)
    # print bottom_left
    top_right = w.wcs_world2pix(end_x, end_y, 1)

    out = fits.PrimaryHDU()
    out.data = data[0].data[int(bottom_left[1]):int(top_right[1]), int(bottom_left[0]):int(top_right[0])]
    out.header = data[0].header
    out.header.update(w[int(bottom_left[1]):int(top_right[1]), int(bottom_left[0]):int(top_right[0])].to_header())
    out.writeto(fits_ou, clobber=1)


def cut_fits(fits_in, fits_ou, x, y, size_x, size_y):
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.nddata import Cutout2D
    from astropy import units as u
    from astropy.wcs import WCS
    """
    Make a cutout of a fits image using Cutout2D and save it to the given file
    :param fits_in: input filename of fits image
    :param fits_ou: output filename of fits image
    :param x: ra or galactic longitude in degrees
    :param y: dec or galactic latitude in degrees
    :param size_x: the width of the image in degrees
    :param size_y: the height of the image in degrees
    """
    # get the coordinate type
    hdu = fits.open(fits_in)
    coord = hdu[0].header['CTYPE1'][:2]
    if coord == 'GL':
        frame = 'galactic'
    elif coord == 'RA':
        frame = 'icrs'
    else:
        print('cut_fits do not understand your coordinate system')

    position = SkyCoord(x*u.degree, y*u.degree, frame=frame)
    wcs = WCS(hdu[0].header)
    size = u.Quantity((size_y, size_x), u.degree)  # note that it is not (size_x, size_y) here
    cutout = Cutout2D(hdu[0].data, position, size, wcs=wcs)
    # hdu[0].header.update(dict(cutout.wcs.to_header().iteritems()))  # python 2.7
    hdu[0].header.update(dict(cutout.wcs.to_header().items()))  # python 3.6
    hdu[0].data = cutout.data
    hdu.writeto(fits_ou, overwrite=True)
    return 0


def convert_spec_to_img(fits_in, fits_ou):
    """
    create img template from spec data for finding pixels in reg
    :param fits_in: the input fits data cube
    :param fits_ou: the outout fits image
    :return: 0
    """
    # change the capital unit
    fits_in_revised = '/home/aquan/tmp/fits_revised.fits'
    s('cp '+fits_in+' '+fits_in_revised)
    fits_in_open = fits.open(fits_in_revised, mode='update')
    fits_in_head = fits_in_open[0].header
    try:
        if fits_in_head['CUNIT1'] == 'DEGREE':
            fits_in_head['CUNIT1'] = 'deg'
    except:
        print('No such a keyword. All is well.')

    try:
        if fits_in_head['CUNIT2'] == 'DEGREE':
            fits_in_head['CUNIT2'] = 'deg'
    except:
        print('No such a keyword. All is well.')

    try:
        if fits_in_head['CUNIT3'] == 'm/S' or 'M/S':
            fits_in_head['CUNIT3'] = 'm/s'
    except:
        print('No such a keyword. All is well.')
    fits_in_open.flush()

    # open fits files
    fits_in_open = fits.open(fits_in_revised)
    fits_in_head = fits_in_open[0].header
    fits_in_data = fits_in_open[0].data

    w = wcs.WCS(fits_in_head)
    out = fits.PrimaryHDU()
    out.header = fits_in_head

    shape = np.shape(fits_in_data)
    dim = len(shape)
    if dim == 2:
        out.data = fits_in_data
    elif dim == 3:
        out.data = fits_in_data[0, :, :]
    elif dim == 4:
        out.data = fits_in_data[0, 0, :, :]

    out.header.update(w.to_header())

    s('rm ' + fits_ou)
    out.writeto(fits_ou)

    # revise header
    fits_ou_open = fits.open(fits_ou, mode='update')
    fits_ou_head = fits_ou_open[0].header
    keyword = np.array(['CRPIX3', 'CTYPE3', 'CRVAL3', 'CDELT3', 'CUNIT3', 'PC01_03', 'PC02_03', 'PC03_03', 'PC03_01', 'PC03_02', 'WCSAXES', 'SPECSYS', 'CDELT4', 'CRPIX4', 'CTYPE4', 'CRVAL4', 'CUNIT4'])
    for key in keyword:
        try:
            fits_ou_head.remove(key)
        except:
            print('No such keyword to remove. All is well.')

    fits_ou_open.flush()
    return 0


def convert_jyperbeam_k(v,bmaj=0.,bmin=0.,fits_img=None):
    """
    ref: http://docs.astropy.org/en/stable/api/astropy.units.brightness_temperature.html#astropy.units.brightness_temperature
    this code give the same results as miriad imstat
    :param v: MHz, the freqency of observation
    :param bmaj: deg, the beam major angle
    :param bmin: deg, the beam minor angle
    :param fits_img: the fits image to provide the bmaj and bmin
    :return: the conversion factor from Jy/beam to K
    """

    if fits_img != None:
        fits_open = fits.open(fits_img)
        head = fits_open[0].header
        data = fits_open[0].data
        bmaj = head['BMAJ']
        bmin = head['BMIN']

    bmaj = bmaj*u.deg
    bmin = bmin*u.deg
    fwhm_to_sigma = 1./(8.*np.log(2.))**0.5
    beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2.)
    freq = v*u.MHz
    equiv = u.brightness_temperature(beam_area, freq)
    return u.Jy.to(u.K, equivalencies=equiv)


def change_fits_with_a_constant(fits_in, constant, operation, fits_ou):
    fits_in_open = fits.open(fits_in)
    fits_in_head = fits_in_open[0].header
    fits_in_data = fits_in_open[0].data
    if operation == '+':
        fits_ou_data = fits_in_data + constant
    elif operation == '-':
        fits_ou_data = fits_in_data - constant
    elif operation == '*':
        fits_ou_data = fits_in_data * constant
    elif operation == '/':
        fits_ou_data = fits_in_data / constant
    else:
        print('which operation: +, -, *, or / ?')
    fits.writeto(fits_ou, fits_ou_data, fits_in_head, clobber=True)


def k_per_pc_to_physical(v, bmaj, bmin):
    """
    convert the K/pc to W m^-3 Hz^-1 sr^-1
    I use a conversion factor from Jy/beam to K, there should be a more physical way to do this, check later
    :param bmaj: the beam major axis length, deg
    :param bmin: the beam monor axis length, deg
    :return: the conversion factor from K/pc to W m^-3 Hz^-1 sr^-1 
    """
    jyperbeam_k = convert_jyperbeam_k(v=v, bmaj=bmaj, bmin=bmin, fits_img=None)

    # convert deg to rad
    bmaj = bmaj * np.pi / 180.
    bmin = bmin * np.pi / 180.

    Jy = 1E-26  # 1 Jy = 1E-26 W m^-2 Hz^-1; ref: Equation 1.3 on page 8 in Wilson, T. L., K. Rohlfs, and S. Huttemeister (2013). Tools of Radio Astronomy
    fwhm_to_sigma = 1. / (8. * np.log(2.)) ** 0.5
    beam_area = 2. * np.pi * (bmaj * bmin * fwhm_to_sigma ** 2.)  # ref http://docs.astropy.org/en/stable/api/astropy.units.brightness_temperature.html#astropy.units.brightness_temperature
    beam_area = beam_area
    Jy_per_beam = Jy / beam_area

    K = Jy_per_beam / jyperbeam_k
    pc = 3.085677581E+16  # 1 pc = 3.085677581E+16 meter, convert parsec to meter
    K_per_pc = K / pc
    return K_per_pc


def change_pix_size(fits_in, fits_tem, fits_ou):
    # note: change_beam_convol_miriad can also change pixel size while changing the beam size
    # change the pixel size, crop to the same size, using miriad regrid
    # casa imregrid: "The imregrid task currently finds the nearest input pixel center and interpolates to the output pixel center. No averaging is done in any direction!"
    # miriad regrid: "REGRID any combination of the first three axes of an image using cubic interpolation, with exclusion of blanked input pixels."
    s('rm -r image.im template.im image_regrid.im')
    s('fits op=xyin in='+fits_in+' out=image.im')
    s('fits op=xyin in='+fits_tem+' out=template.im')
    s('regrid in=image.im tin=template.im axes=1,2 tol=0 out=image_regrid.im')
    s('fits op=xyout in=image_regrid.im out='+fits_ou)
    s('rm -r image.im template.im image_regrid.im')


def join_npy_files(file_ou, *args):
    """
    join any number of numpy array files
    :param file_ou: the output numpy array file
    :param args: the files of numpy arrays those need to be joined
    :return: 0
    """
    array_all = np.load(args[0])
    for array in args[1:]:
        array_one = np.load(array)
        array_all = np.concatenate((array_all, array_one))
    np.save(file_ou, array_all)
    return 0


def pix_val_in_reg_old(reg_ds9, coord, fits_data, fits_head):
    """
    get the pixel values in the region reg_ds9 
    :param reg_ds9: ds9 region in the format of e.g., 'ellipse(a, b, ...', no header
    :param coord: the coordinate system of the region, e.g. 'galactic', 'icrs', etc
    :param fits_img: the fits image
    :return: an numpy array containing the values of pixels within the region
    """
    n_pix_y, n_pix_x = fits_data.shape

    # create a filter using the region
    r = pyregion.parse(coord+';' + reg_ds9).as_imagecoord(fits_head)
    reg_filter = r.get_filter()

    # read the numbers in the region
    splits = reg_ds9.split(',')

    # degree; the Galactic longitude of the center if the region is ellipse, circle or box;
    # one point if the region is polygon
    glon = float(splits[0].split('(')[1])
    glat = float(splits[1])  # degree
    l_box = 1.5  # degree, to make sure all pix in a region are included
    # create a box to search pix in the region
    w = wcs.WCS(fits_head)
    bottom_left = w.wcs_world2pix(glon + l_box, glat - l_box, 1)
    up_right = w.wcs_world2pix(glon - l_box, glat + l_box, 1)

    # get the mean; avoid the pix number out of bounds
    i_min = int(max(bottom_left[0], 0))
    i_max = int(min(up_right[0], n_pix_x))
    j_min = int(max(bottom_left[1], 0))
    j_max = int(min(up_right[1], n_pix_y))

    pix_val = []
    for i in range(i_min, i_max):
        for j in range(j_min, j_max):
            if reg_filter.inside1(i, j):
                pix_val.append(fits_data[j, i])

    return np.asarray(pix_val)


def pix_val_in_reg(reg_ds9, coord, fits_data, fits_head):
    """
    get the pixel values in the region reg_ds9 
    :param reg_ds9: ds9 region in the format of e.g., 'ellipse(a, b, ...', no header
    :param coord: the coordinate system of the region, e.g. 'galactic', 'icrs', etc
    :param fits_img: the fits image
    :return: an numpy array containing the values of pixels within the region
    """
    n_pix_y, n_pix_x = fits_data.shape

    # create a filter using the region
    r = pyregion.parse(coord+';' + reg_ds9).as_imagecoord(fits_head)
    reg_filter = r.get_filter()

    # read the numbers in the region
    splits = reg_ds9.split(',')
    # degree; the Galactic longitude of the center if the region is ellipse, circle or box;
    # one point if the region is polygon
    head = splits[0].split('(') 
    region_type = head[0]
    
    if region_type == 'polygon':
        coords = []
        coords.append(head[1])
        n = len(splits)
        for i in range(1, n-1):
            coords.append(splits[i])
        tail = splits[-1].split(')')[0]
        coords.append(tail)

    coords = [float(i) for i in coords]
    glons = coords[0::2]
    glats = coords[1::2]

    glon = np.mean(glons)
    glat = np.mean(glats)

    glon_delta = (max(glons) - min(glons)) / 2.
    glat_delta = (max(glats) - min(glats)) / 2.

    l_box = max(glon_delta, glat_delta) * 1.2

    # create a box to search pix in the region
    w = wcs.WCS(fits_head)
    bottom_left = w.wcs_world2pix(glon + l_box, glat - l_box, 1)
    up_right = w.wcs_world2pix(glon - l_box, glat + l_box, 1)

    # get the mean; avoid the pix number out of bounds
    i_min = int(max(bottom_left[0], 0))
    i_max = int(min(up_right[0], n_pix_x))
    j_min = int(max(bottom_left[1], 0))
    j_max = int(min(up_right[1], n_pix_y))

    pix_val = []
    for i in range(i_min, i_max):
        for j in range(j_min, j_max):
            if reg_filter.inside1(i, j):
                pix_val.append(fits_data[j, i])

    return np.asarray(pix_val)


def is_nan(num):
    return num != num


def replace_all(file, searchExp, replaceExp):
    """
    replace the searchExp with replace Exp in the file
    :param file: file name
    :param searchExp: old text
    :param replaceExp: new text
    :return: nothing
    """
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)


def download_webpage(url, f_save):
    """
    save a webpage into a file
    :param url: the url of the webpage
    :param f_save: the file to save the webpage
    :return: nothing
    """
    url = url
    response = urlopen(url)
    webContent = response.read()

    f = open(f_save, 'w')
    f.write(webContent)
    f.close


def RA_DEC_to_GLong_GLat(RA, DEC):
     """change the coordinates from (RA, DEC) in degree to (GLong, GLat) in degree"""
     # python 2.7
     """
     a = ICRS(ra=RA, dec=DEC, unit=(u.degree, u.degree))
     b = a.transform_to(Galactic)
     """
     # python 3.6
     c_icrs = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree, frame='icrs')
     c_galactic = c_icrs.transform_to('galactic')  # c_icrs.galactic does the same thing
     # return float(c_galactic.l.degree), float(c_galactic.b.degree)
     return c_galactic.l.degree, c_galactic.b.degree


def GLong_GLat_to_RA_DEC(GLong, GLat):
    """change the coordinates from (GLong, GLat) in degree to (RA, DEC) in degree"""
    # python 2.7
    """
    a = Galactic(l=GLong, b=GLat, unit=(u.degree, u.degree))
    b = a.transform_to(ICRS)
    """

    # python 3.6
    c_galactic = SkyCoord(l=GLong * u.degree, b=GLat * u.degree, frame='galactic')
    c_icrs = c_galactic.transform_to('icrs')  # c_galactic.icrs does the same thing
    return float(c_icrs.ra.degree), float(c_icrs.dec.degree)


def cut_fits_downloaded(filename, xc, yc, xw=1, yw=1, units='pixels', outfile=None, clobber=True, useMontage=False, coordsys='celestial', verbose=False):
    """
    credit: http://code.google.com/p/agpy/source/browse/trunk/agpy/cutout.py
    changed by Hongquan on 25 May 2015: move the imports into the defination part
    Generate a cutout image from a .fits file
    Inputs:
        file  - .fits filename or pyfits HDUList (must be 2D)
        xc,yc - x and y coordinates in the fits files' coordinate system (CTYPE)
        xw,yw - x and y width (pixels or wcs); xw and yw is half width of the output file
        units - specify units to use: either pixels or wcs
        outfile - optional output file
    """     
    try:
        import astropy.io.fits as pyfits
        import astropy.wcs as pywcs
    except ImportError:
        import pyfits
        import pywcs
    import numpy
    try:
        import coords
    except ImportError:
        pass # maybe should do something smarter here, but I want agpy to install...
    try:
        import montage_wrapper as montage
        import os
        CanUseMontage=True
    except ImportError:
        CanUseMontage=False

    class DimensionError(ValueError):
        pass

    if isinstance(filename,str):
        file = pyfits.open(filename)
        opened=True
    elif isinstance(filename,pyfits.HDUList):
        file = filename
        opened=False
    else:
        raise Exception("cutout: Input file is wrong type (string or HDUList are acceptable).")

    head = file[0].header.copy()

    if head['NAXIS'] > 2:
        raise DimensionError("Too many (%i) dimensions!" % head['NAXIS'])
    cd1 = head.get('CDELT1') if head.get('CDELT1') else head.get('CD1_1')
    cd2 = head.get('CDELT2') if head.get('CDELT2') else head.get('CD2_2')
    if cd1 is None or cd2 is None:
        raise Exception("Missing CD or CDELT keywords in header")
    wcs = pywcs.WCS(head)

    if units == 'wcs':
        if coordsys=='celestial' and wcs.wcs.lngtyp=='GLON':
            xc,yc = coords.Position((xc,yc),system=coordsys).galactic()
        elif coordsys=='galactic' and wcs.wcs.lngtyp=='RA':
            xc,yc = coords.Position((xc,yc),system=coordsys).j2000()

    if useMontage and CanUseMontage:
        head['CRVAL1'] = xc
        head['CRVAL2'] = yc
        if units == 'pixels':
            head['CRPIX1'] = xw
            head['CRPIX2'] = yw
            head['NAXIS1'] = int(xw*2)
            head['NAXIS2'] = int(yw*2)
        elif units == 'wcs':
            
            cdelt = numpy.sqrt(cd1**2+cd2**2)
            head['CRPIX1'] = xw   / cdelt
            head['CRPIX2'] = yw   / cdelt
            head['NAXIS1'] = int(xw*2 / cdelt)
            head['NAXIS2'] = int(yw*2 / cdelt)

        head.toTxtFile('temp_montage.hdr',clobber=True)
        newfile = montage.wrappers.reproject_hdu(file[0],header='temp_montage.hdr',exact_size=True)
        os.remove('temp_montage.hdr')
    else:

        xx,yy = wcs.wcs_world2pix(xc,yc,0)

        if units=='pixels':
            xmin,xmax = numpy.max([0,xx-xw]),numpy.min([head['NAXIS1'],xx+xw])
            ymin,ymax = numpy.max([0,yy-yw]),numpy.min([head['NAXIS2'],yy+yw])
        elif units=='wcs':
            xmin,xmax = numpy.max([0,xx-xw/numpy.abs(cd1)]),numpy.min([head['NAXIS1'],xx+xw/numpy.abs(cd1)])
            ymin,ymax = numpy.max([0,yy-yw/numpy.abs(cd2)]),numpy.min([head['NAXIS2'],yy+yw/numpy.abs(cd2)])
        else:
            raise Exception("Can't use units %s." % units)

        if xmax < 0 or ymax < 0:
            raise ValueError("Max Coordinate is outside of map: %f,%f." % (xmax,ymax))
        if ymin >= head.get('NAXIS2') or xmin >= head.get('NAXIS1'):
            raise ValueError("Min Coordinate is outside of map: %f,%f." % (xmin,ymin))

        head['CRPIX1']-=xmin
        head['CRPIX2']-=ymin
        head['NAXIS1']=int(xmax-xmin)
        head['NAXIS2']=int(ymax-ymin)

        if head.get('NAXIS1') == 0 or head.get('NAXIS2') == 0:
            raise ValueError("Map has a 0 dimension: %i,%i." % (head.get('NAXIS1'),head.get('NAXIS2')))

        img = file[0].data[ymin:ymax,xmin:xmax]
        newfile = pyfits.PrimaryHDU(data=img,header=head)
        if verbose: print("Cut image %s with dims %s to %s.  xrange: %f:%f, yrange: %f:%f" % (filename, file[0].data.shape,img.shape,xmin,xmax,ymin,ymax))

    if isinstance(outfile,str):
        newfile.writeto(outfile,clobber=clobber)

    if opened:
        file.close()

    return newfile

     
def move_fits(filename, delta_l, delta_b, out_file='moved.fits'):
    """
    delta_l: increment of Galactic longitude in unit of degree
    delta_b: increment of Galactic latitude in unit of degree
    """

    hdulist = pyfits.open(filename)
    prihdr = hdulist[0].header
    prihdr['CRVAL1']+=delta_l
    prihdr['CRVAL2']+=delta_b
    hdulist.writeto(out_file)
    hdulist.close()


def gcd(ra1,dec1,ra2,dec2):
    """
    # copied from Paul Hancock as I remember
    Great circle distance as calculated by the haversine formula
    ra/dec in degrees
    returns:
    sep in degrees"""
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(np.radians(dlat)/2)**2
    a+= np.cos(np.radians(dec1)) * np.cos(np.radians(dec2)) * np.sin(np.radians(dlon)/2)**2
    sep = np.degrees(2 * np.arcsin(min(1,np.sqrt(a))))
    return sep


def template_generate_ds9_reg(reg_ou):
    """
    generate a ds9 region file
    this is a template only; revise it for using
    :param reg_ou: the output ds9 region file
    :return: 0
    """
    # print(np.arange(2., 10, 3))
    cir_glon = np.arange(25, 40, 2)
    cir_glat = np.zeros(len(cir_glon))
    cir_radius = np.ones(len(cir_glon))*1.5  # deg

    f = open(reg_ou, 'w')
    f.write('# Region file format: DS9 version 4.1\n')
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('galactic\n')

    for i in range(0, len(cir_glon)):
        f.write('circle('+str(cir_glon[i])+','+str(cir_glat[i])+','+str(cir_radius[i]*3600)+'")\n')
    f.close()
    return 0


def cal_dt_from_df_glon(glon, d_sun, gr):
    """
    calculate the distance from the Sun to the Galactic edge with the Galactic longitude of glon, d_sun, and gr
    :param glon: the Galactic longitude of the line of sight, deg
    :param d_sun: the distance from the Sun to the Galactic edge; kpc
    :param gr: the Galactocentric radius of the Milky Way; kpc
    :return: the distance from the Sun to the Galactic edge; kpc
    """
    dt = d_sun * np.cos(glon * np.pi / 180.) + np.sqrt(gr**2 - d_sun**2 * np.sin(glon * np.pi / 180.)**2)
    return dt


def cal_R_gc(glon, df, d_sun):
    """
    cal the distance from an HII region to the Galactic center using Galactic longitude and the distance from the 
    HII region to the Sun
    :param glon: the Galactic longitude, deg
    :param df: the distance from the HII region to the Sun, kpc
    :return: the distance from the from the HII region to the Galactic center
    """
    R_gc_squared = d_sun ** 2 + df ** 2 - 2 * d_sun * df * np.cos(glon * np.pi / 180.)
    R_gc = np.sqrt(R_gc_squared)
    return R_gc


def statistical_te_from_Balsar17(glon, df, d_sun):
    """
    cal the electron temperature of HII regions using the statistical relation in Balser2015ApJ...806..199B.pdf
    :return: the electron temperature of an HII region
    """
    R_gc = cal_R_gc(glon, df, d_sun)
    te = 4928. + 385. * R_gc
    te_err = np.sqrt(277.**2 + (29.*R_gc)**2)
    return te, te_err


def cal_chi2(data, data_err, model):
    """
    calculate the chi^2 of the fitting
    :param data: the observed data
    :param data_err: the error of the observed data
    :param model: the values from the fitting results
    :return: the chi^2
    """
    chi2 = np.sum(((data - model) / data_err) ** 2)
    return chi2


def cal_chi2_pearson(data, model):
    """
    calculate the chi^2 of the fitting using Pearson's chi-squared test: 
    https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test
    :param data: the observed data
    :param data_err: the error of the observed data
    :param model: the values from the fitting results
    :return: the chi^2
    """
    chi2 = np.sum(((data - model) / model) ** 2)
    return chi2


def cal_chi2_no_weight(data, model):
    """
    calculate the chi^2 of the fitting                                                                                       
    :param data: the observed data
    :param model: the values from the fitting results
    :return: the chi^2
    """
    chi2 = np.sum(((data - model) / data) ** 2)
    return chi2


def join_fits_tables(fits_list_in, fits_ou, operation):
    stacked = astropy_table.read(fits_list_in[0], format='fits')
    for fits_file in fits_list_in[1:]:
        new_array = astropy_table.read(fits_file, format='fits')
        if operation == 'v':
            stacked = vstack((stacked, new_array))
        elif operation == 'h':
            stacked = hstack((stacked, new_array))
        else:
            stacked = np.nan
            print('How do you want to operate your fits table? vstack or hstack?')

    stacked.write(fits_ou, format='fits', overwrite=True)
    return 0


def join_votable_vstack(vot_list, vot_ou):
    """
    join the votables in the list together vertically (add rows)
    :param vot_list: the list containing the votables
    :return: a votable with all the tables in the list joined together
    e.g.
    join_votable(vot_list=[path_measure+'emi/emi_072_080_088.vot1', path_measure+'emi/emi_072_080_088.vot2',
                           path_measure+'emi/emi_072_080_088.vot3', path_measure+'emi/emi_072_080_088.vot4',
                           path_measure+'emi/emi_072_080_088.vot5', path_measure+'emi/emi_072_080_088.vot6',
                           path_measure+'emi/emi_072_080_088.vot7'], vot_ou=path_measure+'/emi/test.vot')
    """
    v = []
    for vot in vot_list:
        v.append(astropy_table.read(vot, format='votable'))
    vstack(v).write(vot_ou, format='votable')
    return 0


def cal_angular_resolution(diameter, frequency):
    """
    cal the angular resolution of a single dish radio telescope or an array
    for an array, the 'diameter' is the maximum baseline
    :param diameter: meter, diameter of a single dish or the maximum baseline of an array
    :param frequency: MHz
    :return: angular resolution in unit of arcmin
    """
    diameter = diameter * u.m
    frequency = frequency * u.MHz
    wavelength = const.c / frequency
    resolution = np.arctan(1.21966989 * wavelength / diameter)
    return resolution.to(u.arcmin)


def cal_wavelength(frequency):
    """
    cal the wavelength (meter) at the given frequency (MHz)
    :param frequency: MHz
    :return: the wavelength at the given frequency
    """
    frequency = frequency * u.MHz
    wavelength = const.c / frequency
    return wavelength.to(u.m)


def filter_fits_table_template(cat_in, chi2_max, cat_ou):
    tab = Table.read(cat_in, hdu=1)
    tab = tab[tab['chi2'] < chi2_max]
    tab.write(cat_ou, overwrite = True)
    return 0


def compare_speed_template():
    random.seed('slartibartfast')
    s = [random.random() for i in range(1000)]
    timsort = list.sort
    run_time = min(timeit.Timer('a=s[:]; timsort(a)').repeat(7, 1000))
    print('run time is', run_time)


def convert_csv_to_reg_template(csv_in, reg_ou):
    # read csv
    name, cat = [], []
    glon, glat, radius = [], [], []

    with open(csv_in) as f_in:
        csv_read = csv.reader(f_in)
        for line in csv_read:
            name.append(line[0])
            cat.append(line[1])
            glon.append(line[2])
            glat.append(line[3])
            radius.append(line[4])

    # write header for region file
    with open(reg_ou, 'w') as f_ou:
        f_ou.write('# Region file format: DS9 version 4.1\n')
        f_ou.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" '
                   'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        f_ou.write('galactic\n')

        # write regions
        # starting from 1 to exclude the csv header
        for j in range(1, len(name)):
            region_head = 'circle(' + glon[j] + ',' + glat[j] + ',' + radius[j] + '") # color='
            # region_tail = ' text={' + name[j] + '} \n'
            region_tail = ' \n'
            if cat[j] == 'K':
                color = 'red'
                region = region_head + color + region_tail
            elif cat[j] == 'G':
                color = 'green'
                region = region_head + color + region_tail
            elif cat[j] == 'C':
                color = 'cyan'
                region = region_head + color + region_tail
            elif cat[j] == 'Q':
                region = region_head + color + region_tail
            else:
                print('what color is this region?')
            f_ou.write(region)


def use_multiple_cores_tamplate():
    from multiprocessing import Pool
    def f(x):
        return x * x

    if __name__ == '__main__':
        with Pool(5) as p:
            print(p.map(f, [1, 2, 3]))


def change_names(path, old_name_s, old_name_e, new_name_s, new_name_e):
    for i in range(1, 99999):
        file_old = path + old_name_s + str(i) + old_name_e
        if os.path.isfile(file_old):
            file_new = path + new_name_s + str(i) + new_name_e
            s('cp ' + file_old + ' ' + file_new)
        else:
            break
    return 0


def cal_hii_emission_measure(vt, te):
    """
    Calculate the emission measure of the HII region using the relation
    vt = [0.082 * T^(-1.35) * EM]^0.476
    from Kurtz2005IAUS..227..111K, Mezger1967ApJ...147..471M,
    or Hindson2016PASA...33...20H
    :param vt: the turnover frequency in unit of GHz
    :param te: the electron temperature in unit of Kelvin
    :return: the emission measure in unit of pc cm^-6
    """

    # constants in the equation
    c1 = 0.082
    c2 = -1.35
    c3 = 0.476
    # calculate the EM:
    em = vt ** (1 / c3) / c1 / (te ** c2)
    em = em * u.pc * u.cm ** (-6)
    return em


def cal_radius(ang, dis):
    # calculate the radius of a source in unit of pc
    # ang is the angular size of the source in unit of rad
    # dis is the distance from the source to us in unit of pc
    dis = dis * u.pc
    radius = ang * dis
    return radius


def cal_hii_electron_density(em, radius):
    """
    Calculate the electron density in the HII region from the Emission measure
    :param em: the Emission Measure of the HII region in unit of pc cm^(-6)
    :param radius: the radius of the HII region in unit of pc
    :return: the electron density of the HII region in unit of cm^(-3)
    """
    # give inputs units
    em = em * u.pc * u.cm ** (-6)
    radius = radius * u.pc

    # calculate the electron density
    electron_density = np.sqrt(em / radius)
    return electron_density


def cal_hii_mass(ne, radius):
    """
    Calculate the mass of the ionised gas in the HII region from the electron density
    :param ne: the electron density in the HII region in unit of cm^(-3)
    :param radius: the radius of the HII region in unit of pc
    :return: the mass of the HII region relevant to the Sun
    """
    ne = ne * u.cm ** (-3)
    radius = radius * u.pc

    hii_mass = ne * (const.m_e + const.m_p) * (4. / 3.) * np.pi * radius ** 3
    hii_mass_relevant = hii_mass / const.M_sun
    hii_mass_relevant = hii_mass_relevant.decompose()
    return hii_mass_relevant


def cal_ionising_Lyman_continuum_photon_count(s_int, dis, freq):
    # ionising (Lyman continuum) photon count (N_Ly) in unit of photon s^(-1)
    # s_int: mJy; dis: kpc; freq: GHz
    # ref page 14 in Carpenter1990ApJ...362..147C
    # also ref section 3.4.1 in Hindson2013MNRAS.435.2003H

    s_int = s_int * u.mJy
    dis = dis * u.kpc
    freq = freq * u.GHz
    n_ly = 9e+43 * s_int * dis ** 2 * (freq / 5) ** 0.1

    # make the unit correct
    # (1/10) in n_ly does not affect the final result
    n_ly = n_ly.value * u.s ** (-1)

    return n_ly


def get_def_list():
    file_tool = '/Users/hongquansu/Dropbox/bin/tools.py'
    def_list = []
    with open(file_tool) as f:
        for line in f:
            if line.startswith('def'):
                def_name = line.split()[1].split('(')[0]
                def_list.append(def_name)
    def_list_sorted = sorted(def_list)
    for d in def_list_sorted:
        print(d)


if __name__ == '__main__':
    a = 1
    # get_def_list()
    a = cal_ionising_Lyman_continuum_photon_count(s_int=9800., dis=2.4, freq=0.088)
    print(a)
