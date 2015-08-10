__author__ = 'cmccully'

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astropy import wcs
from scipy.optimize import curve_fit
from astroscrappy import detect_cosmics
import scipy.stats

import utilities
import os

import numpy as np
import pylab as plt
plt.rcParams['figure.facecolor']='white'

def produce_metadata(inputs):
    '''
    Collate image meta data
    This routine takes a list of FITS files, collates the metadata, and outputs them to an astropytable for quick reference.
    Create an astropy table, all entries are saved as strings
    '''
    print 'Initiating metadata table' 
    metadata = Table(names=('filename', 'MJD-OBS', 'FILTER', 'EXPTIME', 'RA0', 'DEC0'), dtype=('S100', 'f8', 'S100', 'S100', 'S100', 'S100'))
    for FITS_file_name in file(inputs['FILE_LIST']):
        fnm = FITS_file_name.split('\n')[0]
        if(os.path.isfile(fnm)):
   	    hdulist = fits.open(fnm) # the split function is used to ignore the text file line breaks
	    row = [] 
	    row.append(FITS_file_name.split('\n')[0])
	    row.append(float(hdulist[0].header['MJD-OBS']))
	    row.append(hdulist[0].header['FILTER'])
	    row.append(hdulist[0].header['EXPTIME'])
	    row.append(hdulist[0].header['RA'])
	    row.append(hdulist[0].header['DEC'])
	    metadata.add_row(row)
	    hdulist.close()
	if(len(fnm)==0):  # ignore blank lines
		continue
	if('#' in fnm[0]): # ignore commented lines
		continue
        if(not os.path.isfile(fnm)):
            print 'WARNING, THE FOLLOWING FILE DOES NOT EXIST:', fnm
    metadata.sort(['MJD-OBS'])
    out_file = '../data/'+'metadata_'+inputs['TAG']+'.tbl'
    metadata.write(out_file, format = 'aastex')
    return metadata

def loop(metadata_table, inputs):
    '''
    Loop through files and update the metadata table.
    This function produces cosmic ray masks and stores the number of masked pixels in metadata.
    Future updates will fit the PSF and store its values.
    '''
    print 'preprocess loop' 
    '''
    # set up directories for plot outputs
    plot_masks_dir = 'plots/out_%s/pre-process/cr_masks/'%inputs['TAG']
    print 'plot_masks_dir', plot_masks_dir
    if(inputs['PLOTMASKS']):
        os.makedirs(plot_masks_dir)
    '''
    # initialize pre-process loop outputs
    metadata_table['num_masked_pixels']=-1
    metadata_table['read_noise']=-1.0
    metadata_table['read_noise'].format = '%1.2f'

    # loop through image files
    for table_index in range(0,len(metadata_table)):
        # open fits file
        hdulist = fits.open(metadata_table['filename'][table_index])
	# get the world coordinate system for this fits file
        bigw = wcs.WCS(hdulist[0].header)
        # create a cosmic ray mask
        mask = produce_cosmic_ray_masks(hdulist, inputs, fnm = metadata_table['filename'][table_index].split('/')[-1])
        # Trim out bad pixel edges. 
	trimsec = utilities.parse_region_keyword(hdulist[0].header['TRIMSEC'])
        # Add a new column to the metadata table to include the number of pixels masked out of the trimmed field.
        metadata_table['num_masked_pixels'][table_index] = np.sum(mask[trimsec])
        # estimate read noise
        read_noise = estimate_read_noise(hdulist, bigw, display=inputs['PLOTREADNOISE'], out = metadata_table['filename'][table_index].split('/')[-1])
	metadata_table['read_noise'][table_index] = read_noise
        # Close the fits file
        hdulist.close()
	# Update the metadata table
        out_file = '../data/'+'metadata_'+inputs['TAG']+'.tbl'
        metadata_table.write(out_file, format = 'aastex')

# Create bad pixel masks
    # cosmic rays
    # Shutter issues
def produce_cosmic_ray_masks(hdulist, inputs, fnm):
    mask, clean = detect_cosmics(hdulist[0].data, sigfrac=0.15, sigclip=4, objlim=4, cleantype='idw')
    if(inputs['PLOTMASKS']):
        plt.rcParams['figure.facecolor']='white'
        plt.rcParams['font.size']=18
        print 'PLOTTING CR MASKS'
        # plots to test out code
        # print 'TRIMSEC', hdulist[0].header['TRIMSEC']
        # print 'trimsec', trimsec
        # print np.sum(mask), np.size(hdulist[0].data[trimsec].ravel())
        plt.figure(figsize=(23,6))
        plt.suptitle(fnm.split('/')[-1].split('.')[0])
	trimsec = utilities.parse_region_keyword(hdulist[0].header['TRIMSEC'])
        ax = plt.subplot(131)
        #print  np.percentile((hdulist[0].data), 99.9)
        #plt.imshow(hdulist[0].data[trimsec], interpolation = 'none', cmap='jet', vmin = np.percentile(hdulist[0].data[trimsec], 0.5), vmax = np.percentile(hdulist[0].data[trimsec], 99.5))
        plt.imshow(np.log10(hdulist[0].data[trimsec]), interpolation = 'none', cmap='jet', vmin = np.log10(np.percentile(hdulist[0].data[trimsec], 0.5)), vmax = np.log10(np.percentile(hdulist[0].data[trimsec], 99.5)))
	plt.ylim(4096,0)
	plt.xlim(0, 4096)
        plt.colorbar()
        plt.subplot(132, sharex=ax, sharey=ax)
        #plt.imshow(np.ma.masked_where(mask[trimsec]==1,hdulist[0].data[trimsec]), interpolation = 'none', cmap='jet', vmin = np.percentile(hdulist[0].data[trimsec], 0.5), vmax = np.percentile(hdulist[0].data[trimsec], 99.5))
        plt.imshow(np.log10(np.ma.masked_where(mask[trimsec]==1,hdulist[0].data[trimsec])), interpolation = 'none', cmap='jet', vmin = np.log10(np.percentile(hdulist[0].data[trimsec], 0.5)), vmax = np.log10(np.percentile(hdulist[0].data[trimsec], 99.5)))
        plt.colorbar()
        plt.xlabel('Pixel (%1.2f arcsec/pixel)'%hdulist[0].header['PIXSCALE'])
        plt.subplot(133, sharex=ax, sharey=ax)
        plt.imshow(mask[trimsec], interpolation = 'none', cmap='gray_r')
        plt.colorbar()
        #plt.subplot(224, sharex=ax, sharey=ax)
        #plt.imshow(clean[trimsec], interpolation = 'none', cmap='jet', vmin = np.percentile(clean, 0.5), vmax = np.percentile(clean, 99.5))
        #plt.colorbar()
	plt.subplots_adjust(left=0.05, right=0.95)
        plt.show()
        #plt.savefig(outfnm)
    return mask


#  Create a noise model
    # Estimate the read noise
def sub_image(ra, dec, bigw, hdulist, pixels, mask='none'):
    x,y = bigw.wcs_world2pix(ra,dec,1)
    if(pixels%2==0):
	    #print 'even'
	    sub_image = hdulist[0].data[y-pixels/2:y+pixels/2,x-pixels/2:x+pixels/2]
    if(pixels%2==1):
	    #print 'odd'
	    sub_image = hdulist[0].data[y-(pixels-1)/2:y+(pixels-1)/2+1,x-(pixels-1)/2:x+(pixels-1)/2+1]
    return sub_image
    #plt.figure()


def gaussian(x, amp, mu, sig):
    return amp*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def rms(array):
	return np.sqrt( np.mean(array**2) - np.mean(array)**2 )

def estimate_read_noise(hdulist, bigw, Npx = 31, display=False, out = 'out'):
    # get valid image limits
    xmin = int((hdulist[0].header['TRIMSEC'])[1:-1].split(',')[0].split(':')[0])
    xmax = int((hdulist[0].header['TRIMSEC'])[1:-1].split(',')[0].split(':')[1])
    ymin = int((hdulist[0].header['TRIMSEC'])[1:-1].split(',')[1].split(':')[0])
    ymax = int((hdulist[0].header['TRIMSEC'])[1:-1].split(',')[1].split(':')[1])
    # initialize the array of read noise estimates
    readnoise_array = []
    # initialize the sample counter
    count = 0
    trials=0
    while(count<1000): # sample 1000 sub-images for the fits file
        # if no acceptable noise estimates are found below, then bail after 20000 trials
        trials += 1
        if(trials>20000):
            return -1
        # sample sub-images of size Npx by Npx for random locations in the fits file image image
        x = np.random.randint(xmin+Npx, xmax-Npx)
        y = np.random.randint(ymin+Npx, ymax-Npx)
        r, d = bigw.wcs_pix2world(x,y, 1)
        img = sub_image(r, d, bigw, hdulist, Npx)

        # Only accept sub-images if they look like Gaussian noise.
        # Cuts on kurtosis, skewness, and maximum value compared to median were determined from simulated Gaussian distributions. The values are selected based on the bounds that reteain 99.9% of the simulated data.
        if(scipy.stats.kurtosis(np.ravel(img))<-0.42 or scipy.stats.kurtosis(np.ravel(img))>0.66 or np.abs(scipy.stats.skew(np.ravel(img)))>0.26 or np.max(img)>np.median(img)+4.5*np.sqrt(np.median(img))):
            continue
        # Next we check for gradient effects. We divide the sub-image into quadrants and estimate their means and rms.
        # We find the maximum difference between means divided by the root-sum-squared of their standard deviations.
        # Finally we cut on the maximum value. This value was determined from a simulation of Guassian noise and finding the bound that retained 99.9% of the data.
        mQ1 = np.mean(img[:(len(img)-1)/2, :(len(img)-1)/2])
        mQ2 = np.mean(img[:(len(img)-1)/2, (len(img)-1)/2:])
        mQ3 = np.mean(img[(len(img)-1)/2:, :(len(img)-1)/2])
        mQ4 = np.mean(img[(len(img)-1)/2:, (len(img)-1)/2:])
        rQ1 = rms(img[:(len(img)-1)/2, :(len(img)-1)/2])
        rQ2 = rms(img[:(len(img)-1)/2, (len(img)-1)/2:])
        rQ3 = rms(img[(len(img)-1)/2:, :(len(img)-1)/2])
        rQ4 = rms(img[(len(img)-1)/2:, (len(img)-1)/2:])
        normed_quad_diffs = [ (mQ1-mQ2)/np.sqrt(rQ1**2 + rQ2**2), (mQ1-mQ3)/np.sqrt(rQ1**2 + rQ3**2), (mQ1-mQ4)/np.sqrt(rQ1**2 + rQ4**2), (mQ2-mQ3)/np.sqrt(rQ2**2 + rQ3**2), (mQ2-mQ4)/np.sqrt(rQ2**2 + rQ4**2), (mQ3-mQ4)/np.sqrt(rQ3**2 + rQ4**2)]
        if( np.max(np.abs(normed_quad_diffs)) > 0.245):
            continue

        # histogram the values of the image
        h,b = np.histogram(img.flat, bins=30)
        v=b[:-1]+(b[1]-b[0])/2.
        # Fit a Gaussian, if the fit fails, then don't accept this trial as a noise sample.

        '''
        plt.figure()
        plt.imshow(hdulist[0].data, interpolation='none')
        plt.colorbar()
        plt.figure()
        plt.subplot(121)
        plt.imshow(img, interpolation  = 'none')
        plt.colorbar()
        plt.subplot(122)
        plt.step(v,h,where='mid')
        plt.show()
        '''
        try:
            p0 = [np.max(h), np.mean(img), np.sqrt( np.mean(img**2)-np.mean(img)**2 )]
            popt, pcov = curve_fit(gaussian, (v), h, p0=p0)
        except:
            continue
        if(popt[0]<=0. or popt[1]<0):
            continue
        if(popt[1]**2 - popt[2] <= 0.):
            continue
        # estimate the readnoise for this subimage based on a Gaussian fit.
        sigma_CCD = popt[2]
        mean_CCD  = popt[1]
        gain = float(hdulist[0].header['GAIN'])

        # If the readnoise calculation is going to fail, then continue
        if(sigma_CCD**2 < mean_CCD/gain):
            continue
        rn = np.sqrt(sigma_CCD**2 - mean_CCD/gain)
        #print rn
        # increment counter only if image is accepted, the gaussian fit is successful, and the readnoise estimate is sensible	
        count += 1
        readnoise_array.append(rn)

    # convert readnoise_array from a list to an array
    readnoise_array = np.array(readnoise_array)
    # estimate the readnoise by using the median of the distribution
    readnoise = np.median(readnoise_array)
    #print np.min(readnoise_array), np.max(readnoise_array)
    #print np.median(readnoise_array), np.mean(readnoise_array), np.sqrt(np.mean(readnoise_array**2) - np.mean(readnoise_array)**2)
    if(display):
        plt.figure(figsize=(9,7))
        plt.subplot(221)
        plt.imshow(img, interpolation='none')
        plt.colorbar()
        plt.title('Random Image Portion')
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.subplot(222)
        plt.step(v,h, lw=2)
        plt.plot(v,gaussian(v,*popt), 'r--', lw=2) 
        plt.title('Distribution of Counts\nw/ Gaussian Fit')
        plt.xlabel('CCD Counts (gain corrected)')
        plt.subplot(212)
        print 'median read noise', np.median(readnoise_array)
        a,b,c = plt.hist(readnoise_array, bins=int(np.max(readnoise_array)-np.min(readnoise_array)))
        plt.plot([np.median(readnoise_array), np.median(readnoise_array)],[0., 1.2*max(a)],'r-', lw=2, label='median')
        plt.legend(loc=1)
        plt.ylim(0.,1.2*max(a))
        plt.xlabel('Read Noise, CCD counts (gain corrected)')
        plt.title('Multiple Read Noise Estimates')
        plt.suptitle(out.split('/')[-1], fontsize=20)
        plt.subplots_adjust(top=0.85, wspace=0.3, hspace=0.4)
        #plt.savefig(out+'.png', dpi=50)
        plt.savefig(out+'.pdf')
        #plt.show()
    return readnoise

    # Use the read noise and the image data to make a noise map

# Solve for WCS tweak
    # fit each star in the catalogue to a Moffat
    # get the image centroids and Moffat parameters.
    # can scale subimages according to the Moffat FWHM
