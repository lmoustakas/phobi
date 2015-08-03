__author__ = 'cmccully'

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astroscrappy import detect_cosmics
import utilities
import os

import numpy as np
import pylab as plt

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

    # loop through image files
    for table_index in range(0,len(metadata_table)):
        # open fits file
        hdulist = fits.open(metadata_table['filename'][table_index])
        # create a cosmic ray mask
        mask = produce_cosmic_ray_masks(hdulist, inputs, fnm = metadata_table['filename'][table_index].split('/')[-1])
        # Trim out bad pixel edges. 
	trimsec = utilities.parse_region_keyword(hdulist[0].header['TRIMSEC'])
        # Add a new column to the metadata table to include the number of pixels masked out of the trimmed field.
        metadata_table['num_masked_pixels'][table_index] = np.sum(mask[trimsec])
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
	    sub_image = image_data[y-pixels/2:y+pixels/2,x-pixels/2:x+pixels/2]
    if(pixels%2==1):
	    #print 'odd'
	    sub_image = image_data[y-(pixels-1)/2:y+(pixels-1)/2+1,x-(pixels-1)/2:x+(pixels-1)/2+1]
    return sub_image
    #plt.figure()

def estimate_read_noise(self, display=0, out = 'out'):
	readnoise=[]
	for k in range(1000):
		#print self.image_data.shape
		x = np.random.randint(100, self.image_data.shape[0]-100)
		y = np.random.randint(100, self.image_data.shape[1]-100)
		r, d = self.bigw.wcs_pix2world(x,y, 1)
		#print 'x,y',x,y
		#print 'ra,dec',r,d
		img = self.sub_image(r, d, 31)
		#img/=float(self.hdulist[0].header['GAIN']) # remove gain correction. read noise is on electrons, not photons ?		
		count = 0
		
		while(np.max(img)>np.median(img)+5.*np.sqrt(np.median(img))):
			count+=1
	  		#print 'trial', count
			x = np.random.randint(100, self.image_data.shape[0]-100)
			y = np.random.randint(100, self.image_data.shape[1]-100)
			r, d = self.bigw.wcs_pix2world(x,y, 1)
			#print 'x,y',x,y
			#print 'ra,dec',r,d
			img = self.sub_image(r, d, 31)
		
		h,b = np.histogram(img.flat, bins=30)
		x=b[:-1]+(b[1]-b[0])/2.
		try:
			popt, pcov = opt.curve_fit(gauss_1d, (x), h, p0=[np.max(h), np.sqrt(np.median(x)), np.median(x)])
		except:
			continue
		if(popt[0]<=0. or popt[2]<0):
			continue
		if(popt[1]**2 - popt[2] <= 0.):
			continue
		rn = np.sqrt(popt[1]**2 - popt[2])
		if(rn == rn):
			readnoise.append(rn)
	 		#print k,'readnoise', popt[1]-np.sqrt(popt[0])
	 		#print 'median', np.median(img)
			#print popt
	self.readnoise = np.median(readnoise)
	print np.median(readnoise)
	if(display>=2):
		plt.figure(figsize=(9,7))
		plt.subplot(221)
		plt.imshow(img, interpolation='none')
		plt.colorbar()
		plt.title('Random Image Portion')
		plt.xlabel('pixel')
		plt.ylabel('pixel')
		plt.subplot(222)
		plt.step(x,h, lw=2)
		plt.plot(x,gauss_1d(x,*popt), 'r--', lw=2) 
		plt.title('Distribution of Counts\nw/ Gaussian Fit')
		plt.xlabel('CCD Counts (gain corrected)')
		plt.subplot(212)
		print 'median read noise', np.median(readnoise)
		a,b,c = plt.hist(readnoise, bins=int(np.max(readnoise)-np.min(readnoise)))
		plt.plot([np.median(readnoise), np.median(readnoise)],[0., 1.2*max(a)],'r-', lw=2, label='median')
		plt.legend(loc=1)
		plt.ylim(0.,1.2*max(a))
		plt.xlabel('Read Noise, CCD counts (gain corrected)')
		plt.title('Multiple Read Noise Estimates')
		plt.suptitle(self.fits_file_name.split('/')[-1], fontsize=20)
		plt.subplots_adjust(top=0.85, wspace=0.3, hspace=0.4)
		plt.savefig(out+'.png', dpi=50)
	#return self.readnoise

    # Use the read noise and the image data to make a noise map

# Solve for WCS tweak
    # fit each star in the catalogue to a Moffat
    # get the image centroids and Moffat parameters.
    # can scale subimages according to the Moffat FWHM
