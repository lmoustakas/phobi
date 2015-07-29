__author__ = 'cmccully'

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astroscrappy import detect_cosmics
import utilities

import numpy as np
import pylab as plt

def produce_meta_data(FITS_file_list, tag):
    '''
    Collate image meta data
    This routine takes a list of FITS files, collates the metadata, and outputs them to an astropytable for quick reference.
    Create an astropy table, all entries are saved as strings
    '''
    meta_data = Table(names=('filename', 'MJD-OBS', 'FILTER', 'EXPTIME', 'RA0', 'DEC0'), dtype=('S100', 'f8', 'S100', 'S100', 'S100', 'S100'))
    for FITS_file_name in file(FITS_file_list):
	hdulist = fits.open(FITS_file_name.split('\n')[0]) # the split function is used to ignore the text file line breaks
	row = [] 
	row.append(FITS_file_name.split('\n')[0])
	row.append(float(hdulist[0].header['MJD-OBS']))
	row.append(hdulist[0].header['FILTER'])
	row.append(hdulist[0].header['EXPTIME'])
	row.append(hdulist[0].header['RA'])
	row.append(hdulist[0].header['DEC'])
	meta_data.add_row(row)
        #for key in  hdulist[0].header.keys():
            #print key,'\t',hdulist[0].header[key]
	#exit()
	hdulist.close()
    meta_data.sort(['MJD-OBS'])
    out_file = '../data/'+'meta_data_'+tag+'.tbl'
    meta_data.write(out_file, format = 'aastex')
    return meta_data

def loop(meta_data_table, tag):
    # Store meta data from the loop
    #print meta_data_table['filename']
    #exit()
    counter = 0
    # want to add columns to the metadata file. overwrite as it runs
    for FITS_file_name in meta_data_table['filename']:
        hdulist = fits.open(FITS_file_name.split('\n')[0])
        mask = produce_cosmic_ray_masks(hdulist, tag)
	trimsec = utilities.parse_region_keyword(hdulist[0].header['TRIMSEC'])
        np.sum(mask), np.sum(mask[trimsec])
        t = Time( meta_data_table['MJD-OBS'][counter], format='mjd' )
        print FITS_file_name, meta_data_table['MJD-OBS'][counter], np.sum(mask), np.sum(mask[trimsec])
        hdulist.close()
	counter+=1

# Create bad pixel masks
    # cosmic rays
    # Shutter issues
def produce_cosmic_ray_masks(hdulist, tag):
    #hdulist = fits.open(FITS_file_name.split('\n')[0]) # the split function is used to ignore the text file line breaks
    mask, clean = detect_cosmics(hdulist[0].data, sigfrac=0.15, sigclip=4, objlim=4, cleantype='idw')
    #if(trim==False):
    #    mask, clean = detect_cosmics(hdulist[0].data, sigfrac=0.15, sigclip=4, objlim=4, cleantype='idw')
    #if(trim==True):
    #    trimsec = utilities.parse_region_keyword(hdulist[0].header['TRIMSEC'])
    #    mask, clean = detect_cosmics(hdulist[0].data[trimsec], sigfrac=0.15, sigclip=4, objlim=4, cleantype='idw')
    '''
    # plots to test out code
    print 'TRIMSEC', hdulist[0].header['TRIMSEC']
    print 'trimsec', trimsec
    print np.sum(mask), np.size(hdulist[0].data[trimsec].ravel())
    plt.figure()
    ax = plt.subplot(221)
    print  np.percentile((hdulist[0].data), 95.)
    plt.imshow(hdulist[0].data[trimsec], interpolation = 'none', cmap='jet', vmin = np.percentile(hdulist[0].data[trimsec], 1.), vmax = np.percentile(hdulist[0].data[trimsec], 95.))
    plt.colorbar()
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(mask, interpolation = 'none', cmap='gray_r')
    plt.colorbar()
    plt.subplot(223, sharex=ax, sharey=ax)
    plt.imshow(clean, interpolation = 'none', cmap='jet', vmin = np.percentile(clean, 1.), vmax = np.percentile(clean, 95.))
    plt.colorbar()
    plt.subplot(224, sharex=ax, sharey=ax)
    plt.imshow(np.ma.masked_where(mask==1,hdulist[0].data[trimsec]), interpolation = 'none', cmap='jet', vmin = np.percentile(hdulist[0].data[trimsec], 1.), vmax = np.percentile(hdulist[0].data[trimsec], 95.))
    plt.colorbar()
    plt.show()
    '''
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
