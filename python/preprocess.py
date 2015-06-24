__author__ = 'cmccully'

# Collate image meta data
    # Filename, MJD, FILTER, EXPTIME, RA0, DEC0

# Create bad pixel masks
    # cosmic rays
    # Shutter issues

# Create a noise model
    # Estimate the read noise
  def image_piece(self,ra, dec, pixels):
    x,y = self.bigw.wcs_world2pix(ra,dec,1)
    if(pixels%2==0):
	    #print 'even'
	    zoom_data = self.image_data[y-pixels/2:y+pixels/2,x-pixels/2:x+pixels/2]
    if(pixels%2==1):
	    #print 'odd'
	    zoom_data = self.image_data[y-(pixels-1)/2:y+(pixels-1)/2+1,x-(pixels-1)/2:x+(pixels-1)/2+1]
    return zoom_data
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
		img = self.image_piece(r, d, 31)
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
			img = self.image_piece(r, d, 31)
		
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
