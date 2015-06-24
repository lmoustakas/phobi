__author__ = 'cmccully'

# For an individual image

# Get the reference catalog (include positional uncertainties)

# Get Photometer targets

# Fit an initial set of parameters to the data.

# Run an emcee loop to get distributions of photometry
def moffat_chi_vals(self,theta,N_pix, flip):
    x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg = theta
    model = self.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip)
    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print np.sum(chi*chi)
    #print x0,y0, N_pix/2
    # NO NEGATIVE AMPLITUDES ALLOWED!
    if(amp0<0 or amp1<0 or amp2<0 or amp3<0):
	return np.inf
    # THE IMAGE CORRECTION HAS TO BE WITHIN THE BOUNDS OF THE IMAGE!
    if(x0>N_pix/2 or x0<-N_pix/2 or y0>N_pix/2 or y0<-N_pix/2):
	return np.inf
    return chi

def moffat_chi_sq(self, theta, N_pix, flip):
	chisq = np.sum((self.moffat_chi_vals(theta, N_pix, flip))**2)
	#print '%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e'%(chisq, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], theta[7], theta[8]) 
	return chisq

def emceeQuadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha, beta, m1, m2, m3, m4, N_px, outputFileTag='out'):
  print '\n\t####################################'
  print   '\t# STARTING EMCEE ###################'
  print   '\t####################################\n'
  obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)
  N_bkg_guess = np.median(obj.image)

  x0,y0, amp1, amp2, amp3, amp4, alpha_qsr, beta_qsr, N_bkg_qsr =FM.qsr_min_parms
  #amp1 = 10**((m1-ZP_mean)/(-2.5))
  #amp2 = 10**((m2-ZP_mean)/(-2.5))
  #amp3 = 10**((m3-ZP_mean)/(-2.5))
  #amp4 = 10**((m4-ZP_mean)/(-2.5))
  #N_bkg = np.median(obj.image)

  #obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)

  def lnprior(_theta):
    _x0, _y0, \
    _amplitude0, \
    _amplitude1, \
    _amplitude2, \
    _amplitude3, \
    _alpha, \
    _beta, \
    _N_bkg = _theta
    #print 'N_bkg_guess', N_bkg_guess
    #exit()
    if  np.abs(_x0)<N_px \
	and np.abs(_y0)<N_px \
	and _amplitude0>0. and _amplitude0<100.*amp1\
	and _amplitude1>0. and _amplitude1<100.*amp2\
	and _amplitude2>0. and _amplitude2<100.*amp3\
	and _amplitude3>0. and _amplitude3<100.*amp4\
	and _alpha>0. and _alpha < 100.*alpha\
	and _beta>0. and  _beta  < 100.*beta \
	and _N_bkg>0. and _N_bkg < 100.*N_bkg_guess:
	  return 0.0  
    else:
        return -np.inf

  def lnprob(theta, obj, N_pix, flip):
    lp = lnprior(theta)
    #print lp, theta
    if not np.isfinite(lp):
        return -np.inf
    mll = obj.moffat_chi_sq(theta, N_pix, flip)
    #print mll
    if not np.isfinite(mll):
        return -np.inf
    #print '%1.2e\t%1.2e'%(lp,mll), theta
    return lp - mll

  # DETERMINE THE IMAGE ORIENTATION
  xv,yv = FM.bigw.wcs_world2pix(ra_qsr,dec_qsr,1)
  r0,d0 = FM.bigw.wcs_pix2world(xv,yv,1)
  r1,d1 = FM.bigw.wcs_pix2world(xv+1,yv+1,1)
  fl = True
  if(r1-r0<0): fl =False

  #################################################
  #################################################
  #################################################

  # RUN MC!
  theta = FM.qsr_min_parms
  #theta = [x0,y0,amp1,amp2,amp3,amp4, alpha, beta, N_bkg]

  print theta

  ndim =len(theta)
  print 'ndim', ndim
  nwalkers = 100
  n_burn_in_iterations = 10
  n_iterations = 10000

  
  prior_vals=theta
  pos = [prior_vals + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  r=np.random.randn(ndim)
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(obj,N_px, fl))

  #print 'burn-in sampler.run_mcmc'
  sampler.run_mcmc(pos, n_burn_in_iterations)

  samples = sampler.chain[:, int(0.5*float(n_burn_in_iterations)):, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  parms = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  #labels = [ r"$x_0$", r"$y_0$", "amp1", "amp2", "amp3", "amp4", r"$\alpha$", "$\beta$", "N$_{bkg}$"]
  labels = [ "x0", "y0", "amp1", "amp2", "amp3", "amp4", "alpha", "beta", "N_bkg"]
  fig= triangle.corner(samples, labels=labels)
  count=0
  for ax in fig.get_axes():
    count+=1
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    if((count-1)%ndim==0 ): 
	ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5, 0.5)

    if(count>=ndim*(ndim-1)): 
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5,  0.5)

  plt.subplots_adjust(bottom=0.075, left=0.075)
  fig.savefig("%s_burn_in_triangle.png"%(outputFileTag))

  for k in range(0,len(parms)):
    print '%s\t%+1.2e\t%+1.2e\t%+1.2e\t%+1.2e'%(labels[k], theta[k], parms[k][0],parms[k][1],parms[k][2])

  #exit()
  
  sampler.reset()
  sampler.run_mcmc(pos, n_iterations)
  samples = sampler.chain[:, int(0.5*float(n_iterations)):, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  parms = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

  #print parms.shape
  #print parms[k][0]

  for k in range(0,len(parms)):
    print '%s\t%+1.2e\t%+1.2e\t%+1.2e\t%+1.2e'%(labels[k], theta[k], parms[k][0],parms[k][1],parms[k][2])

  fig= triangle.corner(samples, labels=labels)
  count=0
  for ax in fig.get_axes():
    count+=1
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    if((count-1)%ndim==0 ): 
	ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5, 0.5)

    if(count>=ndim*(ndim-1)): 
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5,  0.5)

  plt.subplots_adjust(bottom=0.075, left=0.075)
  fig.savefig("%s_triangle.png"%(outputFileTag))

  np.savez(outputFileTag+'_chains.npz', sampler.chain[:, 0:, :].reshape((-1, ndim)))
  print 'chain saving complete'


# Save the light curve.

###########################

# Likelihood function

#
