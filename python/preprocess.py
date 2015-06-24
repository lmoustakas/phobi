__author__ = 'cmccully'

# Collate image meta data
    # Filename, MJD, FILTER, EXPTIME, RA0, DEC0

# Create bad pixel masks
    # cosmic rays
    # Shutter issues

# Create a noise model
    # Estimate the read noise
    # Use the read noise and the image data to make a noise map
    # TODO Flag strong background gradients

# Generate reference catalogs

# Solve for WCS tweak
    # fit each star in the catalogue to a Moffat
    # get the image centroids and Moffat parameters.
    # can scale subimages according to the Moffat FWHM


