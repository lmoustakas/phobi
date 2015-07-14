__author__ = 'cmccully'

import preprocess
import argparse
import datetime


def main():
    parser=argparse.ArgumentParser(description='phobi routine crowded field photometry')
    parser.add_argument("-l","--file_list", default='../data/FITS_file_list.txt', help="directory with FITS file images",type=str)

    args=parser.parse_args() 
    date_string = datetime.datetime.now().strftime("_%Y_%m_%d_%H:%M:%S") # get current time 
    tag = ''.join([(args.file_list).split('/')[-1].split('.')[0], date_string]) # tag data by the list file name with date appended
    # Pre-process images
    preprocess.produce_meta_data(args.file_list, tag)
    # Photometer image
    # Make light curve
    # Infer Time delay
    return

if __name__ == '__main__':
    main()
