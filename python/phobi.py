__author__ = 'cmccully'

import preprocess
import argparse
import datetime


def main():
    # argument parsing
    parser=argparse.ArgumentParser(description='phobi routine crowded field photometry')
    parser.add_argument("-i","--inputs_file", default='./inputs.txt', help="text file with inputs",type=str)
    args=parser.parse_args() 

    # read the input file and return a dictionary
    inputs = eval(open(args.inputs_file).read())
   
    print '========== INPUT FILE:  %s  ==============='%args.inputs_file
    for line in file(args.inputs_file):
        print line.split('\n')[0]
    print '==============================================='
    # if DATE_APPEND flag is True, append the output filename tage with the current date and time.
    # this is useful to prevent overwriting outputs
    if(inputs['DATE_APPEND']): inputs['TAG'] = ''.join([inputs['TAG'], datetime.datetime.now().strftime("_%Y_%m_%d_%H:%M:%S")])

    # Pre-process images
    if(inputs['PREPROCESS']==True):
        metadata = preprocess.produce_metadata(inputs)
        preprocess.loop(metadata, inputs)


    # Make light curve
    # Infer Time delay
    return

if __name__ == '__main__':
    main()
