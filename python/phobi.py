__author__ = 'cmccully'

import preprocess
import argparse
import datetime

def read_inputs(filename):
    # Create inputs dictionary
    inputs = {'FILE_LIST' : '', 'TAG' : '', 'DATE_APPEND' : False, 'PREPROCESS' : False, 'PHOTOMETER' : False, 'INFERDELAY' : False, }
    for line in file(filename):
        if(line[0]!='#'):
            if('FILE_LIST' in line):
                inputs['FILE_LIST'] = line.split()[2] # parse the file list name
            if('TAG' in line):
                inputs['TAG'] = line.split()[2] # parse the tage name
            if('DATE_APPEND' in line):
	        if('True' in line): inputs['DATE_APPEND'] = True # set date append to true
            if('PREPROCESS' in line):
	        if('True' in line): inputs['PREPROCESS'] = True # set date append to true
    if(inputs['DATE_APPEND'] == True):
        date_string = datetime.datetime.now().strftime("_%Y_%m_%d_%H:%M:%S") # get current time 
        inputs['TAG'] = ''.join([inputs['TAG'], date_string]) # tag data by the list file name with date appended
    return inputs

def main():
    # argument parsing
    parser=argparse.ArgumentParser(description='phobi routine crowded field photometry')
    parser.add_argument("-i","--inputs_file", default='./inputs.txt', help="text file with inputs",type=str)
    args=parser.parse_args() 
    # read the input file and return a dictionary
    inputs = read_inputs(args.inputs_file)
    # Pre-process images
    if(inputs['PREPROCESS']==True):
        preprocess.produce_meta_data(inputs['FILE_LIST'], inputs['TAG'])
    # Photometer image
    # Make light curve
    # Infer Time delay
    return

if __name__ == '__main__':
    main()
