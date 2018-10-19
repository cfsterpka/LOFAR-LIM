#!/usr/bin/env python3

#python
import time
from os import mkdir
from os.path import isdir, isfile

##external
import numpy as np
from matplotlib import pyplot as plt

## mine
from utilities import log, processed_data_dir
from read_PSE import read_PSE_timeID
from porta_code import code_logger, pyplot_emulator


if __name__=="__main__":
    ##opening data
    timeID = "D20170929T202255.000Z"
    output_folder = "plot_correlateSSPW_PSE_traces"
    
    PSE_folder = "correlate_SSPW"
    
        
    #### setup directory variables ####
    processed_data_dir = processed_data_dir(timeID)
    
    data_dir = processed_data_dir + "/" + output_folder
    if not isdir(data_dir):
        mkdir(data_dir)
        
    logging_folder = data_dir + '/logs_and_plots'
    if not isdir(logging_folder):
        mkdir(logging_folder)

    #Setup logger and open initial data set
    log.set(logging_folder + "/log_out.txt") ## TODo: save all output to a specific output folder
    log.take_stderr()
    log.take_stdout()
    
    
    #### print details
    print("Time ID:", timeID)
    print("output folder name:", output_folder)
    print("date and time run:", time.strftime("%c") )
    print("PSE folder:", PSE_folder)
    
    
    
    PSE_info = read_PSE_timeID(timeID, PSE_folder)
    PSE_list = PSE_info["PSE_list"]
    ant_locs = PSE_info["ant_locations"]
    
    
    
    
    
    PSE_ID_to_plot = np.arange(24)
    PSE_to_plot = [PSE for PSE in PSE_list if PSE.unique_index in PSE_ID_to_plot]
    
    for PSE in PSE_to_plot:
        fname = data_dir + '/PSE_'+str(PSE.unique_index)
        print(PSE.unique_index, fname)
        PSE.plot_trace_data(ant_locs,  plotter=plt)
                            #pyplot_emulator(code_logger(fname))   )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    