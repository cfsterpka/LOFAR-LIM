#!/usr/bin/env python3

#### just plots a block of data on all antennas
#### usefull for health checks

import numpy as np
import matplotlib.pyplot as plt

from LoLIM.utilities import processed_data_dir, v_air
from LoLIM.IO.raw_tbb_IO import MultiFile_Dal1, filePaths_by_stationName, read_station_delays, read_antenna_pol_flips, read_bad_antennas, read_antenna_delays
from LoLIM.signal_processing import remove_saturation
from LoLIM.findRFI import window_and_filter
from matplotlib.transforms import blended_transform_factory

from scipy.signal import hilbert

##1
from matplotlib.widgets import RadioButtons, Button

tmp_guess_delays = {}
stations_active = {}
default_color = '0.85'
station_buttons = []
station_axes = []
stCBs = []
plots = []

set_ref = False
l = []
m = []
Nplots = 5 #number of simultaneous plots

ref_plot = ()
first_few_names = []
ref_plot_offset = 0.0
ref_plot_peak_location = 0.0
set_stations = []
StationCallBacks = []


ax = 0
##~1

def plot_blocks(timeID, block_size, block_starts, guess_delays, guess_location = None, bad_stations=[], polarization_flips="polarization_flips.txt",
                bad_antennas = "bad_antennas.txt", additional_antenna_delays = "ant_delays.txt", do_remove_saturation = True, do_remove_RFI = True,
                positive_saturation = 2046, negative_saturation = -2047, saturation_post_removal_length = 50, saturation_half_hann_length = 5,
                referance_station = "CS002", amplify = {}, omitPOL=None):
    """plot multiple blocks, for guessing initial delays and finding pulses. If guess_location is None, then guess_delays should be apparent delays,
    if guess_location is a XYZT location, then guess_delays should be real delays. If a station isn't in guess_delays, its' delay is assumed to be zero.
    A station is only not plotted if it is bad_stations. If referance station is in guess_delays, then its delay is subtract from all stations"""

    '''Modified version of plot_multiple_data_all_stations to support usage of clickable GUI for aligning station timings'''


##2
    global first_few_names
    global defualt_color
    global ax

    text = []
##~2

    if referance_station in guess_delays:
        ref_delay = guess_delays[referance_station]
        guess_delays = {sname:delay-ref_delay for sname,delay in guess_delays.items()}

    processed_data_folder = processed_data_dir(timeID)

    polarization_flips = read_antenna_pol_flips( processed_data_folder + '/' + polarization_flips )
    bad_antennas = read_bad_antennas( processed_data_folder + '/' + bad_antennas )
    additional_antenna_delays = read_antenna_delays(  processed_data_folder + '/' + additional_antenna_delays )


    raw_fpaths = filePaths_by_stationName(timeID)
    raw_data_files = {sname:MultiFile_Dal1(raw_fpaths[sname], force_metadata_ant_pos=True, polarization_flips=polarization_flips, bad_antennas=bad_antennas, additional_ant_delays=additional_antenna_delays) \
                      for sname in raw_fpaths.keys() if sname not in bad_stations}

    if guess_location is not None:
        guess_location = np.array(guess_location)

        ref_stat_file = raw_data_files[ referance_station ]
        ant_loc = ref_stat_file.get_LOFAR_centered_positions()[0]
        ref_delay = np.linalg.norm(ant_loc-guess_location[:3])/v_air - ref_stat_file.get_nominal_sample_number()*5.0E-9

        for sname, data_file in raw_data_files.items():
            if sname not in guess_delays:
                guess_delays[sname] = 0.0

            data_file = raw_data_files[sname]
            ant_loc = data_file.get_LOFAR_centered_positions()[0]
            guess_delays[sname] += (np.linalg.norm(ant_loc-guess_location[:3])/v_air - ref_delay)
            guess_delays[sname] -= data_file.get_nominal_sample_number()*5.0E-9

    RFI_filters = {sname:window_and_filter(timeID=timeID,sname=sname) for sname in raw_fpaths.keys() if sname not in bad_stations}
    #RFI_filters = {sname:window_and_filter(timeID=timeID,sname=sname, blocksize=block_size) for sname in raw_fpaths.keys() if sname not in bad_stations}

    data = np.empty(block_size, dtype=np.double)


##3
    #fig = plt.figure(figsize=(10, 12))
    fig, ax = plt.subplots()

    for sname, data_file in raw_data_files.items():
        #print("Loading: "+ sname)
        if referance_station in sname:
           ref_plot = (sname, data_file)
           #ref_plot_offset = guess_delays[sname] #currently unsed
        else:
           plots.append((sname, data_file))


    first_few = [ref_plot]
    for p in range(0,Nplots-1):
        first_few.append(plots.pop(0))


    for sname, data_file in first_few:
        if sname != referance_station:
           first_few_names.append(sname)
##~3

##4
    def more_plots(some_data,First,ax,text):
       global l,m
       ax.clear()

       height = 0
       t0 = np.arange(block_size)*5.0E-9
       transform = blended_transform_factory(plt.gca().transAxes, plt.gca().transData)
       sname_X_loc = 0.0
       n = 0
       for sname, data_file in some_data:

        print("Plotting: "+sname)

        station_delay = -data_file.get_nominal_sample_number()*5.0E-9
        if sname in guess_delays:
            station_delay = guess_delays[sname]

        station_delay_points = int(station_delay/5.0E-9)

        RFI_filter = RFI_filters[sname]
        ant_names = data_file.get_antenna_names()

        num_antenna_pairs = int( len( ant_names )/2 )
        peak_height = 0.0

        print("Block starts at: "+str(block_starts[0]))
        for point in block_starts:
            T = t0 + point*5.0E-9
            for pair in range(num_antenna_pairs):
                data[:] = data_file.get_data(point+station_delay_points, block_size, antenna_index=pair*2)

                if do_remove_saturation:
                    remove_saturation(data, positive_saturation, negative_saturation, saturation_post_removal_length, saturation_half_hann_length)
                if do_remove_RFI:
                    filtered_data = RFI_filter.filter( data )
                else:
                    filtered_data = hilbert(data)
                even_HE = np.abs(filtered_data)

                data[:] = data_file.get_data(point+station_delay_points, block_size, antenna_index=pair*2+1)
                if do_remove_saturation:
                    remove_saturation(data, positive_saturation, negative_saturation, saturation_post_removal_length, saturation_half_hann_length)
                if do_remove_RFI:
                    filtered_data = RFI_filter.filter( data )
                else:
                    filtered_data = hilbert(data)
                odd_HE = np.abs(filtered_data)

                #ax.plot(T, even_HE + height, 'r')
                #ax.plot(T, odd_HE + height, 'g' )

                #for cases where data are hard to align due to signal being tiny:
                found = False
                if len(amplify) > 0:
                    for key in amplify:
                        if sname == key:

                            if omitPOL != 0:
                               ax.plot(T, amplify[key]*even_HE + height, 'r')
                            if omitPOL != 1:
                               ax.plot(T, amplify[key]*odd_HE + height, 'g' )

                            #ax.plot(T, amplify[key]*even_HE + height, 'r') twinx with block_starts && point as x rather than t!  (can I update to latests plot_multiple?)
                            #ax.plot(T, amplify[key]*odd_HE + height, 'g' )

                            found = True
                if not found:
                    if omitPOL != 0:
                       ax.plot(T, even_HE + height, 'r')
                    if omitPOL != 1:
                       ax.plot(T, odd_HE + height, 'g' )
                    #ax2 = ax.twiny()
                    #xmin = min(block_starts)
                    #xmax = max(block_starts)
                    #ax2.set_xlim(xmin, xmax) #can also use ax2.set_xticks(x)
                    #ax2.xaxis.set_ticks_position('both')

                max_even = np.max(even_HE)
                if max_even > peak_height:
                    peak_height = max_even
                max_odd = np.max(odd_HE)
                if max_odd > peak_height:
                    peak_height = max_odd

#        plt.annotate(sname, (points[-1]*5.0E-9+t0[-1], height))
        plt.sca(ax) #somehow this makes it so that the annotations work on all plots except for the first time "Next stations" is clicked
        plt.annotate(sname, (sname_X_loc, height), textcoords=transform, xycoords=transform)

        ax.ticklabel_format(useOffset=False)
        plt.setp(ax.get_xticklabels(), rotation=60, horizontalalignment='right')
        plt.subplots_adjust(bottom=0.1,top=0.98,left=0.21)
        if ref_plot_peak_location != 0.0:
           print("Setting ref peak to " +str(ref_plot_peak_location))
           Cpeak = ref_plot_peak_location #peak center, for plotting
           Dpeak = 0.0001  # +/- added to peak center for plotting
           plt.xlim(Cpeak-Dpeak,Cpeak+Dpeak)
        plt.grid(b=True)

        height += 2*peak_height
        n+=1
##~4





##5
#######################CallBacks############################
#######################CallBacks############################
#######################CallBacks############################
    def onclick(event):
        global stations_active
        global station_buttons
        global ref_plot_peak_location
        global set_ref
        global set_stations
        global default_color
        ix = event.xdata
        iy = event.ydata
        #print("click active stations: ", stations_active)
        falsify = []

        if set_ref:
           ref_plot_peak_location = ix
           print("Set ref peak to: "+str(ref_plot_peak_location))
           set_ref = False
           return #return, don't set anything else!
        n=0
        for s in stations_active.keys():
          if stations_active[s] and s !="":
            falsify.append(s)
            if s not in guess_delays.keys():
                guess_delays[s] = 0.0
            print ("Set " + s + " station delay to: ", guess_delays[s] + (ix-ref_plot_peak_location))
            tmp_guess_delays[s] = guess_delays[s] + (ix - ref_plot_peak_location)

            station_buttons[n].color = default_color
            station_buttons[n].hovercolor = default_color

          n+=1

        if len(falsify)>0:
           for f in falsify:
               stations_active[f] = False
           falsify = []
        plt.draw()


    def nextCB(event):
        global ax
        global plots
        global station_buttons
        global stations_active
        global StationCallBacks
        print("next!")

        stations_active = {}
        next_few = [ref_plot]
        if len(plots) >= Nplots:
         for i in range(0,Nplots-1):
               next_few.append(plots.pop(0))
        elif len(plots) > 0:
         bk = len(plots)
         for i in range(0,len(plots)):
               next_few.append(plots.pop(0))
         for i in range(bk,4):
             StationCallBacks[i].reset(i,"")
             station_buttons[i].label.set_text("")

        else:
          next_few = []
          for i in range(0,Nplots-1):
             StationCallBacks[i].reset(i,"")
             station_buttons[i].label.set_text("")
          print("All stations timings have been set")
          ax.clear()
          ax.annotate("All stations timings have been set", (0.5, 0.5))
          plt.draw()
          return

        height = 0.0
        n = 0
        for sname, data_file in next_few:

            if sname != referance_station:
               first_few_names.append(sname)
               StationCallBacks[n].reset(n,sname)
               station_buttons[n].label.set_text("Set "+sname)
               stations_active[sname] = False
               n+=1

        more_plots(next_few,False,ax,text)
        plt.draw()


    def refpCB(event):
        global ref_plot_peak_location
        global set_ref
        print ("Setting ref peak")# to: ",ix)
        set_ref = True

    def cCB(text):
        Cpeak = float(text) #peak center, for plotting
        Dpeak = 0.0001  # +/- added to peak center for plotting
        plt.xlim(Cpeak+Dpeak,Cpeak-Dpeak)
        plt.draw()

    def writeCB(event):
        print('Writing guess station delays')
        file = open('guess_delays.py','w')
        file.write("\n      self.guess_timing = " +str(ref_plot_peak_location) + "\n")
        #file.write("      self.event_index = ____\n\n")
        file.write('      self.guess_station_delays = {\n')
        for g in guess_delays:
            if g in tmp_guess_delays:
               '''If the delay has been updated, write the update:'''
               file.write('"'+g+'": '+str(tmp_guess_delays[g])+',\n')
            else:
               """Else, just write the old delay:"""
               file.write('"'+g+'": '+str(guess_delays[g])+',\n')
        file.write('}\n')
        print('Done!\n')


    def press(event):
        val = event.key
        print('press', val)
        '''if val == ('0' or  '1' or '2' or '3' or '4'):
           val = int(val)
           station_buttons[val].color = 'blue'
           station_buttons[val] .hovercolor = 'blue'
        '''


    class StationCB:
        global stations_active
        def __init__(self,np,station):
            self.N = np
            self.station = station

        def reset(self,np,station):
            self.N = np
            self.station = station

        def stationCB(self,event):
           print("Station active: "+self.station)
           #print(self.station + ' You pressed: ',event.name)
           #---highlight button with color when pressed state, show buttons in correct order
           station_buttons[self.N].color = 'blue'
           station_buttons[self.N].hovercolor = 'blue'
           stations_active[self.station] = True
######################~CallBacks############################
######################~CallBacks############################
######################~CallBacks############################

    more_plots(first_few,True,ax,text)

    axnext = plt.axes([ 0.05, 0.80,  0.1, 0.07])
    bnext = Button(axnext, 'Next\nStations')
    bnext.on_clicked(nextCB)

    axfile = plt.axes([0.05,0.7,  0.1,  0.07])
    bfile = Button(axfile, 'Write\nGuesses')
    bfile.on_clicked(writeCB)

    axref = plt.axes([0.05, 0.1,  0.13, 0.07])
    bref= Button(axref, 'Set\nRef Peak')
    bref.on_clicked(refpCB)

    #Cbox = plt.axes([0.05, 0.6,  0.1, 0.07])
    #text_box = TextBox(Cbox, 'Peak Center', initial="")
    #text_box.on_submit(cCB)

    start = 0.22
    for f in range(0,Nplots-1):
        fax = plt.axes([ 0.05, start+f*0.1,  0.13, 0.05])
        fbtn = Button(fax, "Set "+first_few_names[f])
        stations_active[first_few_names[f]] = False
        s = StationCB(f,first_few_names[f])
        StationCallBacks.append(s)
        station_axes.append(fax)
        station_buttons.append(fbtn)
        fbtn.on_clicked(s.stationCB)

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('key_press_event', press)

    plt.show()
##~5
