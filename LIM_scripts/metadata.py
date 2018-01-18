#!/usr/bin/env python3

"""This module reads in calibration metadata from file in the early fases of LOFAR. In the future this should be replaced by reading the metadata from the files.

.. moduleauthor:: Sander ter Veen <s.terveen@astro.ru.nl>

Modified by Brian Hare for use with LOFAR for Lightning Imaging
"""

## Imports
import numpy as np
import struct

from utilities import SId_to_Sname, latlonCS002, RTD, MetaData_directory

#### first some simple utility functions ###

def mapAntennasetKeyword(antennaset):
    """Ugly fix to map correct antenna names in input to wrong antenna names
    for metadata module.
    """

    # Strip whitespace
    antennaset = antennaset.strip()
    
    allowed = ["LBA_OUTER", "LBA_INNER", "LBA_X", "LBA_Y", "HBA", "HBA_0", "HBA_1"]
    incorrect = {'LBA_INNER': 'LBA_INNER',
                  'LBA_OUTER': 'LBA_OUTER',
                  'HBA_ZERO': 'HBA_0',
                  'HBA_ONE': 'HBA_1',
                  'HBA_DUAL': 'HBA',
                  'HBA_JOINED': 'HBA',
				  'HBA_ZERO_INNER': 'HBA_0', # Only true for core stations
				  'HBA_ONE_INNER': 'HBA_1',  # Only true for core stations
				  'HBA_DUAL_INNER': 'HBA',   # Only true for core stations
				  'HBA_JOINED_INNER': 'HBA'} # Only true for core stations

    if antennaset in incorrect:
        antennaset = incorrect[antennaset]
    elif antennaset == "HBA_BOTH":
        # This keyword is also wrong but present in file headers
        print( "Keyword " + antennaset + " does not comply with ICD, mapping...")
        antennaset = "HBA"
    
    assert antennaset in allowed

    return antennaset
    
def make_antennaID_filter(antennaIDs):
    """For a list of antennaIDs, return a filter to filter data by antenna.
    example use:
        getStationPhaseCalibration("CS001", "LBA_OUTER")[ make_antennaID_filter(["002000001"]) ]
        
    note: Only works for one station at a time. Assumes that the array you want to filter includes ALL antennas in the appropriate antenna set"""
        
    RCU_id = np.array([int(ID[-3:]) for ID in antennaIDs])
    return RCU_id
    





#### read callibration data ###

def getStationPhaseCalibration(station, antennaset, file_location=None):
    """Read phase calibration data for a station.

    Required arguments:

    ================== ====================================================
    Parameter          Description
    ================== ====================================================
    *station*          station name (as str) or ID.
    *mode*             observation mode.
    ================== ====================================================
    
    returns weights for 512 subbands.

    Examples::

        >>> metadata.getStationPhaseCalibration("CS002","LBA_OUTER")
        array([[ 1.14260161 -6.07397622e-18j,  1.14260161 -6.05283530e-18j,
             1.14260161 -6.03169438e-18j, ...,  1.14260161 +4.68675289e-18j,
             1.14260161 +4.70789381e-18j,  1.14260161 +4.72903474e-18j],
           [ 0.95669876 +2.41800591e-18j,  0.95669876 +2.41278190e-18j,
             0.95669876 +2.40755789e-18j, ...,  0.95669876 -2.41017232e-19j,
             0.95669876 -2.46241246e-19j,  0.95669876 -2.51465260e-19j],
           [ 0.98463207 +6.80081617e-03j,  0.98463138 +6.89975906e-03j,
             0.98463069 +6.99870187e-03j, ...,  0.98299670 +5.71319125e-02j,
             0.98299096 +5.72306908e-02j,  0.98298520 +5.73294686e-02j],
           ...,
           [ 1.03201290 +7.39535744e-02j,  1.03144532 +8.14880844e-02j,
             1.03082273 +8.90182487e-02j, ..., -0.82551740 -6.23731331e-01j,
            -0.82094046 -6.29743206e-01j, -0.81631975 -6.35721497e-01j],
           [ 1.12370332 -1.15296909e-01j,  1.12428451 -1.09484545e-01j,
             1.12483564 -1.03669252e-01j, ..., -0.92476286 +6.48703460e-01j,
            -0.92810503 +6.43912711e-01j, -0.93142239 +6.39104744e-01j],
           [ 1.10043006 -6.18995646e-02j,  1.10075250 -5.58731668e-02j,
             1.10104193 -4.98450938e-02j, ..., -1.01051042 +4.40052904e-01j,
            -1.01290481 +4.34513198e-01j, -1.01526883 +4.28960464e-01j]])

        >>> metadata.getStationPhaseCalibration(122,"LBA_OUTER")
        Calibration data not yet available. Returning 1
        array([[ 1.+0.j,  1.+0.j,  1.+0.j, ...,  1.+0.j,  1.+0.j,  1.+0.j],
           [ 1.+0.j,  1.+0.j,  1.+0.j, ...,  1.+0.j,  1.+0.j,  1.+0.j],
           [ 1.+0.j,  1.+0.j,  1.+0.j, ...,  1.+0.j,  1.+0.j,  1.+0.j],
           ...,
           [ 1.+0.j,  1.+0.j,  1.+0.j, ...,  1.+0.j,  1.+0.j,  1.+0.j],
           [ 1.+0.j,  1.+0.j,  1.+0.j, ...,  1.+0.j,  1.+0.j,  1.+0.j],
           [ 1.+0.j,  1.+0.j,  1.+0.j, ...,  1.+0.j,  1.+0.j,  1.+0.j]])

    """

    # Return mode nr depending on observation mode
    antennasetToMode = {"LBA_OUTER": "1",
                         "LBA_INNER": "3",
                         "HBA": "5",
                         "HBA_0": "5",
                         "HBA_1": "5",
                         }

    antennaset = mapAntennasetKeyword( antennaset )

    if antennaset not in antennasetToMode.keys():
        raise KeyError("Not a valid antennaset " + antennaset)

    modenr = antennasetToMode[antennaset]
    if not isinstance(station, str):
        # Convert a station id to a station name
        station = SId_to_Sname[station]

    stationNr = station[2:]

    # filename
    if file_location is None:
        file_location = MetaData_directory + '/lofar/StaticMetaData/CalTables'

    filename = file_location + '/CalTable_' + stationNr + '_mode' + modenr + '.dat'
    with open(filename, 'rb') as fin:
        # Test for header record above raw data - present in newer caltables (starting 2012)
        line = fin.readline().decode()
        if 'HeaderStart' in line:
            while not 'HeaderStop' in line:
                line = fin.readline().decode()
        else:  # no header present, seek to starting position
            fin.seek(0)
            
        data = np.fromfile(fin, dtype=np.double)
        
    data.resize(512, 96, 2)

    complexdata = np.empty(shape=(512, 96), dtype=complex)
    complexdata.real = data[:, :, 0]
    complexdata.imag = data[:, :, 1]

    return complexdata.transpose()

def convertPhase_to_Timing(phase_calibration, sample_time=5.0e-9):
    """Given the phase calibration of the 512 LOFAR subbands, such as the output of getStationPhaseCalibration, return the timing callibration of each antenna. 
    Not sure how well this works with HBA antennas. Sample time should be seconds per sample. Default is 5 ns"""
    
    phases = np.angle(phase_calibration)
    delays = (phases[:, 1] - phases[:, 0]) * (1024 / (2*np.pi)) * sample_time ## this just finds the slope from the first two points. Are there better methods?
    return delays
    


#def getStationGainCalibration(station, antennaset, file_location=None):
#    """Read phase calibration data for a station.
#
#    Required arguments:
#
#    ================== ====================================================
#    Parameter          Description
#    ================== ====================================================
#    *station*          station name or ID.
#    *mode*             observation mode.
#    ================== ====================================================
#
#    Optional arguments:
#
#    ================== ====================================================
#    Parameter          Description
#    ================== ====================================================
#    *return_as_hArray* Default False
#    ================== ====================================================
#
#    returns one gain per RCU. This gain is calculated using the absolute
#    value from the CalTables assuming these are not frequency dependent.
#    This seems to be true in current (2013-08) tables.
#    """
#
#    cal = getStationPhaseCalibration(station, antennaset, file_location)
#
#    gain = np.abs(cal[:,0])
#    return gain



#### information about cable lengths and delays ####

def getCableDelays(station, antennaset):
    """ Get cable delays in seconds.

    Required arguments:

    ================== ====================================================
    Parameter          Description
    ================== ====================================================
    *station*          Station name or ID e.g. "CS302", 142
    *antennaset*       Antennaset used for this station. Options:

                       * LBA_INNER
                       * LBA_OUTER
                       * LBA_X
                       * LBA_Y
                       * LBA_SPARSE0
                       * LBA_SPARSE1
                       * HBA_0
                       * HBA_1
                       * HBA

    ================== ====================================================

    returns "array of (rcus * cable delays ) for all dipoles in a station"

    """

    # Check station id type
    if not isinstance(station, str):
        # Convert a station id to a station name
        station = SId_to_Sname[station]
        
    antennaset = mapAntennasetKeyword( antennaset ) 
        
    if "LBA_OUTER" == antennaset:
        rcu_connection = "LBL"
    elif "LBA_INNER" == antennaset:
        rcu_connection = "LBH"
    elif antennaset in ['HBA', "HBA_1", "HBA_0"]:
        rcu_connection = "HBA"
    else:
        raise KeyError("Not a valid antennaset " + antennaset)

    cabfilename = MetaData_directory + '/lofar/StaticMetaData/CableDelays/' + station + '-CableDelays.conf'
    cabfile = open(cabfilename)

    cable_delays = np.zeros(96)

    str_line = ''
    while "RCUnr" not in str_line:
        str_line = cabfile.readline()
        if len(str_line) == 0:
            # end of file reached, no data available
            assert False

    str_line = cabfile.readline()
    for i in range(96):
        str_line = cabfile.readline()
        sep_line = str_line.split()
        if rcu_connection == "LBL":
            cable_delays[int(sep_line[0])] = float(sep_line[2]) * 1e-9
        elif rcu_connection == "LBH":
            cable_delays[int(sep_line[0])] = float(sep_line[4]) * 1e-9
        elif rcu_connection == "HBA":
            cable_delays[int(sep_line[0])] = float(sep_line[6]) * 1e-9

    return cable_delays

def getCableLength(station,antennaset):
    
    # Check station id type
    if not isinstance(station, str):
        # Convert a station id to a station name
        station = SId_to_Sname[station]
        
    antennaset = mapAntennasetKeyword( antennaset ) 
        
    if "LBA_OUTER" == antennaset:
        rcu_connection = "LBL"
    elif "LBA_INNER" == antennaset:
        rcu_connection = "LBH"
    elif antennaset in ['HBA', "HBA_1", "HBA_0"]:
        rcu_connection = "HBA"
    else:
        raise KeyError("Not a valid antennaset " + antennaset)

    cabfilename = MetaData_directory + '/lofar/StaticMetaData/CableDelays/' + station + '-CableDelays.conf'
    cabfile = open(cabfilename)

    cable_length = np.zeros(96)

    str_line = ''
    while "RCUnr" not in str_line:
        str_line = cabfile.readline()
        if len(str_line) == 0:
            # end of file reached, no data available
            assert False

    str_line = cabfile.readline()
    for i in range(96):
        str_line = cabfile.readline()
        sep_line = str_line.split()
        if rcu_connection == "LBL":
            cable_length[int(sep_line[0])] = float(sep_line[1])
        elif rcu_connection == "LBH":
            cable_length[int(sep_line[0])] = float(sep_line[3])
        elif rcu_connection == "HBA":
            cable_length[int(sep_line[0])] = float(sep_line[5])

    return cable_length

def antennaset2rcumode(antennaset,filter):
    antennaset = mapAntennasetKeyword( antennaset ) 
    rcumode=dict()
    rcumode[('LBA_INNER','LBA_10_90')]=1
    rcumode[('LBA_OUTER','LBA_10_90')]=2
    rcumode[('LBA_INNER','LBA_30_90')]=3
    rcumode[('LBA_OUTER','LBA_30_90')]=4
    rcumode[('HBA','HBA_110_190')]=5
    rcumode[('HBA','HBA_170_230')]=6
    rcumode[('HBA','HBA_210_250')]=7
    if 'HBA' in antennaset:
        antennaset='HBA'
    return rcumode[(antennaset,filter)]

def getCableAttenuation(station,antennaset,filter=None):

    cable_length = getCableLength(station, antennaset)

    attenuationFactor=dict()
    attenuationFactor[1]=-0.0414#{50:-2.05,80:-3.32,85:-3.53,115:-4.74,130:-5.40}
    attenuationFactor[2]=attenuationFactor[1]
    attenuationFactor[3]=attenuationFactor[1]
    attenuationFactor[4]=attenuationFactor[1]
    attenuationFactor[5]=-0.0734#{50:-3.64,80:-5.87,85:-6.22,115:-8.53,130:-9.52}
    attenuationFactor[6]=-0.0848#{50:-4.24,80:-6.82,85:-7.21,115:-9.70,130:-11.06}
    attenuationFactor[7]=-0.0892#{50:-4.46,80:-7.19,85:-7.58,115:-10.18,130:-11.61}
    if filter==None:
        if 'LBA' in antennaset:
            filter='LBA_30_90'
        else:
            print( "Please specify the filter!")
            filter='HBA_110_190'
    rcumode=antennaset2rcumode(antennaset,filter)
    att=attenuationFactor[rcumode]
    return cable_length*att











#### functions for antenna and station location #####

def getItrfAntennaPosition(station, antennaset):
    """Returns the antenna positions of all the antennas in the station
    in ITRF coordinates for the specified antennaset.
    station can be the name or id of the station. 
    
    Required arguments:

    =================== ==============================================
    Parameter           Description
    =================== ==============================================
    *station*           Name or id of the station. e.g. "CS302" or 142
    *antennaset*        Antennaset used for this station. Options:

                        * LBA_INNER
                        * LBA_OUTER
                        * LBA_X
                        * LBA_Y
                        * LBA_SPARSE0
                        * LBA_SPARSE1
                        * HBA_0
                        * HBA_1
                        * HBA

    =================== ==============================================

    """
    # Check station id type
    if isinstance(station, int):
        # Convert a station id to a station name
        station = SId_to_Sname[station]
        
    antennaset = mapAntennasetKeyword( antennaset ) 
    
    if "LBA" in antennaset:
        antennatype = "LBA"
    elif "HBA" in antennaset:
        antennatype = "HBA"
    

    # Obtain filename of antenna positions
    filename = MetaData_directory + "/lofar/StaticMetaData/AntennaFields/" + station + "-AntennaField.conf"

    # Open file
    f = open(filename, 'r')

    if station[0:2] != "CS":
        if "HBA" in antennaset:
            antennaset = "HBA"

    # Find position of antennaset in file
    str_line = ''
    while antennatype != str_line.strip():
        str_line = f.readline()
        if len(str_line) == 0:
            # end of file reached, no data available
            assert False

    # Find the location of the station. Antenna locations are relative to this
    str_line = f.readline()
    str_split = str_line.split()
    stationX = float(str_split[2])
    stationY = float(str_split[3])
    stationZ = float(str_split[4])

    str_line = f.readline()

    # Get number of antennas and the number of directions
    nrantennas = int(str_line.split()[0])
    nrdir = int(str_line.split()[4])

    antenna_positions = np.empty((2*nrantennas, nrdir), dtype=np.double)
    for i in range(nrantennas):
        line = f.readline().split()
        
        antenna_positions[2*i, 0] = float(line[0]) + stationX
        antenna_positions[2*i, 1] = float(line[1]) + stationY
        antenna_positions[2*i, 2] = float(line[2]) + stationZ
        
        antenna_positions[2*i+1, 0] = float(line[3]) + stationX
        antenna_positions[2*i+1, 1] = float(line[4]) + stationY
        antenna_positions[2*i+1, 2] = float(line[5]) + stationZ

    if antennatype == "LBA":
        # There are three types of feed
        # H for HBA
        # h for lbh
        # l for lbl
        feed = {}
        feed["CS"] = {}
        feed["RS"] = {}
        feed["DE"] = {}
        feed["CS"]["LBA_SPARSE_EVEN"] = "24llhh"
        feed["CS"]["LBA_SPARSE_ODD"] = "24hhll"
        feed["CS"]["LBA_X"] = "48hl"
        feed["CS"]["LBA_Y"] = "48lh"
        feed["CS"]["LBA_INNER"] = "96h"
        feed["CS"]["LBA_OUTER"] = "96l"
        feed["RS"]["LBA_SPARSE_EVEN"] = "24llhh"
        feed["RS"]["LBA_SPARSE_ODD"] = "24hhll"
        feed["RS"]["LBA_X"] = "48hl"
        feed["RS"]["LBA_Y"] = "48lh"
        feed["RS"]["LBA_INNER"] = "96h"
        feed["RS"]["LBA_OUTER"] = "96l"
        feed["DE"]["LBA"] = "192h"
        if station[0:2] == "CS" or "RS":
            feedsel = feed[station[0:2]][antennaset]
            nrset = int(feedsel.split('l')[0].split('h')[0].split('H')[0])
            feeds = ''
            feedsel = feedsel[len(str(nrset)):]
            for i in range(nrset):
                feeds += feedsel

        indexselection = []
        for i in range(len(feeds)):
            if feeds[i] == 'l':
                # The 'l' feeds are the last 96 numbers of the total list
                indexselection.append(i + 96)
            elif feeds[i] == 'h':
                # The 'h' feeds are the first 96 numbers of the total list
                indexselection.append(i)
            else:
                # This selection is not yet supported
                assert False
        antenna_positions = antenna_positions[indexselection]

    return antenna_positions

def getStationPositions(station, antennaset, coordinatesystem):
    """Returns the antenna positions of all the antennas in the station
    relative to the station center for the specified antennaset.
    station can be the name or id of the station. 

    Required arguments:

    ================== ==============================================
    Argument           Description
    ================== ==============================================
    *station*          Name or id of the station. e.g. "CS302" or 142
    *antennaset*       Antennaset used for this station. Options:

                       * LBA_INNER
                       * LBA_OUTER
                       * LBA_X
                       * LBA_Y
                       * LBA_SPARSE0
                       * LBA_SPARSE1
                       * HBA_0
                       * HBA_1
                       * HBA

    *coordinatesystem* WGS84 or ITRF

    output: 
        if coordinatesystem == "WGS84", then return [lat, lon, alt] as a numpy array
        else if coordinatesystem=="ITRF", then return [X, Y, Z] as a numpy array

    """

    # Check if requested antennaset is known
    assert coordinatesystem in ["WGS84", 'ITRF']
    
    # Check station id type
    if isinstance(station, int):
        # Convert a station id to a station name
        station = SId_to_Sname[station]
        
    antennaset = mapAntennasetKeyword( antennaset ) 

    # Obtain filename of antenna positions
    if 'WGS84' in coordinatesystem:
        filename = MetaData_directory + "/lofar/StaticMetaData/AntennaArrays/" + station + "-AntennaArrays.conf"
    else:
        filename = MetaData_directory + "/lofar/StaticMetaData/AntennaFields/" + station + "-AntennaField.conf"

    # Open file
    f = open(filename, 'r')
    
    if "LBA" in antennaset:
        antennaset = "LBA"

    if station[0:2] != "CS":
        if "HBA" in antennaset:
            antennaset = "HBA"

    # Find position of antennaset in file
    str_line = ''
    while antennaset != str_line.strip():
        str_line = f.readline()
        if len(str_line) == 0:
            # end of file reached, no data available
            print( "Antenna set not found in calibration file", filename)
            return None

    # Skip name and station reference position
    str_line = f.readline().split()
    
    A = float(str_line[2]) ## lon in WGS84, X in ITRF
    B = float(str_line[3]) ## lat in WGS84, Y in ITRF
    C = float(str_line[4]) ## alt in WGS84, Z in ITRF

    return np.array([A,B,C])

ITRFCS002 = getStationPositions('CS002', 'LBA_OUTER', coordinatesystem='ITRF')  # ($LOFARSOFT/data/lofar/StaticMetaData/AntennaFields/CS002-AntennaField.conf)

def convertITRFToLocal(itrfpos, phase_center=ITRFCS002, reflatlon=latlonCS002, out=None):
    """
    ================== ==============================================
    Argument           Description
    ================== ==============================================
    *itrfpos*          an ITRF position as 1D numpy array, or list of positions as a 2D array
    *phase_center*     the origin of the coordinate system, in ITRF. Default is CS002.
    *reflatlon*        the rotation of the coordinate system. Is the [lat, lon] (in degrees) on the Earth which defines "UP"
    
    Function returns a 2D numpy array (even if input is 1D).
    Out cannot be same array as itrfpos
    """
    lat = reflatlon[0]/RTD
    lon = reflatlon[1]/RTD
    arg0 = np.array([-np.sin(lon),   -np.sin(lat) * np.cos(lon),   np.cos(lat) * np.cos(lon)])
    arg1 = np.array([np.cos(lon) ,   -np.sin(lat) * np.sin(lon),   np.cos(lat) * np.sin(lon)])
    arg2 = np.array([         0.0,    np.cos(lat),                 np.sin(lat)])
    
    if out is None:
        ret = np.empty(itrfpos.shape, dtype=np.double )
    else:
        ret = out
    
    ret  = np.outer(itrfpos[...,0]-phase_center[0], arg0 )
    ret += np.outer(itrfpos[...,1]-phase_center[1], arg1 )
    ret += np.outer(itrfpos[...,2]-phase_center[2], arg2 )
    
    return ret
    
    