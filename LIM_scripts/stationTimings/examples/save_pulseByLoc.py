#!/usr/bin/env python3

from LoLIM.stationTimings.autoCorrelator3_part1_fromLOC import save_EventByLoc

## these lines are anachronistic and should be fixed at some point
from LoLIM import utilities
utilities.default_raw_data_loc = "/exp_app2/appexp1/public/raw_data"
utilities.default_processed_data_loc = "/home/brian/processed_files"

if __name__=="__main__":
    
    station_delays = {
"CS002": 0.0,
            
"RS208":  6.97951240511e-06 , ##??
"RS210": 0.0, ##??
            
'RS305' :  7.1934989871e-06 , ## diff to guess: 1.37699440122e-09 ## station has issues
            
'CS001' :  2.24401172494e-06 , ## diff to guess: -1.97549250145e-12
'CS003' :  1.4070479021e-06 , ## diff to guess: -2.43902860006e-11
'CS004' :  4.40075444977e-07 , ## diff to guess: -2.42217738996e-11
'CS005' :  -2.17678756484e-07 , ## diff to guess: -2.21583592931e-11
'CS006' :  4.30647454293e-07 , ## diff to guess: -1.65542771957e-11
'CS007' :  3.99377465507e-07 , ## diff to guess: -2.02129216007e-11
'CS011' :  -5.8955385854e-07 , ## diff to guess: -1.23462674702e-12
'CS013' :  -1.81293364466e-06 , ## diff to guess: -2.14074126082e-11
'CS017' :  -8.4551582793e-06 , ## diff to guess: 4.31470974773e-11
'CS021' :  9.23663075135e-07 , ## diff to guess: 0.0
'CS024' :  2.33045304324e-06 , ## diff to guess: 8.82976926023e-12
'CS026' :  -9.26723147045e-06 , ## diff to guess: -9.26723147045e-06
'CS030' :  -2.73965173619e-06 , ## diff to guess: 2.41776414432e-11
'CS031' :  6.28836357272e-07 , ## diff to guess: 6.28836357272e-07
'CS032' :  -1.54223879178e-06 , ## diff to guess: -2.35580980099e-06
'CS101' :  -8.21535661203e-06 , ## diff to guess: 1.21635324342e-10
'CS103' :  -2.85878989548e-05 , ## diff to guess: 7.47105624697e-10
'CS201' :  -1.0510320435e-05 , ## diff to guess: 1.70364032269e-10
'CS301' :  -6.87236667019e-07 , ## diff to guess: -2.36495962731e-11
'CS302' :  -5.26832572532e-06 , ## diff to guess: 3.69871327386e-10
'CS501' :  -9.63031916796e-06 , ## diff to guess: -1.77630731779e-11
'RS106' :  6.78633242596e-06 , ## diff to guess: 1.575478568e-08
'RS205' :  7.09137015323e-06 , ## diff to guess: 9.84835165259e-10 
'RS306' :  7.41317947913e-06 , ## diff to guess: 1.6935786672e-08
'RS307' :  7.80951637749e-06 , ## diff to guess: 5.67451737617e-08
'RS310' :  7.96312229684e-06 , ## diff to guess: 4.82651956563e-07
'RS406' :  6.96541956339e-06 , ## diff to guess: 8.9825482633e-09
'RS407' :  6.83033534647e-06 , ## diff to guess: 3.5607628373e-09
'RS409' :  7.4346353712e-06 , ## diff to guess: 1.55561380515e-07
'RS503' :  6.9238247372e-06 , ## diff to guess: 2.51105167164e-10
'RS508' :  6.36331650134e-06 , ## diff to guess: 2.40127126525e-09
            }
    
    
    station_delays['CS002'] = 0.0 ## need to add referance station
    
    save_EventByLoc(timeID = "D20180813T153001.413Z", 
                    XYZT = [0, 0, 0, 0], 
                    station_timing_delays = station_delays, 
                    pulse_index = 0, 
                    output_folder = "correlate_foundPulses_byLOC", 
                    pulse_width=50,
                    polarization_flips="polarization_flips.txt", 
                    bad_antennas="bad_antennas.txt", 
                    additional_antenna_delays = "ant_delays.txt")
    
    
    