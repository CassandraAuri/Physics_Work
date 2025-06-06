from viresclient import set_token
from viresclient import SwarmRequest
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
from tqdm.notebook import tqdm
import mplcyberpunk
from itertools import chain
from MFA import MFA
from scipy.fft import fft, fftfreq
import asilib
import asilib.asi
from scipy import signal
from scipy.signal import butter, filtfilt, freqz
from numpy.typing import NDArray
import matplotlib.animation as animation
import matplotlib.colors as colors
import cmasher as cmr
import geopack.geopack as gp
from scipy.optimize import curve_fit, fsolve
import streamlit as st
plt.style.use("cyberpunk")
plt.rcParams['figure.facecolor'] =  '121212'
plt.rcParams['axes.facecolor'] =  '121212'
plt.rcParams['savefig.facecolor'] =  '121212'
mplcyberpunk.make_lines_glow()
mplcyberpunk.add_underglow()

from scipy.spatial.transform import Rotation as R
import numpy as np

def find_closest_indices(times1, times2):
    # Convert to numpy arrays
    times1 = np.array(times1)
    times2 = np.array(times2)
    
    # Compute the differences between each time in times1 and all times in times2
    # Resulting in a 2D array where each row contains the absolute differences for one time in times1
    differences = np.abs(times1[:, None] - times2)
    
    # Find the index of the minimum difference for each time in times1
    closest_indices = np.argmin(differences, axis=1)
    
    return closest_indices

def quaternion_inverse_scipy(q):
    # Ensure q is a numpy array
    q = np.asarray(q)
    
    # Create a Rotation object from the quaternion
    rotation = R.from_quat(q)  # Note: scipy uses [x, y, z, w] format
    
    # Compute the inverse rotation
    inverse_rotation = rotation.inv()
    
    
    return inverse_rotation

def long_sep(long1,long2,lag_data,sps):
    """
    Given a lag, find the average longitinduinal seperation
    """
    num_shift = int(lag_data[2] * sps)
    if num_shift > 0:
        adjusted_long1 = long1[num_shift:]
        adjusted_long2 = long2[:-num_shift]
    else:
        adjusted_long1 = long1[:num_shift]
        adjusted_long2 = long2[-num_shift:]
    
    # Ensure both arrays have the same length
    min_length = min(len(adjusted_long1), len(adjusted_long2))
    adjusted_long1 = adjusted_long1[:min_length]
    adjusted_long2 = adjusted_long2[:min_length]
    
    # Calculate the pairwise longitudinal separation
    separations = np.abs(np.array(adjusted_long1) - np.array(adjusted_long2))
    
    # Calculate and return the average separation
    average_separation = np.mean(separations)
    
    return average_separation

def process_string(input_string):
    # Convert the string to lowercase
    lowercased_string = input_string.lower()
    # Remove any spaces
    concatenated_string = lowercased_string.replace(" ", "")
    return concatenated_string

def butter_bandpass(cutoffs, fs, order=4):
    if cutoffs[0] ==0:
        return butter(order, cutoffs[1], fs=fs, btype="low")
    else:
        return butter(order, cutoffs, fs=fs, btype="band")


def butter_bandpass_filter(data, cutoffs, fs, order=4):
    b, a = butter_bandpass(cutoffs, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def create_sampled_datetimes(datetime_tuple, sampling_rate_seconds):
    """
    The create_sampled_datetimes function generates an array of equally spaced datetimes within a specified range. 

    Parameters
        datetime_tuple (tuple): A tuple containing two pandas.Timestamp objects.
        The first element is the start datetime.
        The second element is the end datetime.
        sampling_rate_seconds (int): The interval in seconds between each sampled datetime in the output array.

    Returns
    numpy.ndarray: An array of pandas.Timestamp objects representing the sampled datetimes, starting from the start datetime to the end datetime,
        spaced by the specified sampling rate.
    """
    start_datetime, end_datetime = datetime_tuple
    
    # Calculate the number of samples
    delta = end_datetime - start_datetime
    num_samples = int(delta.total_seconds() // sampling_rate_seconds) + 1
    
    # Create an array of equally spaced datetimes
    sampled_datetimes = [start_datetime + timedelta(seconds=sampling_rate_seconds * i) for i in range(num_samples)]
    
    return np.array(sampled_datetimes)
def offset_series(ts1, ts2, lag):
    if lag > 0:
        offset_ts1 = np.concatenate((np.full(lag, np.nan), ts1))
        offset_ts2 = np.concatenate((ts2, np.full(lag, np.nan)))
    else:
        lag = -lag
        offset_ts1 = np.concatenate((ts1, np.full(lag, np.nan)))
        offset_ts2 = np.concatenate((np.full(lag, np.nan), ts2))
    return offset_ts1, offset_ts2

# Step 4: Pad the Ends
def pad_ends(ts1, ts2):
    max_len = max(len(ts1), len(ts2))
    ts1_padded = np.pad(ts1, (0, max_len - len(ts1)), constant_values=np.nan)
    ts2_padded = np.pad(ts2, (0, max_len - len(ts2)), constant_values=np.nan)
    return ts1_padded, ts2_padded

def sinc_interpolation(x: NDArray, s: NDArray, u: NDArray) -> NDArray:
    """Whittaker–Shannon or sinc or bandlimited geodeticolation.
    Args:
        x (NDArray): signal to be geodeticolated, can be 1D or 2D
        s (NDArray): time points of x (*s* for *samples*)
        u (NDArray): time points of y (*u* for *upsampled*)
    Returns:
        NDArray: geodeticolated signal at time points *u*
    Reference:
        This code is based on https://gist.github.com/endolith/1297227
        and the comments therein.
    TODO:
        * implement FFT based geodeticolation for speed up
    """
    sinc_ = np.sinc((u - s[:, None]) / (s[1] - s[0]))

    return np.dot(x, sinc_)

time_stamper = np.vectorize(lambda x: x.timestamp())

def arrangement(time, array, shape):  # arranges B into a useable format for use later
    barranged = np.zeros((len(time), shape))
    # Re-arranges into proper (n x 3 ) matricies, ugly but works
    for j in range(len(time)):
        for k in range(shape):
            barranged[j][k] = array[j][k]
    return barranged
def find_indices(datetime_array, start_time, end_time):
    """
    Find indices in datetime_array that correspond to dt1 and dt2.

    Parameters:
    datetime_array (np.array): Array of datetime objects.
    dt1 (datetime.datetime): First datetime to find.
    dt2 (datetime.datetime): Second datetime to find.
    sps (int): Samples per second, default is 16.

    Returns:
    int: Index corresponding to dt1.
    int: Index corresponding to dt2.
    """
    # Ensure the array is sorted
    start_index = np.argmin(np.abs(datetime_array - start_time))

    # Find the index of the closest time to the end_time
    end_index = np.argmin(np.abs(datetime_array - end_time))

    return start_index,end_index

def footprint(time, latitude, longitude, altitude, alt,vsw= [-400,0,-40]):
    """
    time, datetime, time for magnetic data for footprint
    vsw velocity of solar wind, tuple of x,y,z
    longitude of satellite in degrees
    latitude of satellite in degrees
    altitude of satellite in km from centre of earth (should be above ~6371)
    THIS CODE ONLY Works in the NOrthern hemisphere

    """
    def cubic(t, a, b, c, d):
            return a*t**3 + b*t**2 + c*t + d
    def x(t, params_x):
        return cubic(t, *params_x)

    def y(t, params_y):
        return cubic(t, *params_y)

    def z(t, params_z):
        return cubic(t, *params_z)
    def radius(t, params_x, params_y, params_z):
        return np.sqrt(x(t, params_x)**2 + y(t, params_y)**2 + z(t, params_z)**2)
    
    def curve_fit_func(xx,yy,zz, differencealt):
        
        r = np.linspace(1, 1.5, 10000)# construct an array of radiuses from 1-1.5

        radius_data=np.sqrt(xx**2+yy**2+zz**2)

        params_x, _ = curve_fit(cubic, radius_data, xx) #Constructs fits on the traces inward since the spatial resolution produced by geopack is limited.
        params_y, _ = curve_fit(cubic, radius_data, yy)
        params_z, _ = curve_fit(cubic, radius_data, zz)

        

        index_closest=np.argmin(np.abs(radius(r, params_x, params_y, params_z)-(alt-differencealt+6371)/6371))#Find the index that produces the closest radius to the altitude

        return x(r[index_closest],params_x ),y(r[index_closest],params_y ),z(r[index_closest],params_z )
    
    t1 = time
    t0 = datetime(1970,1,1) #epoch
    py_dt = pd.to_datetime(t1).to_pydatetime()
    ut = (py_dt-t0).total_seconds()
    lat_sat=np.deg2rad(latitude)
    lon_sat=np.deg2rad(longitude) #converts to radii
    gp.recalc(ut, vsw[0], vsw[1], vsw[2])
    r, theta= gp.geodgeo(altitude,lat_sat,1) #this r accounts for earths oblateness, so we need to find the difference between my 6371 assumption and the real value and account for that
    differencearray= (altitude+6371)-r
    x_gc,y_gc,z_gc = gp.sphcar((r)/6371,theta,lon_sat,1)  #spherical to cartesian
    

    x_gsm, y_gsm, z_gsm = gp.geogsm(x_gc,y_gc,z_gc, 1) #cartesian to gsm

    x_foot,y_foot,z_foot=np.zeros(len(x_gsm)), np.zeros(len(y_gsm)), np.zeros(len(z_gsm)) #initalize an array
    for index in range(len(x_gsm)):
        x_foot_int, y_foot_int, z_foot_int, xx2, yy2,zz2 = gp.trace(x_gsm[index], y_gsm[index], z_gsm[index], dir=1,rlim=3, maxloop=1000 ) #traces each set of lat,lon,alt outward


        x_foot[index],y_foot[index],z_foot[index] = curve_fit_func(xx2,yy2,zz2, differencearray[index])


            

    x_done, y_done, z_done = gp.geogsm(x_foot, y_foot, z_foot, -1)

    alt_sat_done, lat_sat_done,lon_sat_done = np.zeros(len(x_done)), np.zeros(len(x_done)), np.zeros(len(x_done))
    for index in range(len(x_done)):
        
        r_done,theta_done,lon_sat_done[index]= gp.sphcar(x_done[index], y_done[index], z_done[index],-1)

        alt_sat_done[index], lat_sat_done[index]= gp.geodgeo(r_done*6371,theta_done,-1) #TODO check if this is right

    print(alt_sat_done, 'altitude derived from fit')

    if np.any(np.abs(alt_sat_done - alt) > 5):
        raise Exception("One or more values in the footprinting are greater than 5km away from the specified alt. Contact owner for a fix, not your fault")
    print(np.rad2deg(lon_sat_done)-360,np.rad2deg(lat_sat_done) , 'lat and lon' )
    sat_lla=np.array([ np.rad2deg(lat_sat_done), np.rad2deg(lon_sat_done)-360,  alt_sat_done])
    return sat_lla

def unit_array(array):
    arraysum = np.sum(np.abs(array), axis=1)
    # Normalizes and finds unitary
    array_unit = array / arraysum[:, np.newaxis]  # normalizes
    return array_unit

def Graphing_Ratio(space_craft_with_E, efield, bfield, time_E, time_B, user_select, fig, axes):
    """
    Gives the ratio of E/B in the spectral domain, either returning an animation or a plot of the conductivies it derived
    """
    def Logic_for_one_step(index, index_satellite,Bresamp, nperseg, overwrite_indicies=False):
        """
        Finds periodogram, Conductivity and Alfven speed for each window
        """
        if overwrite_indicies == True:
            index_start, index_end= index[0], index[1]
        else:
            window_length=user_select["window_length"]
            index_start,index_end= int(index*16*user_select['sampling_rate']), int(index*16*user_select['sampling_rate'])+16*window_length
        
        
        #TODO do this per space-craft and do it in the ratio desired, currently North over Eas
        #if index==0:
            #print(phase_ENorth_BEast, 'calculated')
        ##print(np.absolute(crossspectral_ENorth_BEast[0]),np.sqrt(np.multiply(powerspec_E_0, powerspec_B_1))[0] )
        #coherence_ENorth_BEast = signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend=None, return_onesided=False,  nperseg=nperseg)
        #coherence_EEast_BNorth = signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend=None,return_onesided=False , nperseg=nperseg)
       
        #coherence_ENorth_BEast, coherence_EEast_BNorth = coherence_ENorth_BEast[sorted_frequencies_indicies], coherence_EEast_BNorth[sorted_frequencies_indicies]

        ##TODO implement a B-lag cross correlation
        if user_select['lag'] == True and index_satellite == 1: 
            #TODO add a condition that if index our of range just put np.Nans
            if np.ceil(lag_data[2])*16<= index_start and np.floor(lag_data[2]*16).astype(int)+index_end< len(Bresamp[index_satellite]):
                idex_lag=-np.round(lag_data[2]*16).astype(int)
                index_start_test,index_end_test =index_start + idex_lag, index_end+idex_lag
                indicies= range(index_start_test,index_end_test)
                frequencies_E_0, powerspec_E_0 = signal.welch(efield[index_satellite][range(index_start_test,index_end_test), 0], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
                _, powerspec_B_0 = signal.welch(Bresamp[index_satellite][range(index_start_test,index_end_test), 0], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False,  nperseg=nperseg)

                _, powerspec_E_1 = signal.welch(efield[index_satellite][range(index_start_test,index_end_test), 1], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
                _, powerspec_B_1 = signal.welch(Bresamp[index_satellite][range(index_start_test,index_end_test),1], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
                sorted_frequencies_indicies= np.argsort(frequencies_E_0)
                frequencies_E_0= frequencies_E_0[sorted_frequencies_indicies]

                powerspec_E_0, powerspec_B_0, powerspec_E_1, powerspec_B_1 = powerspec_E_0[sorted_frequencies_indicies],  powerspec_B_0[sorted_frequencies_indicies], powerspec_E_1[sorted_frequencies_indicies], powerspec_B_1[sorted_frequencies_indicies]
                ratio_EB_01 = np.sqrt((powerspec_E_0/ powerspec_B_1))
                ratio_EB_10 = np.sqrt((powerspec_E_1/ powerspec_B_0))
                _ , crossspectral_ENorth_BEast = signal.csd(efield[index_satellite][range(index_start_test,index_end_test), 0], Bresamp[index_satellite][range(index_start_test,index_end_test),1], 16, window="hann", detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
                _ ,   crossspectral_EEast_BNorth = signal.csd(efield[index_satellite][range(index_start_test,index_end_test), 1], Bresamp[index_satellite][range(index_start_test,index_end_test),0], 16, window="hann", detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)


                crossspectral_ENorth_BEast, crossspectral_EEast_BNorth = crossspectral_ENorth_BEast[sorted_frequencies_indicies], crossspectral_EEast_BNorth[sorted_frequencies_indicies]


                phase_ENorth_BEast= np.angle(crossspectral_ENorth_BEast, deg=True)
                phase_EEast_BNorth = np.angle(crossspectral_EEast_BNorth, deg=True)
                _ ,cross_BB=signal.csd(Bresamp[index_satellite-1][range(index_start, index_end),1], Bresamp[index_satellite][range(index_start_test,index_end_test),1], 16,window="hann", return_onesided=False, detrend=None, nperseg=nperseg )

                phase_BB= np.angle(cross_BB, deg=True)
                cross_BB, phase_BB = cross_BB[sorted_frequencies_indicies], phase_BB[sorted_frequencies_indicies]

                _ ,cross_EE=signal.csd(efield[index_satellite-1][range(index_start, index_end),1], efield[index_satellite][range(index_start_test,index_end_test),1], 16,window="hann", return_onesided=False, detrend=None, nperseg=nperseg )


                phase_EE= np.angle(cross_EE, deg=True)
                cross_EE, phase_BB = cross_EE[sorted_frequencies_indicies], phase_EE[sorted_frequencies_indicies]

               
                return np.array([[frequencies_E_0, [None] * nperseg], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[np.mean([index_start_test/16,index_end_test/16], dtype=int)]*nperseg, [None]*nperseg], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [cross_BB, phase_BB], [cross_EE, phase_EE]]), indicies #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time
            else: 
                indicies=None
                frequencies_E_0, _ =  signal.welch(efield[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
                sorted_frequencies_indicies= np.argsort(frequencies_E_0)
                frequencies_E_0= frequencies_E_0[sorted_frequencies_indicies]

                powerspec_E_0 = [np.nan]*nperseg
                powerspec_B_0 = [np.nan]*nperseg
                powerspec_E_1 = [np.nan]*nperseg
                powerspec_B_1 = [np.nan]*nperseg
                ratio_EB_01 = [np.nan]*nperseg
                ratio_EB_10 = [np.nan]*nperseg
                crossspectral_EEast_BNorth = [np.nan]*nperseg
                crossspectral_ENorth_BEast = [np.nan]*nperseg
                phase_ENorth_BEast = [np.nan]*nperseg
                phase_EEast_BNorth = [np.nan]*nperseg
                cross_BB = [np.nan]*nperseg
                phase_BB= [np.nan]*nperseg
                cross_EE = [np.nan]*nperseg
                phase_EE= [np.nan]*nperseg
                index_start_test,index_end_test = np.nan, np.nan
                
                return np.array([[frequencies_E_0, [None] * nperseg], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[None]*nperseg, [None]*nperseg], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [cross_BB, phase_BB], [cross_EE, phase_EE]]), indicies #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time


        else:
            frequencies_E_0, powerspec_E_0 = signal.welch(efield[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
            _, powerspec_B_0 = signal.welch(Bresamp[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False,  nperseg=nperseg)

            _, powerspec_E_1 = signal.welch(efield[index_satellite][range(index_start, index_end), 1], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
            _, powerspec_B_1 = signal.welch(Bresamp[index_satellite][range(index_start, index_end),1], 16, window="hann",detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
            sorted_frequencies_indicies= np.argsort(frequencies_E_0)
            frequencies_E_0= frequencies_E_0[sorted_frequencies_indicies]

            powerspec_E_0, powerspec_B_0, powerspec_E_1, powerspec_B_1 = powerspec_E_0[sorted_frequencies_indicies],  powerspec_B_0[sorted_frequencies_indicies], powerspec_E_1[sorted_frequencies_indicies], powerspec_B_1[sorted_frequencies_indicies]
            ratio_EB_01 = np.sqrt((powerspec_E_0/ powerspec_B_1))
            ratio_EB_10 = np.sqrt((powerspec_E_1/ powerspec_B_0))
            _ , crossspectral_ENorth_BEast = signal.csd(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16, window="hann", detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)
            _ ,   crossspectral_EEast_BNorth = signal.csd(efield[index_satellite][range(index_start, index_end), 1], Bresamp[index_satellite][range(index_start, index_end),0], 16, window="hann", detrend=None, scaling='spectrum', return_onesided=False, nperseg=nperseg)


            crossspectral_ENorth_BEast, crossspectral_EEast_BNorth = crossspectral_ENorth_BEast[sorted_frequencies_indicies], crossspectral_EEast_BNorth[sorted_frequencies_indicies]


            phase_ENorth_BEast= np.angle(crossspectral_ENorth_BEast, deg=True)
            phase_EEast_BNorth = np.angle(crossspectral_EEast_BNorth, deg=True)
            cross_BB= [None] * nperseg
            phase_BB= [None] * nperseg
            cross_EE =[None] * nperseg
            phase_EE = [None] * nperseg
            indicies=  range(index_start, index_end)
            
            
            return np.array([[frequencies_E_0, [None] * nperseg], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[np.mean([index_start/16,index_end/16], dtype=int)]*nperseg, [None]*nperseg], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [cross_BB, phase_BB], [cross_EE, phase_EE]]), indicies #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time

        

    
    def conductivities(data,times):

        """
        Gives the plot of conductivities and Alfven speed found by looping through window
        """
        def subplot_select():
            #TODO clean up with for loop and apply to other
            length_for_axis = 0
            try:
                length_for_axis += len(user_select["graph_E_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(user_select["graph_B_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(user_select["graph_PF_chosen"])
            except TypeError:
                pass
            if user_select["FAC"] == True:
                length_for_axis += 1
            else:
                pass
            if user_select["Difference"] == True:
                length_for_axis += 1
            else:
                pass
            if user_select["Pixel_intensity"] == True:
                try:
                    length_for_axis += len(user_select["sky_map_values"])
                except TypeError:
                    raise ValueError("Sky Map not Selected")
            if user_select["heatmap"] != None:
                length_for_axis += len(user_select["heatmap"])*len(user_select["satellite_graph"])
            return length_for_axis
        

        length_of_axis=subplot_select()
        for index in range(len(user_select["conductivities"])):
            for k in range(len(user_select["satellite_graph"])):
                if user_select['coordinate_system'] == "North East Centre": 
                    #From North East or East North, are they self similar who knows?
                    if user_select["conductivities"][index] == "ENorth/BEast":
                        conductivies= data[k, :,  3, 0, 0] # we want zeroth frequency term of the EB ratio
                    else: #East over North
                        conductivies= data[k, :,  3, 1, 0]

                else:
                    if index == "Azimuth over Polodial":
                        conductivies= data[k, :,  3, 0, 0] # we want zeroth frequency term of the EB ratio
                    else: #East over North
                        conductivies= data[k, :,  3, 1, 0]
                axes[index + length_of_axis].plot(times[np.array(data[k, :, 4,0, 0], dtype=int)],1/(1.256e-6*conductivies), label= "".join([user_select["satellite_graph"][k], " " ,user_select["conductivities"][index]]))
                axes [index + length_of_axis].set_ylabel("Conductivties (S)")
                axes[ index + length_of_axis].set_title("Conductivity versus time")
                axes[ index + length_of_axis].legend()

                

        


        return 
    
    def Animation_rows():
        #TODO Creates an animation of times series of deviations E and B, plot periodograms, plot E/B ratio for given, and what polarization or coordinate should be graphed for each plot
        axes_length=0
        user_select_selected = [ user_select["E_periodogram"], user_select["B_periodogram"]]
        if user_select["Time_Series"] != None:
            for i in range(len(user_select["satellite_graph"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        else:
            pass
        if user_select["EB_periodogram"] != None:
            for i in range(len(user_select["EB_periodogram"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["EB_cross power"] != None:
            for i in range(len(user_select["EB_cross power"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["EB_cross phase"] != None:
            for i in range(len(user_select["EB_cross power"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["lags_cross_B"] != None:
            for i in range(len(user_select["lags_cross_B"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["lags_cross_E"] != None:
            for i in range(len(user_select["lags_cross_E"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        else:
            pass
        for key in user_select_selected:
            if key != None:
                axes_length += 1
        return axes_length

    def Animation(data,sampled_datetimes, B_sinc, B_resample, efield, indicies, frames, time_E):
        """
        Creates an animation of each window
        """
        fig_ani, axes_ani = plt.subplots(figsize=(10,15),  nrows=Animation_rows())
        axs = np.array(axes_ani)
        if user_select["Time_Series"] != None:
            twin_x_axes=[axes_ani[k].twinx() for k in range(len(user_select["satellite_graph"]))]
        def animate(i):
            """
            Wrapper for animation
            """

            for ax in axs.reshape(-1):
                ax.clear()
            def EB_Periodogram_Plot(axes_used):
                
                try:
                    for l, val in enumerate(user_select["EB_periodogram"]):
                        for k, val_sat in enumerate((user_select["satellite_graph"])):
                            if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                                lag_seconds = int(lag_data[2])
                                lag_fraction = lag_data[2] - lag_seconds
                                lag_nanoseconds = int(lag_fraction * 1e9)
                                # Create timedelta64 objects
                                seconds_delta = np.timedelta64(lag_seconds, 's')
                                nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                                # Add the timedelta to the datetime array
                                delta =seconds_delta + nanoseconds_delta
                            else:
                                delta=np.timedelta64(0,'s')
                            if val == 'ENorth/BEast':
                                index=0
                            else:
                                index=1
                            axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 3,index, :], "-o", label=val) #Frequencies E's
                            axes_ani[axes_used +l ].set_ylabel(f" E over B ratio of {val_sat} m/s ")
                            axes_ani[axes_used + l ].legend()
                            axes_ani[axes_used + l ].set_yscale("log")
                            axes_ani[axes_used + l ].set_xlim(user_select["bandpass"][1])
                except TypeError:
                    for l, val in enumerate(user_select["EB_periodogram"]):
                        for k, val_sat in enumerate((user_select["satellite_graph"])):
                            if val == 'ENorth/BEast':
                                index=0
                            else:
                                index=1
                            axes_ani.plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 3,index, :], "-o", label=val_sat) #Frequencies E's
                        axes_ani.set_ylabel(str(val) + " m/s ")
                        axes_ani.legend()
                        axes_ani.set_yscale("log")
                return axes_used + len(user_select["EB_periodogram"])


            def EB_power_plot(axes_used):
                for l, val in enumerate(user_select["EB_cross power"]):
                    for k, val_sat in enumerate((user_select["satellite_graph"])):
                        if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                            lag_seconds = int(lag_data[2])
                            lag_fraction = lag_data[2] - lag_seconds
                            lag_nanoseconds = int(lag_fraction * 1e9)
                            # Create timedelta64 objects
                            seconds_delta = np.timedelta64(lag_seconds, 's')
                            nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                            # Add the timedelta to the datetime array
                            delta =seconds_delta + nanoseconds_delta
                        else:
                            delta=np.timedelta64(0,'s')
                        if val == 'ENorth/BEast crosspower':
                            index=0
                        else:
                            index=1
                        if user_select['bandpass'][0] == True:
                            bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                        else: bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), data[k, i, 5,index, bandpass], "-o",label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " nT ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_yscale("log")
                    axes_ani[axes_used + l ].set_xlim(user_select["bandpass"][1])
                return  axes_used + len(user_select["EB_cross power"])
            
            def EB_phase_plot(axes_used):
                for l, val in enumerate(user_select["EB_cross phase"]):
                    for k, val_sat in enumerate((user_select["satellite_graph"])):
                        if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                            lag_seconds = int(lag_data[2])
                            lag_fraction = lag_data[2] - lag_seconds
                            lag_nanoseconds = int(lag_fraction * 1e9)
                            # Create timedelta64 objects
                            seconds_delta = np.timedelta64(lag_seconds, 's')
                            nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                            # Add the timedelta to the datetime array
                            delta =seconds_delta + nanoseconds_delta
                        else:
                            delta=np.timedelta64(0,'s')

                        if val == 'ENorth/BEast cross phase':
                            index=0
                        else:
                            index=1
                        if user_select['bandpass'][0] == True:
                            bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                        else: bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), data[k, i, 6,index, bandpass],"-o", label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " nT ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_xlim(user_select["bandpass"][1])
                return  axes_used + len(user_select["EB_cross phase"])
            
            def BB_power_plot(axes_used):
                for l, val in enumerate(user_select["lags_cross_B"]):
                    for k, val_sat in enumerate((user_select["satellite_graph"])):
                        if val == 'B B lag cross power':
                            index=0
                        else:
                            index=1
                        if user_select['bandpass'][0] == True:
                            bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                        else: bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                        if index==0:
                            axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), np.absolute(data[k, i, 7,index, bandpass]), "-o",label=val_sat) #Frequencies E's
                            axes_ani[axes_used + l ].set_yscale("log")
                        else:
                            axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), data[k, i, 7,index, bandpass], "-o",label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val))
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_xlim(user_select["bandpass"][1])
                return  axes_used + len(user_select["lags_cross_B"])
            
            def EE_power_plot(axes_used):
                for l, val in enumerate(user_select["lags_cross_E"]):
                    for k, val_sat in enumerate((user_select["satellite_graph"])):
                        if val == 'E E lag cross power':
                            index=0
                        else:
                            index=1
                        if user_select['bandpass'][0] == True:
                            bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                        else: bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                        if index==0:
                            axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), np.absolute(data[k, i, 7,index, bandpass]), "-o",label=val_sat) #Frequencies E's
                            axes_ani[axes_used + l ].set_yscale("log")
                        else:
                            axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), data[k, i, 7,index, bandpass], "-o",label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val))
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_xlim(user_select["bandpass"][1])
                return  axes_used + len(user_select["lags_cross_E"])
            


            def B_Periodogram_Plot(axes_used):
                
                for l, val in enumerate(user_select["B_periodogram"]):
                    for k, val_sat in enumerate((user_select["satellite_graph"])):
                        if val == 'B_North':
                            index=0
                        else:
                            index=1
                        if user_select['bandpass'][0] == True:
                            bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                        else: bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]

                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), data[k, i, 2,index, bandpass]*1e9, label=val_sat) #Frequencies E's
                    
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " nT ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_yscale("log")
                return  axes_used + len(user_select["B_periodogram"])

            def E_Periodogram_Plot(axes_used):
                for l, val in enumerate(user_select["E_periodogram"]):
                    for k, val_sat in enumerate((user_select["satellite_graph"])):
                        if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                            lag_seconds = int(round(lag_data[2]*16))
                        else:
                            lag_seconds=0
                        if val == 'E_North':
                            index=0
                        else:
                            index=1
                        if user_select['bandpass'][0] == True:
                            bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                        else: bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, bandpass ]), data[k, i, 1,index, bandpass]*1e3, label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " mV/m ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_yscale("log")
                return  axes_used + len(user_select["E_periodogram"])
            

            def Time_Series_plot():
                colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
                #There are n number of axises that are contained for the time series
                for k in range(len(user_select["satellite_graph"])):
                    axes_ani[k].set_title("".join(["Swarm ", user_select["satellite_graph"][k]]))
                    twin_x_axes[k].clear()
                    twin_x_axes[k].yaxis.set_label_position("right")
                    if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                        lag_seconds = int(lag_data[2])
                        lag_fraction = lag_data[2] - lag_seconds
                        lag_nanoseconds = int(lag_fraction * 1e9)
                        # Create timedelta64 objects
                        seconds_delta = np.timedelta64(lag_seconds, 's')
                        nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                        # Add the timedelta to the datetime array
                        delta =seconds_delta + nanoseconds_delta
                    else:
                        delta=np.timedelta64(0,'s')
                    for l in range(len(user_select["Time_Series"])):
                        
                        if user_select['coordinate_system'] == "North East Centre":  #Make one plot over plotted with eachother
                            if indicies[k][i] == None:
                                pass
                            else:
                                if user_select["Time_Series"][l] == 'E_North': #Make Electric field twin x axis.
                                    twin_x_axes[k].plot(time_E[indicies[k][i]] ,  efield[k][indicies[k][i], 0]*1e3, color=colors[k+3], label='E North')
                                    twin_x_axes[k].set_ylabel("E (mV/m)")

                                elif user_select["Time_Series"][l] == "E_East":
                                    
                                    twin_x_axes[k].plot(time_E[indicies[k][i]] ,  efield[k][indicies[k][i], 1]*1e3, color=colors[k+4], label='E East')
                                    twin_x_axes[k].set_ylabel("E (mV/m)")
                                    
                                elif user_select["Time_Series"][l] == "B_East":
                                    axes_ani[k].plot(time_E[indicies[k][i]] ,  B_sinc[k][indicies[k][i], 1]*1e9, color=colors[k], label= 'B East')
                                    axes_ani[k].set_ylabel("B (nT)")
                                    pass
                                    
                                elif user_select["Time_Series"][l] == "B_North":
                                    axes_ani[k].plot(time_E[indicies[k][i]] ,  B_sinc[k][indicies[k][i], 0]*1e9, color=colors[k+1] , label='B North') 
                                    axes_ani[k].set_ylabel("custom_alt=TruenT)")
                    axes_ani[k].legend(loc=2)
                    twin_x_axes[k].legend(loc=1)
                    ##axes_ani[k].set_xlim(time_E[indicies[k][i]].min(),time_E[indicies[k][i]].max() )
                    
            


                
                return len(twin_x_axes)
            axes_used=0
            if user_select["Time_Series"] != None:
                axes_used=Time_Series_plot()
            if user_select["E_periodogram"] != None:
                axes_used = E_Periodogram_Plot(axes_used) 

            if user_select["B_periodogram"] != None:
                axes_used = B_Periodogram_Plot(axes_used)

            if user_select["EB_periodogram"] != None:
                axes_used = EB_Periodogram_Plot(axes_used)

            if user_select["EB_cross power"] != None:
                axes_used  = EB_power_plot(axes_used)
            if user_select["EB_cross phase"] != None:
                axes_used=EB_phase_plot(axes_used)
            if user_select["lags_cross_B"] !=None:
                axes_used=BB_power_plot(axes_used)
            if user_select["lags_cross_E"] !=None:
                axes_used=EE_power_plot(axes_used)
            return

        ani = animation.FuncAnimation(fig=fig_ani, func=animate, frames=frames) #What

        FFwriter = animation.FFMpegWriter(fps=int(6*sampling_rate_seconds))  #TODO given the framerate of 2 frames per second, which is equivalent to 6 seconds imager time = 1 second animation time, calculate the number of windows inbetween given a set step size
        ani.save("animationAlfven.mp4", writer=FFwriter)
        
        return 
    def graph_heatmap(data, datetimes):
        """
        Creates a heatmap of the periodograms in E,B and the ratio of E and B (depending on options selected)
        """
        def subplot_select():
            #TODO clean up with for loop and apply to other
            length_for_axis = 0
            try:
                length_for_axis += len(user_select["graph_E_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(user_select["graph_B_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(user_select["graph_PF_chosen"])
            except TypeError:
                pass
            if user_select["FAC"] == True:
                length_for_axis += 1
            else:
                pass
            if user_select ["E_difference"] != None:
                length_for_axis +=len(user_select["E_difference"])
            if user_select ["B_difference"] != None:
                length_for_axis +=len(user_select["B_difference"])
            if user_select ["PF_difference"] != None:
                length_for_axis +=len(user_select["PF_difference"])
            else:
                pass
            if user_select["Pixel_intensity"] == True:
                try:
                    length_for_axis += len(user_select["sky_map_values"])
                except TypeError:
                    raise ValueError("Sky Map not Selected")
            return length_for_axis
        length_for_axis=subplot_select()
            
        indicies=0
        for i in range(len(user_select["heatmap"])):
            for k in range(len(user_select["satellite_graph"])):
                if user_select['coordinate_system'] == "North East Centre":
                    if user_select["heatmap"][i] == "E_North":
                        index1 = 1
                        index2=  0
                        
                    elif user_select["heatmap"][i] == "E_East":
                        index1 = 1
                        index2=  1
                        scale=1e3
                    elif user_select["heatmap"][i] == "B_North":
                        index1 = 2
                        index2=  0
                        scale=1e3
                    elif user_select["heatmap"][i] == "B_East":
                        index1 = 2
                        index2=  1
                    elif user_select["heatmap"][i] == "ENorth/BEast ratio":
                        index1 = 3
                        index2=  0
                    elif user_select["heatmap"][i] == "EEast/BNorth ratio":
                        index1 = 3
                        index2=  1
                    elif user_select["heatmap"][i] == "ENorth/BEast crosspower":
                        index1 = 5
                        index2=  0
                    elif user_select["heatmap"][i] == "EEast/BNorth crosspower":
                        index1 = 5
                        index2=  1
                    elif user_select["heatmap"][i] == "ENorth/BEast cross phase":
                        index1 = 6
                        index2=  0
                    elif user_select["heatmap"][i] == "EEast/BNorth cross phase":
                        index1 = 6
                        index2=  1
                    elif user_select["heatmap"][i] == 'B B lag cross power':
                        index1 = 7
                        index2=  0
                    elif user_select["heatmap"][i] == 'B B lag cross phase':
                        index1 = 7
                        index2=  1
                    elif user_select["heatmap"][i] == 'E E lag cross power':
                        index1 = 8
                        index2=  0 
                    elif user_select["heatmap"][i] == 'E E lag cross phase':
                        index1 = 8
                        index2=  1
                    
                else:
                    if user_select["heatmap"][i] == "E_Azimuth":
                        index1 = 1
                        index2=  0
                    elif user_select["heatmap"][i] == "E_Polodial":
                        index1 = 1
                        index2=  1
                    elif user_select["heatmap"][i] == "B_Azimuth":
                        index1 = 2
                        index2=  0
                    elif user_select["heatmap"][i] == "B_Polodial":
                        index1 = 2
                        index2=  1
                    elif user_select["heatmap"][i] == "E Azimuth / B Polodial":
                        index1 = 3
                        index2=  0
                    elif user_select["heatmap"][i] == "E Polodial / B Azimuth":
                        index1 = 3
                        index2=  1
                if user_select['bandpass'][0] == True:
                    bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                else: 
                    bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]

                non_nan_mask = ~np.isnan(np.real(data[k,:,4,0,0]))

                # Get the indices of non-NaN elements
                non_nan_indices = np.where(non_nan_mask)[0]

                
                times= datetimes[np.real(np.array(data[k, non_nan_indices, 4,0, 0] , dtype=int))] 
                
                if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                    delta=timedelta(seconds=lag_data[2])
                else:
                    delta=timedelta(seconds=0)

                
                if index1==2: #B's
                    print(np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e9, 25))
                    dev=np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e9, 25)
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta ,
                                np.absolute(data[k, 0, 0, 0, bandpass]), np.real(np.array(data[k, :, index1, index2,bandpass]))[:,  non_nan_indices]*1e9,
                                shading='auto',
                                    norm=colors.LogNorm(vmin=dev), cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif  index1==1: #E's
                    print(np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e3, 25))
                    dev=np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e3, 25)
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta , 
                                    np.absolute(data[k, 0, 0, 0, bandpass]), np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e3 , shading='auto', 
                                    norm=colors.LogNorm(vmin=dev),
                                    cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif index1==3: #E/B ratio
                    print(np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices], 75))
                    dev=np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices], 75)
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta , 
                                    np.absolute(data[k, 0, 0, 0, bandpass]), np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices] , shading='auto', 
                                    norm=colors.LogNorm(vmax=dev),
                                    cmap='winter_r' ) #selects average time, frequencies, and then the periodogram 

                elif index1==5: #cross power
                    print(np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e9*1e3, 25))
                    dev=np.percentile(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e9*1e3, 25)
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta , 
                            np.absolute(data[k, 0, 0, 0, bandpass]), np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]*1e9*1e3 , shading='auto', 
                            norm=colors.LogNorm(vmin=dev),
                            cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif index1 ==6: #phase
                    #print(np.real(data[k, 0, 0, 0, :]), np.real(np.array(data[k, 0, index1, index2, :])).T, 'heatmap')
                    img=axes[length_for_axis+indicies].pcolormesh(times+delta , 
                                np.absolute(data[k, 0, 0,0, bandpass]),
                                (np.around(np.real(np.array(data[k, :, index1, index2, bandpass])[:,  non_nan_indices])/5, decimals=0))*5, shading='auto', cmap=plt.get_cmap('cmr.infinity')) #selects average time, frequencies, and then the periodogram 
                elif index1 ==7 and k==1 or  index1 ==8 and k==1: #B B component only goes once
                    if index2==0: #power
                        masked_data = np.ma.masked_invalid(np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices]) #for lognorm to deal with np.nans
                        img = axes[length_for_axis+indicies].pcolormesh(times + delta, 
                                np.absolute(data[k, 0, 0, 0, bandpass]),masked_data , shading='auto',
                                norm=colors.LogNorm()  ) #selects average time, frequencies, and then the periodogram #selects average time, frequencies, and then the periodogram 
                    if index2 ==1: #phase
                        img = axes[length_for_axis+indicies].pcolormesh(times + delta , 
                                np.absolute(data[k, 0, 0, 0, bandpass]),np.real(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices], shading='auto', cmap=plt.get_cmap('cmr.infinity')) #selects average time, frequencies, and then the periodogram #selects average time, frequencies, and then the periodogram 
                      
                else: #uses absolute value
                     
                     pass  
                              
                if index1 ==7 and k==0 or index1 ==8 and k==0:
                    pass

                else:
                    fig.colorbar(img, ax=axes[length_for_axis+indicies], extend='max', label=user_select["heatmap"][i])
                    axes[length_for_axis+indicies].set_ylabel(
                        "Frequencies"
                    )
                    axes[length_for_axis+indicies].set_xlabel(
                        "Times"
                    )
                    axes[length_for_axis+indicies].set_title(
                        user_select['satellite_graph'][k]
                    )
                    axes[length_for_axis+indicies].set_ylim(user_select["bandpass"][1])
                    axes[length_for_axis+indicies].set_xlim(user_select["time_range"])
                    indicies+=1
                    

        
        return
    
    def graph_singles(data):
        def subplot_select():
            #TODO clean up with for loop and apply to other
            length_for_axis = 0
            try:
                length_for_axis += len(user_select["graph_E_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(user_select["graph_B_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(user_select["graph_PF_chosen"])
            except TypeError:
                pass
            if user_select["FAC"] == True:
                length_for_axis += 1
            else:
                pass
            if user_select ["E_difference"] != None:
                length_for_axis +=len(user_select["E_difference"])
            if user_select ["B_difference"] != None:
                length_for_axis +=len(user_select["B_difference"])
            if user_select ["PF_difference"] != None:
                length_for_axis +=len(user_select["PF_difference"])
            else:
                pass
            if user_select["Pixel_intensity"] == True:
                try:
                    length_for_axis += len(user_select["sky_map_values"])
                except TypeError:
                    raise ValueError("Sky Map not Selected")
            return length_for_axis
        length_for_axis=subplot_select()
            
        indicies=0
        for i in range(len(user_select["singles_graph"])):
            for k in range(len(user_select["satellite_graph"])):
                factor=1
                print('weoweeee')
                if user_select['coordinate_system'] == "North East Centre":
                    if user_select["singles_graph"][i] == "E_North":
                        index1 = 1
                        index2=  0
                        factor=1e3
                    elif user_select["singles_graph"][i] == "E_East":
                        index1 = 1
                        index2=  1
                        factor= 1e3
                    elif user_select["singles_graph"][i] == "B_North":
                        index1 = 2
                        index2=  0
                        factor=1e9
                    elif user_select["singles_graph"][i] == "B_East":
                        index1 = 2
                        index2=  1
                        factor=1e9
                    elif user_select["singles_graph"][i] == "ENorth/BEast ratio":
                        index1 = 3
                        index2=  0
                    elif user_select["singles_graph"][i] == "EEast/BNorth ratio":
                        index1 = 3
                        index2=  1
                    elif user_select["singles_graph"][i] == "ENorth/BEast crosspower":
                        index1 = 5
                        index2=  0
                    elif user_select["singles_graph"][i] == "EEast/BNorth crosspower":
                        index1 = 5
                        index2=  1
                    elif user_select["singles_graph"][i] == "ENorth/BEast cross phase":
                        index1 = 6
                        index2=  0
                    elif user_select["singles_graph"][i] == "EEast/BNorth cross phase":
                        index1 = 6
                        index2=  1
                    elif user_select["singles_graph"][i] == 'B B lag cross power':
                        index1 = 7
                        index2=  0
                    elif user_select["singles_graph"][i] == 'B B lag cross phase':
                        index1 = 7
                        index2=  1
                    elif user_select["singles_graph"][i] == "ENorth/BEast coherence":
                        index1 = 8
                        index2=  0
                    elif user_select["singles_graph"][i] == "EEast/BNorth coherence":
                        index1 = 8
                        index2=  1
                    
                else:
                    if user_select["heatmap"][i] == "E_Azimuth":
                        index1 = 1
                        index2=  0
                    elif user_select["heatmap"][i] == "E_Polodial":
                        index1 = 1
                        index2=  1
                    elif user_select["heatmap"][i] == "B_Azimuth":
                        index1 = 2
                        index2=  0
                    elif user_select["heatmap"][i] == "B_Polodial":
                        index1 = 2
                        index2=  1
                    elif user_select["heatmap"][i] == "E Azimuth / B Polodial":
                        index1 = 3
                        index2=  0
                    elif user_select["heatmap"][i] == "E Polodial / B Azimuth":
                        index1 = 3
                        index2=  1
                if user_select['bandpass'][0] == True:
                    bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                else: 
                    bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                if index1 ==7 and index2==1 or index1==6: #phase
                    axes[length_for_axis+ indicies].plot( np.absolute(data[k, 0, 0, 0, bandpass]),np.real(data[k, 0, index1,index2,bandpass]))#frequency, data
                else:
                    axes[length_for_axis + indicies].plot(np.absolute(data[k, 0, 0,0,bandpass]), np.log10(np.absolute(data[k, 0, index1,index2,bandpass])*factor))
                axes[length_for_axis+indicies].set_title(
                        user_select['satellite_graph'][k]
                    )
                axes[length_for_axis+indicies].set_ylabel(user_select["singles_graph"][i])
                axes[length_for_axis+indicies].set_xlabel("frequencies Hz")
                
            indicies +=1



    time_range = user_select["time_range"]
    len_satellite = len(user_select["satellite_graph"])
    if user_select['heatmap'] !=None:
        sampling_rate_seconds=user_select["sampling_rate"]
        sampled_datetimes = create_sampled_datetimes(time_range, sampling_rate_seconds)

        length_of_windows=int(len(sampled_datetimes) - 1/sampling_rate_seconds *user_select['window_length'])

        window_length=user_select["window_length"]
    else:
        window_length=int((user_select["time_range_single"][1]-user_select["time_range_single"][0]).total_seconds())

    if user_select['nperseg'] == 'window length':
        nperseg=16*window_length
    if user_select['nperseg'] == 'half window length':
        nperseg=8*window_length
    if user_select['nperseg'] == 'quarter window':
        nperseg=4*window_length
    if user_select["heatmap"] !=None:
        data = np.zeros((len_satellite, length_of_windows, 9, 2,  nperseg), dtype=np.complex_) #[satelitte,number of windows, type of data, different polarizations,  data for each window ]
    else:
        data = np.zeros((len_satellite, 1, 9, 2,  nperseg), dtype=np.complex_)

    efield=np.multiply(efield,1e-3)
    
    B_sinc,B_resample = np.zeros(np.shape(efield)),np.zeros(np.shape(efield))
    indicies_total = []
    

    for k in range(len(user_select["satellite_graph"])): #length of satellites

        for i in range(3):

            B_sinc[k][:, i] =sinc_interpolation(bfield[k][:, i]*1e-9, time_B[k],time_E)
            B_resample[k][:,i]= signal.resample(bfield[k][:, i]*1e-9, len(time_E)) #As shown in testing, use sinc in time time domain, resample in spectral domain. Need to resample from 50, to 16 Hz for periodograms

        if user_select['heatmap'] != None:
            indicies_window=[]
            for i in range(length_of_windows): #Loops through each window and at the end stops early so window length doesnt cause error
                data[k, i], indicies = Logic_for_one_step(i, k,  B_resample, nperseg)
                indicies_window.append(indicies) #don't
            indicies_total.append(indicies_window)
        else:
            #User selected two times to be used in periodogram
            indicies_used_for_logic=find_indices(time_E, np.datetime64(user_select["time_range_single"][0]), np.datetime64(user_select["time_range_single"][1]))
            data[k, 0], indicies = Logic_for_one_step(indicies_used_for_logic, k, B_resample, nperseg, True   )

    if user_select["heatmap"] != None:
        graph_heatmap(data, sampled_datetimes)
        if user_select["animation"] != False:
            Animation(data,sampled_datetimes, B_sinc, B_resample, efield, indicies_total, length_of_windows, time_E)
        if user_select["conductivities"] != None:
            conductivities(data, sampled_datetimes)
    else:
        graph_singles(data)


    

    

    return  data
import numpy as np

def remove_spikes_with_interpolation(data: np.ndarray, z_thresh: float = 3.0) -> np.ndarray:
    """
    Removes spikes from an N×3 time series array by detecting large z-scores and interpolating over them.

    Parameters:
        data (np.ndarray): N×3 time series array (time along rows, 3 components as columns).
        z_thresh (float): Z-score threshold for detecting spikes.

    Returns:
        np.ndarray: Cleaned N×3 array with interpolated values.
    """
    cleaned = data.copy()
    N = data.shape[0]

    for i in range(3):  # for each component (column)
        x = cleaned[:, i]
        z = (x - np.mean(x)) / np.std(x)
        spike_idx = np.where(np.abs(z) > z_thresh)[0]

        if len(spike_idx) == 0:
            continue  # nothing to do

        # Indices of non-spike values
        valid_idx = np.setdiff1d(np.arange(N), spike_idx)
        valid_values = x[valid_idx]

        # Linear interpolation over spikes
        x[spike_idx] = np.interp(spike_idx, valid_idx, valid_values)

    return cleaned

def lag(x,y):
    """
    returns an array that is time-corrected to the first
    
    """
    #if len(x) != len(y):
    #    y=np.delete(y,-1, axis=0)
    correlation = signal.correlate(
    x[:,1], y[:,1], mode="full")
    lags = signal.correlation_lags(len(x[:,1]),
                               len(y[:,1]))
    delta=lags[np.argmax(correlation)]/50
    return correlation, lags,delta


def EBplotsNEC(user_select):
    def generalerrors():
        if user_select["Pixel_intensity"] == True:
            try:
                len(user_select["sky_map_values"])
            except TypeError:
                raise Exception("Please select a auroral animation as required")
        differences=[user_select['B_difference'], user_select['E_difference'], user_select['PF_difference']]
        for difference in differences:
            if difference and user_select['lag'] == False:
                raise Exception(f"Need to select lags for swarm A and C before selecting {difference}")
    generalerrors()
    data_returned=None
    set_token(
        "https://vires.services/ows",
        set_default=True,
        token="kmxv5mTTyYwzw4kQ9lsCkGfQHtjjJRVZ",
    )  # key
    plt.style.use("cyberpunk")  # Dark mode!
    print(user_select)
    print(user_select.keys())
    global has_E


    def empharisis_processing(cadence):
        emph = []
        space_craft = []
        data_stuff = []
        for i in range(len(collectionB_01)):
            ds = requester(
                collectionB_01[i],
                measurements[0],
                True,
                asynchronous=False,
                show_progress=False,
            )  # , sampling_step="PT{}S".format(cadence))

            if "".join(("swarm", ds["Spacecraft"][0].lower())) in labels:
                data_stuff.append(ds)
            else:
                pass
        for i in range(len(data_stuff)):
            time_array = data_stuff[i]["Spacecraft"].index.to_numpy()
            # Since time array is one hertz and we want to look for 1/cadence hertz we simply use
            lat_satellite_not_footprint = data_stuff[i]["Latitude"].to_numpy()
            lon_satellite_not_footprint = data_stuff[i]["Longitude"].to_numpy()
            altitude = data_stuff[i]["Radius"].to_numpy() / 1000 - 6378.14
            delete_array = np.linspace(0, cadence - 1, cadence, dtype=int)
            emph.append(
                [
                    np.delete(time_array, delete_array),
                    np.delete(lat_satellite_not_footprint, delete_array),
                    np.delete(lon_satellite_not_footprint, delete_array),
                    np.delete(altitude, delete_array),
                ]
            )
            # emph.append([time_array, lat_satellite_not_footprint,
            # lon_satellite_not_footprint, altitude])
            space_craft.append(data_stuff[i]["Spacecraft"].to_numpy()[0])
        # np.savetxt("emph_0200_0330.csv", emph[0], delimiter=",", fmt="%s")
        return emph, space_craft


    def rows():  # finds the number of rows needed for the produced figure
        length_for_axis = 0
        if user_select['graph_E_chosen'] != None:
            length_for_axis += len(user_select["graph_E_chosen"])

        if user_select['graph_B_chosen'] != None:
            length_for_axis += len(user_select["graph_B_chosen"])

        if user_select['graph_PF_chosen'] != None:
            length_for_axis += len(user_select["graph_PF_chosen"])
        if user_select["heatmap"] != None:
            length_for_axis += len(user_select["heatmap"])*len(user_select["satellite_graph"])
            if 'B B lag cross power' in user_select["heatmap"]:
                length_for_axis -=1
            if 'B B lag cross phase' in user_select["heatmap"]:
                length_for_axis -=1
            if 'E E lag cross power' in user_select["heatmap"]:
                length_for_axis -=1
            if 'E E lag cross phase' in user_select["heatmap"]:
                length_for_axis -=1
        if user_select["singles_graph"] !=None:
            length_for_axis+=len(user_select['singles_graph'])
        if user_select["conductivities"] != None:
            length_for_axis += len(user_select["conductivities"])

        if user_select["FAC"] == True:
            length_for_axis += 1

        if user_select ["E_difference"] != None:
            length_for_axis +=len(user_select["E_difference"])
        if user_select ["B_difference"] != None:
            length_for_axis +=len(user_select["B_difference"])
        if user_select ["PF_difference"] != None:
            length_for_axis +=len(user_select["PF_difference"])
        if user_select ["sky_map_values"] != None and user_select['Pixel_intensity']:
            length_for_axis += len(user_select["sky_map_values"])
        print(length_for_axis, 'rows')

        return length_for_axis


    fig, axes = plt.subplots(
        nrows=rows(),
        figsize=(10, 15),
        constrained_layout=True
    )
    time_range = user_select["time_range"]
    print(time_range)
    has_E = []  # Sets Which space-craft have a corersponding E field
    # Labels of space-craft interested in
    labels = user_select["satellite_graph"]
    # Measurement names from swarm
    measurements = [
        "B_NEC",
        [["VsatN", "VsatE", "VsatC"], ["Vixv", "Viy", "Viz"], "Quality_flags"],
    ]
    measurements_flat = [
        "VsatN",
        "VsatE",
        "VsatC",
        "Vixv",
        "Viy",
        "Viz",
        "Quality_flags",
    ]

    collectionE_16 = [
        "SW_EXPT_EFIA_TCT16",
        "SW_EXPT_EFIB_TCT16",
        "SW_EXPT_EFIC_TCT16",
    ]
    collectionE_02 = [
        "SW_EXPT_EFIA_TCT02",
        "SW_EXPT_EFIB_TCT02",
        "SW_EXPT_EFIC_TCT02",
    ]
    collectionB_50 = ["SW_OPER_MAGA_HR_1B", "SW_OPER_MAGB_HR_1B", "SW_OPER_MAGC_HR_1B"]
    collectionB_01 = ["SW_OPER_MAGA_LR_1B", "SW_OPER_MAGB_LR_1B", "SW_OPER_MAGC_LR_1B"]
    collectionF = [
        "SW_OPER_FACATMS_2F",
        "SW_OPER_FACBTMS_2F",
        "SW_OPER_FACCTMS_2F",
    ]  # Data packages from swarm
    twoHz = ["SW_OPER_EFIA_LP_1B", "SW_OPER_EFIB_LP_1B", "SW_OPER_EFIC_LP_1B"]

    # IF B, not selected but ponyting flux is, we assume 50hz data
    try:
        if user_select["B_frequency"][0] == "1Hz":
            collectionB = collectionB_01
        else:
            collectionB = collectionB_50
    except TypeError:  # B is none
        collectionB = collectionB_50
    try:
        if user_select["E_frequency"][0] == "2Hz":
            collectionE = collectionE_02
        else:
            collectionE = collectionE_16
    except TypeError:  # E is none
        collectionE = collectionE_16

    # Requests data from swarm
    def requester(sc_collection, measurement, residual, sampling_step=None, **kwargs):
        try:
            request = SwarmRequest()
            request.set_collection(sc_collection)
            if residual == True:
                request.set_products(
                    measurements=measurement,
                    models=["CHAOS"],
                    residuals=True,
                    sampling_step=sampling_step,
                )
            else:
                request.set_products(
                    measurements=measurement,
                    models=["CHAOS"],
                    sampling_step=sampling_step,
                )
            data = request.get_between(time_range[0], time_range[1], **kwargs)
            df = data.as_dataframe()
        except:
            df = []
        return df

    def graphingE(label, arrayx, arrayy, arraybandy, satelliteindex, has_twin):
        time=arrayx
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        global  axes_twin_E_y, axes_twin_E_x
        if user_select['lag'] and has_twin==False:
            axes_twin_E_y=[ axes[x].twiny() for x in range(len(user_select["graph_E_chosen"])) ]
        elif user_select['lag'] and rows() ==1 and has_twin==False:
            axes_twin_E_y= axes.twiny()

        if user_select["bandpass"][0] == True and has_twin==False:
            if rows() !=1 :
                axes_twin_E_x=[ axes[x].twinx() for x in range(len(user_select["graph_E_chosen"])) ]
            else:
                axes_twin_E_x= axes.twinx()

        for i in range(len(user_select["graph_E_chosen"])):
            if user_select["coordinate_system"] == "North East Centre":

                if user_select["graph_E_chosen"][i] == "North":

                    index = 0
                elif user_select["graph_E_chosen"][i] == "East":
                    index = 1
                elif user_select["graph_E_chosen"][i] == "Centre":
                    index = 2
            else:
                if user_select["graph_E_chosen"][i] == "Polodial": 
                    index = 0
                elif user_select["graph_E_chosen"][i] == "Azimuthal":
                    index = 1
                elif user_select["graph_E_chosen"][i] == "Mean-field":
                    index = 2
            
            try:
                if user_select['lag'] == True:
                    print(lag_data[3],process_string(label) , 'test')
                    if lag_data[3] == process_string(label):

                        lag_seconds = int(lag_data[2])
                        lag_fraction = lag_data[2] - lag_seconds
                        lag_nanoseconds = int(lag_fraction * 1e9)
                        # Create timedelta64 objects
                        seconds_delta = np.timedelta64(lag_seconds, 's')
                        nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                        
                        # Add the timedelta to the datetime array

                        time = arrayx + seconds_delta + nanoseconds_delta
                        print(arrayx)
                        arrayx_y = arrayx - (seconds_delta + nanoseconds_delta)
                        print(arrayx_y)
                        if rows() ==1:
                            axes_twin_E_y.plot(arrayx_y, [np.nan]*len(arrayx))
                            axes_twin_E_y.tick_params(axis='x', labelcolor=colors[satelliteindex])
                            axes_twin_E_y.set_xlim((min(arrayx_y), max(arrayx_y)))
                        else:
                            axes_twin_E_y[i].plot(arrayx_y, [np.nan]*len(arrayx))
                            axes_twin_E_y[i].tick_params(axis='x', labelcolor=colors[satelliteindex])
                            axes_twin_E_y[i].set_xlim((min(arrayx_y), max(arrayx_y)))
                        print((min(arrayx_y), max(arrayx_y)))

                        
            except NameError:
                    pass


            if user_select["bandpass"][0] == True:
                alpha=0.35
            else:
                alpha=1
            try:
                if has_twin == False and user_select['lag']:
                    axes[i].tick_params(axis='x', labelcolor=colors[satelliteindex])
                axes[i].plot(time, arrayy[:, index], label=label, color=colors[satelliteindex], alpha=alpha)
                axes[i].set_ylabel(
                    r"$V_{{{}}}$ $(m/s)$".format(user_select["graph_E_chosen"][i])
                )
                axes[i].set_title(f"Electric Field in the {user_select['graph_E_chosen'][i]} direction")
                axes[i].legend(loc=2)
                axes[i].set_xlim((time_range[0], time_range[1]))
                if user_select['singles_graph'] != None:
                    axes[i].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                    axes[i].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')

                if user_select["bandpass"][0] == True:
                    axes_twin_E_x[i].plot(time,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=1)
                    axes_twin_E_x[i].set_ylabel(
                    r"$V_{{{}}}$ bandpassed $(m/s)$".format(user_select["graph_E_chosen"][i])
                )
                    axes_twin_E_x[i].legend(loc=1)
            except TypeError:
                axes.plot(time, arrayy[:, index], label=label, color=colors[satelliteindex], alpha=alpha)
                if has_twin == False  and user_select['lag']:
                    axes.tick_params(axis='x', labelcolor=colors[satelliteindex])
                axes.set_ylabel(
                    r"$V_{{{}}}$ $(m/s)$".format(user_select["graph_E_chosen"][i])
                )
                axes.set_title(f"Electric Field in the {user_select['graph_E_chosen'][i]} direction")
                axes.legend(loc=2)
                axes.set_xlim((time_range[0], time_range[1]))
                if user_select['singles_graph'] != None:
                    axes.axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                    axes.axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')

                if user_select["bandpass"][0] == True:
                    axes_twin_E_x.plot(time,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=1)
                    axes_twin_E_x.set_ylabel(
                    r"$V_{{{}}}$ bandpassed $(m/s)$".format(user_select["graph_E_chosen"][i])
                )
                    axes_twin_E_x.legend(loc=1)
                    
        return True

    def graphingB(label, arrayx, arrayy, arraybandy, satelliteindex, has_twin, lag_bool):
        global axes_twin_B_x, axes_twin_B_y
        time=arrayx
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        
        if user_select["graph_B_chosen"] != None:
            length_for_axis = len(user_select["graph_B_chosen"])
        else:
            length_for_axis = 0
        
        if user_select['lag'] and has_twin==False:
            axes_twin_B_y=[ axes[x + length_for_axis].twiny() for x in range(len(user_select["graph_B_chosen"])) ]
        elif user_select['lag'] and rows() ==1 and has_twin==False:
            axes_twin_B_y= axes.twiny()

        if user_select["bandpass"][0] == True and has_twin==False:
            
            axes_twin_B_x=[ axes[x +length_for_axis].twinx() for x in range(len(user_select["graph_B_chosen"])) ] 


        for i in range(len(user_select["graph_B_chosen"])):
            if user_select['coordinate_system'] == "North East Centre":
                if user_select["graph_B_chosen"][i] == "North":
                    index = 0
                elif user_select["graph_B_chosen"][i] == "East":
                    index = 1
                elif user_select["graph_B_chosen"][i] == "Centre":
                    index = 2
            else:
                if user_select["graph_B_chosen"][i] == "Polodial":
                    index = 0
                elif user_select["graph_B_chosen"][i] == "Azimuthal":
                    index = 1
                elif user_select["graph_B_chosen"][i] == "Mean-field":
                    index = 2
            if user_select['lag'] == True and lag_bool:
                print(lag_data[3], label, process_string(label))
                if lag_data[3] == process_string(label) or lag_data[3]==label:
                    lag_seconds = int(lag_data[2])
                    lag_fraction = lag_data[2] - lag_seconds
                    lag_nanoseconds = int(lag_fraction * 1e9)
                    # Create timedelta64 objects
                    seconds_delta = np.timedelta64(lag_seconds, 's')
                    nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                    # Add the timedelta to the datetime array

                    time = arrayx + seconds_delta + nanoseconds_delta
                    print(arrayx)
                    arrayx_y = arrayx - (seconds_delta + nanoseconds_delta)
                    print(arrayx_y)
                    if rows() ==1:
                        axes_twin_B_y.plot(arrayx_y, [np.nan]*len(arrayx))
                        axes_twin_B_y.tick_params(axis='x', labelcolor=colors[satelliteindex])
                        axes_twin_B_y.set_xlim((min(arrayx_y), max(arrayx_y)))
                    else:
                        axes_twin_B_y[i].plot(arrayx_y, [np.nan]*len(arrayx))
                        axes_twin_B_y[i].tick_params(axis='x', labelcolor=colors[satelliteindex])
                        axes_twin_B_y[i].set_xlim((min(arrayx_y), max(arrayx_y)))
                    print((min(arrayx_y), max(arrayx_y)))

            if user_select["bandpass"][0] == True:
                alpha=0.35
            else:
                alpha=1
            if rows() ==1 :
                if has_twin == False and user_select['lag']:
                    axes.tick_params(axis='x', labelcolor=colors[satelliteindex])
                axes.plot(time, arrayy[:, index], label=label, color=colors[satelliteindex], alpha=alpha)
                axes.legend(loc=2)
                axes.set_ylabel(
                    r"$B_{{{}}}$".format(user_select["graph_B_chosen"][i]) + " (nT) "
                )
                axes.set_xlim((time_range[0], time_range[1]))
                axes.set_title(f"Magnetic Field in the {user_select['graph_B_chosen'][i]} direction")


                if user_select["bandpass"][0] == True:
                    axes_twin_B_x[i].plot(time,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=1)
                    axes_twin_B_x[i].set_ylabel(
                    r"$B_{{{}}}$ bandpassed $(nT)$".format(user_select["graph_B_chosen"][i])
                )
                    axes_twin_B_x[i].legend(loc=1)
            else:
                if has_twin == False and user_select['lag']:
                    axes[i + length_for_axis].tick_params(axis='x', labelcolor=colors[satelliteindex])
                    
                axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
                axes[i + length_for_axis].plot(time, arrayy[:, index], label=label, color=colors[satelliteindex], alpha=alpha)
                axes[i + length_for_axis].legend(loc=2)
                axes[i + length_for_axis].set_ylabel(
                    r"$B_{{{}}}$".format(user_select["graph_B_chosen"][i]) + " (nT) "
                )
                axes[i + length_for_axis].set_title(f"Magnetic Field in the {user_select['graph_B_chosen'][i]} direction")

                if user_select["bandpass"][0] == True:
                    axes_twin_B_x[i].plot(time,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=1)
                    axes_twin_B_x[i].set_ylabel(
                    r"$B_{{{}}}$ bandpassed $(nT)$".format(user_select["graph_B_chosen"][i])
                )
                    axes_twin_B_x[i].legend(loc=1)
                if user_select['singles_graph'] != None:
                    axes[i+ length_for_axis].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                    axes[i+ length_for_axis].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')
        return True

    def graphingFlux(label, arrayx, arrayy, index_band, satelliteindex, has_twin):
        global axes_twin_PF,  axes_twin_PF_y
        time=arrayx
        length_for_axis = 0
        # because python starts index at 0
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        try:
            length_for_axis += len(user_select["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_B_chosen"])
        except TypeError:
            pass

        if user_select['lag'] and has_twin==False:
            axes_twin_PF_y=[ axes[x + length_for_axis].twiny() for x in range(len(user_select["graph_PF_chosen"])) ]
        elif user_select['lag'] and rows() ==1 and has_twin==False:
            axes_twin_PF_y= axes.twiny()

        try:
            if user_select["bandpass"][0] == True and has_twin==False:
    
                axes_twin_PF=[ axes[x+length_for_axis].twinx() for x in range(len(user_select["graph_PF_chosen"])) ] 
        except TypeError:
            if user_select["bandpass"][0] == True and has_twin==False:
                axes_twin_PF=axes.twinx() 
        for i in range(len(user_select["graph_PF_chosen"])):
            if user_select['coordinate_system'] == "North East Centre":
                if user_select["graph_PF_chosen"][i] == "North":
                    index = 0
                elif user_select["graph_PF_chosen"][i] == "East":
                    index = 1
                elif user_select["graph_PF_chosen"][i] == "Centre":
                    index = 2
            else:
                if user_select["graph_PF_chosen"][i] == "Polodial":
                    index = 0
                elif user_select["graph_PF_chosen"][i] == "Azimuthal":
                    index = 1
                elif user_select["graph_PF_chosen"][i] == "Mean-field":
                    index = 2
            if user_select['lag'] == True:
                print(lag_data[3], label, process_string(label))
                if lag_data[3] == process_string(label) or lag_data[3]==label:
                    lag_seconds = int(lag_data[2])
                    lag_fraction = lag_data[2] - lag_seconds
                    lag_nanoseconds = int(lag_fraction * 1e9)
                    # Create timedelta64 objects
                    seconds_delta = np.timedelta64(lag_seconds, 's')
                    nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                    # Add the timedelta to the datetime array
                    time = arrayx + seconds_delta + nanoseconds_delta
                    print(arrayx)
                    arrayx_y = arrayx - (seconds_delta + nanoseconds_delta)
                    print(arrayx_y)
                    if rows() ==1:
                        axes_twin_PF_y.plot(arrayx_y, [np.nan]*len(arrayx))
                        axes_twin_PF_y.tick_params(axis='x', labelcolor=colors[satelliteindex])
                        axes_twin_PF_y.set_xlim((min(arrayx_y), max(arrayx_y)))
                    else:
                        axes_twin_PF_y[i].plot(arrayx_y, [np.nan]*len(arrayx))
                        axes_twin_PF_y[i].tick_params(axis='x', labelcolor=colors[satelliteindex])
                        axes_twin_PF_y[i].set_xlim((min(arrayx_y), max(arrayx_y)))
                    print((min(arrayx_y), max(arrayx_y)))



            try:
                if index_band==0:
                    if user_select['lag'] and has_twin==False:
                        axes[i + length_for_axis].tick_params(axis='x', labelcolor=colors[satelliteindex])
                    axes[i + length_for_axis].plot(time, arrayy[index], label=label, color=colors[satelliteindex])

                    axes[i + length_for_axis].set_ylabel(
                        r"$S_{{{}}}$".format(user_select["graph_PF_chosen"][i])
                        + r" $mW m^{-2}$"
                    )
                    axes[i + length_for_axis].set_title(f"Poynting Flux in the {user_select['graph_PF_chosen'][i]} direction")
                    axes[i + length_for_axis].legend(loc=2)
                    axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
                if index_band!=0 and user_select["bandpass"][0] == True:
                    axes_twin_PF[i].plot(time,arrayy[index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                    axes_twin_PF[i].legend(loc=1)
                if user_select['singles_graph'] != None:
                    axes[i+ length_for_axis].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                    axes[i+ length_for_axis].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')


            except TypeError:
                if index_band==0:
                    if user_select['lag']:
                        axes[i + length_for_axis].tick_params(axis='x', labelcolor=colors[satelliteindex])
                    axes.plot(time, arrayy[index], label="".join([label, "bandpassed"]), color=colors[satelliteindex])
                    axes.set_ylabel(
                        r"$S_{{{}}}$".format(user_select["graph_PF_chosen"][i])
                        + r" $mW m^{-2}$"
                    )
                    axes.legend(loc=2)
                    axes.set_xlim((time_range[0], time_range[1]))
                else:
                    axes_twin_PF.plot(time,arrayy[index], label="".join([label, "bandpassed"]),  color=colors[satelliteindex+3], alpha=0.7)
                    axes_twin_PF.legend(loc=1)
                    axes_twin_PF.set_ylabel(
                    r"$S_{{{}}}$ bandpassed".format(user_select["graph_PF_chosen"][i])
                    + r" $mW m^{-2}$"
                )
        return True

    def graphingF(label, arrayx, arrayy):
        length_for_axis = 0

        try:
            length_for_axis += len(user_select["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_B_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_PF_chosen"])
        except TypeError:
            pass

        axes[length_for_axis].plot(arrayx, arrayy, label=label)
        axes[length_for_axis].legend(loc=2)
        axes[length_for_axis].set_ylabel(r"Field Aligned Current $\mu A /m^2$")
        axes[length_for_axis].set_xlim((time_range[0], time_range[1]))
        if user_select['singles_graph'] != None:
                axes[length_for_axis].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                axes[length_for_axis].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')

    def graphingDifference(arrayx, time, label, i, plotted, fieldtype):
        length_for_axis = 0
        try:
            length_for_axis += len(user_select["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_B_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_PF_chosen"])
        except TypeError:
            pass
        if user_select["FAC"] == True:
            length_for_axis += 1
        else:
            pass
        first_non_nan_index = np.where(~np.isnan(arrayx[1]-arrayx[0]))[0][0]
        axes[length_for_axis +i].set_xlim(first_non_nan_index, len(arrayx[0]))
        axes[length_for_axis +i].plot(arrayx[1]-arrayx[0], label='Subtracted')
        max_lim=np.nanmax(np.absolute(arrayx[1]-arrayx[0]))
        axes[length_for_axis +i].set_ylim(-max_lim,max_lim)
        axes_2=axes[length_for_axis +i].twinx()
        for l in range(2):
            axes_2.plot(arrayx[l], label=label[l], alpha=0.2)
        max_lim=np.nanmax(np.absolute(arrayx))
        axes_2.set_ylim(-max_lim,max_lim)
        axes[length_for_axis +i].legend()
        if user_select['singles_graph'] != None:
                axes[length_for_axis +i].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                axes[length_for_axis + i].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')
        axes[length_for_axis + i].set_title(f"Difference of {fieldtype} in the {plotted} direction with appropriate units")

    def Graphing_skymap(pixel, time, spacecraft):
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        length_for_axis = 0

        try:
            length_for_axis += len(user_select["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_B_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(user_select["graph_PF_chosen"])
        except TypeError:
            pass
        if user_select["FAC"] == True:
            length_for_axis += 1
        else:
            pass
        if user_select ["E_difference"] != None:
            length_for_axis +=len(user_select["E_difference"])
        if user_select ["B_difference"] != None:
            length_for_axis +=len(user_select["B_difference"])
        if user_select ["PF_difference"] != None:
            length_for_axis +=len(user_select["PF_difference"])

        if user_select['lag']:
            axes_twin_y=[ axes[x + length_for_axis].twiny() for x in range(len(pixel))]
            
        for i in range(len(pixel)):  # length of platforms basically
            arrayx = time[i]
            for j in range(len(pixel[0])):  # Length of satellites selected
                label="".join(["swarm ", spacecraft[j]])
                for k in range(len(pixel[i][j])):
                    if pixel[i][j][k] == 0:
                        pixel[i][j][k] = np.nan
                if user_select['lag']:
                    if lag_data[3] == process_string(label) or lag_data[3]==label:
                        lag_seconds = int(lag_data[2])
                        lag_fraction = lag_data[2] - lag_seconds
                        lag_nanoseconds = int(lag_fraction * 1e9)
                        # Create timedelta64 objects
                        seconds_delta = np.timedelta64(lag_seconds, 's')
                        nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                        # Add the timedelta to the datetime array
                        arrayx = time[i] + seconds_delta + nanoseconds_delta
                        arrayy= time[i] -( seconds_delta + nanoseconds_delta)
                        print(arrayy)
                        axes_twin_y[i].plot(arrayy, [np.nan]*len(arrayx))
                        axes_twin_y[i].tick_params(axis='x', labelcolor=colors[i])
                        axes_twin_y[i].set_xlim((min(arrayy), max(arrayy)))
                else:
                    arrayx = time[i]
                axes[i + length_for_axis].plot(
                    arrayx, pixel[i][j], label=label)
            axes[i + length_for_axis].legend(loc=2)
            if user_select['lag']:
                if lag_data[3] != process_string(label):
                    axes[i + length_for_axis].tick_params(axis='x', labelcolor=colors[i])
            axes[i + length_for_axis].set_ylabel("Nearest Pixel intensity")
            axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
            # axes[i + length_for_axis].set_xlim(65,70)if user_select['lag'] == True and lag_bool:
                


    def Coordinate_change(lattiude, longitude, radius):  # Coordinate change, from Ivan correspondence
        a, b, e2 = 6378137.0, 6356752.3142, 0.00669437999014  # From DRS80
        lat, lon, h = np.deg2rad(lattiude), np.deg2rad(longitude), radius
        v = a / np.sqrt(1 - e2 * np.sin(lat) * np.sin(lat))  # logic
        x = (v + h) * np.cos(lat) * np.cos(lon)
        y = (v + h) * np.cos(lat) * np.sin(lon)
        z = (v * (1 - e2) + h) * np.sin(lat)
        return np.array([x, y, -1*z])

    def requesterarraylogic():
        nonlocal data_returned
        def V():
            return_data_non_band = []
            return_data_band =[]
            satellites_with_E = []
            times_for_flux = []
            twin_axis=False
            lag_bool=False
            lon_1=None
            Etime=[] #returns empty in case no velocity data found, required
            for i in range(len(collectionE)):
                dsV = requester(  # requests data
                    collectionE[i],
                    measurements_flat,
                    False,
                    asynchronous=False,
                    show_progress=False,
                )
                if len(dsV) != 0:  # checks if empty
                    # Checks if space-craft is selected
                    if "".join(("swarm", dsV["Spacecraft"][0].lower())) in labels:
                        

                        has_E.append(True)
                        satellites_with_E.append(
                            "".join(("swarm", dsV["Spacecraft"][0].lower()))
                        )
                        dsB = requester( 
                        collectionB_50[i], #Mag B, high resolution, 50Hz B (Magnetic field)
                        ["q_NEC_CRF", measurements[0]], #Magnetic field in NEC coordinates
                        True, 
                        asynchronous=False,
                        show_progress=False)
                        indicies=find_closest_indices(dsV.index, dsB.index)
                        quatnecrf=dsB["q_NEC_CRF"].to_numpy()[indicies]
                        Vsat=np.array([dsV["Vixv"] , dsV["Viy"], dsV["Viz"]]).T
                        Etime = dsV.index
                        VNEC=[]
                        for l in range(len(quatnecrf)):
                            inverse_quat = quaternion_inverse_scipy(dsB["q_NEC_CRF"].to_numpy()[indicies][l])
                            rot_NEC_V= inverse_quat.apply(Vsat[l])
                            VNEC.append(rot_NEC_V)
                        VelocityNEC=remove_spikes_with_interpolation(np.array(VNEC))
                        VelocityNECbandpass = np.zeros(np.shape(VelocityNEC))
                        if user_select["bandpass"][0] == True:
                            for l in range(3):  # moving average of bres for all three components
                                VelocityNECbandpass[:, l] = butter_bandpass_filter(data=VelocityNEC[:, l], cutoffs=user_select["bandpass"][1], fs=16)

                        def B_Logic_For_E(lag_bool):
                            nonlocal lon_1
                            # finds the closest index in B for each time
                            def Time_corrections(E_time, B_time):
                                # finds the closest time value for each element in E_time
                                closest_time = np.zeros(len(E_time), dtype=int)
                                # resource intesiive, find something better
                                for i in range(len(E_time)):
                                    closest_time[i] = np.argmin(
                                        np.abs(E_time[i] - B_time)
                                    )
                                return closest_time
                            print(collectionB)
                            dsB = requester(
                                collectionB[i],
                                measurements[0],
                                False,
                                asynchronous=False,
                                show_progress=False,
                            )

                            times_of_b_for_flux = Time_corrections(
                                dsV.index.to_numpy(), dsB.index.to_numpy()
                            )  # Takes times of both E and B and finds the closest values in B to E

                            meanfield=arrangement(dsB.index.to_numpy(),dsB["B_NEC_CHAOS"].to_numpy(),3)[times_of_b_for_flux, :]
                            #TODO if swarm a and c selected
                            global B_lag_earlier 
                            if user_select['lag'] and lag_bool:
                                B_synoptic = arrangement(dsB.index.to_numpy(), dsB["B_NEC"].to_numpy()-dsB["B_NEC_CHAOS"].to_numpy(),3)
                                global lag_data
                                lag_data=list(lag(B_lag_earlier,B_synoptic))
                                lag_data.append("".join(("swarm", dsV["Spacecraft"][0].lower()))) 
                            else:
                                lag_bool=True
                                ##TODO make cleaner
                                B_lag_earlier = arrangement(dsB.index.to_numpy(), dsB["B_NEC"].to_numpy()-dsB["B_NEC_CHAOS"].to_numpy(),3)
                                lon_1=dsB['Longitude'].to_numpy()

                            return times_of_b_for_flux, meanfield,  lag_bool

                        times_for_flux, meanfield,lag_bool = B_Logic_For_E(lag_bool)

                        if user_select['coordinate_system'] == "Mean-field aligned":
                            latitude, longitude, radius = dsV['Latitude'].to_numpy(), dsV['Longitude'].to_numpy(),  dsV["Radius"].to_numpy()  # Gets Emphermis data
                            r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                          
                            VelocityData = MFA(VelocityNEC, meanfield, r_nec.T) #puts data into MFA coordinate system
                            VelocityNECbandpass = MFA(VelocityNECbandpass, meanfield, r_nec.T)
                        else:
                            VelocityData = VelocityNEC

                        # Plots velocity field time seres

                        if user_select["graph_E_chosen"] != None:
                            twin_axis=graphingE(
                                "".join(("Swarm ", dsV["Spacecraft"][0])),
                                dsV.index.to_numpy(),
                                VelocityData,
                                VelocityNECbandpass,
                                sum(has_E)-1, #finds the amounts of trues in has_E then -'s by 1 to make it start at -1, this is so the color is the same as B and flux
                                twin_axis
                            )
                        # For ponyting flux
                        return_data_band.append(VelocityNECbandpass)
                        return_data_non_band.append(VelocityData)  

                    else:
                        # Says theres no E component
                        has_E.append(False)
                else:  # Says theres no E component
                    has_E.append(False)

            return return_data_non_band, return_data_band, times_for_flux, Etime, satellites_with_E

        def B():
            return_data_non_band = []
            return_data_band =[]
            time_array=[]
            twin_axis=False
            blabels=[]
            lag_bool=None
            B_lag_earlier = None
            lon_1=None
            # Goes through every satellite selected
            for i in range(len(user_select["satellite_graph"])):
                if user_select["satellite_graph"][i] == "epop":
                    raise NotImplementedError("Not implemented")
                else:
                    # Data package level, data/modes, **kwargs
                    if user_select["satellite_graph"][i] == "swarma":
                        j = 0
                    if user_select["satellite_graph"][i] == "swarmb":
                        j = 1
                    if user_select["satellite_graph"][i] == "swarmc":
                        j = 2

                    dsmodel_res = requester(
                        collectionB[j],
                        measurements[0],
                        True,
                        asynchronous=False,
                        show_progress=False,
                    )
                    ds = requester(
                        collectionB[j],
                        measurements[0],
                        False,
                        asynchronous=False,
                        show_progress=False,
                    )

                    bmodel_res, bmodel_actual = dsmodel_res["B_NEC_res_CHAOS"],ds["B_NEC_CHAOS"]

                    Bdata = ds["B_NEC"]  # data package
                    time = Bdata.index.to_numpy()
                    Bdata = Bdata.to_numpy()

                    barranged = arrangement(time, Bdata, 3)
                    bmodelarranged = arrangement(time, bmodel_actual, 3)

                    bresarranged = np.subtract(barranged, bmodelarranged)
                    bresarrangedband=np.zeros(np.shape(bresarranged))
                    if user_select["bandpass"][0] == True:
                        for l in range(3):  # moving average of bres for all three components
                            bresarrangedband[:, l] = butter_bandpass_filter(bresarranged[:, l], user_select["bandpass"][1], 50)
                    
                    ##TODO Add MFA
                    if user_select['coordinate_system'] == "Mean-field aligned":
                        latitude, longitude, radius = ds['Latitude'].to_numpy(), ds['Longitude'].to_numpy(),  ds["Radius"].to_numpy()-6.371e6  # Gets Emphermis data
                        
                        r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                        Bdata = MFA(bresarranged, bmodelarranged, r_nec.T) #puts data into MFA coordinate system
                        Bdataband =MFA(bresarrangedband, bmodelarranged,r_nec.T)
                    else:
                        Bdata = bresarranged
                        Bdataband = bresarrangedband
                    if user_select['lag'] and lag_bool == False :
                        B_synoptic = arrangement(ds.index.to_numpy(), ds["B_NEC"].to_numpy()-ds["B_NEC_CHAOS"].to_numpy(),3)
                        global lag_data
                        lag_data=list(lag(B_lag_earlier,B_synoptic))
                        lag_data.append("".join(("swarm", ds["Spacecraft"][0].lower()))) 
                        lag_bool=True
                        print(lag_data)
                        st.write(f"The lag in seconds is {lag_data[2]}")
                        st.write(f"The average longindinual seperation during the pass is {long_sep(lon_1, ds['Longitude'].to_numpy(), lag_data,50)} in degrees")
                    elif user_select['lag'] and lag_bool == None:
                        lag_bool=False
                        ##TODO make cleaner
                        B_lag_earlier = arrangement(ds.index.to_numpy(), ds["B_NEC"].to_numpy()-ds["B_NEC_CHAOS"].to_numpy(),3)
                        lon_1=ds['Longitude'].to_numpy()
                    print(lag_bool)
                    if user_select["graph_B_chosen"] != None: #plots if selected
                        twin_axis=graphingB(
                            "".join(("Swarm ", dsmodel_res["Spacecraft"][0])),
                            time,
                            Bdata,
                            Bdataband,
                            i,
                            twin_axis,
                            lag_bool,
                        )


                    return_data_non_band.append(Bdata)
                    return_data_band.append(Bdataband)
                    time_array.append(time)
                    blabels.append("".join(("Swarm ", dsmodel_res["Spacecraft"][0])))
            else:
                pass
            return return_data_non_band,return_data_band, time_array, blabels

        def F(space_craft_with_E=["swarma", "swarmb", "swarmc"]):
            data_return = []
            data_stuff = []
            for i in range(len(collectionF)):
                ds = requester(
                    collectionF[i],
                    "FAC",
                    False,
                    asynchronous=False,
                    show_progress=False,
                )
                if "".join(("swarm", ds["Spacecraft"][0].lower())) in labels:
                    fac = ds["FAC"]
                    time = fac.index
                    fac = ds["FAC"].to_numpy(fac)

                    graphingF("".join(("Swarm ", ds["Spacecraft"][0])), time, fac)
                    if (
                        "".join(("swarm", ds["Spacecraft"][0].lower()))
                        in space_craft_with_E
                    ):
                        data_return.append(
                            [
                                "".join(("Swarm ", ds["Spacecraft"][0])),
                                time.to_numpy(),
                                fac,
                            ]
                        )
                    data_stuff.append(ds)
                else:
                    pass
            return data_return, data_stuff

        if (
            user_select["graph_E_chosen"] == None
            and user_select["graph_PF_chosen"] == None
            and user_select['E_B_ratio'] == False
        ):
            pass
        else:
            # E field, indexs of B for time, times for plotting flux
            efield_nonband, efield_band, times_for_b, time_E, space_craft_with_E = V() #TODO need efield  here
        if (
            user_select["graph_B_chosen"] == None
            and user_select["graph_PF_chosen"] == None
            and user_select['E_B_ratio'] == False
        ):
            pass
        else:
            bfield_non_band, bfield_band, time_B, blabels = B()

        if user_select["FAC"] == True:
            try:
                FAC_data, Low_hertz_dataframe = F(space_craft_with_E)
            except:
                FAC_data, Low = F()


        def pontying_flux(efieldband,efieldnonband, bfieldband,bfieldnonband):  # Take rough estimate by doing N cross E to get C poynting flux
            return_data = []
            looping_list=[[efieldnonband, bfieldnonband], [efieldband, bfieldband]]
            twin_axis=False
            band_or_no=["nonband", "band"]
            for k in range(2):
                efield=looping_list[k][0]
                bfield=looping_list[k][1]
                for i in range(len(efield)):
                    eflux = efield[i]
                    bflux=np.empty(np.shape(eflux))

                    for l in range(3):
                        bflux[:,l] = sinc_interpolation(x=bfield[i][:,l], s=time_B[i], u=time_E)

                    flux = np.cross(eflux * 1e-3, bflux * 1e-9) * 1.256e6*1e3
                    flux_individual = np.transpose(flux)

                    twin_axis=graphingFlux(space_craft_with_E[i], time_E, flux_individual, k, i, twin_axis)
                    return_data.append(
                        np.transpose(flux_individual)
                    )
            return return_data

        def Difference_plots(bfield,  btime, blabel, efield, etime, elabel, poynting_flux_data,poynting_flux_time, poynting_flux_label):

            data=[[user_select["B_difference"], bfield,btime,blabel], [user_select["E_difference"], efield,etime,elabel], [user_select["PF_difference"],
                    poynting_flux_data, poynting_flux_time, poynting_flux_label]] #creates a list of the different field types the user may want to difference
            field_type=["Magnetic", "Electric", "Poynting Flux"] #list of fields
            indicies_used=0 #holder to say which axes to plot on 
            for i in range(len(data)): #Loops through data list
                if data[i][0] == None: #if not selected pass
                    pass
                else:
                    for k in range(len(data[i][0])): #Loop through different coordinates
                        for index in range(len(data[i][3])): #Loops through satellite

                            if lag_data[3] == process_string(data[i][3][index]):
                                if index==1: #Ask's if the space-craft is the lagged one or fixed one
                                    indexopposite=0
                                else:
                                    indexopposite=1
                                if i==0: #Only magnetic field has an sps 50, rest have 16
                                    lag_used=np.round(lag_data[2]*50).astype(int)
                                else:
                                    lag_used=np.round(lag_data[2]*16).astype(int)
                                if  data[i][0][k] == "North": #Asks the coordinate system #TODO add MFA
                                    index_coordinate=0
                                elif data[i][0][k] == "East":
                                    index_coordinate=1
                                else:
                                    index_coordinate=2

                                if lag_used > 0: #Offsets the time series of the lagged satellite and cuts the ends off to make sure its the same length as the other
                                    offset_ts1 = np.concatenate((np.full(lag_used, np.nan), data[i][1][index][:, index_coordinate]))[:len(data[i][1][indexopposite][:,index_coordinate])]
                                    offset_ts2 = data[i][1][indexopposite][:,index_coordinate]
                                elif lag_used < 0:
                                    offset_ts1 = data[i][1][index][:, index_coordinate]
                                    offset_ts2 = np.concatenate((np.full(-lag_used, np.nan), data[i][1][indexopposite][:,1]))[:len(data[i][1][index][:, index_coordinate])]
                        if len(data[i][2]) ==  len(data[i][1]): #Asks if the length of time series, E's are always synced up but B's differ
                            graphingDifference(
                                [offset_ts1,offset_ts2], #data
                                data[i][2][indexopposite], #time
                                data[i][3], #label
                                indicies_used, #axes 
                                data[i][0][k], # coordinate selected
                                field_type[i] #field type
                            )
                        else:
                            graphingDifference(
                                [offset_ts1,offset_ts2], #data
                                data[i][2], #time
                                data[i][3], #label
                                indicies_used, #axes 
                                data[i][0][k], # coordinate selected
                                field_type[i] #field type
                            )

                        indicies_used+=1 #increments the axes plotted
 
        def skymap():
            pixel_chosen_total = []
            sat_time_each_platform = []
            pixel_chosen_average_total = []
            lat_satellite, lon_satellite = [], []
            for k in range(len(user_select["sky_map_values"])):
                asi_array_code = user_select["sky_map_values"][k][0]
                location_code = user_select["sky_map_values"][k][1]
                alt = int(user_select["sky_map_values"][k][2])


                def ASI_logic():
                    if asi_array_code.lower() == "themis":
                        if alt==90 or alt==110 or alt==150:
                            asi = asilib.asi.themis(
                                location_code,
                                time_range=time_range,
                                alt=alt,
                                custom_alt=False,
                            )
                        else:
                            asi = asilib.asi.themis(
                                location_code,
                                time_range=time_range,
                                alt=alt,
                                custom_alt=True,
                            )
                        cadence = 3
                    elif asi_array_code.lower() == "rego":
                        asi = asilib.asi.rego(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt=True,
                        )
                        cadence = 3
                    elif asi_array_code.lower() == "trex_nir":
                        asi = asilib.asi.trex.trex_nir(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt=True,
                        )
                        cadence = 6
                    elif asi_array_code.lower() == "trex_rgb":
                        cadence = 3
                        asi = asilib.asi.trex.trex_rgb(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            colors="rgb",
                            custom_alt=True,
                        )

                    return asi, cadence

                asi, cadence = ASI_logic()
                emph, space_craft_label = empharisis_processing(cadence)

                # return cloest_time_array

                for i in range(len(emph)):  # Conjunction logic
                    # Trex is 6 second cadence compared to 3 of rego and themos

                    sat_time = np.array(emph[i][0])  # sets timestamp

                    sat_lla = np.array([emph[i][1], emph[i][2], emph[i][3]]).T
    
                    conjunction_obj = asilib.Conjunction(asi, (sat_time, sat_lla))

                    lat_sat=conjunction_obj.sat["lat"].to_numpy()
                    lon_sat=conjunction_obj.sat["lon"].to_numpy()

                    alt_sat=conjunction_obj.sat["alt"].to_numpy()

                    # Converts altitude to assumed auroral height
                    
                    sat_lla= footprint(sat_time[0],lat_sat,lon_sat,alt_sat, alt)
                    conjunction_obj = asilib.Conjunction(asi, (sat_time, sat_lla.T))

                    lat_satellite.append(conjunction_obj.sat["lat"].to_numpy())
                    lon_satellite.append(conjunction_obj.sat["lon"].to_numpy())

                lat, lon = asi.skymap["lat"], asi.skymap["lon"]

                lat[np.isnan(lat)] = np.inf
                lon[np.isnan(lon)] = np.inf

                pixel_chosen = np.zeros((len(emph), len(emph[0][0])))
                values = np.zeros((np.shape(lat)))
                indicies_total = np.zeros((len(emph), len(emph[0][0]), 2), dtype=int)
                index_of_image = 0
                pixel_chosen_average = np.zeros((len(emph), len(emph[0][0])))

                def average(start_index,grid, subregion_size=user_select['pixel_average']):
                    """w
                    Calculate the average of a subregion in a grid.
                    
                    Parameters:
                    grid (np.ndarray): The input 2D grid of values.
                    start_index (tuple): A tuple (row, col) indicating the starting index of the subregion.
                    subregion_size (tuple): A tuple (height, width) indicating the size of the subregion.
                    
                    Returns:
                    float: The average value of the subregion.
                    """
                    row, col = start_index
                    row=int(row)
                    col=int(col)
                    height, width = subregion_size, subregion_size
                    
                    # Extract the subregion
                    ##if (row + height > grid.shape[0]) or (col + width > grid.shape[1]):
                        ##raise ValueError("Subregion exceeds grid boundaries.")
                    subregion = grid[row:row+height, col:col+width]
                    
                    # Calculate the average of the subregion
                    average_value = np.sum(subregion)

                    return average_value

                for i in range(len(emph[0][0])):  # len of time series
                    # Themis starts at time_range[0], rego and trex start a time_range[0] + a cadence
                    if len(asi.data[0]) != len(emph[0][0]) and i == 0:
                        continue
                    elif len(asi.data[0]) != len(emph[0][0]) and i != 0:
                        i -= 1

                    # blind optimzation, loops through time range, satellites,
                    # and lats and lons of the skymap to find the closest pixel to the satellite which can then be used to
                    # find intensities

                    # if(indicies_total[0][0][0] == 0 and indicies_total[0][0][1] == 0):
                    # print(np.shape(lat)[0], np.shape(lat)[1], emph[0][0][i])
                    for satellite in range(len(emph)):
                        for j in range(np.shape(lat)[0]):
                            for k in range(np.shape(lat)[1]):
                                values[j][k] = np.sqrt(
                                    (lat[j][k] - lat_satellite[satellite][i]) ** 2
                                    + (lon[j][k] - lon_satellite[satellite][i]) ** 2
                                )
                        indicies = np.unravel_index(
                            (np.abs(values)).argmin(), values.shape
                        )
                        indicies_total[satellite, i] = indicies
                        try:
                            pixel_chosen[satellite][i] = asi.data[1][index_of_image][
                                indicies[0], indicies[1]
                            ]
                            pixel_chosen_average[satellite][i] = average(
                                indicies, asi.data[1][index_of_image]
                            )
                        except ValueError:  # rgb
                            pixel_chosen[satellite][i] = asi.data[1][index_of_image][
                                indicies[0], indicies[1], 1
                            ]
                            pixel_chosen_average[satellite][i] = average(
                                indicies, asi.data[1][index_of_image][:, :, 1]
                            )
                        except IndexError:
                            pixel_chosen[satellite][i] = 0

                    if np.mod(i, cadence) == 0 and i != 0:
                        index_of_image += 1

                        indicies = np.unravel_index(
                            (np.abs(values)).argmin(), values.shape
                        )
                        indicies_total[satellite, i] = indicies
                        try:
                            pixel_chosen[satellite][i] = asi.data[1][index_of_image][
                                indicies[0], indicies[1]
                            ]
                            pixel_chosen_average[satellite][i] = average(
                                indicies, asi.data[1][index_of_image]
                            )
                        except ValueError:  # rgb
                            pixel_chosen[satellite][i] = asi.data[1][index_of_image][
                                indicies[0], indicies[1], 1
                            ]
                            pixel_chosen_average[satellite][i] = average(
                                indicies, asi.data[1][index_of_image][:, :, 1]
                            )
                        except IndexError:
                            pixel_chosen[satellite][i] = 0
                #TODO implement a GUI option for 5x5 averaged, 3x3 averaged or single
                pixel_chosen_total.append(pixel_chosen)
                sat_time_each_platform.append(sat_time)

                pixel_chosen_average_total.append(pixel_chosen_average)
                #Temporary
                #pixel_chosen_average_total.append()


            return pixel_chosen_average_total, sat_time_each_platform, space_craft_label
        


        if user_select["graph_PF_chosen"] == None:
            pass
        else:
            flux = pontying_flux(efield_band,efield_nonband, bfield_band,bfield_non_band)

        if user_select['E_difference'] !=None or user_select["B_difference"] !=None or user_select["PF_difference"] !=None:
            if user_select['E_difference'] == None:
                efield=None
                etime=None
                elabel=None
            else:
                if user_select["bandpass"][0] == True:
                    efield = efield_band
                else:
                    efield = efield_nonband
                etime=time_E
                elabel=space_craft_with_E
            if user_select["B_difference"] == None:
                bfield=None
                btime=None
                blabel=None
            else:
                if user_select["bandpass"][0] == True:
                    bfield = bfield_band
                else:
                    bfield = bfield_non_band
                btime=time_B
                blabel=blabels
            if user_select["PF_difference"] == None:
                PFfield=None
                PFtime=None
                PFlabel=None
            else:
                PFtime=time_E
                PFlabel=space_craft_with_E
                PFfield=flux

            Difference_plots(bfield, btime, blabel, efield,etime,elabel,PFfield,PFtime,PFlabel)

        ##TODO Implement ratio
        if user_select["E_B_ratio"] == True:
            if user_select["bandpass"][0] == True:
                efield = efield_band
                bfield = bfield_band
            else:
                efield = efield_nonband
                bfield = bfield_non_band
            data_returned=Graphing_Ratio(
                space_craft_with_E, efield, bfield, time_E, time_B, user_select, fig, axes
            )
        if user_select["sky_map_values"] != None and user_select["Pixel_intensity"] == True:

            pixels, time, space_craft = skymap()
            Graphing_skymap(pixels, time, space_craft)

        else:
            pass

    requesterarraylogic()

    fig.suptitle("Time Versus Auroral Parameters From Swarm Spacecraft")

    # for i in range(len(axes)):  # final touches
    # mplcyberpunk.make_lines_glow(axes[i])  # adds glow to plots
    # mplcyberpunk.add_gradient_fill(
    # ax=axes[i], alpha_gradientglow=0.8, gradient_start="zero")

    return fig, axes, data_returned