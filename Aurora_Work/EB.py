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
import pickle
from scipy.fft import fft, fftfreq
import asilib
import asilib.asi
import cdflib
import aacgmv2
from scipy import signal
from scipy.signal import butter, filtfilt, freqz
from numpy.typing import NDArray
import matplotlib.animation as animation
import matplotlib.colors as colors
import cmasher as cmr
def process_string(input_string):
    # Convert the string to lowercase
    lowercased_string = input_string.lower()
    # Remove any spaces
    concatenated_string = lowercased_string.replace(" ", "")
    return concatenated_string

def butter_bandpass(cutoffs, fs, order=4):
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
        #coherence_ENorth_BEast = signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend='linear', return_onesided=False,  nperseg=nperseg)
        #coherence_EEast_BNorth = signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend='linear',return_onesided=False , nperseg=nperseg)
       
        #coherence_ENorth_BEast, coherence_EEast_BNorth = coherence_ENorth_BEast[sorted_frequencies_indicies], coherence_EEast_BNorth[sorted_frequencies_indicies]

        ##TODO implement a B-lag cross correlation
        if user_select['lag'] == True and index_satellite == 1: 
            #TODO add a condition that if index our of range just put np.Nans
            print(np.ceil(lag_data[2])*16, np.floor(lag_data[2]*16).astype(int), index_start, index_end, len(Bresamp[index_satellite]))
            if np.ceil(lag_data[2])*16<= index_start and np.floor(lag_data[2]*16).astype(int)+index_end< len(Bresamp[index_satellite]):
                idex_lag=-np.round(lag_data[2]*16).astype(int)
                index_start_test,index_end_test =index_start + idex_lag, index_end+idex_lag
                indicies= range(index_start_test,index_end_test)
                frequencies_E_0, powerspec_E_0 = signal.welch(efield[index_satellite][range(index_start_test,index_end_test), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
                _, powerspec_B_0 = signal.welch(Bresamp[index_satellite][range(index_start_test,index_end_test), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False,  nperseg=nperseg)

                _, powerspec_E_1 = signal.welch(efield[index_satellite][range(index_start_test,index_end_test), 1], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
                _, powerspec_B_1 = signal.welch(Bresamp[index_satellite][range(index_start_test,index_end_test),1], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
                sorted_frequencies_indicies= np.argsort(frequencies_E_0)
                frequencies_E_0= frequencies_E_0[sorted_frequencies_indicies]

                powerspec_E_0, powerspec_B_0, powerspec_E_1, powerspec_B_1 = powerspec_E_0[sorted_frequencies_indicies],  powerspec_B_0[sorted_frequencies_indicies], powerspec_E_1[sorted_frequencies_indicies], powerspec_B_1[sorted_frequencies_indicies]
                ratio_EB_01 = np.sqrt((powerspec_E_0/ powerspec_B_1))
                ratio_EB_10 = np.sqrt((powerspec_E_1/ powerspec_B_0))
                _ , crossspectral_ENorth_BEast = signal.csd(efield[index_satellite][range(index_start_test,index_end_test), 0], Bresamp[index_satellite][range(index_start_test,index_end_test),1], 16, window="hann", detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
                _ ,   crossspectral_EEast_BNorth = signal.csd(efield[index_satellite][range(index_start_test,index_end_test), 1], Bresamp[index_satellite][range(index_start_test,index_end_test),0], 16, window="hann", detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)


                crossspectral_ENorth_BEast, crossspectral_EEast_BNorth = crossspectral_ENorth_BEast[sorted_frequencies_indicies], crossspectral_EEast_BNorth[sorted_frequencies_indicies]


                phase_ENorth_BEast= np.angle(crossspectral_ENorth_BEast, deg=True)
                phase_EEast_BNorth = np.angle(crossspectral_EEast_BNorth, deg=True)
                _ ,cross_BB=signal.csd(Bresamp[index_satellite-1][range(index_start, index_end),1], Bresamp[index_satellite][range(index_start_test,index_end_test),1], 16,window="hann", return_onesided=False, detrend='linear', nperseg=nperseg )

                phase_BB= np.angle(cross_BB, deg=True)
                cross_BB, phase_BB = cross_BB[sorted_frequencies_indicies], phase_BB[sorted_frequencies_indicies]
                print(frequencies_E_0)
                return np.array([[frequencies_E_0, [None] * nperseg], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[np.mean([index_start_test/16,index_end_test/16], dtype=int)]*nperseg, [None]*nperseg], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [cross_BB, phase_BB]]), indicies #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time
            else: 
                indicies=None
                frequencies_E_0, _ =  signal.welch(efield[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
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
                index_start_test,index_end_test = np.nan, np.nan
                return np.array([[frequencies_E_0, [None] * nperseg], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[None]*nperseg, [None]*nperseg], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [cross_BB, phase_BB]]), indicies #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time


        else:
            frequencies_E_0, powerspec_E_0 = signal.welch(efield[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
            _, powerspec_B_0 = signal.welch(Bresamp[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False,  nperseg=nperseg)

            _, powerspec_E_1 = signal.welch(efield[index_satellite][range(index_start, index_end), 1], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
            _, powerspec_B_1 = signal.welch(Bresamp[index_satellite][range(index_start, index_end),1], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
            sorted_frequencies_indicies= np.argsort(frequencies_E_0)
            frequencies_E_0= frequencies_E_0[sorted_frequencies_indicies]

            powerspec_E_0, powerspec_B_0, powerspec_E_1, powerspec_B_1 = powerspec_E_0[sorted_frequencies_indicies],  powerspec_B_0[sorted_frequencies_indicies], powerspec_E_1[sorted_frequencies_indicies], powerspec_B_1[sorted_frequencies_indicies]
            ratio_EB_01 = np.sqrt((powerspec_E_0/ powerspec_B_1))
            ratio_EB_10 = np.sqrt((powerspec_E_1/ powerspec_B_0))
            _ , crossspectral_ENorth_BEast = signal.csd(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16, window="hann", detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)
            _ ,   crossspectral_EEast_BNorth = signal.csd(efield[index_satellite][range(index_start, index_end), 1], Bresamp[index_satellite][range(index_start, index_end),0], 16, window="hann", detrend='linear', scaling='spectrum', return_onesided=False, nperseg=nperseg)


            crossspectral_ENorth_BEast, crossspectral_EEast_BNorth = crossspectral_ENorth_BEast[sorted_frequencies_indicies], crossspectral_EEast_BNorth[sorted_frequencies_indicies]


            phase_ENorth_BEast= np.angle(crossspectral_ENorth_BEast, deg=True)
            phase_EEast_BNorth = np.angle(crossspectral_EEast_BNorth, deg=True)
            cross_BB= [None] * nperseg
            phase_BB= [None] * nperseg
            indicies=  range(index_start, index_end)
            
        
        return np.array([[frequencies_E_0, [None] * nperseg], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[np.mean([index_start/16,index_end/16], dtype=int)]*nperseg, [None]*nperseg], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [cross_BB, phase_BB]]), indicies #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time

        

    
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
                if user_select['coordinate_system'][0] == "North East Centre": 
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
        print(user_select_selected)
        if user_select["Time_Series"] != None:
            for i in range(len(user_select["satellite_graph"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        else:
            pass
        if user_select["EB_periodogram"] != None:
            for i in range(len(user_select["EB_cross power"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["EB_cross power"] != None:
            for i in range(len(user_select["EB_cross power"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["EB_cross phase"] != None:
            for i in range(len(user_select["EB_cross power"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if user_select["lags_cross"] != None:
            for i in range(len(user_select["lags_cross"])):
                user_select_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        else:
            pass
        print(user_select_selected)
        for key in user_select_selected:
            print(key)
            if key != None:
                axes_length += 1
        print(axes_length, 'axes_length')
        return axes_length

    def Animation(data,sampled_datetimes, B_sinc, B_resample, efield, indicies, frames, time_E):
        """
        Creates an animation of each window
        """
        print(Animation_rows())
        fig_ani, axes_ani = plt.subplots(figsize=(10,15),  nrows=Animation_rows(), constrained_layout=True)
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
                            axes_ani[axes_used + k ].plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 3,index, :], "-o", label=val) #Frequencies E's
                            axes_ani[axes_used + k ].set_ylabel(f" E over B ratio of {val_sat} m/s ")
                            axes_ani[axes_used + k ].legend()
                            axes_ani[axes_used + k ].set_yscale("log")
                            axes_ani[axes_used + k ].set_xlim(user_select["bandpass"][1])
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
                    if i==0:
                        print(data[k, 0, 6,index, :], 'plotted')
                    axes_ani[axes_used + l ].set_xlim(user_select["bandpass"][1])
                return  axes_used + len(user_select["EB_cross phase"])
            
            def BB_power_plot(axes_used):
                for l, val in enumerate(user_select["lags_cross"]):
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
                return  axes_used + len(user_select["lags_cross"])
            


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

                        print((np.real(data[k, 0, 0, 0, :]), val_sat))
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
                        
                        if user_select['coordinate_system'][0] == "North East Centre":  #Make one plot over plotted with eachother
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
                                    axes_ani[k].set_ylabel("B(nT)")
                    axes_ani[k].legend(loc=2)
                    twin_x_axes[k].legend(loc=1)
                    ##axes_ani[k].set_xlim(time_E[indicies[k][i]].min(),time_E[indicies[k][i]].max() )
                    
            


                
                return len(twin_x_axes)
            axes_used=0
            if user_select["Time_Series"] != None:
                axes_used=Time_Series_plot()
                print(axes_used, 'axestimeree')
            if user_select["E_periodogram"] != None:
                axes_used = E_Periodogram_Plot(axes_used) 
                print(axes_used, 'axesE')
            if user_select["B_periodogram"] != None:
                axes_used = B_Periodogram_Plot(axes_used)
                print(axes_used, 'axesB')
            if user_select["EB_periodogram"] != None:
                axes_used = EB_Periodogram_Plot(axes_used)
                print('EB_Periodogram_Plot')
            if user_select["EB_cross power"] != None:
                axes_used  = EB_power_plot(axes_used)
            if user_select["EB_cross phase"] != None:
                axes_used=EB_phase_plot(axes_used)
            if user_select["lags_cross"] !=None:
                axes_used=BB_power_plot(axes_used)
            return
        print('ani start', frames)
        ani = animation.FuncAnimation(fig=fig_ani, func=animate, frames=frames) #What
        print('ani end')
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
            if user_select["Difference"] == True:
                length_for_axis += 1
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
                if user_select['coordinate_system'][0] == "North East Centre":
                    if user_select["heatmap"][i] == "E_North":
                        index1 = 1
                        index2=  0
                        
                    elif user_select["heatmap"][i] == "E_East":
                        index1 = 1
                        index2=  1
                    elif user_select["heatmap"][i] == "B_North":
                        index1 = 2
                        index2=  0
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
                    elif user_select["heatmap"][i] == "ENorth/BEast coherence":
                        index1 = 8
                        index2=  0
                    elif user_select["heatmap"][i] == "EEast/BNorth coherence":
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
                print(non_nan_indices, 'indiciess')
                print(data[k, non_nan_indices, 4,0, 0])
                
                times= datetimes[np.real(np.array(data[k, non_nan_indices, 4,0, 0] , dtype=int))] 
                
                if user_select['lag']  == True and lag_data[3] == user_select['satellite_graph'][k]: ##TODO implement
                    delta=timedelta(seconds=lag_data[2])
                else:
                    delta=timedelta(seconds=0)

                
                if index1==2: #B's
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta ,
                                np.absolute(data[k, 0, 0, 0, bandpass]), np.real(np.array(data[k, :, index1, index2,bandpass]))[:,  non_nan_indices],
                                shading='auto',
                                    norm=colors.LogNorm(), cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif  index1==1: #E's
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta , 
                                    np.absolute(data[k, 0, 0, 0, bandpass]), np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices] , shading='auto', 
                                    norm=colors.LogNorm(),
                                    cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif index1==3: #E/B ratio
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta , 
                                    np.absolute(data[k, 0, 0, 0, bandpass]), np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices] , shading='auto', 
                                    norm=colors.LogNorm(),
                                    cmap='winter' ) #selects average time, frequencies, and then the periodogram 

                elif index1==5: #cross power
                    img = axes[length_for_axis+indicies].pcolormesh(times+delta , 
                            np.absolute(data[k, 0, 0, 0, bandpass]), np.absolute(np.array(data[k, :, index1, index2, bandpass]))[:,  non_nan_indices] , shading='auto', 
                            norm=colors.LogNorm(),
                            cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif index1 ==6: #phase
                    #print(np.real(data[k, 0, 0, 0, :]), np.real(np.array(data[k, 0, index1, index2, :])).T, 'heatmap')
                    img=axes[length_for_axis+indicies].pcolormesh(times+delta , 
                                np.absolute(data[k, 0, 0,0, bandpass]),
                                (np.around(np.real(np.array(data[k, :, index1, index2, bandpass])[:,  non_nan_indices])/5, decimals=0))*5, shading='auto', cmap=plt.get_cmap('cmr.infinity')) #selects average time, frequencies, and then the periodogram 
                elif index1 ==7 and k==1: #B B component only goes once
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
                              
                if index1 ==7 and k==0:
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
            if user_select["Difference"] == True:
                length_for_axis += 1
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
                if user_select['coordinate_system'][0] == "North East Centre":
                    if user_select["singles_graph"][i] == "E_North":
                        index1 = 1
                        index2=  0
                        
                    elif user_select["singles_graph"][i] == "E_East":
                        index1 = 1
                        index2=  1
                    elif user_select["singles_graph"][i] == "B_North":
                        index1 = 2
                        index2=  0
                    elif user_select["singles_graph"][i] == "B_East":
                        index1 = 2
                        index2=  1
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
                print(np.shape(data))
                if user_select['bandpass'][0] == True:
                    bandpass=np.where((np.real(data[k, 0, 0, 0, :]) >= user_select['bandpass'][1][0]) & (np.real(data[k, 0, 0, 0, :]) <= user_select['bandpass'][1][1]))[0]
                else: 
                    bandpass =np.where(np.real((data[k, 0, 0, 0, :]) >= 0))[0]
                if index1 ==7 and index2==1 or index1==6: #phase
                    axes[length_for_axis+ indicies].plot( np.absolute(data[k, 0, 0, 0, bandpass]),np.real(data[k, 0, index1,index2,bandpass]))#frequency, data
                else:
                    print(np.absolute(data[k, 0, 0,0,bandpass]))
                    print(np.absolute(data[k, 0, index1,index2,bandpass]), 'test1234')
                    axes[length_for_axis + indicies].plot(np.absolute(data[k, 0, 0,0,bandpass]), np.absolute(data[k, 0, index1,index2,bandpass]))
                    axes[length_for_axis + indicies].set_yscale("log")
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
        print(len(sampled_datetimes), 1/sampling_rate_seconds *user_select['window_length'])
        length_of_windows=int(len(sampled_datetimes) - 1/sampling_rate_seconds *user_select['window_length'])
        print(length_of_windows)
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
        data = np.zeros((len_satellite, length_of_windows, 8, 2,  nperseg), dtype=np.complex_) #[satelitte,number of windows, type of data, different polarizations,  data for each window ]
    else:
        data = np.zeros((len_satellite, 1, 8, 2,  nperseg), dtype=np.complex_)

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


    

    

    return 

def lag(x,y):
    """
    returns an array that is time-corrected to the first
    
    """
    #if len(x) != len(y):
    #    y=np.delete(y,-1, axis=0)
    print(np.shape(x), np.shape(y))
    correlation = signal.correlate(
    x[:,1], y[:,1], mode="full")
    lags = signal.correlation_lags(len(x[:,1]),
                               len(y[:,1]))
    delta=lags[np.argmax(correlation)]/50
    return correlation, lags,delta


def EBplotsNEC(user_select):
    
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
        if user_select["singles_graph"] !=None:
            length_for_axis+=len(user_select['singles_graph'])
        if user_select["conductivities"] != None:
            length_for_axis += len(user_select["conductivities"])

        if user_select["FAC"] == True:
            length_for_axis += 1

        if user_select["Difference"] == True:
            length_for_axis += 1
        if user_select["Pixel_intensity"] == True:
            try:
                length_for_axis += len(user_select["sky_map_values"])
            except TypeError:
                raise ValueError("Sky Map not Selected")
        print(length_for_axis, 'rows')
        return length_for_axis

    fig, axes = plt.subplots(
        nrows=rows(),
        figsize=(10, 15),
        constrained_layout=True,
    )
    time_range = user_select["time_range"]
    has_E = []  # Sets Which space-craft have a corersponding E field
    # Labels of space-craft interested in
    labels = user_select["satellite_graph"]
    # Measurement names from swarm
    measurements = [
        "B_NEC",
        [["VsatN", "VsatE", "VsatC"], ["Evx", "Evy", "Evz"], "Quality_flags"],
    ]
    measurements_flat = [
        "VsatN",
        "VsatE",
        "VsatC",
        "Evx",
        "Evy",
        "Evz",
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
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        if user_select["bandpass"][0] == True and has_twin==False:
            global axes_twin_E
            axes_twin_E=[ axes[x].twinx() for x in range(len(user_select["graph_E_chosen"])) ] 
        for i in range(len(user_select["graph_E_chosen"])):
            if user_select["coordinate_system"][0] == "North East Centre":

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
                    if lag_data[3] == process_string(label):

                        lag_seconds = int(lag_data[2])
                        lag_fraction = lag_data[2] - lag_seconds
                        lag_nanoseconds = int(lag_fraction * 1e9)
                        # Create timedelta64 objects
                        seconds_delta = np.timedelta64(lag_seconds, 's')
                        nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                        # Add the timedelta to the datetime array
                        arrayx = arrayx + seconds_delta + nanoseconds_delta

            except NameError:
                pass
            axes[i].plot(arrayx, arrayy[:, index], label=label, color=colors[satelliteindex])
            axes[i].set_ylabel(
                r"$E_{{{}}}$ $(mV/m)$".format(user_select["graph_E_chosen"][i])
            )
            axes[i].legend(loc=2)
            axes[i].set_xlim((time_range[0], time_range[1]))
            if user_select['singles_graph'] != None:
                axes[i].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                axes[i].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')

            if user_select["bandpass"][0] == True:
                axes_twin_E[i].plot(arrayx,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                axes_twin_E[i].set_ylabel(
                r"$E_{{{}}}$ bandpassed $(mV/m)$".format(user_select["graph_E_chosen"][i])
            )
                axes_twin_E[i].legend(loc=1)
        return True

    def graphingB(label, arrayx, arrayy, arraybandy, satelliteindex, has_twin):
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        
        if user_select["graph_E_chosen"] != None:
            length_for_axis = len(user_select["graph_E_chosen"])
        else:
            length_for_axis = 0
        
        if user_select["bandpass"][0] == True and has_twin==False:
            global axes_twin_B
            axes_twin_B=[ axes[x +length_for_axis].twinx() for x in range(len(user_select["graph_B_chosen"])) ] 
        for i in range(len(user_select["graph_B_chosen"])):
            if user_select['coordinate_system'][0] == "North East Centre":
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
            if user_select['lag'] == True:
                print(lag_data[3], process_string(label))
                if lag_data[3] == process_string(label):
                    lag_seconds = int(lag_data[2])
                    lag_fraction = lag_data[2] - lag_seconds
                    lag_nanoseconds = int(lag_fraction * 1e9)
                    # Create timedelta64 objects
                    seconds_delta = np.timedelta64(lag_seconds, 's')
                    nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                    # Add the timedelta to the datetime array
                    arrayx = arrayx + seconds_delta + nanoseconds_delta

            axes[i + length_for_axis].plot(arrayx, arrayy[:, index], label=label, color=colors[satelliteindex])
            axes[i + length_for_axis].legend(loc=2)
            axes[i + length_for_axis].set_ylabel(
                r"$B_{{{}}}$".format(user_select["graph_B_chosen"][i]) + " (nT) "
            )
            axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))

            if user_select["bandpass"][0] == True:
                axes_twin_B[i].plot(arrayx,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                axes_twin_B[i].set_ylabel(
                r"$B_{{{}}}$ bandpassed $(nT)$".format(user_select["graph_B_chosen"][i])
            )
                axes_twin_B[i].legend(loc=1)
            if user_select['singles_graph'] != None:
                axes[i+ length_for_axis].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                axes[i+ length_for_axis].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')
        return True

    def graphingFlux(label, arrayx, arrayy, index_band, satelliteindex, has_twin):
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
        try:
            if user_select["bandpass"][0] == True and has_twin==False:
                global axes_twin_PF
                axes_twin_PF=[ axes[x+length_for_axis].twinx() for x in range(len(user_select["graph_PF_chosen"])) ] 
        except TypeError:
            print('singleaxes')
            if user_select["bandpass"][0] == True and has_twin==False:
                axes_twin_PF=axes.twinx() 
        for i in range(len(user_select["graph_PF_chosen"])):
            if user_select['coordinate_system'][0] == "North East Centre":
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
                if lag_data[3] == process_string(label):
                    lag_seconds = int(lag_data[2])
                    lag_fraction = lag_data[2] - lag_seconds
                    lag_nanoseconds = int(lag_fraction * 1e9)
                    # Create timedelta64 objects
                    seconds_delta = np.timedelta64(lag_seconds, 's')
                    nanoseconds_delta = np.timedelta64(lag_nanoseconds, 'ns')

                    # Add the timedelta to the datetime array
                    arrayx = arrayx + seconds_delta + nanoseconds_delta

            try:
                if index_band==0:
                    p1,=axes[i + length_for_axis].plot(arrayx, arrayy[index], label=label, color=colors[satelliteindex])

                    axes[i + length_for_axis].set_ylabel(
                        r"$S_{{{}}}$".format(user_select["graph_PF_chosen"][i])
                        + r" $mW m^{-2}$"
                    )
                    axes[i + length_for_axis].legend(loc=2)
                    axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
                if index_band!=0 and user_select["bandpass"][0] == True:
                    axes_twin_PF[i].plot(arrayx,arrayy[index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                    axes_twin_PF[i].legend(loc=1)
                if user_select['singles_graph'] != None:
                    axes[i+ length_for_axis].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                    axes[i+ length_for_axis].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')


            except TypeError:
                if index_band==0:
                    axes.plot(arrayx, arrayy[index], label="".join([label, "bandpassed"]), color=colors[satelliteindex])
                    axes.set_ylabel(
                        r"$S_{{{}}}$".format(user_select["graph_PF_chosen"][i])
                        + r" $mW m^{-2}$"
                    )
                    axes.legend(loc=2)
                    axes.set_xlim((time_range[0], time_range[1]))
                else:
                    axes_twin_PF.plot(arrayx,arrayy[index], label="".join([label, "bandpassed"]),  color=colors[satelliteindex+3], alpha=0.7)
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

    def graphingDifference(arrayx, time, label):
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

        axes[length_for_axis].plot(time, arrayx[1]-arrayx[0], label='Subtracted')
        max_lim=np.nanmax(np.absolute(arrayx[1]-arrayx[0]))
        axes[length_for_axis].set_ylim(-max_lim,max_lim)
        axes_2=axes[length_for_axis].twinx()
        for i in range(2):
            axes_2.plot(time, arrayx[i], label=label[i], alpha=0.2)
        max_lim=np.nanmax(np.absolute(arrayx))
        print(max_lim, 'maxes')
        axes_2.set_ylim(-max_lim,max_lim)
        axes[length_for_axis].legend()
        if user_select['singles_graph'] != None:
                axes[length_for_axis].axvline(user_select["time_range_single"][0], color='orchid', linestyle='dashed')
                axes[length_for_axis].axvline(user_select["time_range_single"][1], color='orchid', linestyle='dashed')

    def Graphing_skymap(pixel, time, spacecraft):
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
        for i in range(len(pixel)):  # length of platforms basically
            for j in range(len(pixel[0])):  # Length of satellites selected
                for k in range(len(pixel[i][j])):
                    if pixel[i][j][k] == 0:
                        pixel[i][j][k] = np.nan
                print(i + length_for_axis, i, length_for_axis, 'skymap')
                axes[i + length_for_axis].plot(
                    time[i], pixel[i][j], label="".join(["swarm ", spacecraft[j]])
                )
            axes[i + length_for_axis].legend(loc=2)
            axes[i + length_for_axis].set_ylabel("Nearest Pixel intensity")
            axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
            # axes[i + length_for_axis].set_xlim(65,70)

    def Coordinate_change(lattiude, longitude, radius):  # Coordinate change, from Ivan correspondence
        a, b, e2 = 6378137.0, 6356752.3142, 0.00669437999014  # From DRS80
        lat, lon, h = np.deg2rad(lattiude), np.deg2rad(longitude), radius
        v = a / np.sqrt(1 - e2 * np.sin(lat) * np.sin(lat))  # logic
        x = (v + h) * np.cos(lat) * np.cos(lon)
        y = (v + h) * np.cos(lat) * np.sin(lon)
        z = (v * (1 - e2) + h) * np.sin(lat)
        return np.array([x, y, -1*z])

    def requesterarraylogic():
        def E():
            return_data_non_band = []
            return_data_band =[]
            returned_times = []
            satellites_with_E = []
            times_for_flux = []
            twin_axis=False
            lag_bool=False
            B_lag_earlier=[]
            for i in range(len(collectionE)):
                dsE = requester(  # requests data
                    collectionE[i],
                    measurements_flat,
                    False,
                    asynchronous=False,
                    show_progress=False,
                )
                if len(dsE) != 0:  # checks if empty
                    # Checks if space-craft is selected
                    if "".join(("swarm", dsE["Spacecraft"][0].lower())) in labels:
                        
                        print(has_E, dsE["Spacecraft"][0].lower())
                        has_E.append(True)
                        satellites_with_E.append(
                            "".join(("swarm", dsE["Spacecraft"][0].lower()))
                        )
                        velocity = dsE[measurements[1][0]].to_numpy()  # Velocity of satellite in NEC
                        velocity_unit = unit_array(velocity)

                        Electric = dsE[measurements[1][1]].to_numpy()  # Electric field in satellite coordinates
                        ElectricNEC = np.multiply(velocity_unit, Electric)  # transforms satellite coordinates into NEC

                        ElectricNECbandpass=np.zeros(np.shape(ElectricNEC))
                        if user_select["bandpass"][0] == True:
                            for l in range(3):  # moving average of bres for all three components
                                ElectricNECbandpass[:, l] = butter_bandpass_filter(data=ElectricNEC[:, l], cutoffs=user_select["bandpass"][1], fs=16)


                            


                        def B_Logic_For_E(lag_bool):
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

                            dsB = requester(
                                collectionB[i],
                                measurements[0],
                                False,
                                asynchronous=False,
                                show_progress=False,
                            )

                            times_of_b_for_flux = Time_corrections(
                                dsE.index.to_numpy(), dsB.index.to_numpy()
                            )  # Takes times of both E and B and finds the closest values in B to E

                            meanfield=arrangement(dsB.index.to_numpy(),dsB["B_NEC_CHAOS"].to_numpy(),3)[times_of_b_for_flux, :]
                            #TODO if swarm a and c selected
                            global B_lag_earlier 
                            if user_select['lag'] and lag_bool:
                                B_synoptic = arrangement(dsB.index.to_numpy(), dsB["B_NEC"].to_numpy()-dsB["B_NEC_CHAOS"].to_numpy(),3)
                                global lag_data
                                lag_data=list(lag(B_lag_earlier,B_synoptic))
                                lag_data.append("".join(("swarm", dsE["Spacecraft"][0].lower()))) 
                            else:
                                lag_bool=True
                                ##TODO make cleaner
                                B_lag_earlier = arrangement(dsB.index.to_numpy(), dsB["B_NEC"].to_numpy()-dsB["B_NEC_CHAOS"].to_numpy(),3)
                            

                            return times_of_b_for_flux, meanfield,  lag_bool

                        times_for_flux, meanfield,lag_bool = B_Logic_For_E(lag_bool)

                        if user_select['coordinate_system'][0] == "Mean-field aligned":
                            latitude, longitude, radius = dsE['Latitude'].to_numpy(), dsE['Longitude'].to_numpy(),  dsE["Radius"].to_numpy()  # Gets Emphermis data
                            r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                          
                            ElectricData = MFA(ElectricNEC, meanfield, r_nec.T) #puts data into MFA coordinate system
                            ElectricNECbandpass = MFA(ElectricNECbandpass, meanfield, r_nec.T)
                        else:
                            ElectricData = ElectricNEC

                        # Plots electric field time seres
                        print(sum(has_E)-1, 'E length')
                        if user_select["graph_E_chosen"] != None:
                            twin_axis=graphingE(
                                "".join(("Swarm ", dsE["Spacecraft"][0])),
                                dsE.index.to_numpy(),
                                ElectricData,
                                ElectricNECbandpass,
                                sum(has_E)-1, #finds the amounts of trues in has_E then -'s by 1 to make it start at -1, this is so the color is the same as B and flux
                                twin_axis
                            )
                        # For ponyting flux
                        return_data_band.append(ElectricNECbandpass)
                        return_data_non_band.append(ElectricData)
                        returned_times = dsE.index.to_numpy()  # Times for plot

                    else:
                        # Says theres no E component
                        has_E.append(False)
                else:  # Says theres no E component
                    has_E.append(False)

            return return_data_non_band, return_data_band, times_for_flux, returned_times, satellites_with_E

        def B():
            return_data_non_band = []
            return_data_band =[]
            time_array=[]
            twin_axis=False
            blabels=[]
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
                    if user_select['coordinate_system'][0] == "Mean-field aligned":
                        latitude, longitude, radius = ds['Latitude'].to_numpy(), ds['Longitude'].to_numpy(),  ds["Radius"].to_numpy()  # Gets Emphermis data
                        r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                        Bdata = MFA(bresarranged, bmodelarranged, r_nec.T) #puts data into MFA coordinate system
                        Bdataband =MFA(bresarrangedband, bmodelarranged,r_nec.T)
                    else:
                        Bdata = bresarranged
                        Bdataband = bresarrangedband
                    print(i, "bcall")
                    if user_select["graph_B_chosen"] != None: #plots if selected
                        twin_axis=graphingB(
                            "".join(("Swarm ", dsmodel_res["Spacecraft"][0])),
                            time,
                            Bdata,
                            Bdataband,
                            i,
                            twin_axis
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
                    fac = pd.Series.to_numpy(fac)

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
            efield_nonband, efield_band, times_for_b, time_E, space_craft_with_E = E()
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
                        [space_craft_with_E[i], time_E, flux_individual[2]]
                    )
            return return_data

        def Difference_plots(bfield,btime, blabel):
            # flux is 16Hz, change to 2Hz to match FAC,
            for index in range(len(bfield)):
                print(lag_data[3], blabel, process_string(blabel[index]))
                if lag_data[3] == process_string(blabel[index]):
                    if index==1:
                        indexopposite=0
                    else:
                        indexopposite=1
                    lag_used=np.round(lag_data[2]*50).astype(int)
                    print(lag_used)
                    if lag_used > 0:
                        offset_ts1 = np.concatenate((np.full(lag_used, np.nan), bfield[index][:, 1]))[:len(bfield[indexopposite][:,1])]
                        offset_ts2 = bfield[indexopposite][:,1]
                    elif lag_used < 0:
                        offset_ts1 = bfield[index][:, 1]
                        offset_ts2 = np.concatenate((np.full(-lag_used, np.nan), bfield[indexopposite][:,1]))[:len(bfield[index][:, 1])]
            graphingDifference(
                [offset_ts1,offset_ts2],
                btime[indexopposite],
                blabel
            )

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
                        asi = asilib.asi.themis(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt="interp",
                        )
                        cadence = 3
                    elif asi_array_code.lower() == "rego":
                        asi = asilib.asi.rego(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt="interp",
                        )
                        cadence = 3
                    elif asi_array_code.lower() == "trex_nir":
                        asi = asilib.asi.trex.trex_nir(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt="interp",
                        )
                        cadence = 6
                    elif asi_array_code.lower() == "trex_rgb":
                        cadence = 3
                        asi = asilib.asi.trex.trex_rgb(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            colors="rgb",
                            custom_alt="interp",
                        )

                    return asi, cadence

                asi, cadence = ASI_logic()
                emph, space_craft_label = empharisis_processing(cadence)

                # return cloest_time_array

                for i in range(len(emph)):  # Conjunction logic
                    # Trex is 6 second cadence compared to 3 of rego and themos

                    sat_time = np.array(emph[i][0])  # sets timestamp
                    print(sat_time, 'sat_time')
                    sat_lla = np.array([emph[i][1], emph[i][2], emph[i][3]]).T

                    conjunction_obj = asilib.Conjunction(asi, (sat_time, sat_lla))

                    # Converts altitude to assumed auroral height

                    conjunction_obj.lla_footprint(alt=alt)

                    lat_satellite.append(conjunction_obj.sat["lat"].to_numpy())
                    lon_satellite.append(conjunction_obj.sat["lon"].to_numpy())
                print(lat_satellite)
                print(lon_satellite)
                print('EB, lat, lon')

                lat, lon = asi.skymap["lat"], asi.skymap["lon"]

                lat[np.isnan(lat)] = np.inf
                lon[np.isnan(lon)] = np.inf

                pixel_chosen = np.zeros((len(emph), len(emph[0][0])))
                values = np.zeros((np.shape(lat)))
                indicies_total = np.zeros((len(emph), len(emph[0][0]), 2), dtype=int)
                index_of_image = 0
                pixel_chosen_average = np.zeros((len(emph), len(emph[0][0])))

                def average(start_index,grid, subregion_size=user_select['pixel_average'][0]):
                    """
                    Calculate the average of a subregion in a grid.
                    
                    Parameters:
                    grid (np.ndarray): The input 2D grid of values.
                    start_index (tuple): A tuple (row, col) indicating the starting index of the subregion.
                    subregion_size (tuple): A tuple (height, width) indicating the size of the subregion.
                    
                    Returns:
                    float: The average value of the subregion.
                    """
                    print(start_index)
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
                    print(average_value)
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

        if user_select["Difference"] == True:
            Difference_plots(bfield_band, time_B, blabels)

        ##TODO Implement ratio
        if user_select["E_B_ratio"] == True:
            if user_select["bandpass"][0] == True:
                efield = efield_band
                bfield = bfield_band
            else:
                efield = efield_nonband
                bfield = bfield_non_band
            Graphing_Ratio(
                space_craft_with_E, efield, bfield, time_E, time_B, user_select, fig, axes
            )
        if user_select["sky_map_values"] != None and user_select["Pixel_intensity"] == True:

            pixels, time, space_craft = skymap()
            Graphing_skymap(pixels, time, space_craft)

        else:
            pass

    requesterarraylogic()

    fig.suptitle("Time Versus Auroral Parameters From Swarm Spacecraft", y=1)

    # for i in range(len(axes)):  # final touches
    # mplcyberpunk.make_lines_glow(axes[i])  # adds glow to plots
    # mplcyberpunk.add_gradient_fill(
    # ax=axes[i], alpha_gradientglow=0.8, gradient_start="zero")

    return fig, axes, 