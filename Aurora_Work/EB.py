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


def unit_array(array):
    arraysum = np.sum(np.abs(array), axis=1)
    # Normalizes and finds unitary
    array_unit = array / arraysum[:, np.newaxis]  # normalizes
    return array_unit

def Graphing_Ratio(space_craft_with_E, efield, bfield, time_E, time_B, query_dict, fig, axes):
    """
    Gives the ratio of E/B in the spectral domain, either returning an animation or a plot of the conductivies it derived
    """
    def Logic_for_one_step(index, index_satellite,Bresamp):
        """
        Finds periodogram, Conductivity and Alfven speed for each window
        """
        window_length=query_dict["window_length"]
        index_start,index_end= int(index*16*query_dict['sampling_rate']), int(index*16*query_dict['sampling_rate'])+16*window_length
        print(index_start,index_end, index)
        print(len(time_E))
        
        #TODO do this per space-craft and do it in the ratio desired, currently North over East
        frequencies_E_0, powerspec_E_0 = signal.welch(efield[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=8*window_length)
        frequencies_B_0, powerspec_B_0 = signal.welch(Bresamp[index_satellite][range(index_start, index_end), 0], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=8*window_length)

        frequencies_E_1, powerspec_E_1 = signal.welch(efield[index_satellite][range(index_start, index_end), 1], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=8*window_length)
        frequencies_B_1, powerspec_B_1 = signal.welch(Bresamp[index_satellite][range(index_start, index_end),1], 16, window="hann",detrend='linear', scaling='spectrum', return_onesided=False, nperseg=8*window_length)
        sorted_frequencies_indicies= np.argsort(frequencies_E_0)
        frequencies_E_0, frequencies_B_0 = frequencies_E_0[sorted_frequencies_indicies], frequencies_B_0[sorted_frequencies_indicies]

        powerspec_E_0, powerspec_B_0, powerspec_E_1, powerspec_B_1 = powerspec_E_0[sorted_frequencies_indicies],  powerspec_B_0[sorted_frequencies_indicies], powerspec_E_1[sorted_frequencies_indicies], powerspec_B_1[sorted_frequencies_indicies]
        ratio_EB_01 = np.sqrt((powerspec_E_0/ powerspec_B_1))
        ratio_EB_10 = np.sqrt((powerspec_E_1/ powerspec_B_0))
        freq , crossspectral_ENorth_BEast = signal.csd(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16, window="hann", detrend='linear', scaling='spectrum', return_onesided=False, nperseg=8*window_length)
        freq,  crossspectral_EEast_BNorth = signal.csd(efield[index_satellite][range(index_start, index_end), 1], Bresamp[index_satellite][range(index_start, index_end),0], 16, window="hann", detrend='linear', scaling='spectrum', return_onesided=False, nperseg=8*window_length)


        crossspectral_ENorth_BEast, crossspectral_EEast_BNorth = crossspectral_ENorth_BEast[sorted_frequencies_indicies], crossspectral_EEast_BNorth[sorted_frequencies_indicies]


        phase_ENorth_BEast= np.angle(crossspectral_ENorth_BEast, deg=True)
        phase_EEast_BNorth = np.angle(crossspectral_EEast_BNorth, deg=True)
        #if index==0:
            #print(phase_ENorth_BEast, 'calculated')
        ##print(np.absolute(crossspectral_ENorth_BEast[0]),np.sqrt(np.multiply(powerspec_E_0, powerspec_B_1))[0] )
        coherence_ENorth_BEast = signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend='linear', nperseg=8*window_length)
        coherence_EEast_BNorth = signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend='linear', nperseg=8*window_length)
        #print(efield[index_satellite][range(index_start, index_end), 0], 'E')
        #print(Bresamp[index_satellite][range(index_start, index_end),1],'B')
        #print(signal.coherence(efield[index_satellite][range(index_start, index_end), 0], Bresamp[index_satellite][range(index_start, index_end),1], 16,window="hann", detrend='linear', nperseg=8), 'scipy')
        ##print(np.shape(crossspectral_ENorth_BEast),np.shape(np.angle(crossspectral_ENorth_BEast)), np.shape(powerspec_B_1), np.shape(coherence_EEast_BNorth))

        welch_phase_north_east=np.rad2deg(np.angle(powerspec_E_0) - np.angle(powerspec_B_1))
        welch_phase_east_north=np.angle(powerspec_E_1, deg=True) - np.angle(powerspec_B_0, deg=True)
        return np.array([[frequencies_E_0, [None] * len(frequencies_E_0)], [np.sqrt(powerspec_E_0), np.sqrt(powerspec_E_1)], [np.sqrt(powerspec_B_0), np.sqrt(powerspec_B_1)], [ratio_EB_01, ratio_EB_10], [[np.mean([index_start/16,index_end/16], dtype=int)]*len(frequencies_B_0), [None]*len(frequencies_B_0)], [np.absolute(crossspectral_ENorth_BEast), np.absolute(crossspectral_EEast_BNorth)], [phase_ENorth_BEast, phase_EEast_BNorth], [welch_phase_north_east, welch_phase_east_north]]), range(index_start, index_end)  #all frequencies are the same ##TODO times need to be same length for numpy to work, create array of average time

        

    
    def conductivities(data,times):

        """
        Gives the plot of conductivities and Alfven speed found by looping through window
        """
        def subplot_select():
            #TODO clean up with for loop and apply to other
            length_for_axis = 0
            try:
                length_for_axis += len(query_dict["graph_E_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(query_dict["graph_B_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(query_dict["graph_PF_chosen"])
            except TypeError:
                pass
            if query_dict["FAC"] == True:
                length_for_axis += 1
            else:
                pass
            if query_dict["Difference"] == True:
                length_for_axis += 1
            else:
                pass
            if query_dict["Pixel_intensity"] == True:
                try:
                    length_for_axis += len(query_dict["sky_map_values"])
                except TypeError:
                    raise ValueError("Sky Map not Selected")
            if query_dict["heatmap"] != None:
                length_for_axis += len(query_dict["heatmap"])*len(query_dict["satellite_graph"])
            return length_for_axis
        

        length_of_axis=subplot_select()
        for index in range(len(query_dict["conductivities"])):
            for k in range(len(query_dict["satellite_graph"])):
                if query_dict['coordinate_system'][0] == "North East Centre": 
                    #From North East or East North, are they self similar who knows?
                    if query_dict["conductivities"][index] == "ENorth/BEast":
                        conductivies= data[k, :,  3, 0, 0] # we want zeroth frequency term of the EB ratio
                    else: #East over North
                        conductivies= data[k, :,  3, 1, 0]

                else:
                    if index == "Azimuth over Polodial":
                        conductivies= data[k, :,  3, 0, 0] # we want zeroth frequency term of the EB ratio
                    else: #East over North
                        conductivies= data[k, :,  3, 1, 0]
                axes[index + length_of_axis].plot(times[np.array(data[k, :, 4,0, 0], dtype=int)],1/(1.256e-6*conductivies), label= "".join([query_dict["satellite_graph"][k], " " ,query_dict["conductivities"][index]]))
                axes [index + length_of_axis].set_ylabel("Conductivties (S)")
                axes[ index + length_of_axis].set_title("Conductivity versus time")
                axes[ index + length_of_axis].legend()

                

        


        return
    
    def Animation_rows():
        #TODO Creates an animation of times series of deviations E and B, plot periodograms, plot E/B ratio for given, and what polarization or coordinate should be graphed for each plot
        axes_length=0
        query_dict_selected = [ query_dict["E_periodogram"], query_dict["B_periodogram"]]
        print(query_dict_selected)
        if query_dict["Time_Series"] != None:
            for i in range(len(query_dict["satellite_graph"])):
                query_dict_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        else:
            pass
        if query_dict["EB_periodogram"] != None:
            for i in range(len(query_dict["EB_cross power"])):
                query_dict_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if query_dict["EB_cross power"] != None:
            for i in range(len(query_dict["EB_cross power"])):
                query_dict_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        if query_dict["EB_cross phase"] != None:
            for i in range(len(query_dict["EB_cross power"])):
                query_dict_selected.append([True]) #need a graph for each satellite, creates an array of length satellite that will register for the following for loop
        else:
            pass
        print(query_dict_selected)
        for key in query_dict_selected:
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
        if query_dict["Time_Series"] != None:
            twin_x_axes=[axes_ani[k].twinx() for k in range(len(query_dict["satellite_graph"]))]
        def animate(i):
            """
            Wrapper for animation
            """

            for ax in axs.reshape(-1):
                ax.clear()
            def EB_Periodogram_Plot(axes_used):
                try:
                    for l, val in enumerate(query_dict["EB_periodogram"]):
                        for k, val_sat in enumerate((query_dict["satellite_graph"])):
                            if val == 'ENorth/BEast':
                                index=0
                            else:
                                index=1
                            axes_ani[axes_used + k ].plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 3,index, :], "-o", label=val) #Frequencies E's
                            axes_ani[axes_used + k ].set_ylabel(f" E over B ratio of {val_sat} m/s ")
                            axes_ani[axes_used + k ].legend()
                            axes_ani[axes_used + k ].set_yscale("log")
                            axes_ani[axes_used + k ].set_xlim(query_dict["bandpass"][1])
                except TypeError:
                    for l, val in enumerate(query_dict["EB_periodogram"]):
                        for k, val_sat in enumerate((query_dict["satellite_graph"])):
                            if val == 'ENorth/BEast':
                                index=0
                            else:
                                index=1
                            axes_ani.plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 3,index, :], "-o", label=val_sat) #Frequencies E's
                        axes_ani.set_ylabel(str(val) + " m/s ")
                        axes_ani.legend()
                        axes_ani.set_yscale("log")
                return axes_used + len(query_dict["EB_periodogram"])


            def EB_power_plot(axes_used):
                for l, val in enumerate(query_dict["EB_cross power"]):
                    for k, val_sat in enumerate((query_dict["satellite_graph"])):
                        if val == 'ENorth/BEast crosspower':
                            index=0
                        else:
                            index=1
                        #print(data[k, i, 5,index, :], 'power')
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 5,index, :], "-o",label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " nT ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_yscale("log")
                    axes_ani[axes_used + l ].set_xlim(query_dict["bandpass"][1])
                return  axes_used + len(query_dict["EB_cross power"])
            
            def EB_phase_plot(axes_used):
                for l, val in enumerate(query_dict["EB_cross phase"]):
                    for k, val_sat in enumerate((query_dict["satellite_graph"])):
                        if val == 'ENorth/BEast cross phase':
                            index=0
                        else:
                            index=1
                        positive_frequencies=np.where(data[k, i, 0, 0, :] >= 0)[0]
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, positive_frequencies ]), data[k, i, 6,index, positive_frequencies],"-o", label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " nT ")
                    axes_ani[axes_used + l ].legend()
                    if i==0:
                        print(data[k, 0, 6,index, :], 'plotted')
                    axes_ani[axes_used + l ].set_xlim(query_dict["bandpass"][1])
                return  axes_used + len(query_dict["EB_cross phase"])

            def B_Periodogram_Plot(axes_used):
                for l, val in enumerate(query_dict["B_periodogram"]):
                    for k, val_sat in enumerate((query_dict["satellite_graph"])):
                        if val == 'B_North':
                            index=0
                        else:
                            index=1
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 2,index, :]*1e9, label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " nT ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_yscale("log")
                    axes_ani[axes_used + l ].set_xlim(query_dict["bandpass"][1])
                return  axes_used + len(query_dict["B_periodogram"])

            def E_Periodogram_Plot(axes_used):
                for l, val in enumerate(query_dict["E_periodogram"]):
                    for k, val_sat in enumerate((query_dict["satellite_graph"])):
                        if val == 'E_North':
                            index=0
                        else:
                            index=1
                        axes_ani[axes_used + l ].plot(np.absolute(data[k, i, 0, 0, : ]), data[k, i, 1,index, :]*1e3, label=val_sat) #Frequencies E's
                    axes_ani[axes_used + l ].set_ylabel(str(val) + " mV/m ")
                    axes_ani[axes_used + l ].legend()
                    axes_ani[axes_used + l ].set_yscale("log")
                    axes_ani[axes_used + l ].set_xlim(query_dict["bandpass"][1])
                return  axes_used + len(query_dict["E_periodogram"])
            

            def Time_Series_plot():
                colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
                #There are n number of axises that are contained for the time series
                for k in range(len(query_dict["satellite_graph"])):
                    axes_ani[k].set_title("".join(["Swarm ", query_dict["satellite_graph"][k]]))
                    twin_x_axes[k].clear()
                    twin_x_axes[k].yaxis.set_label_position("right")
                    for l in range(len(query_dict["Time_Series"])):
                        if query_dict['coordinate_system'][0] == "North East Centre":  #Make one plot over plotted with eachother
                            if query_dict["Time_Series"][l] == 'E_North': #Make Electric field twin x axis.
                                twin_x_axes[k].plot(time_E[indicies[k][i]], efield[k][indicies[k][i], 0], color=colors[k+3], label='E North')
                                twin_x_axes[k].set_ylabel("E (mV/m)")

                            elif query_dict["Time_Series"][l] == "E_East":
                                
                                twin_x_axes[k].plot(time_E[indicies[k][i]], efield[k][indicies[k][i], 1], color=colors[k+4], label='E East')
                                twin_x_axes[k].set_ylabel("E (mV/m)")
                                
                            elif query_dict["Time_Series"][l] == "B_East":
                                axes_ani[k].plot(time_E[indicies[k][i]] , B_sinc[k][indicies[k][i], 1], color=colors[k], label= 'B East')
                                axes_ani[k].set_ylabel("B (nT)")
                                pass
                                
                            elif query_dict["Time_Series"][l] == "B_North":
                                axes_ani[k].plot(time_E[indicies[k][i]] , B_sinc[k][indicies[k][i], 0], color=colors[k+1] , label='B North') 
                                axes_ani[k].set_ylabel("B(nT)")
                    axes_ani[k].legend(loc=2)
                    twin_x_axes[k].legend(loc=1)
                    
            


                
                return len(twin_x_axes)
            axes_used=0
            if query_dict["Time_Series"] != None:
                axes_used=Time_Series_plot()
                print(axes_used, 'axestime')
            if query_dict["E_periodogram"] != None:
                axes_used = E_Periodogram_Plot(axes_used) 
                print(axes_used, 'axesE')
            if query_dict["B_periodogram"] != None:
                axes_used = B_Periodogram_Plot(axes_used)
                print(axes_used, 'axesB')
            if query_dict["EB_periodogram"] != None:
                axes_used = EB_Periodogram_Plot(axes_used)
                print('EB_Periodogram_Plot')
            if query_dict["EB_cross power"] != None:
                axes_used  = EB_power_plot(axes_used)
            if query_dict["EB_cross phase"] != None:
                EB_phase_plot(axes_used)
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
                length_for_axis += len(query_dict["graph_E_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(query_dict["graph_B_chosen"])
            except TypeError:
                pass
            try:
                length_for_axis += len(query_dict["graph_PF_chosen"])
            except TypeError:
                pass
            if query_dict["FAC"] == True:
                length_for_axis += 1
            else:
                pass
            if query_dict["Difference"] == True:
                length_for_axis += 1
            else:
                pass
            if query_dict["Pixel_intensity"] == True:
                try:
                    length_for_axis += len(query_dict["sky_map_values"])
                except TypeError:
                    raise ValueError("Sky Map not Selected")
            return length_for_axis
        length_for_axis=subplot_select()
            

        for i in range(len(query_dict["heatmap"])):
            for k in range(len(query_dict["satellite_graph"])):
                if query_dict['coordinate_system'][0] == "North East Centre":
                    if query_dict["heatmap"][i] == "E_North":
                        index1 = 1
                        index2=  0
                    elif query_dict["heatmap"][i] == "E_East":
                        index1 = 1
                        index2=  1
                    elif query_dict["heatmap"][i] == "B_North":
                        index1 = 2
                        index2=  0
                    elif query_dict["heatmap"][i] == "B_East":
                        index1 = 2
                        index2=  1
                    elif query_dict["heatmap"][i] == "ENorth/BEast ratio":
                        index1 = 3
                        index2=  0
                    elif query_dict["heatmap"][i] == "EEast/BNorth ratio":
                        index1 = 3
                        index2=  1
                    elif query_dict["heatmap"][i] == "ENorth/BEast crosspower":
                        index1 = 5
                        index2=  0
                    elif query_dict["heatmap"][i] == "EEast/BNorth crosspower":
                        index1 = 5
                        index2=  1
                    elif query_dict["heatmap"][i] == "ENorth/BEast cross phase":
                        index1 = 6
                        index2=  0
                    elif query_dict["heatmap"][i] == "EEast/BNorth cross phase":
                        index1 = 6
                        index2=  1
                    elif query_dict["heatmap"][i] == "ENorth/BEast coherence":
                        index1 = 7
                        index2=  0
                    elif query_dict["heatmap"][i] == "EEast/BNorth coherence":
                        index1 = 7
                        index2=  1
                else:
                    if query_dict["heatmap"][i] == "E_Azimuth":
                        index1 = 1
                        index2=  0
                    elif query_dict["heatmap"][i] == "E_Polodial":
                        index1 = 1
                        index2=  1
                    elif query_dict["heatmap"][i] == "B_Azimuth":
                        index1 = 2
                        index2=  0
                    elif query_dict["heatmap"][i] == "B_Polodial":
                        index1 = 2
                        index2=  1
                    elif query_dict["heatmap"][i] == "E Azimuth / B Polodial":
                        index1 = 3
                        index2=  0
                    elif query_dict["heatmap"][i] == "E Polodial / B Azimuth":
                        index1 = 3
                        index2=  1
                if index1==2:
                    img = axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].pcolormesh(datetimes[np.array(data[k, :, 4,0, 0], dtype=int)] ,
                                np.absolute(data[k, 0, 0, 0, :]), np.real(np.array(data[k, :, index1, index2, :])).T ,
                                shading='auto',
                                    norm=colors.LogNorm(vmin=np.real(np.array(data[k, :, index1, index2, :])).T.min(),
                                                         vmax=np.real(np.array(data[k, :, index1, index2, :])).T.max()), cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif index1==3: #Change scaling
                    img = axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].pcolormesh(datetimes[np.array(data[k, :, 4,0, 0], dtype=int)] , 
                                    np.absolute(data[k, 0, 0, 0, :]), np.absolute(np.array(data[k, :, index1, index2, :])).T , shading='auto', 
                                    norm=colors.LogNorm(vmin=np.absolute(np.array(data[k, :, index1, index2, :])).T.min(), vmax=np.absolute(np.array(data[k, :, index1, index2, :])).T.max()*1e-1),
                                    cmap='winter_r' ) #selects average time, frequencies, and then the periodogram 

                elif index1==5:
                    img = axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].pcolormesh(datetimes[np.array(data[k, :, 4,0, 0], dtype=int)] , 
                            np.absolute(data[k, 0, 0, 0, :]), np.absolute(np.array(data[k, :, index1, index2, :])).T , shading='auto', 
                            norm=colors.LogNorm(vmin=np.absolute(np.array(data[k, :, index1, index2, :])).T.min(), vmax=np.absolute(np.array(data[k, :, index1, index2, :])).T.max()),
                            cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                elif index1 ==6: #phase
                    #print(np.real(data[k, 0, 0, 0, :]), np.real(np.array(data[k, 0, index1, index2, :])).T, 'heatmap')
                    positive_frequencies=np.where(data[k, 0, 0, 0, :] >= 0)[0]
                    img=axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].pcolormesh(datetimes[np.array(data[k, :, 4,0, 0], dtype=int)] , 
                                np.absolute(data[k, 0, 0,0, positive_frequencies]),
                                (np.around(np.real(np.array(data[k, :, index1, index2, positive_frequencies]))/5, decimals=0))*5, shading='auto', cmap=plt.get_cmap('cmr.infinity')) #selects average time, frequencies, and then the periodogram 
                elif index1 ==7: #not complex, allow non magnitude
                   
                    img = axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].pcolormesh(datetimes[np.array(data[k, :, 4,0, 0], dtype=int)] , 
                            np.absolute(data[k, 0, 0, 0, :]), np.real(np.array(data[k, :, index1, index2, :]).T) , shading='auto', cmap='winter' ) #selects average time, frequencies, and then the periodogram #selects average time, frequencies, and then the periodogram 
                else: #uses absolute value
                     #print(np.real(np.array(data[k, :, index1, index2, :])).T.min(), np.real(np.array(data[k, :, index1, index2, :])).T.max()) #selects average time, frequencies, and then the periodogram )
                     img = axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].pcolormesh(datetimes[np.array(data[k, :, 4,0, 0], dtype=int)] , np.absolute(data[k, 0, 0, 0, :]), np.real(np.array(data[k, :, index1, index2, :])).T , shading='auto', norm=colors.LogNorm(vmin=np.real(np.array(data[k, :, index1, index2, :])).T.min(), vmax=np.real(np.array(data[k, :, index1, index2, :])).T.max()), cmap='winter' ) #selects average time, frequencies, and then the periodogram 
                #print(index1)
                fig.colorbar(img, ax=axes[len(query_dict["satellite_graph"])*i+k+length_for_axis], extend='max', label=query_dict["heatmap"][i])
                axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].set_ylabel(
                    "Frequencies"
                )
                axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].set_xlabel(
                    "Times"
                )
                axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].set_title(
                   query_dict['satellite_graph'][k]
                )
                axes[len(query_dict["satellite_graph"])*i+k + length_for_axis].set_ylim(query_dict["bandpass"][1])
        
        return
    
    time_range = query_dict["time_range"]
    sampling_rate_seconds=query_dict["sampling_rate"]
    sampled_datetimes = create_sampled_datetimes(time_range, sampling_rate_seconds)
    len_satellite = len(query_dict["satellite_graph"])
    print(len(sampled_datetimes), 1/sampling_rate_seconds *query_dict['window_length'])
    length_of_windows=int(len(sampled_datetimes) - 1/sampling_rate_seconds *query_dict['window_length'])
    print(length_of_windows)

    efield=np.multiply(efield,1e-3)
    data = np.zeros((len_satellite, length_of_windows, 8, 2,  (8*query_dict["window_length"])), dtype=np.complex_) #[satelitte,number of windows, type of data, different polarizations,  data for each window ]
    B_sinc,B_resample = np.zeros(np.shape(efield)),np.zeros(np.shape(efield))
    indicies_total = []
    for k in range(len(query_dict["satellite_graph"])): #length of satellites
        for i in range(3):
            B_sinc[k][:, i] =sinc_interpolation(bfield[k][:, i]*1e-9, time_B[k],time_E)
            B_resample[k][:,i]= signal.resample(bfield[k][:, i]*1e-9, len(time_E)) #As shown in testing, use sinc in time time domain, resample in spectral domain. Need to resample from 50, to 16 Hz for periodograms
        indicies_window=[]
        for i in range(length_of_windows): #Loops through each window and at the end stops early so window length doesnt cause error
            data[k, i], indicies = Logic_for_one_step(i, k,  B_resample)
            indicies_window.append(indicies) #don't
        indicies_total.append(indicies_window)
    if query_dict["heatmap"] != None:
        graph_heatmap(data, sampled_datetimes)
    if query_dict["animation"] != False:
        Animation(data,sampled_datetimes, B_sinc, B_resample, efield, indicies_total, length_of_windows, time_E)
    if query_dict["conductivities"] != None:
        conductivities(data, sampled_datetimes)

    

    return 


def EBplotsNEC(query_dict):
    
    set_token(
        "https://vires.services/ows",
        set_default=True,
        token="kmxv5mTTyYwzw4kQ9lsCkGfQHtjjJRVZ",
    )  # key
    plt.style.use("cyberpunk")  # Dark mode!
    print(query_dict)
    print(query_dict.keys())
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
            altitude = data_stuff[i]["Radius"].to_numpy() / 1000 - 6378.1370
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
        if query_dict['graph_E_chosen'] != None:
            length_for_axis += len(query_dict["graph_E_chosen"])

        if query_dict['graph_B_chosen'] != None:
            length_for_axis += len(query_dict["graph_B_chosen"])

        if query_dict['graph_PF_chosen'] != None:
            length_for_axis += len(query_dict["graph_PF_chosen"])
        if query_dict["heatmap"] != None:
            length_for_axis += len(query_dict["heatmap"])*len(query_dict["satellite_graph"])
        if query_dict["conductivities"] != None:
            length_for_axis += len(query_dict["conductivities"])

        if query_dict["FAC"] == True:
            length_for_axis += 1

        if query_dict["Difference"] == True:
            length_for_axis += 1
        if query_dict["Pixel_intensity"] == True:
            try:
                length_for_axis += len(query_dict["sky_map_values"])
            except TypeError:
                raise ValueError("Sky Map not Selected")

        return length_for_axis

    fig, axes = plt.subplots(
        nrows=rows(),
        figsize=(10, 15),
        constrained_layout=True,
    )
    time_range = query_dict["time_range"]
    has_E = []  # Sets Which space-craft have a corersponding E field
    # Labels of space-craft interested in
    labels = query_dict["satellite_graph"]
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
        if query_dict["B_frequency"][0] == "1Hz":
            collectionB = collectionB_01
        else:
            collectionB = collectionB_50
    except TypeError:  # B is none
        collectionB = collectionB_50
    try:
        if query_dict["E_frequency"][0] == "2Hz":
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
        if query_dict["bandpass"][0] == True and has_twin==False:
            global axes_twin_E
            axes_twin_E=[ axes[x].twinx() for x in range(len(query_dict["graph_E_chosen"])) ] 
        for i in range(len(query_dict["graph_E_chosen"])):
            if query_dict["coordinate_system"][0] == "North East Centre":

                if query_dict["graph_E_chosen"][i] == "North":

                    index = 0
                elif query_dict["graph_E_chosen"][i] == "East":
                    index = 1
                elif query_dict["graph_E_chosen"][i] == "Centre":
                    index = 2
            else:
                if query_dict["graph_E_chosen"][i] == "Polodial": 
                    index = 0
                elif query_dict["graph_E_chosen"][i] == "Azimuthal":
                    index = 1
                elif query_dict["graph_E_chosen"][i] == "Mean-field":
                    index = 2
            axes[i].plot(arrayx, arrayy[:, index], label=label, color=colors[satelliteindex])
            axes[i].set_ylabel(
                r"$E_{{{}}}$ $(mV/m)$".format(query_dict["graph_E_chosen"][i])
            )
            axes[i].legend(loc=2)
            axes[i].set_xlim((time_range[0], time_range[1]))

            if query_dict["bandpass"][0] == True:
                axes_twin_E[i].plot(arrayx,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                axes_twin_E[i].set_ylabel(
                r"$E_{{{}}}$ bandpassed $(mV/m)$".format(query_dict["graph_E_chosen"][i])
            )
                axes_twin_E[i].legend(loc=1)
        return True

    def graphingB(label, arrayx, arrayy, arraybandy, satelliteindex, has_twin):
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        
        if query_dict["graph_E_chosen"] != None:
            length_for_axis = len(query_dict["graph_E_chosen"])
        else:
            length_for_axis = 0
        
        if query_dict["bandpass"][0] == True and has_twin==False:
            global axes_twin_B
            axes_twin_B=[ axes[x +length_for_axis].twinx() for x in range(len(query_dict["graph_B_chosen"])) ] 
        for i in range(len(query_dict["graph_B_chosen"])):
            if query_dict['coordinate_system'][0] == "North East Centre":
                if query_dict["graph_B_chosen"][i] == "North":
                    index = 0
                elif query_dict["graph_B_chosen"][i] == "East":
                    index = 1
                elif query_dict["graph_B_chosen"][i] == "Centre":
                    index = 2
            else:
                if query_dict["graph_B_chosen"][i] == "Polodial":
                    index = 0
                elif query_dict["graph_B_chosen"][i] == "Azimuthal":
                    index = 1
                elif query_dict["graph_B_chosen"][i] == "Mean-field":
                    index = 2
            axes[i + length_for_axis].plot(arrayx, arrayy[:, index], label=label, color=colors[satelliteindex])
            axes[i + length_for_axis].legend(loc=2)
            axes[i + length_for_axis].set_ylabel(
                r"$B_{{{}}}$".format(query_dict["graph_B_chosen"][i]) + " (nT) "
            )
            axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))

            if query_dict["bandpass"][0] == True:
                axes_twin_B[i].plot(arrayx,arraybandy[:, index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                axes_twin_B[i].set_ylabel(
                r"$B_{{{}}}$ bandpassed $(nT)$".format(query_dict["graph_B_chosen"][i])
            )
                axes_twin_B[i].legend(loc=1)
        return True

    def graphingFlux(label, arrayx, arrayy, index_band, satelliteindex, has_twin):
        length_for_axis = 0
        # because python starts index at 0
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
        try:
            length_for_axis += len(query_dict["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_B_chosen"])
        except TypeError:
            pass
        try:
            if query_dict["bandpass"][0] == True and has_twin==False:
                global axes_twin_PF
                axes_twin_PF=[ axes[x+length_for_axis].twinx() for x in range(len(query_dict["graph_PF_chosen"])) ] 
        except TypeError:
            print('singleaxes')
            if query_dict["bandpass"][0] == True and has_twin==False:
                axes_twin_PF=axes.twinx() 
        for i in range(len(query_dict["graph_PF_chosen"])):
            if query_dict['coordinate_system'][0] == "North East Centre":
                if query_dict["graph_PF_chosen"][i] == "North":
                    index = 0
                elif query_dict["graph_PF_chosen"][i] == "East":
                    index = 1
                elif query_dict["graph_PF_chosen"][i] == "Centre":
                    index = 2
            else:
                if query_dict["graph_PF_chosen"][i] == "Polodial":
                    index = 0
                elif query_dict["graph_PF_chosen"][i] == "Azimuthal":
                    index = 1
                elif query_dict["graph_PF_chosen"][i] == "Mean-field":
                    index = 2
            try:
                if index_band==0:
                    p1,=axes[i + length_for_axis].plot(arrayx, arrayy[index], label=label, color=colors[satelliteindex])

                    axes[i + length_for_axis].set_ylabel(
                        r"$S_{{{}}}$".format(query_dict["graph_PF_chosen"][i])
                        + r" $mW m^{-2}$"
                    )
                    axes[i + length_for_axis].legend(loc=2)
                    axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
                if index_band!=0 and query_dict["bandpass"][0] == True:
                    axes_twin_PF[i].plot(arrayx,arrayy[index], label="".join([label, "bandpassed"]), color=colors[satelliteindex+3], alpha=0.7)
                    axes_twin_PF[i].legend(loc=1)


            except TypeError:
                if index_band==0:
                    axes.plot(arrayx, arrayy[index], label="".join([label, "bandpassed"]), color=colors[satelliteindex])
                    axes.set_ylabel(
                        r"$S_{{{}}}$".format(query_dict["graph_PF_chosen"][i])
                        + r" $mW m^{-2}$"
                    )
                    axes.legend(loc=2)
                    axes.set_xlim((time_range[0], time_range[1]))
                else:
                    axes_twin_PF.plot(arrayx,arrayy[index], label="".join([label, "bandpassed"]),  color=colors[satelliteindex+3], alpha=0.7)
                    axes_twin_PF.legend(loc=1)
                    axes_twin_PF.set_ylabel(
                    r"$S_{{{}}}$ bandpassed".format(query_dict["graph_PF_chosen"][i])
                    + r" $mW m^{-2}$"
                )
        return True

    def graphingF(label, arrayx, arrayy):
        length_for_axis = 0

        try:
            length_for_axis += len(query_dict["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_B_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_PF_chosen"])
        except TypeError:
            pass

        axes[length_for_axis].plot(arrayx, arrayy, label=label)
        axes[length_for_axis].legend(loc=2)
        axes[length_for_axis].set_ylabel(r"Field Aligned Current $\mu A /m^2$")
        axes[length_for_axis].set_xlim((time_range[0], time_range[1]))

    def graphingDifference(label, arrayx, arrayy, i):
        length_for_axis = 0
        try:
            length_for_axis += len(query_dict["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_B_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_PF_chosen"])
        except TypeError:
            pass
        if query_dict["FAC"] == True:
            length_for_axis += 1
        else:
            pass

        if label == "swarma":
            axes[length_for_axis].plot(arrayx, -1 * arrayy[0], label=label)
        else:
            axes[length_for_axis].plot(arrayx, arrayy[0], label=label)
        axes[length_for_axis].set_ylabel(r"$S_{{{}}}$ (solid)".format("centre"))
        axes[length_for_axis].set_ylim(-1, 1)
        global ax2
        if i == 0:
            ax2 = axes[length_for_axis].twinx()
            ax2.set_ylabel("FAC (DOTTED)")
            ax2.set_ylim(-1, 1)
        ax2.plot(arrayx, arrayy[1], linestyle="dotted")
        axes[length_for_axis].legend(loc=2)
        axes[length_for_axis].set_xlim((time_range[0], time_range[1]))
        # axes[length_for_axis].set_xlim(axes[length_for_axis].get_xlim()[::-1])

    def Graphing_skymap(pixel, time, spacecraft):
        length_for_axis = 0
        try:
            length_for_axis += len(query_dict["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_B_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_PF_chosen"])
        except TypeError:
            pass
        if query_dict["FAC"] == True:
            length_for_axis += 1
        else:
            pass
        if query_dict["Difference"] == True:
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
                        if query_dict["bandpass"][0] == True:
                            for l in range(3):  # moving average of bres for all three components
                                ElectricNECbandpass[:, l] = butter_bandpass_filter(data=ElectricNEC[:, l], cutoffs=query_dict["bandpass"][1], fs=16)


                            


                        def B_Logic_For_E():
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

                            return times_of_b_for_flux, meanfield

                        times_for_flux, meanfield = B_Logic_For_E()

                        if query_dict['coordinate_system'][0] == "Mean-field aligned":
                            latitude, longitude, radius = dsE['Latitude'].to_numpy(), dsE['Longitude'].to_numpy(),  dsE["Radius"].to_numpy()  # Gets Emphermis data
                            r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                          
                            ElectricData = MFA(ElectricNEC, meanfield, r_nec.T) #puts data into MFA coordinate system
                            ElectricNECbandpass = MFA(ElectricNECbandpass, meanfield, r_nec.T)
                        else:
                            ElectricData = ElectricNEC

                        # Plots electric field time seres
                        print(sum(has_E)-1, 'E length')
                        if query_dict["graph_E_chosen"] != None:
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
            # Goes through every satellite selected
            for i in range(len(query_dict["satellite_graph"])):
                if query_dict["satellite_graph"][i] == "epop":
                    raise NotImplementedError("Not implemented")
                else:
                    # Data package level, data/modes, **kwargs
                    if query_dict["satellite_graph"][i] == "swarma":
                        j = 0
                    if query_dict["satellite_graph"][i] == "swarmb":
                        j = 1
                    if query_dict["satellite_graph"][i] == "swarmc":
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
                    if query_dict["bandpass"][0] == True:
                        for l in range(3):  # moving average of bres for all three components
                            bresarrangedband[:, l] = butter_bandpass_filter(bresarranged[:, l], query_dict["bandpass"][1], 50)
                    

                    ##TODO Add MFA
                    if query_dict['coordinate_system'][0] == "Mean-field aligned":
                        latitude, longitude, radius = ds['Latitude'].to_numpy(), ds['Longitude'].to_numpy(),  ds["Radius"].to_numpy()  # Gets Emphermis data
                        r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                        Bdata = MFA(bresarranged, bmodelarranged, r_nec.T) #puts data into MFA coordinate system
                        Bdataband =MFA(bresarrangedband, bmodelarranged,r_nec.T)
                    else:
                        Bdata = bresarranged
                        Bdataband = bresarrangedband
                    print(i, "bcall")
                    if query_dict["graph_B_chosen"] != None: #plots if selected
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
            else:
                pass
            return return_data_non_band,return_data_band, time_array

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
            query_dict["graph_E_chosen"] == None
            and query_dict["graph_PF_chosen"] == None
            and query_dict['E_B_ratio'] == False
        ):
            pass
        else:
            # E field, indexs of B for time, times for plotting flux
            efield_nonband, efield_band, times_for_b, time_E, space_craft_with_E = E()
        if (
            query_dict["graph_B_chosen"] == None
            and query_dict["graph_PF_chosen"] == None
            and query_dict['E_B_ratio'] == False
        ):
            pass
        else:
            bfield_non_band, bfield_band, time_B = B()

        if query_dict["FAC"] == True:
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

        def Difference_plots(flux, fac):
            # flux is 16Hz, change to 2Hz to match FAC,
            closest_time = np.zeros(len(fac[0][1]), dtype=int)
            # resource intesiive, find something better
            for i in range(len(fac[0][1])):
                closest_time[i] = np.argmin(
                    np.abs(fac[0][1][i] - flux[0][1])
                )  # finds the closest value of time_flux for each fac_time
            for i in range(len(flux)):
                flux_time_corrected = flux[i][2][closest_time]
                graphingDifference(
                    flux[i][0],
                    fac[i][1],
                    [
                        flux_time_corrected / np.max(np.abs(flux_time_corrected)),
                        np.array(fac[i][2]) / np.max(np.abs(fac[i][2])),
                    ],
                    i,
                )

        def skymap():
            pixel_chosen_total = []
            sat_time_each_platform = []
            pixel_chosen_average_total = []
            lat_satellite, lon_satellite = [], []
            for k in range(len(query_dict["sky_map_values"])):
                asi_array_code = query_dict["sky_map_values"][k][0]
                location_code = query_dict["sky_map_values"][k][1]
                alt = int(query_dict["sky_map_values"][k][2])

                def ASI_logic():
                    if asi_array_code.lower() == "themis":
                        asi = asilib.asi.themis(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt="geodetic",
                        )
                        cadence = 3
                    elif asi_array_code.lower() == "rego":
                        asi = asilib.asi.rego(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt="geodetic",
                        )
                        cadence = 3
                    elif asi_array_code.lower() == "trex_nir":
                        asi = asilib.asi.trex.trex_nir(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            custom_alt="geodetic",
                        )
                        cadence = 6
                    elif asi_array_code.lower() == "trex_rgb":
                        cadence = 3
                        asi = asilib.asi.trex.trex_rgb(
                            location_code,
                            time_range=time_range,
                            alt=alt,
                            colors="rgb",
                            custom_alt="geodetic",
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

                    # Converts altitude to assumed auroral height

                    conjunction_obj.lla_footprint(alt=alt)
                    lat_sat, lon_sat, ignored = (
                        aacgmv2.convert_latlon_arr(  # Converts to magnetic coordinates
                            in_lat=conjunction_obj.sat["lat"].to_numpy(),
                            in_lon=conjunction_obj.sat["lon"].to_numpy(),
                            height=alt,
                            dtime=datetime(2021, 3, 18, 8, 0),
                            method_code="G2A",
                        )
                    )

                    lat_satellite.append(lat_sat)
                    lon_satellite.append(lon_sat)

                lat, lon = asi.skymap["lat"], asi.skymap["lon"]
                for indi in range(len(lat)):
                    lat[indi], lon[indi], ignored = (
                        aacgmv2.convert_latlon_arr(  # Converts to magnetic coordinates
                            in_lat=lat[indi],
                            in_lon=lon[indi],
                            height=alt,
                            dtime=datetime(2021, 3, 18, 8, 0),
                            method_code="G2A",
                        )
                    )
                lat[np.isnan(lat)] = np.inf
                lon[np.isnan(lon)] = np.inf

                pixel_chosen = np.zeros((len(emph), len(emph[0][0])))
                values = np.zeros((np.shape(lat)))
                indicies_total = np.zeros((len(emph), len(emph[0][0]), 2), dtype=int)
                index_of_image = 0
                pixel_chosen_average = np.zeros((len(emph), len(emph[0][0])))

                def average(index, grid, average_pixels=5):
                    value = 0
                    bounds = (average_pixels - 1) / 2
                    for row in range(
                        average_pixels,
                        int(index[0] - bounds),
                    ):
                        for col in range(
                            average_pixels,
                            int(index[1] - bounds),
                        ):
                            value = value + grid[row][col]
                    return value / average_pixels**2

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

                pixel_chosen_total.append(pixel_chosen)
                sat_time_each_platform.append(sat_time)

                pixel_chosen_average_total.append(pixel_chosen)

            return pixel_chosen_average_total, sat_time_each_platform, space_craft_label
        


        if query_dict["graph_PF_chosen"] == None:
            pass
        else:
            flux = pontying_flux(efield_band,efield_nonband, bfield_band,bfield_non_band)
        try:
            if query_dict["Difference"] == True:
                Difference_plots(flux, FAC_data)
        except TypeError:
            pass

        ##TODO Implement ratio
        if query_dict["E_B_ratio"] == True:
            if query_dict["bandpass"][0] == True:
                efield = efield_band
                bfield = bfield_band
            else:
                efield = efield_nonband
                bfield = bfield_non_band
            Graphing_Ratio(
                space_craft_with_E, efield, bfield, time_E, time_B, query_dict, fig, axes
            )
        if query_dict["sky_map_values"] != None and query_dict["Pixel_intensity"] == True:

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