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
from scipy.signal import butter, lfilter, freqz
from numpy.typing import NDArray

def butter_lowpass(cutoff, fs, order=25):
    return butter(order, cutoff, fs=fs, btype="low", analog=False)


def butter_lowpass_filter(data, cutoff, fs, order=25):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def sinc_interpolation(x: NDArray, s: NDArray, u: NDArray) -> NDArray:
    """Whittaker–Shannon or sinc or bandlimited interpolation.
    Args:
        x (NDArray): signal to be interpolated, can be 1D or 2D
        s (NDArray): time points of x (*s* for *samples*)
        u (NDArray): time points of y (*u* for *upsampled*)
    Returns:
        NDArray: interpolated signal at time points *u*
    Reference:
        This code is based on https://gist.github.com/endolith/1297227
        and the comments therein.
    TODO:
        * implement FFT based interpolation for speed up
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

def period_graph(ax, arg, *kwargs):
    return

def time_series_graph(ax,arg,*kwargs):
    return

def Graphing_Ratio(space_craft_with_E, efield, bfield, time_E, time_B, query_dict):
    """
    Gives the ratio of E/B in the spectral domain, either returning an animation or a plot of the conductivies it derived
    """
    def conductivities():

        return
    def Animation():

        return 
    

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

    def moving_average(a):
        # gives 20 second averaging
        print((time_range[1] - time_range[0]).total_seconds())
        n = int(len(a) / ((time_range[1] - time_range[0]).total_seconds() / 20))
        y_padded = np.pad(a, (n // 2, n - 1 - n // 2), mode="edge")
        y_smooth = np.convolve(y_padded, np.ones((n,)) / n, mode="valid")
        return y_smooth

    def rows():  # finds the number of rows needed for the produced figure
        length_for_axis = 0
        if query_dict['graph_E_chosen'] != None:
            length_for_axis += len(query_dict["graph_E_chosen"])

        if query_dict['graph_B_chosen'] != None:
            length_for_axis += len(query_dict["graph_B_chosen"])

        if query_dict['graph_PF_chosen'] != None:
            length_for_axis += len(query_dict["graph_PF_chosen"])

        if query_dict["FAC"] == True:
            length_for_axis += 1

        if query_dict["Difference"] == True:
            length_for_axis += 1

        if query_dict["E_B_ratio"] == True:
            length_for_axis += 1
        if query_dict["Pixel_intensity"] == True:
            try:
                length_for_axis += len(query_dict["sky_map_values"])
            except TypeError:
                raise ValueError("Sky Map not Selected")

        return length_for_axis

    fig, axes = plt.subplots(
        nrows=rows(),
        figsize=(15, 10),
        sharex=False,
        sharey=False,
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

    def graphingE(label, arrayx, arrayy):
        for i in range(len(query_dict["graph_E_chosen"])):
            print(query_dict["coordinate_system"][0] == "North East Centre")
            if query_dict["coordinate_system"][0] == "North East Centre":
                print("test1")
                if query_dict["graph_E_chosen"][i] == "North":
                    print("test2")
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

            axes[i].plot(arrayx, arrayy[:, index], label=label)
            axes[i].set_ylabel(
                r"$E_{{{}}}$ $(mV/m)$".format(query_dict["graph_E_chosen"][i])
            )
            axes[i].legend(loc=2)
            axes[i].set_xlim((time_range[0], time_range[1]))

    def graphingB(label, arrayx, arrayy):
        if query_dict["graph_E_chosen"] != None:
            length_for_axis = len(query_dict["graph_E_chosen"])
        else:
            length_for_axis = 0
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
            axes[i + length_for_axis].plot(arrayx, arrayy[:, index], label=label)
            axes[i + length_for_axis].legend(loc=2)
            axes[i + length_for_axis].set_ylabel(
                r"$B_{{{}}}$".format(query_dict["graph_B_chosen"][i]) + " (nT) "
            )
            axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))

    def graphingFlux(label, arrayx, arrayy):
        length_for_axis = 0
        # because python starts index at 0
        try:
            length_for_axis += len(query_dict["graph_E_chosen"])
        except TypeError:
            pass
        try:
            length_for_axis += len(query_dict["graph_B_chosen"])
        except TypeError:
            pass
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
                axes[i + length_for_axis].plot(arrayx, arrayy[index], label=label)
                print(i,label)
                print(np.shape(arrayy[index]))
                axes[i + length_for_axis].set_ylabel(
                    r"$S_{{{}}}$".format(query_dict["graph_PF_chosen"][i])
                    + r" $mW m^{-2}$"
                )
                axes[i + length_for_axis].legend(loc=2)
                axes[i + length_for_axis].set_xlim((time_range[0], time_range[1]))
            except TypeError:
                print("The hell")
                axes.plot(arrayx, arrayy[index], label=label)
                axes.set_ylabel(
                    r"$S_{{{}}}$".format(query_dict["graph_PF_chosen"][i])
                    + r" $mW m^{-2}$"
                )
                axes.legend(loc=2)
                axes.set_xlim((time_range[0], time_range[1]))

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
            return_data = []
            returned_times = []
            satellites_with_E = []
            times_for_flux = []
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
                        has_E.append(True)
                        satellites_with_E.append(
                            "".join(("swarm", dsE["Spacecraft"][0].lower()))
                        )
                        velocity = dsE[measurements[1][0]].to_numpy()  # Velocity of satellite in NEC
                        velocity_unit = unit_array(velocity)

                        Electric = dsE[measurements[1][1]].to_numpy()  # Electric field in satellite coordinates
                        ElectricNEC = np.multiply(velocity_unit, Electric)  # transforms satellite coordinates into NEC

                        for l in range(3):  # moving average of bres for all three components
                            ElectricNEC[:, l] = ElectricNEC[:, l] - moving_average(
                                ElectricNEC[:, l]
                            )

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
                            print(np.shape(ElectricNEC), np.shape(meanfield), np.shape(r_nec))
                            ElectricData = MFA(ElectricNEC, meanfield, r_nec.T) #puts data into MFA coordinate system
                        else:
                            ElectricData = ElectricNEC

                        # Plots electric field time seres
                        if query_dict["graph_E_chosen"] != None:
                            graphingE(
                                "".join(("Swarm ", dsE["Spacecraft"][0])),
                                dsE.index.to_numpy(),
                                ElectricData,
                            )
                        # For ponyting flux
                        return_data.append(ElectricData)
                        returned_times = dsE.index.to_numpy()  # Times for plot

                    else:
                        # Says theres no E component
                        has_E.append(False)
                else:  # Says theres no E component
                    has_E.append(False)

            return return_data, times_for_flux, returned_times, satellites_with_E

        def B():
            return_data = []
            time_array=[]
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
                    for l in range(3):  # moving average of bres for all three components
                        bresarranged[:, l] = bresarranged[:, l] - moving_average(
                            bresarranged[:, l]
                        )

                    ##TODO Add MFA
                    if query_dict['coordinate_system'][0] == "Mean-field aligned":
                        latitude, longitude, radius = ds['Latitude'].to_numpy(), ds['Longitude'].to_numpy(),  ds["Radius"].to_numpy()  # Gets Emphermis data
                        r_nec = Coordinate_change(latitude, longitude, radius) #Changes coordinates of r to NEC
                        Bdata = MFA(bresarranged, bmodelarranged, r_nec.T) #puts data into MFA coordinate system
                    else:
                        Bdata = bresarranged

                    if query_dict["graph_B_chosen"] != None: #plots if selected
                        graphingB(
                            "".join(("Swarm ", dsmodel_res["Spacecraft"][0])),
                            time,
                            Bdata,
                        )

                    if (query_dict["graph_PF_chosen"] != None and has_E[j] == True):  # only needs to pass data back if we need to look at pyonting flux E cross B
                        return_data.append(Bdata)
                    time_array.append(time)
            else:
                pass
            return return_data, time_array

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
        ):
            pass
        else:
            # E field, indexs of B for time, times for plotting flux
            efield, times_for_b, time_E, space_craft_with_E = E()
        if (
            query_dict["graph_B_chosen"] == None
            and query_dict["graph_PF_chosen"] == None
        ):
            pass
        else:
            bfield,time_B = B()

        if query_dict["FAC"] == True:
            try:
                FAC_data, Low_hertz_dataframe = F(space_craft_with_E)
            except:
                FAC_data, Low = F()
                print("done F")
        print("test1")

        def pontying_flux():  # Take rough estimate by doing N cross E to get C poynting flux
            return_data = []
            
            for i in range(len(efield)):
                eflux = efield[i]
                bflux=np.empty(np.shape(eflux))

                for l in range(3):
                    print(np.shape(bfield[i][:,l]))
                    bflux[:,l] = sinc_interpolation(x=bfield[i][:,l], s=time_B[i], u=time_E)

                print(np.shape(eflux), np.shape(bflux))
                flux = np.cross(eflux * 1e-3, bflux * 1e-9) * 1.256e6*1e3
                flux_individual = np.transpose(flux)
                print(space_craft_with_E, 'craft')
                print(np.shape(flux_individual))
                graphingFlux(space_craft_with_E[i], time_E, flux_individual)
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
                print(query_dict["sky_map_values"])
                asi_array_code = query_dict["sky_map_values"][k][0]
                location_code = query_dict["sky_map_values"][k][1]
                alt = int(query_dict["sky_map_values"][k][2])

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
                        print(np.shape(lat)[0],np.shape(lat)[1])
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
                        print(i ,"thingy")
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
                np.savetxt(f"sattime.csv", sat_time)
                np.savetxt(f"pixel_inten.csv", pixel_chosen)

                pixel_chosen_average_total.append(pixel_chosen)

            return pixel_chosen_average_total, sat_time_each_platform, space_craft_label
        


        if query_dict["graph_PF_chosen"] == None:
            pass
        else:
            flux = pontying_flux()
        try:
            if query_dict["Difference"] == True:
                Difference_plots(flux, FAC_data)
        except TypeError:
            pass

        ##TODO Implement ratio
        if query_dict["E_B_ratio"] == True:
            pass
            graphing_ratio(
                space_craft_with_E, efield, bfield, time_E, time_B, query_dict
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
