import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from EB import EBplotsNEC
import matplotlib.animation as animation
from moviepy.editor import VideoFileClip, concatenate_videoclips, clips_array
import moviepy
from datetime import datetime, timedelta, date, time
import asilib
import asilib.asi
import mplcyberpunk
import aacgmv2
from scipy import spatial
from scipy.optimize import minimize
from viresclient import set_token
from viresclient import SwarmRequest
import asilib.skymap

plt.style.use("cyberpunk")

st.title("Date and conjunction time")


# https://stackoverflow.com/questions/47792242/rounding-time-off-to-the-nearest-second-python
def round_seconds(obj: datetime) -> datetime:
    if obj.microsecond >= 500_000:
        obj += timedelta(seconds=1)
    return obj.replace(microsecond=0)


def FixNaNs(arr):
    if len(arr.shape) > 1:
        raise Exception("Only 1D arrays are supported.")
    idxs = np.nonzero(arr == arr)[0]

    if len(idxs) == 0:
        return None

    ret = arr

    for i in range(0, idxs[0]):
        ret[i] = ret[idxs[0]]

    for i in range(idxs[-1] + 1, ret.size):
        ret[i] = ret[idxs[-1]]

    return ret


collectionB_50 = ["SW_OPER_MAGA_HR_1B", "SW_OPER_MAGB_HR_1B", "SW_OPER_MAGC_HR_1B"]
measurements = [
    "B_NEC",
    [["VsatN", "VsatE", "VsatC"], ["Ehx", "Ehy", "Ehz"], "Quality_flags"],
]


def requester(
    sc_collection,
    measurement,
    residual,
    sampling_step=None,
    dict=None,
    cadence=0,
    **kwargs,
):
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
                measurements=measurement, models=["CHAOS"], sampling_step=sampling_step
            )
        data = request.get_between(
            dict["time_range"][0],
            dict["time_range"][1] + timedelta(seconds=cadence),
            **kwargs,
        )
        df = data.as_dataframe()
    except:
        df = []
    return df


def empharisis_processing(dict, cadence):
    emph = []
    space_craft = []
    data_stuff = []
    for i in range(len(collectionB_50)):
        ds = requester(
            collectionB_50[i],
            measurements[0],
            True,
            asynchronous=False,
            show_progress=False,
            sampling_step="PT{}S".format(cadence),
            dict=dict,
            cadence=cadence,
        )
        if "".join(("swarm", ds["Spacecraft"][0].lower())) in dict["satellite_graph"]:
            data_stuff.append(ds)
        else:
            pass
    for i in range(len(data_stuff)):
        time_array = data_stuff[i]["Spacecraft"].index
        # Since time array is one hertz and we want to look for 1/cadence hertz we simply use
        lat_satellite_not_footprint = data_stuff[i]["Latitude"].to_numpy()
        lon_satellite_not_footprint = data_stuff[i]["Longitude"].to_numpy()
        altitude = data_stuff[i]["Radius"].to_numpy() / 1000 - 6378.1370
        # delete_array = np.linspace(0, cadence-1, cadence, dtype=int)
        delete_array = 0
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
    return emph


def Animation(timerange):
    """
    Animation is required embedded within the GUI as this is what cooperates with streamlit, should be its own python file. However it is self-contained within
    """
    def graphing_animation(dict):
        data_3, data_6 = empharisis_processing(dict, 3), empharisis_processing(dict, 6)
        time_range = dict["time_range"]
        platforms = dict["satellite_graph"]
        save_file = []

        # Finds the footprint of selected region with each selected spacecraft

        def animator_map():
            print(movie_generator, "movie_gen")
            for i, (time, image, axes, _) in enumerate(movie_generator):
                # clears time series
                try:
                    # Gets rid of satellite position (scatter plot) and cleans the satellite-track plot
                    for j in range(len(to_remove_scatter)):
                        to_remove_scatter[j].set_visible(False)

                    for j in range(len(to_remove_plots)):
                        to_remove_plots[j].pop()

                except NameError:  # If not initialized, pass
                    pass

                ax.set_title(time)
                ax.set_ylabel("Magnetic latitude")
                ax.set_xlabel("Magnetic Longitude")

                ax.set_xlim(
                    axes.get_xlim()[0], axes.get_xlim()[1]
                )  # sets the limits to be firm, they are the from the image generated by the asi
                ax.set_ylim(axes.get_ylim()[0], axes.get_ylim()[1])

                to_remove_scatter = []
                to_remove_plots = []
                for satellite in range(len(platforms)):  # Loop through satellites
                    to_remove_scatter.append(
                        ax.scatter(# Scatter plot for satellite position from empherasis data
                            lon[ indicies_total[satellite][i][0], indicies_total[satellite][i][1]], 
                            lat[indicies_total[satellite][i][0], indicies_total[satellite][i][1]],
                            s=100,
                            color="red",
                        )
                    )
                    to_remove_plots.append(
                        ax.plot(
                            lon_satellite[satellite],
                            lat_satellite[satellite],
                            color="blue",
                        )
                    )

        def animator_fisheye():

            for i, (time, image, _, im) in enumerate(movie_generator):
                # Plot the entire satellite track, its current location, and a 20x20 km box
                # around its location.
                for j in range(len(sat_azel_pixels_total)):
                    ax.plot(sat_azel_pixels_total[j][:, 0],
                            sat_azel_pixels_total[j][:, 1])
                    ax.scatter(sat_azel_pixels_total[j][i, 0], sat_azel_pixels_total[j][i, 1],
                                marker='o', s=50)
        fig, ax = plt.subplots(  # intializes plots
            1,
            1,
            figsize=(15, 10),
            constrained_layout=True,
        )
        # Loops through the number of stations selected by the user
        for k in range(len(dict["sky_map_values"])):
            asi_array_code = dict["sky_map_values"][k][0]
            location_code = dict["sky_map_values"][k][1]
            alt = int(dict["sky_map_values"][k][2])
            plt.cla()

            def ASI_logic():
                if asi_array_code.lower() == "themis":
                    frame_rate = 2
                    asi = asilib.asi.themis(
                        location_code, time_range=time_range, alt=alt, custom_alt="interp"
                    )
                elif asi_array_code.lower() == "rego":
                    frame_rate = 2
                    asi = asilib.asi.rego(
                        location_code, time_range=time_range, alt=alt, custom_alt="interp"
                    )
                elif asi_array_code.lower() == "trex_nir":
                    frame_rate = 1
                    asi = asilib.asi.trex.trex_nir(
                        location_code, time_range=time_range, alt=alt, custom_alt="interp"
                    )

                elif asi_array_code.lower() == "trex_rgb":
                    frame_rate = 2
                    asi = asilib.asi.trex.trex_rgb(
                        location_code,
                        time_range=time_range,
                        alt=alt,
                        colors="rgb",
                        custom_alt="interp"
                    )
                else:
                    raise NotImplementedError("How did you get this to happen")
                print(dict["sky_map_values"][k][3], "map test")
                if dict["sky_map_values"][k][3] == 'Map':
                    movie_generator = asi.animate_map_gen(  # initaliziation
                        ax=ax, overwrite=True, ffmpeg_params={"framerate": frame_rate}
                    )
                else:
                    movie_generator = asi.animate_fisheye_gen(  # initaliziation
                        ax=ax, overwrite=True, ffmpeg_params={"framerate": frame_rate}
                    )

                return asi, movie_generator

            asi, movie_generator = ASI_logic()

            lat_satellite, lon_satellite, conj_obj_array,sat_azel_pixels_total = [], [], [], []  
            # Creates empty arrays for data from satellites to go, which is there rendered in animator
            for i in range(len(platforms)):  # length of spacecraft REFACTOR
                # Trex is 6 second cadence compared to 3 of rego and themos
                if asi_array_code.lower() == "trex_nir":
                    data = data_6
                else:
                    data = data_3
                sat_time = np.array(data[i][0])  # sets timestamp
                sat_lla = np.array([data[i][1], data[i][2], data[i][3]]).T

                conjunction_obj = asilib.Conjunction(asi, (sat_time, sat_lla))

                # Converts altitude to assumed auroral height

                conjunction_obj.lla_footprint(alt=alt)
                if dict["sky_map_values"][k][3] == 'Map':
                    lat_satellite.append(conjunction_obj.sat["lat"].to_numpy())
                    lon_satellite.append(conjunction_obj.sat["lon"].to_numpy())
                else:
                    conj_obj_array.append(conjunction_obj)
                    sat_azel, sat_azel_pixels = conjunction_obj.map_azel()
                    sat_azel_pixels_total.append(sat_azel_pixels)

            # from asilib
            if dict["sky_map_values"][k][3] == 'Map': #creates a map in lat lon
                pixel_chosen = np.zeros((len(platforms), len(data[0][0])))
                lat, lon = asi.skymap["lat"], asi.skymap["lon"]

                values = np.zeros((np.shape(lat)[0], np.shape(lat)[1]))
                non_blind_search = 20
                indicies_total = np.zeros((len(platforms), len(data[0][0]), 2), dtype=int)

                for satellite in range(len(platforms)):
                    lat_satellite[satellite], lon_satellite[satellite], ignored = (
                        aacgmv2.convert_latlon_arr(  # Converts to magnetic coordinates
                            in_lat=lat_satellite[satellite],
                            in_lon=lon_satellite[satellite],
                            height=alt,
                            dtime=datetime(2021, 3, 18, 8, 0),
                            method_code="G2A",
                        )
                    )
                for i in range(len(lon)):
                    lat[i], lon[i], ignored = (
                        aacgmv2.convert_latlon_arr(  # Converts to magnetic coordinates
                            in_lat=lat[i],
                            in_lon=lon[i],
                            height=alt,
                            dtime=datetime(2021, 3, 18, 8, 0),
                            method_code="G2A",
                        )
                    )
                lat[np.isnan(lat)] = np.inf
                lon[np.isnan(lon)] = np.inf
                for i in range(len(data[0][0])):  # len of time series
                    # Themis starts at time_range[0], rego and trex start a time_range[0] + a cadence
                    if len(asi.data[0]) != len(data[0][0]) and i == 0:
                        continue
                    elif len(asi.data[0]) != len(data[0][0]) and i != 0:
                        i -= 1

                    # blind optimzation, loops through time range, satellites,
                    # and lats and lons of the skymap to find the closest pixel to the satellite which can then be used to
                    # find intensities

                    # if(i == 0):  # intialization
                    for satellite in range(len(platforms)):  # blind search
                        for j in range(np.shape(lat)[0]):
                            for k in range(np.shape(lat)[1]):
                                values[j][k] = np.sqrt(
                                    (lat[j][k] - lat_satellite[satellite][i]) ** 2
                                    + (lon[j][k] - lon_satellite[satellite][i]) ** 2
                                )

                        indicies = np.unravel_index((np.abs(values)).argmin(), values.shape)
                        indicies_total[satellite, i] = indicies
                        print(indicies)
                    try:
                        pixel_chosen[satellite][i] = asi.data[1][i][
                            indicies[0], indicies[1]
                        ]
                    except ValueError:  # rgb
                        pixel_chosen[satellite][i] = asi.data[1][i][
                            indicies[0], indicies[1], 0
                        ]
                    except IndexError:
                        pixel_chosen[satellite][i] = 0
                animator_map()
                movie_container = "mp4"
                movie_address = (
                    f'{time_range[0].strftime("%Y%m%d_%H%M%S")}_'
                    f'{time_range[1].strftime("%H%M%S")}_'
                    f"{asi_array_code.lower()}_{location_code.lower()}_{str(alt)}km_map.{movie_container}"
                )  # file address of movie saved by asilib

                movie_address_total = (
                    asilib.config["ASI_DATA_DIR"] / "animations" / movie_address
                )  # full address from C:

                # Saves address so movie.py can load it in the GUI
                save_file.append(movie_address_total)
            else:
                animator_fisheye()

                movie_container = "mp4"
                movie_address = (
                    f'{time_range[0].strftime("%Y%m%d_%H%M%S")}_'
                    f'{time_range[1].strftime("%H%M%S")}_'
                    f"{asi_array_code.lower()}_{location_code.lower()}_{str(alt)}km_fisheye.{movie_container}"
                )  # file address of movie saved by asilib

                movie_address_total = (
                    asilib.config["ASI_DATA_DIR"] / "animations" / movie_address
                )  # full address from C:

                # Saves address so movie.py can load it in the GUI
                save_file.append(movie_address_total)
        return save_file

    def Animation_GUI():
        st.title("Animation Interface:")
        if "station_count" not in st.session_state:
            st.session_state["station_count"] = 1
        st.number_input(
            label="Number of Stations to animate",
            min_value=1,
            max_value=4,
            value=st.session_state["station_count"],
            key="station_count",
        )

        def station_logic():
            count = st.session_state["station_count"]
            col = st.columns(count)  # number of columns selected
            # numpy doesnt have string arrays, only character arrays
            # initializes empty array for columns selected
            values = [[None, None]] * count

            def station_GUI():
                for i in range(count):
                    # initalizies our site state (project already initialized in line 41)
                    if "".join([str(i), "site"]) not in st.session_state:
                        st.session_state["".join([str(i), "site"])] = []

                    if "".join([str(i), "height"]) not in st.session_state:
                        st.session_state["".join([str(i), "site"])] = []

                    with col[i]:  # each column
                        st.multiselect("Name of Project", ["REGO", "THEMIS", "trex_nir", "trex_rgb"],key="".join([str(i), "project"]),max_selections=1)  # Sets the project
                        # if no project, do not execute site parameters
                        if st.session_state["".join([str(i), "project"])] != []:
                            # if rego project, select site
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "REGO"
                            ):
                                st.multiselect(
                                    "Name of Site",
                                    ["FSMI", "GILL", "RESU", "TALO", "rank", "fsim"],
                                    key="".join([str(i), "site"]),
                                    max_selections=1,
                                )

                            # if project is themis, select themis site
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "THEMIS"
                            ):
                                st.multiselect(
                                    "Name of Site",
                                    [
                                        "FSMI",
                                        "GILL",
                                        "FSIM",
                                        "rank",
                                        "talo",
                                        "atha",
                                        "tpas",
                                        "inuv",
                                    ],
                                    key="".join([str(i), "site"]),
                                    max_selections=1,
                                )
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "trex_nir"
                            ):
                                st.multiselect(
                                    "Name of Site",
                                    ["rabb", "gill"],
                                    key="".join([str(i), "site"]),
                                    max_selections=1,
                                )
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "trex_rgb"
                            ):
                                st.multiselect(
                                    "Name of Site",
                                    ["rabb", "gill", "fsmi", "pina", "yknf"],
                                    key="".join([str(i), "site"]),
                                    max_selections=1,
                                )


                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "REGO"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=230,
                                    step=10,
                                    key="".join([str(i), "height"]),
                                )
                                st.multiselect(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    "Fisheye",
                                    key="".join([str(i), "map"]),
                                    max_selections=1,
                                )
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "THEMIS"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=170,
                                    step=10,
                                    key="".join([str(i), "height"])
                                )
                                st.multiselect(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    "Fisheye",
                                    key="".join([str(i), "map"]),
                                    max_selections=1,
                                )
                                
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "trex_nir"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=230,
                                    step=10,
                                    key="".join([str(i), "height"]),
                                )
                                st.multiselect(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    "Fisheye",
                                    key="".join([str(i), "map"]),
                                    max_selections=1,
                                )
                            if (
                                st.session_state["".join([str(i), "project"])][0]
                                == "trex_rgb"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=170,
                                    step=10,
                                    key="".join([str(i), "height"]),
                                )
                                st.multiselect(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    "Fisheye",
                                    key="".join([str(i), "map"]),
                                    max_selections=1,
                                )


                        else:
                            pass

            def value_setter():  # sets values of selected sites and projects
                nonlocal values
                for i in range(count):  # goes through all columns
                    # doesn't add values if empty
                    if st.session_state["".join([str(i), "project"])] == []:
                        pass

                    else:
                        # doesn't add values if empty
                        if (
                            st.session_state["".join([str(i), "site"])] != []
                            and st.session_state["".join([str(i), "height"])]
                            and st.session_state["".join([str(i), "site"])] != [] #????
                        ):
                            project = st.session_state["".join([str(i), "project"])][0]

                            sites = st.session_state["".join([str(i), "site"])][0]

                            heights = st.session_state["".join([str(i), "height"])]
                            try:
                                maps = st.session_state["".join([str(i), "map"])][0]
                            except IndexError:
                                maps=None

                            values[i] = [project, sites, heights, maps]
                        else:
                            pass

            station_GUI()
            value_setter()

            return values

        global skymap_values
        skymap_values = station_logic()

        Animation_dict = {
            "time_range": timerange,
            "satellite_graph": st.session_state["Satellite_Graph"],
            "sky_map_values": skymap_values,
        }

        def Animation_function_caller():
            def Animate_graph():
                # Gets the figures and axes from cache
                fig, axes = st.session_state["Graph"]
                n = int(
                    (
                        Animation_dict["time_range"][1]
                        - (Animation_dict["time_range"][0] + timedelta(seconds=3))
                    ).total_seconds()
                    / 3
                )  # Finds the number of frames needed
                time = np.array(
                    [
                        (Animation_dict["time_range"][0])
                        + timedelta(seconds=i * 3)  # Creates array for x-axis (time)
                        for i in range(n)
                    ]
                )
                axes_changed = axes

                def Update(i):
                    # Goes through each axes and draw a vertical line at selected time to show where animation is
                    try:
                        for j in range(len(axes)):
                            lines = axes_changed[j].axvline(
                                time[i], linewidth=1, linestyle="dashed", color="red"
                            )
                    except TypeError:
                        lines = axes_changed.axvline(
                            time[i], linewidth=1, linestyle="dashed", color="red"
                        )  # Move to other monitor

                lin_ani = animation.FuncAnimation(
                    fig, Update, frames=n
                )  # Creates animation
                FFwriter = animation.FFMpegWriter(fps=2)  # Writes to mp4
                lin_ani.save("animationgraph.mp4", writer=FFwriter)

            animation_strings = graphing_animation(Animation_dict)

            try:
                clip1 = VideoFileClip(r"{}".format(animation_strings[0]))
            except IndexError:
                clip1 = None
            try:
                clip2 = VideoFileClip(r"{}".format(animation_strings[1]))
            except IndexError:
                clip2 = None
            try:
                clip3 = VideoFileClip(r"{}".format(animation_strings[2]))
            except IndexError:
                clip3 = None
            try:
                clip4 = VideoFileClip(r"{}".format(animation_strings[3]))
            except IndexError:
                clip4 = None

            if clip3 == None and clip2 == None and clip4 == None:
                combined = clips_array([[clip1]])
            elif clip4 == None and clip3 == None:
                combined = clips_array([[clip1, clip2]])
            elif clip4 == None:
                combined = clips_array([[clip1, clip2, clip3]])
            else:
                combined = clips_array([[clip1, clip2], [clip3, clip4]])
            if "Graph" in st.session_state and st.session_state["animation"] != True:
                Animate_graph()
                clip_graph = VideoFileClip("animationgraph.mp4")
                if clip3 == None and clip2 == None and clip4 == None:
                    combined = clips_array([[clip1, clip_graph]])
                elif clip4 == None and clip3 == None:
                    combined = clips_array([[clip1, clip2, clip_graph]])
                elif clip4 == None:
                    combined = clips_array([[clip1, clip2], [clip3, clip_graph]])
                else:
                    combined = clips_array(
                        [[clip1, clip2], [clip3, clip4], [clip_graph, clip_graph]]
                    )
            if "Graph" in st.session_state and st.session_state["animation"] == True:
                Animate_graph()
                clip_graph = VideoFileClip("animationgraph.mp4")
                alfven_graph = VideoFileClip("animationAlfven.mp4")
                if clip3 == None and clip2 == None and clip4 == None:
                    combined = clips_array([[clip1, clip_graph, alfven_graph]])
                elif clip4 == None and clip3 == None:
                    combined = clips_array([[clip1, clip2] [alfven_graph, clip_graph]])
                elif clip4 == None:
                    combined = clips_array([[clip1, clip2], [clip3, clip_graph], [clip_graph, alfven_graph]])
                else:
                    combined = clips_array(
                        [[clip1, clip2], [clip3, clip4], [alfven_graph, clip_graph]]
                    )

                

            combined.write_videofile("animation_display.mp4")
            st.video("animation_display.mp4")
            st.session_state["Animation_logic_completed"] = True

        # calls to make one column but should dynamiically update
        # station_logic(1)
        button_for_animation = st.button(
            label="Render graphs", key="Animation_executer"
        )
        if button_for_animation == True:
            Animation_function_caller()

    Animation_GUI()


def Graph():
    st.title("Graph Interface:")

    def Graph_options_B(coord_options):
        st.multiselect(
            label="What directions of B would you like to graph",
            options=coord_options,
            key="B_options_to_use",
        )
        st.multiselect(
            label="What frequency would you like to use",
            key="Frequency_B",
            options=["1Hz", "50Hz"],
            max_selections=1,
        )

    def Graph_options_E(coord_options):
        st.multiselect(
            label="What directions of E would you like to graph",
            options=coord_options,
            key="E_options_to_use",
        )
        st.multiselect(
            label="What frequency would you like to use",
            key="Frequency_E",
            options=["2Hz", "16Hz"],
            max_selections=1,
        )

    def Graph_options_F(ignored):
        pass

    def Graph_options_PF(coord_options):
        st.multiselect(
            label="What directions of Ponyting Flux would you like to graph",
            options=coord_options,
            key="PF_options_to_use",
        )

    options_for_graphs = ["B", "E", "FAC", "Poynting flux"]
    Graph_functions = [
        Graph_options_B,
        Graph_options_E,
        Graph_options_F,
        Graph_options_PF,
    ]

    def GUI_interface():
        graphs = st.multiselect(
            label="What would you like to graph",
            options=options_for_graphs,
            key="Graph_select",
            default=None,
        )

        coordinate_system = st.multiselect(
            label="What coordinate system would you like it in",
            options=["North East Centre", "Mean-field aligned"],
            max_selections=1,
            key="Coordinate_system",
        )
        return coordinate_system, graphs

    def Drop_down_menus(coordinate_system, graphs):
        try:  # tries to generate column but doesn't work if the index of the coordinate system is 0
            # sets coordinate system to give to the functions in Graph_functions
            if coordinate_system[0] == "North East Centre":
                coord_options = ["North", "East", "Centre"]
            elif coordinate_system[0] == "Mean-field aligned":
                coord_options = ["Mean-field", "Azimuthal", "Polodial"]
            try:  # st.columns doesn't like len(0)
                # creates columns for each variable to graph ie: B, E etc
                col = st.columns(len(graphs))
                for i in range(len(graphs)):  # Initializes columns
                    # 2D for loop to see if selected graph option corresponds to which Graph_function
                    for j in range(len(options_for_graphs)):
                        with col[i]:  # GUI columns
                            # if graphs options is the same index as the Graph_functions, calls graph functions
                            if graphs[i] == options_for_graphs[j]:
                                # calls index in function array, thus calling the corresponding function with *args of the coordinate system
                                Graph_functions[j](coord_options)
            except st.errors.StreamlitAPIException:
                pass
        except IndexError:
            pass

    coordinate_system_selected, graphs_selected = GUI_interface()
    Drop_down_menus(coordinate_system_selected, graphs_selected)
    if "Difference" not in st.session_state:
        st.session_state["Difference"] = False
    if "E_B_ratio" not in st.session_state:
        st.session_state["E_B_ratio"] = False
    if "Pixel_intensity" not in st.session_state:
        st.session_state["Pixel_intensity"] = False
    if "Filtering" not in st.session_state:
        st.session_state["Filtering"] = False
    if "Heatmap" not in st.session_state:
        st.session_state["Heatmap"] = False
    if "Conductivies" not in st.session_state:
        st.session_state["Conductivities"] = False
    if "Alfven_Animation" not in st.session_state:
        st.session_state["Alfven_Animation"] = False


    st.checkbox(
        label="would you like to find the normalized difference between FAC and Pyonting flux (must select both FAC and pyonting flux centre)",
        key="Difference",
        value=st.session_state["Difference"],
    )
    st.checkbox(
        label=r"Would you like to find the ratio of $E/B$",
        value=st.session_state["E_B_ratio"],
        key="E_B_ratio",
    )
    # (need to select from "Animation Auroral)
    st.checkbox(
        label=r"Would you like to graph the pixel intesity",
        value=st.session_state["Pixel_intensity"],
        key="Pixel_intensity",
    )
    st.checkbox(
        label=r"Would you like to bandpass filter for Alfven (50km, 1km)",
        value=st.session_state["Filtering"],
        key="Filtering",
    )

    if "low_pass" not in st.session_state:
        st.session_state["low_pass"] = None
    if "high_pass" not in st.session_state:
        st.session_state["high_pass"] = None

    if "sampling_rate" not in st.session_state:
        st.session_state["sampling_rate"] = None
    
    if "Window_Length" not in st.session_state:
        st.session_state["Window_Length"] = None


    if "heatmap_graphs" not in st.session_state:
        st.session_state["heatmap_graphs"] = None
    if "conductivity_graphs" not in st.session_state:
        st.session_state["conductivity_graphs"] = None
    if "Alfven_Animation" not in st.session_state:
        st.session_state["Alfven_Animation"] = None
        
    if "Time_Series_Graph" not in st.session_state:
        st.session_state["Time_Series_Graph"] = None
    if "E_Peridogram_Graph" not in st.session_state:
        st.session_state["E_Peridogram_Graph"] = None
    if "B_Peridogram_Graph" not in st.session_state:
        st.session_state["B_Peridogram_Graph"] = None
    if "E/B_Periodogram_Graph" not in st.session_state:
        st.session_state["E/B_Periodogram_Graph"] = None

    if st.session_state["Filtering"] == True:
        st.select_slider(label="Please select the low pass in Hz", value=0.2, options=[0.1,0.2,0.5, 1, 2, 4], key="low_pass")
        st.select_slider(label="Please select the low pass in Hz", value=7, options=[0.5, 1, 2, 4, 6, 7, 7.9], key="high_pass")
    
    if st.session_state["E_B_ratio"] == True:
        st.select_slider(label="Please select the running window interval (sampling rate)", value=1, options=[0.1,0.2,0.5, 1, 2], key="sampling_rate")
        st.select_slider(label="Please select the running window length", value=4, options=[2,3,4,5,6,7,7.9], key="Window_Length")


        if st.session_state["Coordinate_system"][0] == "North East Centre":
            options=[["E_North", "E_East"], ["B_North", "B_East"], ["ENorth/BEast", "EEast/BNorth"]]
        else:
            options=[["E_Azimuthal", "E_Polodial"], ["B_Azimuthal", "B_Polodial"], ["EAzimuthal/BPolodial", "EPolodial/BAzimuthal"]]


        st.checkbox(label="Would you like to graph a heatmap of the event through timed (frequency versus time versus ampltiude spectra)", key="Heatmap")
        st.checkbox(label="Would you like to graph the conductivies derived the EB ratio versus time", key="Conductivies")
        st.checkbox(label="Would you like to create an animation of the peridograms through time", key="Alfven_Animation")

        if st.session_state["Heatmap"] == True:
            st.multiselect(label="Choose the variable in the heatmap(s)", options=sum(options, []), key="heatmap_graphs")

        if st.session_state["Conductivies"] == True:
            st.multiselect(label="Choose the polarization used for conductivies", options=options[2], key="conductivity_graphs")

        if st.session_state["Alfven_Animation"] == True:

            st.multiselect(label="Please select the plots in the Alfven animation", options=["Time Series", "E Periodogram", "B Periodogram", "E/B Periodogram"], key="Alfven_graphs")

            if "Time Series" in st.session_state["Alfven_graphs"]:
                print(options[:-1])
                st.multiselect(label="Please select the time series you want to plot through", options=sum(options[:-1], []), key="Time_Series_Graph")

            if "E Periodogram" in st.session_state["Alfven_graphs"]:
                st.multiselect(label="Please select the E polarization you want to plot", options=options[0], key="E_Peridogram_Graph")

            if "B Periodogram" in st.session_state["Alfven_graphs"]:
                st.multiselect(label="Please select the B polarization you want to plot", options=options[1], key="B_Peridogram_Graph")

            if "E/B Periodogram" in st.session_state["Alfven_graphs"]:
                st.multiselect(label="Please select the E/B polarization you would like to plot", options=options[2], key="E/B_Periodogram_Graph")
            
        


        


def Render_Graph(timerange):
    parameters = [ 
        "Satellite_Graph",
        "Satellite_Graph",
        "B_options_to_use",
        "Frequency_B",
        "E_options_to_use",
        "PF_options_to_use",
        "Frequency_E",
    ]
    for i in range(len(parameters)):
        # Initializes all the parameters
        if parameters[i] not in st.session_state:
            st.session_state[parameters[i]] = None
    # Index error because if Graph hasn't been selected, [0] doesn't work as stated in many comments
    try:
        np.reshape([np.where(np.array(st.session_state["Graph_select"]) == "FAC")], -1)[
            0
        ]  # Finds if any index exists named FAC in the graphs selected
        FAC_boolean = True
    except IndexError:
        FAC_boolean = False
    # try:  # [0] doesnt work
    try:
        count = st.session_state["station_count"]
        values = [[None, None]] * count

        def value_setter():  # sets values of selected sites and projects
            nonlocal values
            for i in range(count):  # goes through all columns
                # doesn't add values if empty
                if st.session_state["".join([str(i), "project"])] == []:
                    pass

                else:
                    # doesn't add values if empty
                    if st.session_state["".join([str(i), "site"])] != []:
                        project = st.session_state["".join([str(i), "project"])][0]

                        sites = st.session_state["".join([str(i), "site"])][0]

                        heights = st.session_state["".join([str(i), "height"])]
                        
                        values[i] = [project, sites, heights]
                    else:
                        pass
            return values

        skymap_values = value_setter()
    except KeyError:
        skymap_values = None

    dict = {
        "time_range": timerange,
        "satellite_graph": st.session_state["Satellite_Graph"],
        "coordinate_system": st.session_state["Coordinate_system"],
        "graph_B_chosen": st.session_state["B_options_to_use"],
        "B_frequency": st.session_state["Frequency_B"],
        "E_frequency": st.session_state["Frequency_E"],
        "graph_E_chosen": st.session_state["E_options_to_use"],
        "graph_PF_chosen": st.session_state["PF_options_to_use"],
        "FAC": FAC_boolean,
        "Difference": st.session_state["Difference"],
        "E_B_ratio": st.session_state["E_B_ratio"],
        "Pixel_intensity": st.session_state["Pixel_intensity"],
        "sky_map_values": skymap_values,
        "bandpass": [st.session_state["Filtering"] ,[st.session_state["low_pass"], st.session_state["high_pass"]]],
        "heatmap": st.session_state["heatmap_graphs"],
        "conductivities": st.session_state["conductivity_graphs"],
        "animation": st.session_state["Alfven_Animation"],
        "Time_Series":st.session_state["Time_Series_Graph"],
        "E_periodogram": st.session_state["E_Peridogram_Graph"],
        "B_periodogram": st.session_state["B_Peridogram_Graph"],
        "EB_periodogram": st.session_state["E/B_Periodogram_Graph"],
        "sampling_rate": st.session_state["sampling_rate"],
        "window_length":st.session_state["Window_Length"]

    }

    if dict["coordinate_system"][0] == "North East Centre":
        figaxes = EBplotsNEC(dict)
        st.session_state["Graph"] = figaxes
    print("axes")
    print(dict["coordinate_system"][0])
    if dict["coordinate_system"][0] == "Mean-field aligned":
        figaxes = EBplotsNEC(dict)
        print("called")
        st.session_state["Graph"] = figaxes
    return


def Main():
    if "time_start" not in st.session_state:
        st.session_state["time_start"] = time(8, 0)
    if "time_end" not in st.session_state:
        st.session_state["time_end"] = time(8, 30)
    if "date" not in st.session_state:
        st.session_state["date"] = date(2021, 3, 18)
    if "second_start" not in st.session_state:
        st.session_state["second_start"] = "0"
    if "second_end" not in st.session_state:
        st.session_state["second_end"] = "0"
    Animation_value = st.sidebar.checkbox(
        label="Would you like to display an auroral animation",
        value=False,
        help=None,
        on_change=None,
        args=None,
        kwargs=None,
        disabled=False,
        label_visibility="visible",
    )
    Graph_value = st.sidebar.checkbox(
        label="Would you like to a graph",
        value=False,
        help=None,
        on_change=None,
        args=None,
        kwargs=None,
        disabled=False,
        label_visibility="visible",
    )
    date_range = st.date_input(label="Date of conjunction", key="date")
    time_rangestart = st.time_input(
        label="Time to the start of the conjunction",
        step=60,
        key="time_start",
        value=st.session_state["time_start"],
    )

    time_rangeend = st.time_input(
        label="Time to the end of the conjunction",
        step=60,
        key="time_end",
        value=st.session_state["time_end"],
    )
    seconds_start = st.selectbox(
        "Seconds_start",
        ("0", "6", "8", "10", "24", "35", "36", "42", "48", "55"),
        key="second_start",
    )
    seconds_end = st.selectbox(
        "Seconds_end", ("0", "12", "24", "30", "36", "42", "48", "55"), key="second_end"
    )
    global timerange
    timerange = (
        datetime.combine(date_range, time_rangestart)
        + timedelta(seconds=int(seconds_start)),
        datetime.combine(date_range, time_rangeend)
        + timedelta(seconds=int(seconds_end)),
    )

    st.multiselect(
        label="What satellites would you like to Graph",
        options=[
            "swarma",
            "swarmb",
            "swarmc",
            "epop",
            "dsmp16",
        ],
        key="Satellite_Graph",
    )
    if Graph_value == True:
        Graph()
        button = st.button(label="Render graphs", key="Graph_executer")
        print(button)
        if st.session_state["Graph_executer"] == True:
            Render_Graph(timerange)
            print("done")

    if "Graph" in st.session_state:
        st.pyplot(st.session_state["Graph"][0])

    if Animation_value == True:
        Animation(timerange)


if __name__ == "__main__":
    Main()
