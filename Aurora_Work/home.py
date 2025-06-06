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
import pandas as pd
from scipy import spatial
from scipy.optimize import minimize
from viresclient import set_token
from viresclient import SwarmRequest
import asilib.skymap
import warnings
import geopack.geopack as gp
warnings.filterwarnings('ignore')
plt.style.use("cyberpunk")
from scipy.optimize import curve_fit, fsolve
import os
from cycler import cycler
plt.rcParams['figure.facecolor'] =  '121212'
plt.rcParams['axes.facecolor'] =  '121212'
plt.rcParams['savefig.facecolor'] =  '121212'
mplcyberpunk.make_lines_glow()
mplcyberpunk.add_underglow()
asilib.config["ASI_DATA_DIR"] = os.path.dirname(os.path.abspath(__file__))
print(asilib.__version__)


from scipy.optimize import curve_fit
import geopack.geopack as gp
def footprint(time, latitude, longitude, altitude, alt,vsw= [-400,0,0]):
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
    
    t1 = time[0]
    t0 = datetime(1970,1,1) #epoch
    py_dt = pd.to_datetime(t1).to_pydatetime()
    ut = (py_dt-t0).total_seconds()
    lat_sat=np.deg2rad(latitude)
    lon_sat=np.deg2rad(longitude) #converts to radii
    gp.recalc(ut)
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

st.title("Cassandra Litwinowich's Auroral Website!")
st.header("How to Use: (Buttons for drop down menus on left, you will need to scroll down)")
st.divider()
st.write("First we need to find suitable conjunctions. Well first you may ask what a conjunction is. \
            A Conjunction is a closer approach between a satellite with either \
            another satellite or a ground based camera. Using Swarm Aurora \
            [link]('https://swarm-aurora.com/conjunctionFinder/') \
            we can find Swarm satellite conjunctions with itself or ground-based auroral cameras. \
            For a more automated way to look for auroras look at Aurorax [link]('https://aurorax.space']) .\
            A good event has already been preselect for swarm A and C over THEMIS (camera constellation) Inuvik. \
              This event has good Electric Field Data and demonstrates a poleward boundary intensification [link]('https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023143')")
st.divider()
st.write("Great we now understand where we can get conjunctions, we the next question is what \
             relevant information for studying the aurora can we get out these satellite-satellite or \
             satellite-camera conjunctions. First lets tackle the satellites. Satellites give us an in-situ measurement of the Electric Fields and Magnetic fields. \
             Subtracting the mean field generated by the Earth's dynamo, which is done with the CHAOS models [link]('https://earth-planets-space.springeropen.com/articles/10.1186/s40623-016-0486-1'), \
             we are left the pertubation caused the ions entering the ionosphere.The Electric Field is derived from the velocities and is as \
              important as the Magnetic field for understanding the dynamics of the aurora. Additionally, the Poynting Flux and Field Aligned current, \
              are secondary tools derived from these quantities that give us a further view into the dynamics of the aurora. \
             The EB ratio is for more advanced users who are encouraged to contact me at cmckenna@ualberta.ca for information on its workings and are not meant for general viewing at this time. Use with caution!")
st.divider()
st.write("Lastly, lets discuss the ground based cameras. There are four constellations implemented in this GUI. See below for their specifics. Additionally another thing to note is altitudes. \
              The Swarm Satellites are layed out in configurations that change over time but generally swarm A and C are linked together at an altitude of 420km while swarm B is higher at 500km. \
            Most Auroral activity occurs at much lower altitudes  \
            This means that the activity seen in the aurora travells along the magnetic field lines and is then seen by the satellites. Careful care has been taken to ensure this mapping is accurate. Lastly the fisheye camera option is generally preferred over the map function \
            ")
st.divider()
st.write("THEMIS: Oldest camera set, panchromatic, altitude should be 110km unless otherwise specified. Additionally, resolution and quality varies greatly, sampling rate 1 image every 3seconds")
st.write("REGO : Old camera, Focuses on Nitrogen Emissions  \
             [link]('https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023758#:~:text=Despite%20these%20intrinsic%20differences%2C%20the,energy%20emission%20red-line%20aurora.'),\
              altitude should be 150-230, generally 190km works well. Additionally, resolution and quality greatly deteriotes with the age of the camera (2013 launch date), 1 image every 3 seconds")
st.write("TREx-NIR: New camera, focuses on Near Infared, altitude should always be 150km, lower resolution than other cameras and longer integration time: 1 image every 6 seconds")
st.write("TREx-RGB:  New camera, RECOMMENDED,  focuses on Optical Emission, altitude be 110km unless otherwise specified, works very well, 1 image every 3 seconds")




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
        try: #If one of the magnetometers on swarm is not working
            print(ds["Spacecraft"][0].lower() , "test")
            if "".join(("swarm", ds["Spacecraft"][0].lower())) in dict["satellite_graph"]:
                data_stuff.append(ds)
            else:
                pass
        except KeyError: 
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


def graphing_animation(dict):
    data_3, data_6 = empharisis_processing(dict, 3), empharisis_processing(dict, 6)
    time_range = dict["time_range"]
    platforms = dict["satellite_graph"]
    save_file = []
    fig, ax = plt.subplots(  # intializes plots
        1,
        1,
        figsize=(10, 15),
        constrained_layout=True,
    )

    # Finds the footprint of selected region with each selected spacecraft

    def animator_map():

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
    
    # Loops through the number of stations selected by the user
    for k in range(len(dict["sky_map_values"])):

        asi_array_code = dict["sky_map_values"][k][0]
        location_code = dict["sky_map_values"][k][1]
        alt = int(dict["sky_map_values"][k][2])
        plt.cla()

        def ASI_logic():
            if asi_array_code.lower() == "themis":
                frame_rate = 2
                if alt==90 or alt==110 or alt==150:
                    asi = asilib.asi.themis(
                        location_code, time_range=time_range, alt=alt, custom_alt=False 
                    )
                else:
                    asi = asilib.asi.themis(
                        location_code, time_range=time_range, alt=alt, custom_alt=True
                    )
            elif asi_array_code.lower() == "rego":
                frame_rate = 2
                asi = asilib.asi.rego(
                    location_code, time_range=time_range, alt=alt, custom_alt=True
                )
            elif asi_array_code.lower() == "trex_nir":
                frame_rate = 1
                asi = asilib.asi.trex.trex_nir(
                    location_code, time_range=time_range, alt=alt, custom_alt=True
                )

            elif asi_array_code.lower() == "trex_rgb":
                frame_rate = 2
                asi = asilib.asi.trex.trex_rgb(
                    location_code,
                    time_range=time_range,
                    alt=alt,
                    colors="rgb",
                    custom_alt=True,
                    acknowledge=True
                )
            else:
                raise NotImplementedError("How did you get this to happen")

            if dict["sky_map_values"][k][3] == 'Map':
                movie_generator = asi.animate_map_gen(  # initaliziation
                    ax=ax, overwrite=True, ffmpeg_params={"framerate": frame_rate}
                )
            else:
                movie_generator = asi.animate_fisheye_gen(  # initaliziation
                    ax=ax, overwrite=True, ffmpeg_params={"framerate": frame_rate}, color_norm='lin'
                )

            return asi, movie_generator

        asi, movie_generator = ASI_logic()


        lat_satellite, lon_satellite, conj_obj_array,sat_azel_pixels_total = [], [], [], []  
        # Creates empty arrays for data from satellites to go, which is there rendered in animator
        for i in range(len(platforms)):  # length of spacecraft REFACTOR TODO broken here
            print(len(platforms))
            # Trex is 6 second cadence compared to 3 of rego and themos
            if asi_array_code.lower() == "trex_nir":
                data = data_6
            else:
                data = data_3
            sat_time = np.array(data[i][0])  # sets timestamp
            sat_lla = np.array([data[i][1], data[i][2], data[i][3]]).T

            conjunction_obj = asilib.Conjunction(asi, (sat_time, sat_lla))

            # Converts altitude to assumed auroral height
        # Converts altitude to assumed auroral height
            lat_sat=conjunction_obj.sat["lat"].to_numpy()
            lon_sat=conjunction_obj.sat["lon"].to_numpy()

            alt_sat=conjunction_obj.sat["alt"].to_numpy()

            sat_lla = footprint(sat_time,lat_sat,lon_sat,alt_sat, alt)
            print("footdone")
            conjunction_obj = asilib.Conjunction(asi, (sat_time, sat_lla.T))

            lat_satellite.append(conjunction_obj.sat["lat"].to_numpy())
            lon_satellite.append(conjunction_obj.sat["lon"].to_numpy())



            alt_sat=conjunction_obj.sat["alt"].to_numpy()

            if dict["sky_map_values"][k][3] != 'Map':
                conj_obj_array.append(conjunction_obj)
                sat_azel, sat_azel_pixels = conjunction_obj.map_azel()
                sat_azel_pixels_total.append(sat_azel_pixels)
        print("question")
        print(dict["sky_map_values"][k][3])
        # from asilib
        if dict["sky_map_values"][k][3] == 'Map': #creates a map in lat lon
            pixel_chosen = np.zeros((len(platforms), len(data[0][0])))
            lat, lon = asi.skymap["lat"], asi.skymap["lon"]

            values = np.zeros((np.shape(lat)[0], np.shape(lat)[1]))
            non_blind_search = 20
            indicies_total = np.zeros((len(platforms), len(data[0][0]), 2), dtype=int)
            print("startomg aacg,v2")
            lat[np.isnan(lat)] = np.inf
            lon[np.isnan(lon)] = np.inf
            print("step 2")
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
            print("animator map")
            animator_map()
            movie_container = "mp4"
            movie_address = (
                f'{time_range[0].strftime("%Y%m%d_%H%M%S")}_'
                f'{time_range[1].strftime("%H%M%S")}_'
                f"{asi_array_code.lower()}_{location_code.lower()}_map.{movie_container}"
            )  # file address of movie saved by asilib

            # Saves address so movie.py can load it in the GUI
            st.write(r"Aurora_Work/animations{}".format( movie_address))
            save_file.append( r"Aurora_Work/animations/{}".format( movie_address))
        else:
            animator_fisheye()

            movie_container = "mp4"
            movie_address = (
                f'{time_range[0].strftime("%Y%m%d_%H%M%S")}_'
                f'{time_range[1].strftime("%H%M%S")}_'
                f"{asi_array_code.lower()}_{location_code.lower()}_fisheye.{movie_container}"
            )  # file address of movie saved by asilib

            # Saves address so movie.py can load it in the GUI
        
            save_file.append( r"Aurora_Work/animations/{}".format( movie_address))

    return save_file

def Animation_function_caller(Animation_dict):
    def Animate_graph():
        # Gets the figures and axes from cache
        fig, axes,data = st.session_state["Graph"]
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
    if st.session_state['Animation_asi'] == True:
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
    else:
        clip1,clip2,clip3,clip4=None,None,None,None

    if st.session_state['Animation_asi'] == True and "Graph" not in st.session_state:
        if clip3 == None and clip2 == None and clip4 == None:
            combined = clips_array([[clip1]])
        elif clip4 == None and clip3 == None:
            combined = clips_array([[clip1, clip2]])
        elif clip4 == None:
            combined = clips_array([[clip1, clip2, clip3]])
        else:
            combined = clips_array([[clip1, clip2], [clip3, clip4]])

    if "Graph" in st.session_state and 'Alfven_graphs' not in st.session_state:
        Animate_graph()
        clip_graph = VideoFileClip("animationgraph.mp4")
        if clip1 == None:
                combined = clips_array([[clip_graph, alfven_graph]])
        elif clip3 == None and clip2 == None and clip4 == None:
            combined = clips_array([[clip1, clip_graph]])
        elif clip4 == None and clip3 == None:
            combined = clips_array([[clip1, clip2, clip_graph]])
        elif clip4 == None:
            combined = clips_array([[clip1, clip2], [clip3, clip_graph]])
        else:
            combined = clips_array(
                [[clip1, clip2], [clip3, clip4], [clip_graph, clip_graph]]
            )
    if "Graph" in st.session_state and 'Alfven_graphs' in st.session_state:

        Render_Graph(timerange) #gets a new copy
        Animate_graph()
        clip_graph = VideoFileClip("animationgraph.mp4")
        alfven_graph = VideoFileClip("animationAlfven.mp4")
        if clip1 == None:
            combined = clips_array([[clip_graph, alfven_graph]])
        elif clip3 == None and clip2 == None and clip4 == None:
            combined = clips_array([[clip1, clip_graph, alfven_graph]])
        elif clip4 == None and clip3 == None:
            combined = clips_array([[clip1, clip2], [alfven_graph, clip_graph]])
        elif clip4 == None:
            combined = clips_array([[clip1, clip2], [clip3, clip_graph], [clip_graph, alfven_graph]])
        else:
            combined = clips_array(
                [[clip1, clip2], [clip3, clip4], [alfven_graph, clip_graph]]
            )



        

    combined.write_videofile("animation_display.mp4")
    st.video("animation_display.mp4")
    st.session_state["Animation_logic_completed"] = True

def Animation(timerange):
    """
    Animation is required embedded within the GUI as this is what cooperates with streamlit, should be its own python file. However it is self-contained within
    """
    

    def Animation_GUI():
        st.title("Animation Interface:")
        if "station_count" not in st.session_state:
            st.session_state["station_count"] = 1
        if st.session_state["station_count"] == None:
            st.session_state["station_count"] = 1
        st.number_input(
            label="Number of Stations to animate",
            min_value=1,
            max_value=4,
            key="station_count",
        )

        def station_logic():
            count = st.session_state["station_count"]
            print(count)
            print(st.session_state["station_count"])
            col = st.columns(count)  # number of columns selected
            # numpy doesnt have string arrays, only character arrays
            # initializes empty array for columns selected
            values = [[None, None]] * count

            def station_GUI():
                for i in range(count):
                    # initalizies our site state (project already initialized in line 41)
                    if "".join([str(i), "sites"]) not in st.session_state:
                        st.session_state["".join([str(i), "sites"])] = []

                    if "".join([str(i), "heights"]) not in st.session_state:
                        st.session_state["".join([str(i), "sites"])] = []

                    with col[i]:  # each column
                        st.selectbox("Name of Project", [ "THEMIS", "REGO", "trex_nir", "trex_rgb"],key="".join([str(i), "project"]))  # Sets the project
                        # if no project, do not execute site parameters
                        if st.session_state["".join([str(i), "project"])] != []:
                            # if rego project, select site
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "REGO"
                            ):
                                st.selectbox(
                                    "Name of Site",
                                    ["FSMI", "GILL", "RESU", "TALO", "rank", "fsim", "sach"],
                                    key="".join([str(i), "s"]),
                                )

                            # if project is themis, select themis site
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "THEMIS"
                            ):

                                st.selectbox(
                                    "Name of Site",
                                    options=(
                                        "inuv",
                                        "FSMI",
                                        "GILL",
                                        "FSIM",
                                        "rank",
                                        "talo",
                                        "atha",
                                        "tpas",
                                        "gako",
                                        "mcgr",
                                        "kuuj"
                                    ),
                                    key="".join([str(i), "s"]),
                                )
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "trex_nir"
                            ):
                                st.selectbox(
                                    "Name of Site",
                                    ["rabb", "gill"],
                                    key="".join([str(i), "s"]),
                                )
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "trex_rgb"
                            ):
                                st.selectbox(
                                    "Name of Site",
                                    ["rabb", "gill", "fsmi", "pina", "yknf", 'atha'],
                                    key="".join([str(i), "s"]),
                                )


                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "REGO"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=230,
                                    value=190,
                                    step=10,
                                    key="".join([str(i), "heights"]),
                                )
                                st.selectbox(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    key="".join([str(i), "imagetype"]),
                                )
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "THEMIS"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=170,
                                    value=110,
                                    step=10,
                                    key="".join([str(i), "heights"])
                                )
                                st.selectbox(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    key="".join([str(i), "imagetype"]),
                                )
                                
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "trex_nir"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=230,
                                    value=150,
                                    step=10,
                                    key="".join([str(i), "heights"]),
                                )
                                st.selectbox(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    key="".join([str(i), "imagetype"]),
                                )
                            if (
                                st.session_state["".join([str(i), "project"])]
                                == "trex_rgb"
                            ):
                                st.number_input(
                                    "Height of Skymap",
                                    min_value=70,
                                    max_value=170,
                                    value=110,
                                    step=10,
                                    key="".join([str(i), "heights"]),
                                )
                                st.selectbox(
                                    "Type of Image",
                                    ["Fisheye", "Map"],
                                    key="".join([str(i), "imagetype"]),
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
                            st.session_state["".join([str(i), "s"])] != []
                            and st.session_state["".join([str(i), "heights"])] 
                        ):
                            project = st.session_state["".join([str(i), "project"])]

                            sites = st.session_state["".join([str(i), "s"])]

                            heights = st.session_state["".join([str(i), "heights"])]
                            try:
                                maps = st.session_state["".join([str(i), "imagetype"])]
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

        

        # calls to make one column but should dynamiically update
        # station_logic(1)
        button_for_animation = st.button(
            label="Render Animations", key="Animation_executer_asi"
        )
        if button_for_animation == True:
            Animation_function_caller(Animation_dict)

    Animation_GUI()


def Graph():
    st.title("Graph Interface:")
    

    def Graph_options_B(coord_options):
        st.multiselect(
            label="What directions of B would you like to graph",
            options=coord_options,
            key="B_options_to_use",
            help="Subtracted from the CHAOS model to get rid of the mean field, the prodominent component due to the aurora is the East or Polodial Component"
        )
        st.selectbox(
            label="What frequency would you like to use",
            key="Frequency_B",
            options=[ "50Hz", "1Hz"],
        )
        st.checkbox(
            label="Would you like to difference theses V's versus each satellite (requires the LAG option to selected which requires swarm A and C to be selected)",
            key="B_difference",
        )

    def Graph_options_E(coord_options):

        st.multiselect(
            label="What directions of V would you like to graph",
            options=coord_options,
            key="E_options_to_use",
            help="There is no electric field in a static ionosphere, main component is North or Azimuthal"
        )
        st.selectbox(
            label="What frequency would you like to use",
            key="Frequency_E",
            options=["2Hz", "16Hz"],
        )
        st.checkbox(
            label="Would you like to difference theses V's versus each satellite (requires the LAG option to selected which requires swarm A and C to be selected)",
            key="E_difference",
        )

    def Graph_options_F(ignored):
        pass

    def Graph_options_PF(coord_options):

        st.multiselect(
            label="What directions of Ponyting Flux would you like to graph",
            options=coord_options,
            key="PF_options_to_use",
            help="Poynting flux is S=E cross B so main component is centre, E and B not required to be selected"
        )
        st.checkbox(
            label="Would you like to difference theses Ponyting Flux's versus each satellite (requires the LAG option to selected which requires swarm A and C to be selected)",
            key="PF_difference",
        )
    if "B_difference" not in st.session_state:
        st.session_state["B_difference"] = False
    if "E_difference" not in st.session_state:
        st.session_state["E_difference"] = False
    if "PF_difference" not in st.session_state:
        st.session_state["PF_difference"] = False
    options_for_graphs = ["B", "V", "Field Aligned Current", "Poynting flux"]
    Graph_functions = [
        Graph_options_B,
        Graph_options_E,
        Graph_options_F,
        Graph_options_PF,
    ]

    def GUI_interface():
        graphs = st.multiselect(
            label="What would you like to graph the satellite-related measurements (eg: Electric field, Pixel intensity etc)",
            options=options_for_graphs,
            key="Graph_select",
            default=None,
        )

        coordinate_system = st.selectbox(
            label="What coordinate system would you like it in",
            options=["North East Centre", "Mean-field aligned"],
            help="North East Centre is more common for quick study, mean field untested!!!!",
            key="Coordinate_system",
        )
        return coordinate_system, graphs

    def Drop_down_menus(coordinate_system, graphs):
        try:  # tries to generate column but doesn't work if the index of the coordinate system is 0
            # sets coordinate system to give to the functions in Graph_functions
            if coordinate_system == "North East Centre":
                coord_options = ["North", "East", "Centre"]
            elif coordinate_system == "Mean-field aligned":
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
    if "pixel_average" not in st.session_state:
        st.session_state["pixel_average"] = None   
    if st.session_state['Pixel_intensity'] == True:
        st.selectbox(label="Number of pixels to be averaged", options=[1,3,5,7],  key="pixel_average")
    if "Filtering" not in st.session_state:
        st.session_state["Filtering"] = False
    if "Heatmap" not in st.session_state:
        st.session_state["Heatmap"] = False
    if "Conductivies" not in st.session_state:
        st.session_state["Conductivities"] = False
    if "Alfven_Animation" not in st.session_state:
        st.session_state["Alfven_Animation"] = False


    st.checkbox(
        label="would you like to find the normalized difference between Field Aligned Current and Pyonting flux (must select both FAC and pyonting flux centre)",
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
        help="Graphs the pixel intensity from auroral animation, this is extremely expensive computationally (be warned), both finding a footprint for the satellite as well as using a minimization function find the closest pixel"
    )
    st.checkbox(
        label=r"Would you like to bandpass filter your satellite measurements",
        value=st.session_state["Filtering"],
        key="Filtering",
        help="The swarm satellites are travelling at approximately 10km/s so a 0.1Hz and above gets rid of ionospheric signatures greater than 100km in azimuthal length"
    )
    if "swarma" in st.session_state["Satellite_Graph"] and "swarmc" in st.session_state["Satellite_Graph"]:
        st.checkbox(label='Would you like to lag the plots by the synoptic B scale for offsetting as well as using magnitude latitude for the basic plots', value=False, key='lag')

    if "low_pass" not in st.session_state:
        st.session_state["low_pass"] = None
    if "high_pass" not in st.session_state:
        st.session_state["high_pass"] = None

    if "sampling_rate" not in st.session_state:
        st.session_state["sampling_rate"] = None
    
    if "Window_Length" not in st.session_state:
        st.session_state["Window_Length"] = None
    if "nperseg" not in st.session_state:
        st.session_state["nperseg"] = None


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
    if "EB_cross power" not in st.session_state:
        st.session_state["EB_cross power"] = None
    if "EB_cross phase" not in st.session_state:
        st.session_state["EB_cross phase"] = None
    if "lags_cross_B" not in st.session_state:
        st.session_state["lags_cross_B"] = None
    if "lags_cross_E" not in st.session_state:
        st.session_state["lags_cross_E"] = None
    if 'lag' not in st.session_state:
        st.session_state['lag'] = None
    if "graph_type" not in st.session_state:
        st.session_state['graph_type'] = ["heatmap"]
    if "singles_graphs" not in st.session_state:
        st.session_state["singles_graphs"] = None

    global timerange_singles
    timerange_singles=None

    if st.session_state["Filtering"] == True:
        st.select_slider(label="Please select the low pass in Hz", value=0.2, options=[0,0.05,0.1,0.2,0.5, 1, 2, 4], key="low_pass")
        st.select_slider(label="Please select the low pass in Hz", value=7, options=[0.5, 1, 2, 4, 6, 7, 7.9], key="high_pass")
    
    if st.session_state["E_B_ratio"] == True:
        st.multiselect(label="Would you like a heatmap (running windows) or a single interval", options=("heatmap", "single value"), key="graph_type", default=["heatmap"])
        if st.session_state["Coordinate_system"] == "North East Centre":
            options=[["E_North", "E_East"], ["B_North", "B_East"], ["ENorth/BEast ratio", "EEast/BNorth ratio"], ["ENorth/BEast crosspower", "EEast/BNorth crosspower"], ["ENorth/BEast cross phase", "EEast/BNorth cross phase"], ['B B lag cross power', 'B B lag cross phase'], ['E E lag cross power', 'E E lag cross phase']]
        else:
            options=[["E_Azimuthal", "E_Polodial"], ["B_Azimuthal", "B_Polodial"], ["EAzimuthal/BPolodial", "EPolodial/BAzimuthal"]]
        st.write(st.session_state['graph_type'])
        print(st.session_state['graph_type'])
        if st.session_state['graph_type'] == 'heatmap':
            st.select_slider(label="Please select the running window interval (sampling rate)", value=1, options=[0.1,0.2,0.5, 1, 2], key="sampling_rate")
            st.select_slider(label="Please select the running window length", value=4, options=[2,3,4,5,6,8,10, 20,30,40,60,119], key="Window_Length")
            st.select_slider(label="Please select the number of samples per segment (note 16sps)", value='window length', options=['window length', 'half window length', 'quarter window'], key="nperseg")

            st.checkbox(label="Would you like to graph a heatmap of the event through timed (frequency versus time versus ampltiude spectra)", key="Heatmap")
            st.checkbox(label="Would you like to graph the conductivies derived the EB ratio versus time", key="Conductivies")
            st.checkbox(label="Would you like to create an animation of the peridograms through time", key="Alfven_Animation")

            if st.session_state["Heatmap"] == True:
                if st.session_state['lag'] == True:
                    st.multiselect(label="Choose the variable in the heatmap(s)", options=sum(options, []), key="heatmap_graphs")
                else:
                    st.multiselect(label="Choose the variable in the heatmap(s)", options=sum(options[:-2], []), key="heatmap_graphs")

            if st.session_state["Conductivies"] == True:
                st.multiselect(label="Choose the polarization used for conductivies", options=options[2], key="conductivity_graphs")

            if st.session_state["Alfven_Animation"] == True:

                st.multiselect(label="Please select the plots in the Alfven animation", options=["Time Series", "E Periodogram", "B Periodogram", "E/B Periodogram", "Cross Power Spectrum", "Cross phase", 'B B lag'], key="Alfven_graphs")

                if "Time Series" in st.session_state["Alfven_graphs"]:

                    st.multiselect(label="Please select the time series you want to plot through", options=sum(options[:-5], []), key="Time_Series_Graph")

                if "E Periodogram" in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the E polarization you want to plot", options=options[0], key="E_Peridogram_Graph")

                if "B Periodogram" in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the B polarization you want to plot", options=options[1], key="B_Peridogram_Graph")

                if "E/B Periodogram" in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the E/B polarization you would like to plot", options=options[2], key="E/B_Periodogram_Graph")
                
                if "Cross Power Spectrum" in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the E/B cross power you would like to plot", options=options[3], key="EB_cross power")

                if "Cross phase" in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the E/B cross phase you would like to plot", options=options[4], key="EB_cross phase")
                if 'B B lag' in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the B lags you would like to plot", options=options[5], key="lags_cross_B")
                if 'E E lag' in st.session_state["Alfven_graphs"]:
                    st.multiselect(label="Please select the B lags you would like to plot", options=options[6], key="lags_cross_E")
        elif graph_type[0] == "single value":
            time_rangestart = st.time_input(
            label="Start time of interval for a single EB ratio plot",
            step=60,
            key="time_start_single",
            value=st.session_state["time_start"],
            )

            time_rangeend = st.time_input(
                label="end time of interval for a single EB ratio plot",
                step=60,
                key="time_end_single",
                value=st.session_state["time_end"],
            )
            seconds_start = st.slider(
                "Seconds start",
                min_value=0,
                max_value=60,
                step=1,
                value=0,
                key="second_start_single",
            )
            seconds_end = st.slider(
                "Seconds ends",
                min_value=0,
                max_value=60,
                step=1,
                value=60,
                key="second_end_single",
            )
            
            timerange_singles = (
                datetime.combine(st.session_state['date'], st.session_state['time_start_single'])
                + timedelta(seconds=int(seconds_start)),
                datetime.combine(st.session_state['date'], st.session_state['time_end_single'])
                + timedelta(seconds=int(seconds_end)),
            )
            if st.session_state['lag'] == True:
                st.multiselect(label="Choose the variable in the single graphs", options=sum(options, []), key="singles_graphs")
            else:
                st.multiselect(label="Choose the variable in the singles", options=sum(options[:-1], []), key="singles_graphs")
            st.select_slider(label="Please select the number of samples per segment (note 16sps)", value='window length', options=['window length', 'half window length', 'quarter window'], key="nperseg")



            ##TODO IMPLEMENT THE Cross PHASE and POWEr
            
        


        


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
        np.reshape([np.where(np.array(st.session_state["Graph_select"]) == "Field Aligned Current")], -1)[
            0
        ]  # Finds if any index exists named FAC in the graphs selected
        FAC_boolean = True
    except IndexError:
        FAC_boolean = False
    # try:  # [0] doesnt work


    def value_setter(count):  # sets values of selected sites and projects
        nonlocal values
        for i in range(count):  # goes through all columns
            # doesn't add values if empty
            if st.session_state["".join([str(i), "project"])] == []:
                pass

            else:
                # doesn't add values if empty
                if st.session_state["".join([str(i), "s"])] != []:
                    project = st.session_state["".join([str(i), "project"])]

                    sites = st.session_state["".join([str(i), "s"])]

                    heights = st.session_state["".join([str(i), "heights"])]
                    
                    values[i] = [project, sites, heights]
                else:
                    pass
        return values

    if st.session_state["station_count"] !=None:
        count = st.session_state["station_count"]
        values = [[None, None]] * count
        skymap_values = value_setter(count)
    else:
        skymap_values = None
        print("Make none")



    if st.session_state['B_difference'] == True:
        B_difference=st.session_state['B_options_to_use']
    else:
        B_difference = None
    if st.session_state['E_difference'] == True:
        E_difference=st.session_state['E_options_to_use']
    else:
        E_difference = None
    if st.session_state['PF_difference'] == True:
        PF_difference=st.session_state['PF_options_to_use']
    else:
        PF_difference = None


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
        "window_length":st.session_state["Window_Length"],
        "EB_cross power": st.session_state["EB_cross power"],
        "EB_cross phase": st.session_state["EB_cross phase"],
        "lags_cross_B": st.session_state["lags_cross_B"],
        "lags_cross_E": st.session_state["lags_cross_E"],
        "nperseg": st.session_state["nperseg"],
        "lag": st.session_state["lag"],
        "time_range_single": timerange_singles,
        "singles_graph": st.session_state["singles_graphs"],
        "pixel_average": st.session_state["pixel_average"],
        "B_difference": B_difference,
        "E_difference": E_difference,
        "PF_difference": PF_difference,

    }

    if dict["coordinate_system"] == "North East Centre":
        figaxes = EBplotsNEC(dict)
        st.session_state["Graph"] = figaxes

    if dict["coordinate_system"] == "Mean-field aligned":
        fig,axes,data = EBplotsNEC(dict)

        st.session_state["Graph"] = [fig,axes]
    return


def Main():
    if "time_start" not in st.session_state:
        st.session_state["time_start"] = time(7, 38)
    if "time_end" not in st.session_state:
        st.session_state["time_end"] = time(7, 40)
    if "date" not in st.session_state:
        st.session_state["date"] = date(2020, 3, 29)
    if "second_start" not in st.session_state:
        st.session_state["second_start"] = 0
    if "second_end" not in st.session_state:
        st.session_state["second_end"] = 0
    if "Animation_asi" not in st.session_state:
        st.session_state['Animation_asi'] = False
    Animation_value = st.sidebar.checkbox(
        label="Would you like to display an auroral animation",
        label_visibility="visible",
        key='Animation_asi'
    )
    Graph_value = st.sidebar.checkbox(
        label="Would you like to use graph in-situ satellite measurements ",
    )
    date_range = st.date_input(label="Date of conjunction", key="date")
    time_rangestart = st.time_input(
        label="Start time of Interval",
        step=60,
        key="time_start",
        value=st.session_state["time_start"],
    )

    time_rangeend = st.time_input(
        label="End time of interval",
        step=60,
        key="time_end",
        value=st.session_state["time_end"],
    )
    seconds_start = st.slider(
        "Seconds start",
        min_value=0,
        max_value=60,
        step=1,
        value=0,
        key="second_start",
        help="If you would like to change the seconds of the start ie, to be start time + x sec"
    )
    seconds_end = st.slider(
        "Seconds ends",
        min_value=0,
        max_value=60,
        step=1,
        value=0,
        key="second_end",
        help="If you would like to change the seconds of the end ie, to be end time + x sec"
    )
    global timerange
    timerange = (
        datetime.combine(date_range, time_rangestart)
        + timedelta(seconds=int(seconds_start)),
        datetime.combine(date_range, time_rangeend)
        + timedelta(seconds=int(seconds_end)),
    )
    

    st.multiselect(
        label="What satellites would you like to use in your graph",
        options=[
            "swarma",
            "swarmb",
            "swarmc",
        ],
        default=["swarma","swarmc"],
        key="Satellite_Graph",
    )
    if "Animation_executer" not in st.session_state:
        st.session_state["Animation_executer"] =False
    if Graph_value == True:
        Graph()
        button = st.button(label="Render graphs", key="Graph_executer")
        if st.session_state["Alfven_Animation"] == True:
            button = st.button(label="Animate Alfven Graphs", key="Animation_executer")

        if st.session_state["Graph_executer"] == True:
            Render_Graph(timerange)
        if st.session_state["Animation_executer"] == True:
            if st.session_state["Alfven_Animation"] == True:
                Animation_dict = {
                "time_range": timerange,
                "satellite_graph": st.session_state["Satellite_Graph"],
                }
                Animation_function_caller(Animation_dict)



    if "Graph" in st.session_state:
        st.pyplot(st.session_state["Graph"][0])

    if Animation_value == True:
        Animation(timerange)


if __name__ == "__main__":
    Main()