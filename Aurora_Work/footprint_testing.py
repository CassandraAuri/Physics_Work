from datetime import datetime, timedelta
import string
import matplotlib.pyplot as plt
import matplotlib.dates
import matplotlib.patches
import matplotlib.gridspec
import numpy as np
import pyaurorax.ephemeris as ae
import asilib
import asilib.asi
time_range = (datetime(2021, 4, 12, 6, 24),
              datetime(2021, 4, 12, 6, 26))

area_box_km = (20, 20)
#skymap_dict = asilib.load_skymap(asi_array_code, location_code, time_range[0])
print(asilib.config['ASI_DATA_DIR'])
global platforms
platforms = ["swarmb"]
sky_map_values = ['trex_rgb', 'rabb']
rearth = 6378.1370


def emph(cadence):
    print(platforms, "test")

    def altitude(interval):

        altitude_total = []

        print(len(em1.data))
        for i in range(len(platforms)):
            alt = []
            alt_total_not_arranged = []
            for k in range(len(em1.data)):
                alt_total_not_arranged = []
                alt.append(
                    em1.data[k].metadata["radial_distance"]-rearth)
            for j in range(len(alt)-1):
                alt_total_not_arranged.append(np.linspace(
                    start=alt[j], stop=alt[j+1], num=interval))
            alt_total_not_arranged = np.reshape(
                alt_total_not_arranged, -1)
            altitude_total.append(alt_total_not_arranged)
        return altitude_total

    def coordinates_extrapolation(array):
        n = int((time_range[1] - time_range[0]
                 ).total_seconds() / cadence)  # 3 second␣
        latvelocitytotal, lonvelocitytotal = [], []
        interval = int(60/cadence)
        # every 3 seconds theres a time component
        time = np.array([time_range[0] + timedelta(seconds=i*cadence)
                         for i in range(n)])
        for i in range(len(platforms)):
            latvelocity = []
            lonvelocity = []
            for j in range(len(array[0][0])-1):
                latvelocity.append(np.linspace(
                    array[0][i][j], array[0][i][j+1], interval))
                lonvelocity.append(np.linspace(
                    array[1][i][j], array[1][i][j+1], interval))
            latvelocitytotal.append(np.reshape(latvelocity, -1))
            lonvelocitytotal.append(np.reshape(lonvelocity, -1))

        # for i in range(len(platforms)):
        # for j in range(len(platforms)):
        altitude_array = altitude(interval)
        return [time, latvelocitytotal, lonvelocitytotal, altitude_array]

    # image_bounds = imager_bounds()
    global em1
    em1 = ae.search(time_range[0], time_range[1],
                    platforms=platforms)
    lattotal, lontotal = [], []
    for k in range(len(platforms)):
        latstarts, lonstarts = [], [],
        for i in range(len(em1.data)):  # 3 spacecraft
            # sees if the data corresponds to the correct space-craft
            if(em1.data[i].data_source.platform == platforms[k]):
                latstarts.append(em1.data[i].location_geo.lat)  # starting
                lonstarts.append((em1.data[i].location_geo.lon))  # ending
        lattotal.append(latstarts)
        lontotal.append(lonstarts)
    return coordinates_extrapolation(
        np.array([lattotal, lontotal]))

 # [time, satellite_lattitute,.] #add everything in main


def main():
    print(__name__)

    save_file = []
    asi_array_code = sky_map_values[0]
    location_code = sky_map_values[1]
    fig, ax = plt.subplots(
        3, 1, figsize=(7, 10), gridspec_kw={'height_ratios': [4, 1, 1]}, constrained_layout=True
    )

    # retrieves emphermaris data in the format [time, [sattelite_longitudes], [satellite_lattitudes], [satelliite_alitiudes]]
    data = emph(3)
    sat_time = data[0]

    # Finds the footprint of selected region with each selected spacecraft

    def animator():

        for i, (time, image, _, im) in enumerate(movie_generator):
            # Plot the entire satellite track, its current location, and a 20x20 km box
            # around its location.
            print(np.shape(np.array(sat_azel_pixels_total)))
            for j in range(len(sat_azel_pixels_total)):
                ax[0].plot(sat_azel_pixels_total[j][:, 0],
                           sat_azel_pixels_total[j][:, 1])
                ax[0].scatter(sat_azel_pixels_total[j][i, 0], sat_azel_pixels_total[j][i, 1],
                              marker='o', s=50)

            # Annotate the location_code and satellite info in the top-left corner.

    def ASI_logic():

        alt = 110  # footprint value
        if(asi_array_code.lower() == 'themis'):
            asi = asilib.asi.themis(
                location_code, time_range=time_range, alt=alt)
        elif(asi_array_code.lower() == 'rego'):
            asi = asilib.asi.rego(
                location_code, time_range=time_range, alt=alt)
        elif(asi_array_code.lower() == 'trex_nir'):
            asi = asilib.asi.trex.trex_nir(
                location_code, time_range=time_range, alt=alt)
        elif(asi_array_code.lower() == 'trex_rgb'):
            asi = asilib.asi.trex.trex_rgb(
                location_code, time_range=time_range, alt=alt)
        return asi
    asi = ASI_logic()
    # Initiate the movie generator function. Any errors with the data will be␣

    movie_generator = asi.animate_fisheye_gen(  # initaliziation
        ax=ax[0], overwrite=True, cardinal_directions='NE', azel_contours=True
    )
    # Use the generator to get the images and time stamps to estimate mean the ASI
    # brightness along the satellite path and in a (20x20 km) box.
    sat_azel_pixels_total, nearest_pixel_intensity_total, area_intensity_total, area_mask_total = [], [
    ], [], []  # Creates empty arrays for data from satellites to go, which is there rendered in animator

    def Animation_logic_generator():
        nonlocal sat_azel_pixels_total, nearest_pixel_intensity_total, area_intensity_total, area_mask_total

        for i in range(len(platforms)):  # length of spacecraft
            conjunction_obj = asilib.Conjunction(asi, (sat_time, np.array(  # lat_sattelite, long, alt
                [data[1][i], data[2][i], data[3][i]]).T))

            # Normally the satellite time stamps are not the same as the ASI.
            # You may need to call Conjunction.interp_sat() to find the LLA coordinates
            # at the ASI timestamps.
            # Map the satellite track to the imager's azimuth and elevation coordinates and
            # image pixels. NOTE: the mapping is not along the magnetic field lines! You need
            # to install IRBEM and then use conjunction.lla_footprint() before
            # calling conjunction_obj.map_azel.
            sat_azel, sat_azel_pixels = conjunction_obj.map_azel()

            # Need to change masked NaNs to 0s so we can plot the rectangular area contours.
            sat_azel_pixels_total.append(sat_azel_pixels)
    # sat_azel_pixels, area_box_mask_2, asi_brightness_2 for each satellite
    Animation_logic_generator()
    animator()

    print(
        f'Movie saved in {asilib.config["ASI_DATA_DIR"] / "animations"}', )

    movie_container = 'mp4'
    movie_address = f'{time_range[0].strftime("%Y%m%d_%H%M%S")}_' \
        f'{time_range[1].strftime("%H%M%S")}_' \
        f'{asi_array_code.lower()}_{location_code.lower()}_fisheye.{movie_container}'

    movie_address_total = asilib.config["ASI_DATA_DIR"] / \
        'animations'/movie_address
    print(movie_address_total)
    save_file.append(movie_address_total)
    return save_file


if __name__ == '__main__':  # arg parse look up
    test = main()
