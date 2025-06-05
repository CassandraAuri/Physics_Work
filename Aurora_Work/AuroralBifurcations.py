from datetime import datetime, timedelta
from IPython.display import Video
import numpy as np

import matplotlib.pyplot as plt
import asilib
import asilib.asi
import asilib.map
plt.style.use('dark_background')
location_code = 'fsmi'
time_range = [datetime(2022, 11,11, 5, 23, 0), datetime(2022, 11,11, 5, 26, 0)]

# loglevel is to suppress the verbose ffmpeg output.
asi = asilib.asi.trex_rgb(location_code, time_range=time_range, alt=110)

lat_bounds = (asi.meta['lat']-4, asi.meta['lat']+4)
lon_bounds = (asi.meta['lon']-7, asi.meta['lon'])

ex = asilib.map.create_simple_map(lon_bounds=lon_bounds, lat_bounds=lat_bounds)
plt.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99)
gen = asi.animate_map_gen(overwrite=True, ax=ex, ffmpeg_params={'loglevel':'quiet'}, asi_label=False, color_bounds=[10, 100])

for time, image, ax, im in gen:
    if 'time_label' in locals():
        # This is one way I found to clean up an added plotting object.
        time_label.remove()
    time_label = ex.text(0.99, 0.99, f'location: {location_code} | time: {time}',
                         ha='right', va='top', transform=ex.transAxes, fontsize=15)

