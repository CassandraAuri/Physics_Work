import streamlit as st

import asilib
import asilib.asi
from datetime import datetime
import os
import matplotlib.pyplot as plt
datetimetuple=(datetime(2021,4,12,6,20), datetime(2021,4,12,6,30))
asi= asilib.asi.trex_rgb(time_range=datetimetuple, alt=110, location_code='rabb')
lat_bounds = (asi.meta['lat']-5, asi.meta['lat']+5)
lon_bounds = (asi.meta['lon']-8, asi.meta['lon']+8)

ex= asilib.map.create_simple_map(lon_bounds=lon_bounds, lat_bounds=lat_bounds)
plt.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99)

gen = asi.animate_map_gen(overwrite=True, ax=ex, ffmpeg_params={'loglevel':'quiet'}, asi_label=False)

for time, image, ax, im in gen:
    if 'time_label' in locals():
        # This is one way I found to clean up an added plotting object.
        time_label.remove()
    time_label = ex.text(0.99, 0.99, f'location: RABB | time: {time}',
                         ha='right', va='top', transform=ex.transAxes, fontsize=15)
