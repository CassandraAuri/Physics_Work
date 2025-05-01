from datetime import datetime, timedelta
from IPython.display import Video
import numpy as np

import matplotlib.pyplot as plt
import asilib
import asilib.asi
import asilib.map

plt.style.use('dark_background')

print(f'asilib version: {asilib.__version__}')
location_code = 'GAKO'
time_range = [datetime(2009, 2, 27, 10, 3, 0), datetime(2009, 2, 27, 10, 10, 0)]

# loglevel is to suppress the verbose ffmpeg output.
asi = asilib.asi.themis(location_code, time_range=time_range)
asi.animate_fisheye(overwrite=True, ffmpeg_params={'loglevel':'quiet'})
plt.close()  # To show a clean output in this tutorial---it is often unnecessary.
#9 February 2007