import numpy as np
import asilib
import asilib.asi
from datetime import datetime
x_1 = np.linspace(-123, -105, 200)
y_1= np.linspace(62.3, 62.45, 200)

x_2 = np.linspace(-123, -105, 200)
y_2= np.linspace(62.5, 63.6, 200)
time_array=(datetime(2022,12,19,14,4), datetime(2022,12,19,14,6))
asi=asilib.asi.trex_rgb(location_code='yknf', alt=110, time_range=(time_array[0], time_array[-1]), redownload=True)


movie_generator = asi.animate_map_gen()

for i, (time, image, _, im) in enumerate(movie_generator):
    print( image[200,200])