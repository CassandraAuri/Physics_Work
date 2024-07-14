import streamlit as st

import asilib
import asilib.asi
from datetime import datetime
import os
asilib.config["ASI_DATA_DIR"] = os.path.dirname(os.path.abspath(__file__))
asilib.config['ASILIB_DIR'] = os.path.dirname(os.path.abspath(__file__))
time_range=[datetime(2021,4,16,6,39), datetime(2021,4,16,7,42)]
asi = asilib.asi.trex.trex_rgb(
    'fsmi',
    time_range=time_range,
    alt=110,
    colors="rgb",
    custom_alt=True
)
