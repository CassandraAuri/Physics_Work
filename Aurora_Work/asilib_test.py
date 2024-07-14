import streamlit as st

import asilib
import asilib.asi
from datetime import datetime
import os
asilib.config["ASI_DATA_DIR"] = os.path.dirname(os.path.abspath(__file__)) #Changes directory to the src of so writing can happen
asilib.config['ASILIB_DIR'] = os.path.dirname(os.path.abspath(__file__)) #Changes directory to the src of so writing can happen
asilib.config["HERE"] = os.path.dirname(os.path.abspath(__file__))
asilib.config["acknowledged_asis"] = []
asilib.acknowledge["CONFIG_PATH"] =  os.path.dirname(os.path.abspath(__file__))
st.write(asilib.acknowledge["CONFIG_PATH"])
time_range=[datetime(2021,4,16,6,39), datetime(2021,4,16,7,42)]
asi = asilib.asi.trex.trex_rgb(
    'fsmi',
    time_range=time_range,
    alt=110,
    colors="rgb",
    custom_alt=True
)
