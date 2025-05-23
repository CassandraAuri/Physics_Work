import streamlit as st

import asilib
import asilib.asi
from datetime import datetime
import os
st.write(os.path.dirname(os.path.abspath(__file__)))
asilib.config["ASI_DATA_DIR"] = os.path.dirname(os.path.abspath(__file__)) #Changes directory to the src of so writing can happen
asilib.config['ASILIB_DIR'] = os.path.dirname(os.path.abspath(__file__)) #Changes directory to the src of so writing can happen
asilib.config["HERE"] = os.path.dirname(os.path.abspath(__file__))


time_range=datetime(2021,4,16,6,39, 0)
try:
    st.write(
          'start'
    )
    asi = asilib.asi.trex.trex_rgb(
        'fsmi',
        time=time_range,
        alt=110,
        colors="rgb",
        custom_alt=True
    )
except PermissionError:
        st.write('erro1')
        asi = asilib.asi.trex.trex_rgb(
        'fsmi',
        time=time_range,
        alt=110,
        colors="rgb",
        custom_alt=True
    )

