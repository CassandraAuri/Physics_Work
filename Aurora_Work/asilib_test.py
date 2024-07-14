import streamlit as st

import asilib
import asilib.asi
from datetime import datetime
time_range=[datetime(2021,4,16,6,39), datetime(2021,4,16,7,42)]
asi = asilib.asi.trex.trex_rgb(
    'fsmi',
    time_range=time_range,
    alt=110,
    colors="rgb",
    custom_alt=True
)