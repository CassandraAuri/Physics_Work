import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from apexpy import Apex
import lompe
from lompe.model.cmodel import Cmodel
from lompe.data_tools import dataloader
event = '2022-12-19'
hour = 14
minute = 5
t0 = dt.datetime(int(event[0:4]), int(event[5:7]), int(event[8:10]), hour, minute)
DT = dt.timedelta(seconds = 2 * 60) # will select data from time +/- DT

# apex object for magnetic coordinates
a = Apex(t0.year)

# locations of data
tempfile_path = './sample_dataset/' # to put the processed data
basepath = tempfile_path + 'raw/' # unprocessed downloads

example_event = '2022-12-19'
# get the processed data set
ssies = dataloader.read_ssusi(event, 'north', basepath, tempfile_path)
iridium = pd.read_hdf(dataloader.read_iridium(example_event, basepath, tempfile_path))