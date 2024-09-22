import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from apexpy import Apex
import lompe
from lompe.model.cmodel import Cmodel
from lompe.data_tools import dataloader
from viresclient import set_token
from scipy.signal import butter, filtfilt, freqz
from viresclient import SwarmRequest
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
import geopack.geopack as gp
from lompe.utils.conductance import hardy_EUV
from scipy.spatial.transform import Rotation as R
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
def quaternion_inverse_scipy(q):
    # Ensure q is a numpy array
    q = np.asarray(q)
    
    # Create a Rotation object from the quaternion
    rotation = R.from_quat(q)  # Note: scipy uses [x, y, z, w] format
    
    # Compute the inverse rotation
    inverse_rotation = rotation.inv()
    
    
    return inverse_rotation
def unit_array(array):
    arraysum = np.sum(np.abs(array), axis=1)
    # Normalizes and finds unitary
    array_unit = array / arraysum[:, np.newaxis]  # normalizes
    return array_unit

def find_closest_indices(times1, times2):
    # Convert to numpy arrays
    times1 = np.array(times1)
    times2 = np.array(times2)
    
    # Compute the differences between each time in times1 and all times in times2
    # Resulting in a 2D array where each row contains the absolute differences for one time in times1
    differences = np.abs(times1[:, None] - times2)
    
    # Find the index of the minimum difference for each time in times1
    closest_indices = np.argmin(differences, axis=1)
    
    return closest_indices

def butter_bandpass(cutoffs, fs, order=4):
    if cutoffs[0] ==0:
        return butter(order, cutoffs[1], fs=fs, btype="low")
    else:
        return butter(order, cutoffs, fs=fs, btype="band")


def butter_bandpass_filter(data, cutoffs, fs, order=4):
    b, a = butter_bandpass(cutoffs, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def requester(sc_collection, measurement, residual, t0, DT,sampling_step=None, **kwargs):
    try:
        request = SwarmRequest()
        request.set_collection(sc_collection)
        if residual == True:
            request.set_products(
                measurements=measurement,
                models=["CHAOS"],
                residuals=True,
                sampling_step=sampling_step,
            )
        else:
            request.set_products(
                measurements=measurement,
                models=["CHAOS"],
                sampling_step=sampling_step,
            )
        data = request.get_between((t0 - DT),(t0 + DT), **kwargs)
        df = data.as_dataframe()
    except:
        df = []
    return df

def arrangement(time, array, shape):  # arranges B into a useable format for use later
    barranged = np.zeros((len(time), shape))
    # Re-arranges into proper (n x 3 ) matricies, ugly but works
    for j in range(len(time)):
        for k in range(shape):
            barranged[j][k] = array[j][k]
    return barranged
event = '2022-12-19'
file_name = 'all_stations_all2022.netcdf'
# locations of data
tempfile_path = './sample_dataset/' # to put the processed data
basepath = tempfile_path + 'raw/' # unprocessed downloads
supermag = pd.read_hdf(dataloader.read_smag(event, basepath, tempfile_path, file_name))
iridium = pd.read_hdf(dataloader.read_iridium(event, basepath, tempfile_path))
superdarn = pd.read_hdf(dataloader.read_sdarn(event, basepath, tempfile_path))

def footprint(time, latitude, longitude, altitude, alt,vsw= [-400,0,0]):
    """
    time, datetime, time for magnetic data for footprint
    vsw velocity of solar wind, tuple of x,y,z
    longitude of satellite in degrees
    latitude of satellite in degrees
    altitude of satellite in km from centre of earth (should be above ~6371)
    THIS CODE ONLY Works in the NOrthern hemisphere

    """
    def cubic(t, a, b, c, d):
            return a*t**3 + b*t**2 + c*t + d
    def x(t, params_x):
        return cubic(t, *params_x)

    def y(t, params_y):
        return cubic(t, *params_y)

    def z(t, params_z):
        return cubic(t, *params_z)
    def radius(t, params_x, params_y, params_z):
        return np.sqrt(x(t, params_x)**2 + y(t, params_y)**2 + z(t, params_z)**2)
    
    def curve_fit_func(xx,yy,zz, differencealt):
        
        r = np.linspace(1, 1.5, 100000)# construct an array of radiuses from 1-1.5

        radius_data=np.sqrt(xx**2+yy**2+zz**2)

        params_x, _ = curve_fit(cubic, radius_data, xx) #Constructs fits on the traces inward since the spatial resolution produced by geopack is limited.
        params_y, _ = curve_fit(cubic, radius_data, yy)
        params_z, _ = curve_fit(cubic, radius_data, zz)

        

        index_closest=np.argmin(np.abs(radius(r, params_x, params_y, params_z)-(alt-differencealt+6371)/6371))#Find the index that produces the closest radius to the altitude

        return x(r[index_closest],params_x ),y(r[index_closest],params_y ),z(r[index_closest],params_z )
    
    t1 = time
    tepoch = datetime(1970,1,1) #epoch
    ut = (t1-tepoch).total_seconds()
    lat_sat=np.deg2rad(latitude)
    lon_sat=np.deg2rad(longitude) #converts to radii
    gp.recalc(ut)
    r, theta= gp.geodgeo(altitude,lat_sat,1) #this r accounts for earths oblateness, so we need to find the difference between my 6371 assumption and the real value and account for that
    differencearray= (altitude+6371)-r
    x_gc,y_gc,z_gc = gp.sphcar((r)/6371,theta,lon_sat,1)  #spherical to cartesian
    

    x_gsm, y_gsm, z_gsm = gp.geogsm(x_gc,y_gc,z_gc, 1) #cartesian to gsm

    x_foot,y_foot,z_foot=np.zeros(len(x_gsm)), np.zeros(len(y_gsm)), np.zeros(len(z_gsm)) #initalize an array
    for index in range(len(x_gsm)):

        x_foot_int, y_foot_int, z_foot_int, xx, _,zz = gp.trace(x_gsm[index], y_gsm[index], z_gsm[index], dir=1,rlim=1.6, maxloop=300 ) #traces each set of lat,lon,alt outward
        _, _, _, xx2,yy2,zz2 = gp.trace(x_foot_int, y_foot_int, z_foot_int, dir=-1,rlim=100, maxloop=1000 )#Traces inward
        x_foot[index],y_foot[index],z_foot[index] = curve_fit_func(xx2,yy2,zz2, differencearray[index])




            

    x_done, y_done, z_done = gp.geogsm(x_foot, y_foot, z_foot, -1)

    alt_sat_done, lat_sat_done,lon_sat_done = np.zeros(len(x_done)), np.zeros(len(x_done)), np.zeros(len(x_done))
    for index in range(len(x_done)):
        
        r_done,theta_done,lon_sat_done[index]= gp.sphcar(x_done[index], y_done[index], z_done[index],-1)

        alt_sat_done[index], lat_sat_done[index]= gp.geodgeo(r_done*6371,theta_done,-1) #TODO check if this is right



    if np.any(np.abs(alt_sat_done - alt) > 5):
        raise Exception("One or more values in the footprinting are greater than 5km away from the specified alt. Contact owner for a fix, not your fault")

    sat_lla=np.array([ np.rad2deg(lon_sat_done)-360,np.rad2deg(lat_sat_done), alt_sat_done])
    return sat_lla

def Lompe(t0, DT):
    def Grid():
        #grid
        position = (-90.3718,80.454) # lon, lat
        orientation = 90 # east, north
        L, W = 9000e3, 9000e3 # dimensions of grid (L is along orientation vector)
        Lres, Wres = 200.e3, 200.e3 # resolution of grid (Lres is along orientation vector)
        grid = lompe.cs.CSgrid(lompe.cs.CSprojection(position, orientation), L, W, Lres, Wres, R = 6481.2e3)
        return grid

    def Condutance(grid):
        Kp = 2 # for Hardy conductance model
        SH = lambda lon = grid.lon, lat = grid.lat: hardy_EUV(lon, lat, Kp, t0, 'hall', F107=152    )
        SP = lambda lon = grid.lon, lat = grid.lat: hardy_EUV(lon, lat, Kp, t0, 'pedersen', F107=152)
        return SH, SP
    
    def Supermag():
        sm = supermag[t0 - DT : t0 + DT].dropna()

        # make the Lompe data object
        B = np.vstack((sm.Be.values, sm.Bn.values, sm.Bu.values))
        coords = np.vstack((sm.lon.values, sm.lat.values))
        # lompe.Data is the Lompe data format
        #sm_data = lompe.Data(B * 1e-9, coords, datatype = 'ground_mag', scale = 100e-9)

        sm_data = lompe.Data(B * 1e-9, coords, datatype = 'ground_mag', error = 10e-9, iweight = 1)
        return sm_data
    
    def Iridium():
        

        # select the time interval of interest
        irid = iridium[(iridium.time >= t0) & (iridium.time <= t0 + DT)]

        # make the Lompe data object
        irid_B = np.vstack((irid.B_e.values, irid.B_n.values, irid.B_r.values))
        irid_coords = np.vstack((irid.lon.values, irid.lat.values, irid.r.values))
        # lompe.Data is the Lompe data format

        iridium_data = lompe.Data(irid_B * 1e-9, irid_coords, datatype = 'space_mag_fac', error = 70e-9, iweight=0.7)
        return iridium_data
    
    def SwarmB():
        collectionB_01 = ["SW_OPER_MAGA_LR_1B", "SW_OPER_MAGB_LR_1B", "SW_OPER_MAGC_LR_1B"]
        swarmmags=[]
        alt=110
        BUsed_swarmb=[]
        coords_swarmb=[]
        for i in range(3):
            dsB = requester(
                collectionB_01[i],
                "B_NEC",
                True,
                t0,
                DT,
                asynchronous=False,
                show_progress=False,
            )
            

            Bdata = dsB["B_NEC_res_CHAOS"]  # data package
            time = Bdata.index.to_numpy()
            Bdata = Bdata.to_numpy()

            barranged = arrangement(time, Bdata, 3)

            if( np.mean(dsB['Latitude'].to_numpy())>55):
                print(np.mean(dsB['Latitude'].to_numpy()))
                coords= footprint( t0,dsB['Latitude'].to_numpy(), dsB['Longitude'].to_numpy(), (dsB["Radius"].to_numpy()-6.371e6)/1e3, alt)
                coords[2] = dsB["Radius"].to_numpy() #After talking to Karl, keep altitude but change lat and lon appropriately

            else:
                coords=np.array([dsB['Longitude'].to_numpy(), dsB['Latitude'].to_numpy(), dsB["Radius"].to_numpy() ])
            Bused=np.array([ barranged[:,1],barranged[:,0], -1*barranged[:,2]])

            if i==1:
                BUsed_swarmb.append(Bused)
                coords_swarmb.append(coords)

            swarmmags.append(lompe.Data(Bused*1e-9,coords, datatype='space_mag_full', error=5e-9,iweight=1))

        return swarmmags

        
    def SSUIS():
        username = 'Cassandra Mckenna'
        email = 'cmckenna@ualberta.ca'
        aff = 'University of Alberta'

        sat = [16,17,18] # DMSP satellite number
        ssies_data_arr=[]
        madrigal_kwargs = {'user_fullname' : username, 'user_email' : email, 'user_affiliation' : aff}
        for i in range(len(sat)):
            ssies = pd.read_hdf(dataloader.read_ssies(event, sat[i], basepath, tempfile_path, **madrigal_kwargs))

            # select the time interval of interest
            ssies = ssies[t0 - DT : t0 + DT]

            # make the Lompe data object
            v_crosstrack = np.abs(ssies.hor_ion_v).values
            coords = np.vstack((ssies.glon.values, ssies.gdlat.values))
            los  = np.vstack((ssies['le'].values, ssies['ln'].values))
            # lompe.Data is the Lompe data format
            #ssies_data = lompe.Data(np.abs(ssies.hor_ion_v).values, coords, datatype = 'convection', scale = 500, LOS = los)

            #'scale' kw deprecated in favor of 'error' and 'iweight' (importance weight) keywords

            ssies_data_arr.append(lompe.Data(np.abs(ssies.hor_ion_v).values, coords, datatype = 'convection', error = 50, iweight=1, LOS = los))
        return ssies_data_arr
    


    def Superdarn():
        # select the time interval of interest (also remove very high speeds)
        sd = superdarn.loc[(superdarn.index >= t0 - DT) & (superdarn.index <= t0 + DT) & (superdarn.vlos < 2000)].dropna()

        # make the Lompe data object
        vlos = sd['vlos'].values
        coords = np.vstack((sd['glon'].values, sd['glat'].values))
        los  = np.vstack((sd['le'].values, sd['ln'].values))
        # lompe.Data is the Lompe data format
        #sd_data = lompe.Data(vlos, coordinates = coords, LOS = los, datatype = 'convection', scale = 500 ) 

        sd_data = lompe.Data(vlos, coordinates = coords, LOS = los, datatype = 'convection', error = 50, iweight=0.6 ) 
        return sd_data
    


    def SwarmE():
        #Swarm A's magnetics are not great but swarm B's are okay, we use this to our advantage to increase our resolution in E

        measurements_E = [
                "VsatN",
                "VsatE",
                "VsatC",
                "Vixv",
                "Viy",
                "Viz",
                "Quality_flags",
                'Bx',
                'By',
                'Bz'
            ]
        ds = requester( 
            "SW_EXPT_EFIB_TCT02", #Mag B, high resolution, 50Hz B (Magnetic field)
            measurements_E, #Magnetic field in NEC coordinates
            True, 
            t0,
            DT,
            asynchronous=False,
            show_progress=False)
        dsB = requester( 
            'SW_OPER_MAGB_HR_1B', #Mag B, high resolution, 50Hz B (Magnetic field)
            ["q_NEC_CRF"], #Magnetic field in NEC coordinates
            False, 
            t0,
            DT,
            asynchronous=False,
            show_progress=False)
        Etime=ds.index
        indicies=find_closest_indices(ds.index, dsB.index)
        quatnecrf=dsB["q_NEC_CRF"].to_numpy()[indicies]
        quaternions = []
        vsat=np.array([ds["Vixv"] , ds["Viy"], ds["Viz"]]).T
        vnec=[]
        for i in range(len(quatnecrf)):
            inverse_quat = quaternion_inverse_scipy(dsB["q_NEC_CRF"].to_numpy()[indicies][i])
            rot_NEC_V= inverse_quat.apply(vsat[i])
            vnec.append(rot_NEC_V)

        vnec=np.array(vnec)
        if np.mean(ds['Latitude'].to_numpy())>55:
            coords_E= footprint( t0,ds['Latitude'].to_numpy(), ds['Longitude'].to_numpy(), (ds["Radius"].to_numpy()-6.371e6)/1e3, 110 )[:-1]
        else:
            coords_E = np.array([ds['Latitude'].to_numpy(), ds['Longitude'].to_numpy()])
        vused=np.array([vnec[:,1], vnec[:,0]])
        vused = np.where((vused > 4000) | (vused < -4000), np.nan, vused)

        Eswarm = lompe.Data(vused, coordinates = coords_E, datatype = 'convection', error = 10, iweight=1.0 ) 

        return Eswarm
    

    grid_L = Grid()
    SH, SP = Condutance(grid_L)
    supermag_L = Supermag()
    irdium_L = Iridium()
    swarmb_L = SwarmB()
    superdarn_L = Superdarn()
    swarme_L = SwarmE()
    SSIES_L = SSUIS()

    model = lompe.Emodel(grid_L, Hall_Pedersen_conductance = (SH, SP))

    # add data
    model.add_data(irdium_L)
    model.add_data(supermag_L)
    model.add_data(swarmb_L[0],swarmb_L[1], swarmb_L[2] )
    model.add_data(SSIES_L[0],SSIES_L[1], SSIES_L[2] )
    model.add_data(superdarn_L)
    model.add_data(swarme_L)


    # run inversion
    model.run_inversion(l1 = 10, l2 = 1)

    return model
    



def Animation():
    def animate(i):
        print(i)
        plt.clf()

        #Given 50 frames, with 30 seconds between each frame, that comes 25 minutes, we want the integration time to be 3 minute on each side so 6 min
        t0 = datetime(2022,12,19,13,50) + timedelta(seconds=30*i)
        DT = dt.timedelta(seconds = 5 * 60)
        a = Apex(t0.year)



        model = Lompe(t0, DT)
        def Convection_plot():
            lompe.visualization.plot_potential(ax[0], model)

            lompe.visualization.plot_quiver(ax[0], model, 'convection')
            lompe.visualization.plot_datasets(ax[0], model, 'convection')

            lon = np.linspace(-123, -105, 200)
            lat= np.linspace(62.3, 62.45, 200)



            x, y= model.grid_J.projection.geo2cube(lon, lat)
            scat = ax[0].plot(x, y, color='purple', label='Location of Arc Based on Ewogram')
            x, y= model.grid_J.projection.geo2cube( -114.3718, 62.4540)
            lon = np.linspace(-123, -105, 200)
            lat= np.linspace(62.5, 63.6, 200)
            x, y= model.grid_J.projection.geo2cube(lon, lat)
            scat = ax[0].plot(x, y, color='purple')

            plt.legend()
            lompe.visualization.plot_mlt(ax[0], model, t0, a)
            return
        
        def FAC_plot():
            lompe.visualization.plot_quiver( ax[1], model, 'space_mag_fac')
            lompe.visualization.plot_contour( ax[1], model, 'fac')


            lompe.visualization.plot_datasets(ax[1], model, 'space_mag_fac')
            lompe.visualization.plot_datasets(ax[1], model, 'space_mag_full', np.linspace(-0.8, 0.8, 40) * 1e-6 * 2)

            lon = np.linspace(-123, -105, 200)
            lat= np.linspace(62.3, 62.45, 200)



            x, y= model.grid_J.projection.geo2cube(lon, lat)
            scat = ax[1].plot(x, y, color='purple', label='Location of Arc Based on Ewogram')
            x, y= model.grid_J.projection.geo2cube( -114.3718, 62.4540)
            lon = np.linspace(-123, -105, 200)
            lat= np.linspace(62.5, 63.6, 200)
            x, y= model.grid_J.projection.geo2cube(lon, lat)
            scat = ax[1].plot(x, y, color='purple')
            lompe.visualization.plot_mlt(ax[1], model, t0, a)
            return
        
        #Convection_plot()
        #FAC_plot()
        fig = lompe.lompeplot(model, include_data = True, time = t0, apex = a, 
                      colorscales = {'fac'        : np.linspace(-0.8, 0.8, 40) * 1e-6 * 2,
                                     'ground_mag' : np.linspace(-380, 380, 50) * 1e-9 / 3, # upward component
                                     'hall'       : np.linspace(0, 6, 32), # mho
                                     'pedersen'   : np.linspace(0, 6, 32)}, # mho
                      quiverscales={'ground_mag'      : 600*1e-9, 
                                    'space_mag_fac'   : 1200*1e-9, 
                                    'space_mag_full'  : 1200*1e-9, 
                                    'electric_current': 1000* 1e-3,
                                    'convection': 7000,
                                    })
        plt.suptitle(f"centred at time {t0} with +_ {DT} seconds")

    fig, ax = plt.subplots(figsize=(10,15), nrows=2)  #Plot convection and FAC 
    animate(25)
    plt.show()
    #anim = FuncAnimation(
    #    fig,
    #    animate,
    #    frames = 50,
    #    interval = 1000/30,
    #)
    


    
    FFwriter = animation.FFMpegWriter(fps=2)  # Writes to mp4
    #anim.save("lompe25minlompeplot.mp4", writer=FFwriter)


        
if __name__ == "__main__":

    Animation()

