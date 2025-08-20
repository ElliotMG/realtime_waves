"""
Adapted from Steve Woolnough's realtime_waves_v1.py
found on /ncas_climate_vol2/public/seasia/waves/
Uses the real-time methodology of Yang et al. (2021),
based on the original methodology of Yang et al. (2003).
"""
def read_data(filename, latmax=None):
    import xarray as xr
    import numpy as np
    if latmax is None:
        latmax = 24
    print('reading in data')
    fields = xr.open_dataset(filename)
    u = fields['ua'].sel(lat=slice(latmax, -latmax))
    v = fields['va'].sel(lat=slice(latmax, -latmax))
    z = fields['zg'].sel(lat=slice(latmax, -latmax))
    u = u.where(~np.isnan(u), 0)
    v = v.where(~np.isnan(v), 0)
    z = z.where(~np.isnan(z), 0)
    lons = fields['lon'].values
    lats = u['lat'].values
    press = fields['pressure'].values
    times = fields['time'].values
    print('u shape:', u.shape, 'NaNs:', np.isnan(u).sum())
    print('v shape:', v.shape, 'NaNs:', np.isnan(v).sum())
    print('z shape:', z.shape, 'NaNs:', np.isnan(z).sum())
    print('lats:', lats)
    return u, v, z, lons, lats, press, times

def uz_to_qr(u,z,g_on_c) :
    """"
    transform u,z to q, r using q=z*(g/c) + u; r=z*(g/c) - u 
    """
    print('transforming u,z to q,r')
    q=z*g_on_c+u
    r=z*g_on_c-u

    return q,r

def filt_project (qf,rf,vf,lats,y0,waves,pmin,pmax,kmin,kmax,c_on_g)  :
    
    import numpy as np
    print('projecting PCFs')
# find size of arrays 
    nf=qf.shape[0]  
    nz=qf.shape[1]
    nlats=lats.size
    nk=qf.shape[3]
# Find frequencies and wavenumbers corresponding to pmin,pmax and kmin,kmax in coeff matrices
    f=np.fft.fftfreq(nf,0.25)
    fmin=np.where((f >= 1./pmax))[0][0]
    fmax=(np.where((f > 1./pmin))[0][0])-1
    f1p=fmin
    f2p=fmax+1
    f1n=nf-fmax 
    f2n=nf-fmin+1
    k1p=kmin
    k2p=kmax+1
    k1n=nk-kmax
    k2n=nk-kmin+1
 
### note that need to adjust for the fact that array referencing doesn't include last point!!!!!!

# Define the parabolic cylinder functions
    spi2=np.sqrt(2*np.pi)
    dsq=np.array([spi2,spi2,2*spi2,6*spi2]) # Normalization for the 1st 4 parabolic CF
    d=np.zeros([dsq.size,nlats])
    y=lats[:]/y0
    ysq=y**2
    d[0,:]=np.exp(-ysq/4.0)
    d[1,:]=y*d[0,:]
    d[2,:]=(ysq-1.0)*d[0,:]
    d[3,:]=y*(ysq-3.0)*d[0,:]
    dlat=np.abs(lats[0]-lats[1])
    
    qf_Kel=np.zeros([nf,nz,nk],dtype=complex)
    qf_mode=np.zeros([dsq.size,nf,nz,nk],dtype=complex)
    rf_mode=np.zeros([dsq.size,nf,nz,nk],dtype=complex)
    vf_mode=np.zeros([dsq.size,nf,nz,nk],dtype=complex)

# reorder the spectral coefficents to make the latitudes the last dimension
    qf=np.transpose(qf,[0,1,3,2])
    rf=np.transpose(rf,[0,1,3,2])
    vf=np.transpose(vf,[0,1,3,2])

    for m in np.arange(dsq.size) :
        if m == 0:   # For eastward moving Kelvin Waves
            qf_Kel[f1n:f2n,:,k1p:k2p]=np.sum(qf[f1n:f2n,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
            qf_Kel[f1p:f2p,:,k1n:k2n]=np.sum(qf[f1p:f2p,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        # For westward moving waves
        qf_mode[m,f1n:f2n,:,k1n:k2n]=np.sum(qf[f1n:f2n,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        qf_mode[m,f1p:f2p,:,k1p:k2p]=np.sum(qf[f1p:f2p,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        rf_mode[m,f1n:f2n,:,k1n:k2n]=np.sum(rf[f1n:f2n,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        rf_mode[m,f1p:f2p,:,k1p:k2p]=np.sum(rf[f1p:f2p,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        vf_mode[m,f1n:f2n,:,k1n:k2n]=np.sum(vf[f1n:f2n,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        vf_mode[m,f1p:f2p,:,k1p:k2p]=np.sum(vf[f1p:f2p,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
    
    uf_wave=np.zeros([waves.size,nf,nz,nlats,nk],dtype=complex)
    zf_wave=np.zeros([waves.size,nf,nz,nlats,nk],dtype=complex)
    vf_wave=np.zeros([waves.size,nf,nz,nlats,nk],dtype=complex)

    for w in np.arange(waves.size) :
        if waves[w] == 'Kelvin' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*qf_Kel[:,:,:]*d[0,j]
                zf_wave[w,:,:,j,:]=\
                    0.5*qf_Kel[:,:,:]*d[0,j]*c_on_g

        if waves[w] == 'WMRG' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*qf_mode[1,:,:,:]*d[1,j]
                zf_wave[w,:,:,j,:]=\
                    0.5*qf_mode[1,:,:,:]*d[1,j]*c_on_g
                vf_wave[w,:,:,j,:]=\
                    vf_mode[0,:,:,:]*d[0,j]

        if waves[w] == 'R1' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[2,:,:,:]*d[2,j]-rf_mode[0,:,:,:]*d[0,j])
                zf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[2,:,:,:]*d[2,j]+rf_mode[0,:,:,:]*d[0,j])*c_on_g
                vf_wave[w,:,:,j,:]=\
                    vf_mode[1,:,:,:]*d[1,j]
        if waves[w] == 'R2' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[3,:,:,:]*d[3,j]-rf_mode[1,:,:,:]*d[1,j])
                zf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[3,:,:,:]*d[3,j]+rf_mode[1,:,:,:]*d[1,j])*c_on_g
                vf_wave[w,:,:,j,:]=\
                    vf_mode[2,:,:,:]*d[2,j]

    return uf_wave,zf_wave,vf_wave

def write_data(u_wave, z_wave, v_wave, lons, lats, press, times, waves):
    import xarray as xr
    import numpy as np
    outfile = f'/gws/nopw/j04/kscale/USERS/emg/wave_data/DYAMOND3/waves_D3_{sim.name}.nc'
    print('writing data to file')
    # Create a dataset
    ds = xr.Dataset(
        {
            "u_wave": (["wave_type", "time", "pressure", "latitude", "longitude"], u_wave, {
                "units": "m s-1",
                "long_name": "Wave Zonal Wind",
                "standard_name": "eastward_wind"
            }),
            "v_wave": (["wave_type", "time", "pressure", "latitude", "longitude"], v_wave, {
                "units": "m s-1",
                "long_name": "Wave Meridional Wind",
                "standard_name": "northward_wind"
            }),
            "z_wave": (["wave_type", "time", "pressure", "latitude", "longitude"], z_wave, {
                "units": "m",
                "long_name": "Wave Geopotential Height",
                "standard_name": "geopotential_height"
            }),
        },
        coords={
            "longitude": (["longitude"], lons, {
                "units": "degrees_east",
                "long_name": "longitude"
            }),
            "latitude": (["latitude"], lats, {
                "units": "degrees_north",
                "long_name": "latitude"
            }),
            "pressure": (["pressure"], press, {
                "units": "hPa",
                "long_name": "pressure"
            }),
            "time": (["time"], times, {
                "long_name": "time",
                # Uncomment and adjust the units if needed
                # "units": "days since 2016-06-22 00:00:00"
            }),
            "wave_type": (["wave_type"], waves, {
                "long_name": "wave_type"
            }),
        },
        attrs={
            "description": "Wave field anomalies using realtime methodology of Yang et al. (2021)",
            "history": "Created by realtime waves.py"
        }
    )

    # Save the dataset to a NetCDF file
    ds.to_netcdf(outfile)

    print(f"File saved to {outfile}")


    
#---------------------------------------------------------
# start of main routine

# Load packages
import time
import numpy as np
import math as maths
import cartopy.crs as ccrs
import intake
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import easygems.healpix as egh
import healpy as hp
import sys

#---------------------------------------------------------
# Set up specifics for the hackthon catalogue data
import warnings

# Suppress specific FutureWarnings matching the message pattern when using cat[...].to_dask()
warnings.filterwarnings(
    "ignore",
    message=".*The return type of `Dataset.dims` will be changed.*",
    category=FutureWarning,
)

cat = intake.open_catalog('https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml')['online']

sim = cat['um_CTC_km4p4_RAL3P3_n1280_GAL9_nest']

# # =========================================================
# # PREPROCESSING STEP FOR DYAMOND3 (HACKATHON) DATA 
# # =========================================================

ds = (sim(zoom=8, time='PT3H', chunks="auto").to_dask().pipe(egh.attach_coords)) # Zoom=6 is about 1deg, PT3H = 3-hourly resolution
uhp=ds.ua

lon = np.arange(0, 360, 1)
lat = np.arange(90, -91, -1)
lon2d,lat2d = np.meshgrid(lon,lat)

pix = hp.ang2pix(ds.crs.healpix_nside, 
               lon2d,
               lat2d,
               nest=True,
               lonlat=True)

u_ll = uhp.isel(cell=(("lat","lon"), pix))
u_ll_sized = u_ll.sel(pressure=[100,200,300,400,500,700,850,925])

vhp = ds.va
v_ll = vhp.isel(cell=(("lat","lon"), pix))
v_ll_sized = v_ll.sel(pressure=[100,200,300,400,500,700,850,925])

Zhp = ds.zg
Z_ll = Zhp.isel(cell=(("lat","lon"), pix))
Z_ll_sized = Z_ll.sel(pressure=[100,200,300,400,500,700,850,925])


# try:
print('about to create dataset')
ds_uvz = xr.Dataset(
    {
        'ua': (('time', 'pressure', 'lat', 'lon'), u_ll_sized.data),
        'va': (('time', 'pressure', 'lat', 'lon'), v_ll_sized.data),
        'zg': (('time', 'pressure', 'lat', 'lon'), Z_ll_sized.data),
    },
    coords={
        'pressure': u_ll_sized.pressure.data,
        'lat': lat,
        'lon': lon,
    }
)
print('writing...')
ds_uvz.to_netcdf(f'/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND3/uvz_D3_{sim.name}.nc')
print('written')
# except Exception as e:
    # print(f"Error: {e}")

# sys.exit(0)
# -----------------------------------------------------
# MAIN PROCESSING STEP
# -----------------------------------------------------
start_time = time.time()

# define some phyiscal parameters
g=9.8
beta=2.3e-11
radea=6.371e6
spd=86400.
ww=2*np.pi/spd


# Define some parameters spefic to the methodology
latmax=24.  #   +/- latitude range over which to process data.
kmin=2      # minimum zonal wavenumber
kmax=40     # maximum zonal wavenumber
pmin=2.0      # minimum period (in days)
pmax=30.0   # maximum period (in days)
y0=6.0      # meridional trapping scale (degrees)
waves=np.array(['Kelvin','WMRG','R1','R2']) # List of wave types to output


y0real= 2*np.pi*radea*y0/360.0   # convert trapping scale to metres
ce=2*y0real**2*beta
g_on_c=g/ce
c_on_g=ce/g


#read data
u,v,z, lons,lats,press,times = read_data(f'/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND3/uvz_D3_{sim.name}.nc',latmax=24)


#convert u,z to q,r
q,r = uz_to_qr(u,z,g_on_c)


# Fourier transform in time and longitude
qf=np.fft.fft2 (q,axes=(0,-1))
rf=np.fft.fft2 (r,axes=(0,-1))
vf=np.fft.fft2 (v,axes=(0,-1))


# Project onto individual wave modes
uf_wave,zf_wave,vf_wave=filt_project(qf,rf,vf,lats,y0,waves,pmin,pmax,kmin,kmax,c_on_g)


# Inverse Fourier transform in time and longitude
u_wave=np.real(np.fft.ifft2 (uf_wave,axes=(1,-1)))
z_wave=np.real(np.fft.ifft2 (zf_wave,axes=(1,-1)))
v_wave=np.real(np.fft.ifft2 (vf_wave,axes=(1,-1)))


print(u_wave[:,0,0,0,0])
print(z_wave[:,0,0,0,0])
print(v_wave[:,0,0,0,0])


write_data(u_wave,z_wave,v_wave,lons,lats,press,times,waves)
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Waves calculated in {elapsed_time:.2f} seconds.")

# End of main routine
# ----------------------------------------------------------------



##################################################################
######### FOR LOOPING OVER ENSEMBLE MEMBERS ######################
##################################################################

# import time
# import numpy as np

# start_time = time.time()

# # Define the range of members to loop over
# members = [f'em0{i}' for i in range(1, 6)]  # Generates ['em01', 'em02', 'em03', 'em04', 'em05']
# season = '20200120'

# # Define some physical parameters
# g = 9.8
# beta = 2.3e-11
# radea = 6.371e6
# spd = 86400.
# ww = 2 * np.pi / spd

# # Define some parameters specific to the methodology
# latmax = 24.  # +/- latitude range over which to process data.
# kmin = 2      # minimum zonal wavenumber
# kmax = 40     # maximum zonal wavenumber
# pmin = 2.0    # minimum period (in days)
# pmax = 30.0   # maximum period (in days)
# y0 = 6.0      # meridional trapping scale (degrees)
# waves = np.array(['Kelvin', 'WMRG', 'R1', 'R2'])  # List of wave types to output

# y0real = 2 * np.pi * radea * y0 / 360.0  # Convert trapping scale to metres
# ce = 2 * y0real**2 * beta
# g_on_c = g / ce
# c_on_g = ce / g

# # Loop over each member
# for member in members:
#     model = f'CTC_N2560_GAL9_{member}'
#     print(f"Processing model: {model}")

#     # Read data
#     input_file = f'/gws/nopw/j04/kscale/USERS/emg/data/wave_data/KSE/KS_DW/{model}_an40d_{season}.nc'
#     u, v, z, lons, lats, press, times = read_data(input_file, latmax=24)

#     # Convert u, z to q, r
#     q, r = uz_to_qr(u, z, g_on_c)

#     # Fourier transform in time and longitude
#     qf = np.fft.fft2(q, axes=(0, -1))
#     rf = np.fft.fft2(r, axes=(0, -1))
#     vf = np.fft.fft2(v, axes=(0, -1))

#     # Project onto individual wave modes
#     uf_wave, zf_wave, vf_wave = filt_project(qf, rf, vf, lats, y0, waves, pmin, pmax, kmin, kmax, c_on_g)

#     # Inverse Fourier transform in time and longitude
#     u_wave = np.real(np.fft.ifft2(uf_wave, axes=(1, -1)))
#     z_wave = np.real(np.fft.ifft2(zf_wave, axes=(1, -1)))
#     v_wave = np.real(np.fft.ifft2(vf_wave, axes=(1, -1)))

#     print(u_wave[:, 0, 0, 0, 0])
#     print(z_wave[:, 0, 0, 0, 0])
#     print(v_wave[:, 0, 0, 0, 0])

#     # Write data to file
#     write_data(u_wave, z_wave, v_wave, lons, lats, press, times, waves)

# end_time = time.time()
# elapsed_time = end_time - start_time
# print(f"Waves calculated for all members in {elapsed_time:.2f} seconds.")
# #########------------------------------------------------###################