"""
Adapted from Steve Woolnough's realtime_waves_v1.py
found on /ncas_climate_vol2/public/seasia/waves/
Uses the real-time methodology of Yang et al. (2021),
based on the original methodology of Yang et al. (2003).
"""
def read_data(filename, latmax=None):
    import xarray as xr
    import numpy as np
    """
    Read in u, v, z data from a NetCDF file using xarray.
    Here is where we define maximum latitude being 24N-S.
    """
    if latmax is None:
        latmax = 24
    print('reading in data')
    # Open the dataset using xarray
    fields = xr.open_dataset(filename)

    def get_var(ds, possible_names, label):
            for name in possible_names:
                if name in ds:
                    print(f"Found {label}: '{name}'")
                    return ds[name]
            print(f"None of {possible_names} found for {label}")
            raise KeyError(f"None of {possible_names} found in dataset for {label}")

    # Map variable aliases
    u = get_var(fields, ['u', 'x_wind'], 'u')
    v = get_var(fields, ['v', 'y_wind'], 'v')
    z = get_var(fields, ['z', 'ht', 'geopotential_height'], 'z')
    time = get_var(fields, ['t', 'time'], 'time')
    pressure = get_var(fields, ['pressure', 'level', 'pressure_level', 'p'], 'pressure')
    latitude = get_var(fields, ['latitude', 'lat'], 'latitude')
    longitude = get_var(fields, ['longitude', 'lon'], 'longitude')
    
    # Filter latitude range where abs(lat) <= latmax
    latrange =latitude.where(abs(latitude) <= latmax, drop=True)

    # Select the data within the latitude range
    u = u.sel(latitude=latrange)
    v = v.sel(latitude=latrange)
    z = z.sel(latitude=latrange)

    # Extract coordinates
    lons = longitude.values
    lats = latrange.values
    press = pressure.values
    times = time.values
    
    return u, v, z, lons, lats, press, times
    print(f'data {filename} successfully read')

def uz_to_qr(u,z,g_on_c) :
    """"
    transform u,z to q, r using q=z*(g/c) + u; r=z*(g/c) - u 
    """
    print('transforming u,z to q,r')
    q=z*g_on_c+u
    r=z*g_on_c-u

    return q,r
    print('transformation complete')

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
    print('PCF projection complete')

def write_data(u_wave, z_wave, v_wave, lons, lats, press, times, waves):
    import xarray as xr
    import numpy as np
    outfile = f'/gws/nopw/j04/kscale/USERS/emg/data/wave_data/KSE/KS_DW/waves_{model}_{season}.nc'
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

import time
import numpy as np

start_time = time.time()
model = 'analysis'
season = '20200120'
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
#waves=np.array(['Kelvin','WMRG']) # List of wave types to output

y0real= 2*np.pi*radea*y0/360.0   # convert trapping scale to metres
ce=2*y0real**2*beta
g_on_c=g/ce
c_on_g=ce/g


### read in 90 days of u,v,z data at 6 hourly time resolution and equatorial grid point
###  Methodology test and applied with 1 degree resolution data, but not dependent on it
### For Met Office operational forecasts the data would be 83 days of analysis and 7 days of forecast

#read data
u,v,z, lons,lats,press,times = read_data(f'/gws/nopw/j04/kscale/USERS/emg/data/wave_data/KSE/KS_DW/MOA_20191211_20190228.nc',latmax=24)

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
###-----------------------------------------------------------------

# ######### FOR LOOPING OVER ENSEMBLE MEMBERS ######################
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