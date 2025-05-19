# realtime_waves

##### Adapted from Steve Woolnough's original `realtime_waves` python program.
* Streamlined and optimised for `xarray` by Elliot McKinnon-Gray. Uses the real-time methodology of Yang et al. (2021), based on the original methodology of Yang et al. (2003).

Steps are as follows:

1. Reads in data using the <strong>`read_data`</strong> function. Can be analysis, obs, model (or a combination of both). Data provided needs to contain $(u,v,Z)(t,p,\lambda,\phi)$ data. Amend lines 24-26 with the name of these fields in your data. UM (Met Office) default for these is `x_wind`, `y_wind`, `geopotential_height`.

2. Transforms $(u,z) \rightarrow (q,r)$ through the simple transformation $q=(Z g/c_e)+u,$ $r=(Z g/c_e)-u$ in function <strong>`uz_to_qr`</strong>. Here, $c_e$ is the constant of separation derived from the trapping scale $y_0$ calculated in the notebook below. This transformation is necessary for neat solutions to obtain wave-filtered $(u,v,Z)$. See Gill & Clark (1974), and Gill (1980) for details.

3. Filters out equatorial wave time and spatial scales using FFTs.

4. Projects the PCFs (Parabolic Cylinder Functions - see the Yang papers) onto $q$ and $r$, and then reverses the transformation back to $(u,Z)$ in <strong>`filt_project`</strong>. $v$ doesn't require undergoing any transformation for reasons explained in Yang et al. 2003.

5. Writes the output data to a NetCDF file in <strong>`write_data`</strong>. 

* The output is stored as `['u_wave','v_wave','z_wave'](wave_type,longitude,latitude,pressure)`. So for example to extract just Kelvin wave $u$ at the equator at 850hPa for all longitudes and times (e.g. to plot a Hovmoller), this would be extracted by the command in `xarray`: `uKW0N = ds['u_wave'].sel(wave_type='Kelvin',latitude=0,pressure=850)`.


<strong>References</strong><br>
<sub>
Gill, A. E. (1980) Some simple solutions for heat-induced tropical circulation. <em>QJRMS</em>, <strong>106</strong>, 447-462.<br>
Gill, A. E., and A. J. Clarke (1974) Wind-induced upwelling, coastal currents and sea-level changes. <em>Deep Sea Research</em>, <strong>21</strong>, 325-345.<br>
Yang, G.-Y., B. Hoskins, and J. Slingo (2003) Convectively coupled equatorial waves: A new methodology for identifying wave structures in observational data. <em>J. Atmos. Sci.</em>, <strong>60</strong>, 1637-1654.<br>
Yang, G.-Y., S. Ferrett, S. Woolnough, J. Methven, and C. Holloway (2021) Real-time identification of equatorial waves and evaluation of waves in global forecasts. <em>Weather and Forecasting</em>, <strong>36</strong>, 171-193.
</sub>
