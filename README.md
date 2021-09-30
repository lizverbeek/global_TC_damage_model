# Global Tropical Cyclone Risk model

Python implementation of the global TC risk model.

**Requirements:** [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/), [pandas](https://pandas.pydata.org/), [geopandas](http://geopandas.org/), [GDAL](https://gdal.org/), [netCDF4](https://unidata.github.io/netcdf4-python/), [xarray](http://xarray.pydata.org/en/stable/index.html#), [Fiona](https://pypi.org/project/Fiona/), [rasterio](https://rasterio.readthedocs.io/en/latest/#), [matplotlib](https://matplotlib.org/)
\
\
This project contains the following modules:

### **Storm track to windfield conversion**
* create_2D_windfields.py
* ibtracs_preprocessing
* holland_model.py
* windfield_functions.py

These modules are used to convert storm tracks (as provided by the [IBTrACS dataset](https://www.ncdc.noaa.gov/ibtracs/)) to 2D windfields, including preprocessing steps to select relevant tracks, filter out missing values etc. 
After preprocessing, storm tracks are saved as NumPy array files (.npy) per variable needed for windfield conversion (storm name, latitude, longitude, basin, pressure, Rmax, timestep, windspeed and year). After conversion, the 2D windfields are saved as GeoTIFF files per storm.
\
\
**Data conversion**
1. convert_main.py
2. convert_functions.py

These modules can be used to convert storm hazard files given in NetCDF format to GeoTIFF files and to match the projection and resolution of these files to the projection and resolution of the exposure files, which should be provided as GeoTIFF files.
\
\
**Damage computation**
1. damage_functions.py
2. damage_main.py

These modules form the actual damage computation per country for a given vulnerability function and climate scenario.
Input arguments are: 

