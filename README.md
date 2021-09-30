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

These modules are used to convert storm tracks (as provided by the [IBTrACS dataset](https://www.ncdc.noaa.gov/ibtracs/)) to 2D windfields, including preprocessing steps to select relevant tracks, filter out missing values etc. \
After preprocessing, storm tracks are saved as NumPy array files (.npy) per variable needed for windfield conversion (storm name, latitude, longitude, basin, pressure, Rmax, timestep, windspeed and year). After conversion, the 2D windfields are saved as GeoTIFF files per storm. \

### **Data conversion**
* convert_main.py
* convert_functions.py

These modules can be used to convert storm hazard files given in NetCDF format to GeoTIFF files and to match the projection and resolution of these files to the projection and resolution of the exposure files, which should be provided as GeoTIFF files. \

### **Damage computation**
* damage_main.py
* damage_functions.py

These modules form the actual damage computation per country for a given vulnerability function and climate scenario.\
Input arguments are: *vulnerability function*, *country*, *climate scenario*. \
There are two vulnerability functions implemented: "sigmoid" or "power". The implementation of these functions can be found in damage_functions.py . \
The second argument specifies the country for which damage is computed. The current implementation allows for damages to be aggregated on national or regional levels, which should be specified in the damage_main.py file. \
Finally, the third input argument of this file specifies the climate scenario for which damages are calculated: "current" or "future". \


## Full model
In order to run the entire model for historical storm tracks, scripts should be executed in the following order:
1. create_2D_windfields.py
2. convert_main.py
3. damage_main.py

To run the model for hazards given as NetCDF files, only the second and third step should be executed.


### Results
The output of the damage computation is given as CSV files with estimated damage per region (or country if damage is computed on national level) per return period. From this, risk can be computed as the weighted sum of damage for all storms, where each storm is weighted with a probability of 1/RP (return period).
