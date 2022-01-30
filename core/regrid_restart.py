import numpy as np
import xarray as xr
import xesmf as xe
import numpy as np

filename="/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF/input_data/GEOSChem.BoundaryConditions.20180101_0000z_2x2p5.nc4"

ds = xr.open_dataset(filename)

ds_out_dim = xr.Dataset(
	{
		"lat": (["lat"], np.concatenate([[-89.5],np.arange(-88.0,89.0, 2.0), [89.5]])),
		"lon": (["lon"], np.arange(-180.0,178.0, 2.5)),
	}
)

regridder = xe.Regridder(ds, ds_out_dim, "bilinear")

ds_out = regridder(ds)

ds_out.to_netcdf(filename)