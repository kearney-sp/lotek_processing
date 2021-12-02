"""
Lotek processing helper functions for converting datasets
"""
import xarray as xr
from rasterio import features
import json
import pandas as pd
import geopandas as gpd


def shp2mask(shp, xr_object, fill=0, dtype='int16', **kwargs):
    """
    Convert shapefile to dask array

    Parameters
    ----------
    shp :           (GeoDataFrame) shapefile to be converted
    xr_object :     (2D DataArray) dask array to be used as a template.
                    Should have coordinates 'y' and 'x' representing spatial coordinates
    fill :          (integer, float) value to fill grid cells without a value inherited from 'shp' (value type must match 'dtype')
    dtype :         (string) data type for output DataArray
    **kwargs :      (dictionary) additional arguments to be passed to rasterio.features.rasterize()

    Returns
    -------
    DataArray (2D) with coordinates matching 'xr_object'
    """
    raster = features.rasterize(shp, fill=fill, transform=xr_object.transform,
                                out_shape=xr_object.shape, dtype=dtype, **kwargs)
    return xr.DataArray(raster,
                        coords=(xr_object.coords['y'].values, xr_object.coords['x']),
                        dims=('y', 'x'))


def getFeatures(gdf):
    """
    Parse features from GeoDataFrame in such a manner that rasterio wants them

    Parameters
    ----------
    gdf :    (GeoDataFrame) shapefile/vector

    Returns
    -------
    List of coordinates?
    """
    return [i['geometry'] for i in json.loads(gdf.to_json())['features']]


def multi2single(gpdf):
    """
    Convert multi-polygon GeoDataFrame geometries into single-polygon geometries

    Parameters
    ----------
    gpdf :    (GeoDataFrame) shapefile/vector with multi-feature geometries

    Returns
    -------
    GeoDataFrame with single-feature geometries
    """
    gpdf_singlepoly = gpdf[gpdf.geometry.type == 'Polygon']
    gpdf_multipoly = gpdf[gpdf.geometry.type == 'MultiPolygon']

    for i, row in gpdf_multipoly.iterrows():
        Series_geometries = pd.Series(row.geometry)
        df = pd.concat([gpd.GeoDataFrame(row, crs=gpdf_multipoly.crs).T]*len(Series_geometries), ignore_index=True)
        df['geometry']  = Series_geometries
        gpdf_singlepoly = pd.concat([gpdf_singlepoly, df])

    gpdf_singlepoly.reset_index(inplace=True, drop=True)
    return gpdf_singlepoly