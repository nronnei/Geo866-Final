SELECT write_raster(
  ST_AsGDALRaster(ST_Union(rast), 'Gtiff'),
  '/home/ronneini/sql_transfer/wy_srtm_raw.tif'
)
FROM srtm_global_100
WHERE ST_Intersects(
  rast,
  ST_MakeEnvelope(-111.00, 44.00, -109.00, 45.00, 4326)
);
