SELECT write_raster(
  ST_AsGDALRaster(ST_Union(rast), 'Gtiff'),
  '/home/ronneini/sql_transfer/ned_1as_nad83_raw.tif'
)
FROM ned1arcsec
WHERE ST_Intersects(
  rast,
  ST_MakeEnvelope(-111.00, 43.00, -108.00, 45.00, 4269)
);
