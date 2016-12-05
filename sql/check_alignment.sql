WITH base AS (
SELECT rid, rast
FROM public.ned1arcsec
WHERE ST_Intersects(
ST_MakeEnvelope(-111.00, 43.00, -108.00, 45.00, 4269),
ST_Envelope(ned1arcsec.rast) )
)
SELECT base.rid, ST_SameAlignment(base.rast, top.rast) AS same,
ST_NotSameAlignmentReason(base.rast, top.rast) AS reason
FROM base, (SELECT rast FROM base LIMIT 1) As top
WHERE NOT ST_SameAlignment(base.rast, top.rast);
