# astro-ephemeris

Free daily planetary ephemeris calculations for the Three Keys Telegram forecast channel.

## What this repository does

- Uses **Skyfield** with the JPL **DE440s** planetary/lunar ephemeris.
- Calculates geocentric apparent tropical ecliptic longitudes for the Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune and Pluto.
- Calculates daily motion and a retrograde/direct flag.
- Searches the Moscow calendar day for exact conjunction, sextile, square, trine and opposition times.
- Writes the result to `data/daily.json`.

## Schedule

GitHub Actions runs daily at 16:00 UTC (19:00 Moscow), preparing the next Moscow calendar day's data ahead of the Telegram post scheduled for 21:00 Moscow.

A manual run is also available from **Actions → Build daily ephemeris → Run workflow**. Optionally enter a target date in `YYYY-MM-DD`; leave it blank to calculate tomorrow in Moscow.

## Data endpoint

After the first successful run, the latest JSON will be available at:

`https://raw.githubusercontent.com/doitreallyperfect-web/astro-ephemeris/main/data/daily.json`

This public raw URL can be read by Make's HTTP module. It is not yet connected to Make.

## Important limitations and validation

- The script reports astronomical calculations; it does not generate astrological interpretations.
- The Moon phase label is a broad phase category derived from Sun–Moon elongation, not a precise phase-event timestamp.
- Planet positions are sampled at 12:00 Moscow time; exact aspect events are searched throughout the Moscow calendar day.
- Before using this feed for public forecasts, compare selected dates and aspect times with an independent ephemeris such as Swiss Ephemeris/Astrodienst. Do not treat the workflow as production-validated until that comparison and a successful Actions run are complete.

## Run locally

```bash
python -m pip install -r requirements.txt
python ephemeris.py
```

To calculate a specific date:

```bash
TARGET_DATE=2026-11-20 python ephemeris.py
```
