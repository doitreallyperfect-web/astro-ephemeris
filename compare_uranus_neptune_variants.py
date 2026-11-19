#!/usr/bin/env python3
"""Four-way like-for-like diagnostic of the 2027-01-15 Uranus/Neptune sextile.

Variants:
1) DE440s planetary-system barycenters (Skyfield)
2) Swiss Ephemeris Moshier planet centers
3) JPL Horizons planet centers (799, 899)
4) JPL Horizons planetary-system barycenters (7, 8)

All Horizons requests use geocentric apparent observer ecliptic longitude/latitude
(quantities 31), airless apparent coordinates, and UTC sampling.
"""
import csv, io, json, math, time, re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zoneinfo import ZoneInfo

import swisseph as swe
from scipy.optimize import brentq
from skyfield.api import Loader
from skyfield.framelib import ecliptic_frame

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"data/uranus_neptune_four_way.json"
MOSCOW=ZoneInfo("Europe/Moscow")
UTC=timezone.utc
FLAGS=swe.FLG_MOSEPH|swe.FLG_SPEED
START=datetime(2027,1,15,6,0,tzinfo=UTC)
STOP=datetime(2027,1,15,18,0,tzinfo=UTC)
STEP=timedelta(minutes=1)

def wrap(x): return (x+180.0)%360.0-180.0
def jd(dt):
    u=dt.astimezone(UTC)
    return swe.julday(u.year,u.month,u.day,u.hour+u.minute/60+u.second/3600)

def moshier_lon(body,dt):
    x,_=swe.calc_ut(jd(dt),body,FLAGS)
    return x[0]%360.0

def horizons_series(target):
    params={
      "format":"json","COMMAND":f"'{target}'","CENTER":"'500@399'",
      "MAKE_EPHEM":"'YES'","EPHEM_TYPE":"'OBSERVER'",
      "START_TIME":"'2027-01-15 06:00'","STOP_TIME":"'2027-01-15 18:00'",
      "STEP_SIZE":"'1 m'","QUANTITIES":"'31'","ANG_FORMAT":"'DEG'",
      "CSV_FORMAT":"'YES'","APPARENT":"'AIRLESS'","TIME_TYPE":"'UT'",
      "CAL_FORMAT":"'CAL'","EXTRA_PREC":"'YES'"
    }
    url="https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(params)
    with urlopen(url,timeout=60) as response:
        payload=json.loads(response.read().decode("utf-8"))
    if payload.get("error"):
        raise RuntimeError(f"Horizons target {target}: {payload['error']}")
    text=payload.get("result","")
    if "$$SOE" not in text or "$$EOE" not in text:
        raise RuntimeError(f"Horizons target {target}: missing ephemeris delimiters; response: {text[:1200]}")
    table=text.split("$$SOE",1)[1].split("$$EOE",1)[0]
    rows=[]
    date_re=re.compile(r"(\\d{4}-[A-Za-z]{3}-\\d{2}\\s+\\d{2}:\\d{2}(?::\\d{2}(?:\\.\\d*)?)?)")
    number_re=re.compile(r"[-+]?(?:\\d+\\.\\d*|\\.\\d+|\\d+)(?:[Ee][-+]?\\d+)?")
    for line in table.splitlines():
        match=date_re.search(line)
        if not match:
            continue
        stamp=match.group(1)
        dt=None
        for fmt in ("%Y-%b-%d %H:%M:%S.%f","%Y-%b-%d %H:%M:%S","%Y-%b-%d %H:%M"):
            try:
                dt=datetime.strptime(stamp,fmt).replace(tzinfo=UTC)
                break
            except ValueError:
                pass
        if dt is None:
            continue
        nums=number_re.findall(line[match.end():])
        if len(nums)<2:
            continue
        rows.append((dt,float(nums[0])%360.0))
    rows=sorted({dt:lon for dt,lon in rows}.items())
    rows=[(dt,lon) for dt,lon in rows]
    if len(rows)<100:
        preview="\\n".join(table.splitlines()[:8])
        raise RuntimeError(f"Horizons target {target}: parsed only {len(rows)} rows; raw table preview: {preview}")
    return rows

def interpolated_root(u_rows,n_rows):
    nmap={t:lon for t,lon in n_rows}
    common=[(t,u,nmap[t]) for t,u in u_rows if t in nmap]
    if len(common)<100: raise RuntimeError("Horizons Uranus/Neptune time grids do not align")
    values=[wrap(n-u+60.0) for _,u,n in common] # target relative longitude -60
    for i in range(len(common)-1):
        if values[i]==0 or values[i]*values[i+1]<0:
            t0,u0,n0=common[i]; t1,u1,n1=common[i+1]
            # unwrap the signed angular residual locally and linearly interpolate
            y0=values[i]; y1=y0+wrap(values[i+1]-y0)
            frac=-y0/(y1-y0)
            return t0+(t1-t0)*frac
    raise RuntimeError("No -60 degree sextile root found in Horizons interval")

def root_from_callable(lon_u,lon_n):
    grid=[START+STEP*i for i in range(int((STOP-START)/STEP)+1)]
    vals=[wrap(lon_n(t)-lon_u(t)+60.0) for t in grid]
    for i in range(len(grid)-1):
        if vals[i]*vals[i+1]<0:
            sec0=(grid[i]-START).total_seconds()
            sec1=(grid[i+1]-START).total_seconds()
            root=brentq(lambda sec: wrap(lon_n(START+timedelta(seconds=sec))-lon_u(START+timedelta(seconds=sec))+60.0),sec0,sec1,xtol=.01)
            return START+timedelta(seconds=root)
    raise RuntimeError("No -60 degree sextile root found")

def main():
    swe.set_ephe_path(None)
    cache=Path(".skyfield-cache"); cache.mkdir(exist_ok=True)
    load=Loader(str(cache)); ts=load.timescale(); eph=load("de440s.bsp")
    earth=eph["earth"]; ubar=eph["uranus barycenter"]; nbar=eph["neptune barycenter"]
    def sky_lon(body,dt):
        return float(earth.at(ts.from_datetime(dt)).observe(body).apparent().frame_latlon(ecliptic_frame)[1].degrees)%360
    variants={}
    variants["DE440s_barycenters"]={
      "uranus":"uranus barycenter","neptune":"neptune barycenter",
      "root_utc":root_from_callable(lambda t:sky_lon(ubar,t),lambda t:sky_lon(nbar,t)).isoformat()}
    variants["Moshier_planet_centers"]={
      "uranus":"Swiss Ephemeris Uranus center","neptune":"Swiss Ephemeris Neptune center",
      "root_utc":root_from_callable(lambda t:moshier_lon(swe.URANUS,t),lambda t:moshier_lon(swe.NEPTUNE,t)).isoformat()}
    # JPL Horizons: 799/899 are planet centers; 7/8 are system barycenters.
    hu=horizons_series("799"); hn=horizons_series("899")
    hbar_u=horizons_series("7"); hbar_n=horizons_series("8")
    variants["Horizons_planet_centers"]={
      "uranus":"799 Uranus center","neptune":"899 Neptune center",
      "root_utc":interpolated_root(hu,hn).isoformat(),
      "sampling":"1 minute; linear interpolation between adjacent samples"}
    variants["Horizons_system_barycenters"]={
      "uranus":"7 Uranus-system barycenter","neptune":"8 Neptune-system barycenter",
      "root_utc":interpolated_root(hbar_u,hbar_n).isoformat(),
      "sampling":"1 minute; linear interpolation between adjacent samples"}
    # Convert roots to Moscow and compute pairwise offsets.
    for item in variants.values():
        dt=datetime.fromisoformat(item["root_utc"])
        item["root_msk"]=dt.astimezone(MOSCOW).isoformat(timespec="milliseconds")
    roots={k:datetime.fromisoformat(v["root_utc"]) for k,v in variants.items()}
    ref=roots["Horizons_planet_centers"]
    for k,item in variants.items():
        item["difference_from_Horizons_planet_centers_seconds"]=round((roots[k]-ref).total_seconds(),3)
    result={
      "date":"2027-01-15","aspect":"Uranus-Neptune sextile","signed_relative_angle_deg":-60,
      "coordinate_convention":"geocentric apparent tropical ecliptic-of-date longitude",
      "center":"Earth geocenter (500@399) for Horizons; Earth observer for Skyfield and Moshier",
      "horizons_settings":{"quantities":"31 (observer ecliptic longitude/latitude)","apparent":"airless","step":"1 minute","resolved_targets":HORIZONS_META},
      "variants":variants,
      "note":"DE440s contains planetary-system barycenter targets for these outer planets; it does not provide separate Uranus/Neptune center targets in this kernel. Horizons supplies center and barycenter targets independently. Moshier provides planet centers, not a separate system-barycenter solution."
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__": main()
