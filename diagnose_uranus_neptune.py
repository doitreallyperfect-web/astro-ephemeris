#!/usr/bin/env python3
"""Diagnose the Uranus-Neptune sextile timing discrepancy on 2027-01-15."""
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import swisseph as swe
from scipy.optimize import brentq
from skyfield.api import Loader
from skyfield.framelib import ecliptic_frame

MOSCOW=ZoneInfo("Europe/Moscow")
UTC=timezone.utc
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"data/uranus_neptune_diagnostic.json"
FLAGS=swe.FLG_MOSEPH|swe.FLG_SPEED

def wrap(x): return (x+180.0)%360.0-180.0
def jd(dt):
    u=dt.astimezone(UTC)
    hour=u.hour+u.minute/60+u.second/3600+u.microsecond/3.6e9
    return swe.julday(u.year,u.month,u.day,hour)

def main():
    target_date=datetime(2027,1,15,tzinfo=MOSCOW)
    start=target_date.astimezone(UTC)
    cache=Path(".skyfield-cache"); cache.mkdir(exist_ok=True)
    load=Loader(str(cache)); ts=load.timescale(); eph=load("de440s.bsp")
    earth=eph["earth"]; uranus=eph["uranus barycenter"]; neptune=eph["neptune barycenter"]

    def sky_lon(body,dt):
        t=ts.from_datetime(dt)
        return float(earth.at(t).observe(body).apparent().frame_latlon(ecliptic_frame)[1].degrees)%360
    def sky_state(dt):
        u=sky_lon(uranus,dt); n=sky_lon(neptune,dt)
        h=timedelta(minutes=30)
        u0=sky_lon(uranus,dt-h); u1=sky_lon(uranus,dt+h)
        n0=sky_lon(neptune,dt-h); n1=sky_lon(neptune,dt+h)
        rel=wrap(n-u)
        speed=wrap((n1-u1)-(n0-u0))*24.0 # degrees/day; samples are 1 hour apart
        return {"uranus_longitude_deg":u,"neptune_longitude_deg":n,
          "relative_longitude_deg":rel,"sextile_orb_deg":wrap(rel-60.0),
          "relative_speed_deg_per_day":speed}

    def swe_lon(body,dt):
        x,_=swe.calc_ut(jd(dt),body,FLAGS)
        return x[0]%360
    def swe_state(dt):
        u=swe_lon(swe.URANUS,dt); n=swe_lon(swe.NEPTUNE,dt)
        xu,_=swe.calc_ut(jd(dt),swe.URANUS,FLAGS)
        xn,_=swe.calc_ut(jd(dt),swe.NEPTUNE,FLAGS)
        return {"uranus_longitude_deg":u,"neptune_longitude_deg":n,
          "relative_longitude_deg":wrap(n-u),"sextile_orb_deg":wrap(n-u-60.0),
          "relative_speed_deg_per_day":xn[3]-xu[3]}

    # Search the entire Moscow day for sign-changing brackets separately
    # in each ephemeris and each sextile branch; never assume a fixed bracket.
    day_end=start+timedelta(days=1)
    def find_roots(func):
        grid=[start+timedelta(minutes=10*i) for i in range(145)]
        vals=[func(dt) for dt in grid]
        found=[]
        for i in range(len(grid)-1):
            if vals[i]==0:
                found.append(grid[i])
            elif vals[i]*vals[i+1] < 0:
                a=(grid[i]-start).total_seconds(); b=(grid[i+1]-start).total_seconds()
                root=brentq(lambda sec: func(start+timedelta(seconds=sec)),a,b,xtol=0.01)
                found.append(start+timedelta(seconds=root))
        return found
    target_date_report=target_date+timedelta(hours=12,minutes=26)
    candidates=[]
    for target in (60.0,-60.0):
        sky_roots=find_roots(lambda dt: wrap(sky_lon(neptune,dt)-sky_lon(uranus,dt)-target))
        swe_roots=find_roots(lambda dt: wrap(swe_lon(swe.NEPTUNE,dt)-swe_lon(swe.URANUS,dt)-target))
        for sky_dt in sky_roots:
            if swe_roots:
                swe_dt=min(swe_roots,key=lambda dt:abs((dt-sky_dt).total_seconds()))
                candidates.append((abs((sky_dt-target_date_report).total_seconds()),target,sky_dt,swe_dt))
    if not candidates:
        raise RuntimeError("No Uranus-Neptune sextile root found in the Moscow day")
    _,target,sky_dt,swe_dt=min(candidates,key=lambda x:x[0])
    reported_sky=target_date+timedelta(hours=12,minutes=26)
    reported_swe=target_date+timedelta(hours=12,minutes=30,seconds=36)
    points=[]
    for label,dt in [("JPL exact root",sky_dt),("Moshier exact root",swe_dt),
                     ("reported JPL minute",reported_sky),("reported Moshier time",reported_swe)]:
        points.append({"label":label,"time_msk":dt.isoformat(timespec="milliseconds"),
          "jpl_de440s":sky_state(dt),"moshier":swe_state(dt)})
    result={"date":"2027-01-15","aspect":"Uranus-Neptune sextile","signed_target_angle_deg":target,
      "primary_ephemeris_targets":{"Uranus":"uranus barycenter","Neptune":"neptune barycenter"},
      "independent_ephemeris":"Swiss Ephemeris Moshier planet centers",
      "jpl_exact_root_msk":sky_dt.isoformat(timespec="milliseconds"),
      "moshier_exact_root_msk":swe_dt.isoformat(timespec="milliseconds"),
      "root_difference_seconds":round(abs((sky_dt-swe_dt).total_seconds()),3),
      "points":points,
      "interpretation_note":"Compare the relative longitudes and relative speeds at identical timestamps. The primary uses planetary-system barycenters while Moshier calculates planet centers; this is a potential non-like-for-like comparison, not yet a proven cause."}
    OUT.parent.mkdir(exist_ok=True); OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__": main()
