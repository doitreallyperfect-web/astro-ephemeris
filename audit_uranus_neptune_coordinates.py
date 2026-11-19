#!/usr/bin/env python3
"""Same-instant coordinate audit for Uranus and Neptune, 2027-01-15."""
import json, re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zoneinfo import ZoneInfo

import swisseph as swe
from skyfield.api import Loader
from skyfield.framelib import ecliptic_frame

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"data/uranus_neptune_coordinate_audit.json"
UTC=timezone.utc
MOSCOW=ZoneInfo("Europe/Moscow")
DAY=datetime(2027,1,15,tzinfo=MOSCOW)
TIMES=[DAY+timedelta(hours=12,minutes=26,seconds=29),
       DAY+timedelta(hours=12,minutes=27,seconds=7),
       DAY+timedelta(hours=12,minutes=30,seconds=37)]
BODIES={"Uranus":(swe.URANUS,"799","7"),"Neptune":(swe.NEPTUNE,"899","8")}
FLAGS={
 "Moshier_default_apparent":swe.FLG_MOSEPH|swe.FLG_SPEED,
 "Moshier_no_nutation":swe.FLG_MOSEPH|swe.FLG_SPEED|swe.FLG_NONUT,
 "Moshier_J2000":swe.FLG_MOSEPH|swe.FLG_SPEED|swe.FLG_J2000,
 "Moshier_true_position":swe.FLG_MOSEPH|swe.FLG_SPEED|swe.FLG_TRUEPOS,
 "Moshier_no_aberration":swe.FLG_MOSEPH|swe.FLG_SPEED|swe.FLG_NOABERR,
}
def jd(dt):
 u=dt.astimezone(UTC)
 return swe.julday(u.year,u.month,u.day,u.hour+u.minute/60+u.second/3600+u.microsecond/3.6e9)
def wrap(x): return (x+180)%360-180
def interpolate(rows, instant):
 keys=sorted(rows)
 before=max(k for k in keys if k<=instant)
 after=min(k for k in keys if k>=instant)
 if before==after: return rows[before]
 fraction=(instant-before).total_seconds()/(after-before).total_seconds()
 lon1,lat1=rows[before]; lon2,lat2=rows[after]
 lon=(lon1+fraction*wrap(lon2-lon1))%360
 lat=lat1+fraction*(lat2-lat1)
 return lon,lat
def horizons(target):
 p={"format":"json","COMMAND":f"'{target}'","CENTER":"'500@399'","MAKE_EPHEM":"'YES'","EPHEM_TYPE":"'OBSERVER'",
 "START_TIME":"'2027-01-15 06:00'","STOP_TIME":"'2027-01-15 18:00'","STEP_SIZE":"'1 m'","QUANTITIES":"'31'",
 "ANG_FORMAT":"'DEG'","CSV_FORMAT":"'YES'","APPARENT":"'AIRLESS'","TIME_TYPE":"'UT'","CAL_FORMAT":"'CAL'","EXTRA_PREC":"'YES'"}
 with urlopen("https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(p),timeout=60) as r: data=json.loads(r.read().decode())
 if data.get("error"): raise RuntimeError(f"Horizons {target}: {data['error']}")
 txt=data["result"]
 if "$$SOE" not in txt: raise RuntimeError(f"Horizons {target}: no SOE: {txt[:800]}")
 table=txt.split("$$SOE",1)[1].split("$$EOE",1)[0]
 rows={}
 dr=re.compile(r"(\d{4}-[A-Za-z]{3}-\d{2}\s+\d{2}:\d{2}(?::\d{2}(?:\.\d*)?)?)")
 nr=re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")
 for line in table.splitlines():
  m=dr.search(line)
  if not m: continue
  stamp=m.group(1); dt=None
  for fmt in ("%Y-%b-%d %H:%M:%S.%f","%Y-%b-%d %H:%M:%S","%Y-%b-%d %H:%M"):
   try: dt=datetime.strptime(stamp,fmt).replace(tzinfo=UTC); break
   except ValueError: pass
  nums=nr.findall(line[m.end():])
  if dt and len(nums)>=2: rows[dt]=(float(nums[0])%360,float(nums[1]))
 if len(rows)<100: raise RuntimeError(f"Horizons {target}: parsed {len(rows)} rows")
 return rows,{"targetname":data.get("targetname"),"signature":data.get("signature")}
def main():
 swe.set_ephe_path(None)
 cache=Path(".skyfield-cache"); cache.mkdir(exist_ok=True)
 load=Loader(str(cache)); ts=load.timescale(); eph=load("de440s.bsp")
 earth=eph["earth"]
 hdata={}; metadata={}
 for name,(_,center,bary) in BODIES.items():
  for label,target in [("center",center),("barycenter",bary)]:
   hdata[(name,label)],metadata[f"{name}_{label}"]=horizons(target)
 results=[]
 for dt in TIMES:
  instant=dt.astimezone(UTC)
  t=ts.from_datetime(instant); entry={"time_msk":dt.isoformat(),"time_utc":instant.isoformat(),"bodies":{}}
  for name,(swe_body,_,_) in BODIES.items():
   sky_body=eph[f"{name.lower()} barycenter"]
   sky=earth.at(t).observe(sky_body).apparent().frame_latlon(ecliptic_frame)
   sky_lon=float(sky[1].degrees)%360; sky_lat=float(sky[0].degrees)
   h={}
   for label in ("center","barycenter"):
    lon,lat=interpolate(hdata[(name,label)],instant)
    h[label]={"longitude_deg":lon,"latitude_deg":lat}
   swiss={}
   for flag_name,flag in FLAGS.items():
    x,_=swe.calc_ut(jd(dt),swe_body,flag)
    swiss[flag_name]={"longitude_deg":x[0]%360,"latitude_deg":x[1],
      "delta_lon_vs_horizons_center_arcsec":wrap(x[0]-h["center"]["longitude_deg"])*3600,
      "delta_lon_vs_horizons_barycenter_arcsec":wrap(x[0]-h["barycenter"]["longitude_deg"])*3600}
   entry["bodies"][name]={
    "skyfield_de440s_barycenter":{"longitude_deg":sky_lon,"latitude_deg":sky_lat},
    "horizons":h,"swiss_ephemeris":swiss,
    "skyfield_minus_horizons_barycenter_arcsec":wrap(sky_lon-h["barycenter"]["longitude_deg"])*3600,
    "horizons_center_minus_barycenter_arcsec":wrap(h["center"]["longitude_deg"]-h["barycenter"]["longitude_deg"])*3600}
  results.append(entry)
 out={"date":"2027-01-15","coordinate_convention":"geocentric apparent tropical ecliptic-of-date; interpolated to identical UTC instants",
 "times":results,"horizons_target_metadata":metadata,
 "method_notes":["Horizons quantities 31; 1-minute sampling with linear interpolation to each requested second.",
 "Skyfield uses ecliptic_frame (true ecliptic and equinox of date).",
 "Swiss Ephemeris flag variants isolate nutation, J2000 frame, true-position and aberration effects.",
 "This diagnostic compares coordinates, not aspect-root times."]}
 OUT.parent.mkdir(exist_ok=True); OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
 print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=="__main__": main()
