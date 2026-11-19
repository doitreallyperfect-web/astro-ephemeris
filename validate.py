#!/usr/bin/env python3
"""Cross-check JPL DE440s/Skyfield longitudes against independent Moshier."""
import json, os, subprocess, sys
from pathlib import Path
import swisseph as swe

ROOT=Path(__file__).resolve().parent
DAILY=ROOT/"data/daily.json"
OUT=ROOT/"data/validation.json"
DATES=["2026-11-20","2026-12-15","2027-01-15"]
BODIES={"Sun":swe.SUN,"Moon":swe.MOON,"Mercury":swe.MERCURY,"Venus":swe.VENUS,
"Mars":swe.MARS,"Jupiter":swe.JUPITER,"Saturn":swe.SATURN,"Uranus":swe.URANUS,
"Neptune":swe.NEPTUNE,"Pluto":swe.PLUTO}
SIGNS=["Овен","Телец","Близнецы","Рак","Лев","Дева","Весы","Скорпион","Стрелец","Козерог","Водолей","Рыбы"]
def diffsec(a,b): return (((a-b+180)%360)-180)*3600
def main():
    swe.set_ephe_path(None)
    flags=swe.FLG_MOSEPH|swe.FLG_SPEED
    reports=[]
    try:
        for day in DATES:
            env=os.environ.copy(); env["TARGET_DATE"]=day
            subprocess.run([sys.executable,"ephemeris.py"],cwd=ROOT,env=env,check=True,capture_output=True,text=True)
            data=json.loads(DAILY.read_text(encoding="utf-8"))
            y,m,d=map(int,day.split("-")); jd=swe.julday(y,m,d,9.0) # 12:00 Moscow
            comp={}
            for name,bid in BODIES.items():
                xx,ret=swe.calc_ut(jd,bid,flags)
                lon=data["positions"][name]["longitude"]
                comp[name]={"skyfield_jpl_longitude":round(lon,8),"moshier_longitude":round(xx[0]%360,8),
                "difference_arcseconds":round(diffsec(lon,xx[0]),3),"moshier_speed_deg_per_day":round(xx[3],8),
                "sign_moshier":SIGNS[int((xx[0]%360)//30)],"ephemeris_flags_returned":int(ret)}
            reports.append({"date":day,"comparisons":comp,
              "max_absolute_difference_arcseconds":round(max(abs(x["difference_arcseconds"]) for x in comp.values()),3)})
    finally:
        env=os.environ.copy(); env.pop("TARGET_DATE",None)
        subprocess.run([sys.executable,"ephemeris.py"],cwd=ROOT,env=env,check=True,capture_output=True,text=True)
    warnings=[]; failures=[]
    for r in reports:
        for n,x in r["comparisons"].items():
            d=abs(x["difference_arcseconds"])
            if d>120: failures.append(f'{r["date"]} {n}: {d:.1f} arcsec')
            elif d>30: warnings.append(f'{r["date"]} {n}: {d:.1f} arcsec')
    result={"validation_source":"Swiss Ephemeris Moshier (explicit MOSEPH flag)",
    "primary_source":"JPL DE440s via Skyfield","coordinate_comparison":"geocentric apparent tropical ecliptic longitude",
    "comparison_time":"12:00 Europe/Moscow","dates":reports,
    "max_absolute_difference_arcseconds":round(max(abs(x["difference_arcseconds"]) for r in reports for x in r["comparisons"].values()),3),
    "warning_threshold_arcseconds":30.0,"failure_threshold_arcseconds":120.0,
    "warnings":warnings,"failures":failures,"status":"FAIL" if failures else ("PASS_WITH_WARNINGS" if warnings else "PASS"),
    "limitations":["Checks planetary longitudes only, not exact aspect event times.",
    "DE440s uses planetary-system barycenters for several outer planets while Moshier returns planet centers, which can cause small physical offsets.",
    "Moshier is an independent analytical ephemeris, not Swiss compressed DE431 files."]}
    OUT.parent.mkdir(exist_ok=True); OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"max_abs_arcsec":result["max_absolute_difference_arcseconds"],"warnings":warnings,"failures":failures},ensure_ascii=False,indent=2))
    if failures: raise SystemExit(1)
if __name__=="__main__": main()
