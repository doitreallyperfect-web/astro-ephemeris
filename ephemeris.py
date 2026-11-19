#!/usr/bin/env python3
"""Build a tropical, geocentric daily ephemeris using Skyfield + JPL DE440s.

The report is for the next calendar day in Europe/Moscow. Longitudes are
apparent geocentric ecliptic longitudes of date. Exact major aspects are
located numerically inside that Moscow calendar day.
"""
from __future__ import annotations

import json
import math
import os
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
from scipy.optimize import brentq
from skyfield.api import Loader
from skyfield.framelib import ecliptic_frame

MOSCOW = ZoneInfo("Europe/Moscow")
UTC = timezone.utc
OUT = Path("data/daily.json")
PLANETS = {
    "Sun": "sun", "Moon": "moon", "Mercury": "mercury",
    "Venus": "venus", "Mars": "mars barycenter", "Jupiter": "jupiter barycenter",
    "Saturn": "saturn barycenter", "Uranus": "uranus barycenter",
    "Neptune": "neptune barycenter", "Pluto": "pluto barycenter",
}
SIGNS = ["Овен", "Телец", "Близнецы", "Рак", "Лев", "Дева",
         "Весы", "Скорпион", "Стрелец", "Козерог", "Водолей", "Рыбы"]
ASPECTS = [(0, "соединение"), (60, "секстиль"), (90, "квадрат"),
           (120, "трин"), (180, "оппозиция")]


def wrap_delta(value: float) -> float:
    return (value + 180.0) % 360.0 - 180.0


def zodiac(longitude: float) -> tuple[str, float]:
    lon = longitude % 360.0
    idx = int(lon // 30)
    return SIGNS[idx], lon - idx * 30.0


def main() -> None:
    now = datetime.now(MOSCOW)
    target_date = date.fromisoformat(os.environ["TARGET_DATE"]) if os.getenv("TARGET_DATE") else now.date() + timedelta(days=1)
    start_local = datetime.combine(target_date, time.min, tzinfo=MOSCOW)
    end_local = start_local + timedelta(days=1)
    start_utc, end_utc = start_local.astimezone(UTC), end_local.astimezone(UTC)

    cache = Path(".skyfield-cache")
    cache.mkdir(exist_ok=True)
    load = Loader(str(cache))
    ts = load.timescale()
    eph = load("de440s.bsp")
    earth = eph["earth"]
    bodies = {name: eph[key] for name, key in PLANETS.items()}
    ecliptic = ecliptic_frame

    def instant(dt: datetime):
        return ts.from_datetime(dt)

    def longitude(body, dt: datetime) -> float:
        t = instant(dt)
        apparent = earth.at(t).observe(body).apparent()
        return float(apparent.frame_latlon(ecliptic)[1].degrees) % 360.0

    # A half-hour grid is refined with Brent's method to locate exact aspect times.
    sample_count = 49
    sample_dts = [start_utc + (end_utc - start_utc) * (i / (sample_count - 1))
                  for i in range(sample_count)]
    sample_hours = np.array([(dt - start_utc).total_seconds() / 3600 for dt in sample_dts])
    all_lons = {name: np.array([longitude(body, dt) for dt in sample_dts])
                for name, body in bodies.items()}

    noon_local = datetime.combine(target_date, time(12, 0), tzinfo=MOSCOW)
    positions = {}
    for name, body in bodies.items():
        lon = longitude(body, noon_local)
        before = longitude(body, noon_local - timedelta(hours=1))
        after = longitude(body, noon_local + timedelta(hours=1))
        speed = wrap_delta(after - before) / 2.0
        sign, degree = zodiac(lon)
        positions[name] = {
            "longitude": round(lon, 6), "sign": sign,
            "degree_in_sign": round(degree, 4),
            "daily_motion_deg_per_day": round(speed * 24, 6),
            "retrograde": speed < -0.0005,
        }

    aspects = []
    names = list(bodies)
    for i, first in enumerate(names):
        for second in names[i + 1:]:
            # Unwrap the relative longitude so it is continuous over the day.
            raw = (all_lons[second] - all_lons[first]) * math.pi / 180.0
            relative = np.unwrap(raw) * 180.0 / math.pi
            for angle, label in ASPECTS:
                targets = set()
                for val in relative:
                    center = round((val - angle) / 360)
                    for k in range(center - 1, center + 2):
                        targets.add(angle + 360 * k)
                    if angle not in (0, 180):
                        center2 = round((val + angle) / 360)
                        for k in range(center2 - 1, center2 + 2):
                            targets.add(-angle + 360 * k)
                    elif angle == 180:
                        center2 = round((val + 180) / 360)
                        for k in range(center2 - 1, center2 + 2):
                            targets.add(-180 + 360 * k)
                for target in sorted(targets):
                    fvals = relative - target
                    for j in range(len(fvals) - 1):
                        if fvals[j] == 0:
                            root_hour = sample_hours[j]
                        elif fvals[j] * fvals[j + 1] > 0:
                            continue
                        else:
                            def f(hour):
                                dt = start_utc + timedelta(hours=float(hour))
                                return wrap_delta(longitude(bodies[second], dt) -
                                                  longitude(bodies[first], dt) - target)
                            # Use the unwrapped target branch locally; wrap_delta is
                            # valid because each bracket is at most 30 minutes wide.
                            root_hour = brentq(f, sample_hours[j], sample_hours[j + 1])
                        event_dt = start_utc + timedelta(hours=float(root_hour))
                        if start_utc <= event_dt < end_utc:
                            # Avoid duplicates caused by equivalent 360-degree branches.
                            key = (first, second, label, round(event_dt.timestamp() / 60))
                            if not any(a["_key"] == key for a in aspects):
                                aspects.append({
                                    "_key": key, "planet_1": first, "planet_2": second,
                                    "aspect": label,
                                    "exact_time_msk": event_dt.astimezone(MOSCOW).isoformat(timespec="minutes"),
                                })

    aspects.sort(key=lambda x: x["exact_time_msk"])
    for item in aspects:
        item.pop("_key", None)

    moon_lon = positions["Moon"]["longitude"]
    sun_lon = positions["Sun"]["longitude"]
    elongation = (moon_lon - sun_lon) % 360.0
    moon_phase = "новолуние" if elongation < 45 or elongation >= 315 else (
        "первая четверть" if 45 <= elongation < 135 else (
        "полнолуние" if 135 <= elongation < 225 else "последняя четверть"))

    result = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "timezone": "Europe/Moscow",
        "date": target_date.isoformat(),
        "weekday": ["понедельник", "вторник", "среда", "четверг", "пятница", "суббота", "воскресенье"][target_date.weekday()],
        "calculation": {
            "source": "JPL DE440s via Skyfield",
            "coordinates": "geocentric apparent tropical ecliptic-of-date longitudes",
            "position_timestamp_msk": noon_local.isoformat(timespec="minutes"),
            "aspect_window_msk": [start_local.isoformat(), end_local.isoformat()],
            "aspect_search_step_minutes": 30,
            "note": "Planetary positions are calculated values; interpretive text is not generated here."
        },
        "moon": {"phase_approximation": moon_phase, "elongation_from_sun_deg": round(elongation, 4)},
        "positions": positions,
        "exact_aspects": aspects,
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT} for {target_date}; {len(aspects)} exact aspects.")


if __name__ == "__main__":
    main()
