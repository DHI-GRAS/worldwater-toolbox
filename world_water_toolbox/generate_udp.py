import json
from pathlib import Path

import openeo
from world_water_toolbox.wwt import generate_water_extent_udp, generate_water_indicators_udp

# Choose a single backend 
#backend = "openeo.dataspace.copernicus.eu"
backend = "openeofed.dataspace.copernicus.eu"
# backend = "openeo-dev.vito.be"
#backend = "openeo.vito.be"
connection = openeo.connect(backend)

udp = generate_water_extent_udp(connection)

with open(Path(__file__).parent / "single_month_extent_udp.json","w+") as f:
    json.dump(udp,f,indent=2)

indicators_udp = generate_water_indicators_udp(connection)
with open(Path(__file__).parent / "water_indicators.json", "w+") as f:
    json.dump(indicators_udp, f, indent=2)