import json
from pathlib import Path

import openeo
from wwt import generate_water_extent_udp


udp = generate_water_extent_udp(openeo.connect("openeo.vito.be"))

with open(Path(__file__).parent / "single_month_extent_udp.json","w+") as f:

    json.dump(udp,f,indent=2)