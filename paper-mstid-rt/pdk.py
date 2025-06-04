import pandas as pd
import pydeck as pdk

# Sample ray data: Replace with your actual PHaRLAP or simulation output
# Format: a list of rays, each ray is a list of (lon, lat, alt_km) points
hf_rays = [
    [(-95, 39, 0), (-94, 40, 100), (-93, 42, 250), (-92, 45, 0)],
    [(-97, 38, 0), (-96, 39, 80), (-95, 41, 220), (-94, 43, 0)],
    [(-93, 37, 0), (-92, 38, 90), (-91, 40, 240), (-90, 42, 0)],
]

# Flatten into a dataframe for pydeck
ray_segments = []
for ray in hf_rays:
    for i in range(len(ray) - 1):
        start = ray[i]
        end = ray[i + 1]
        ray_segments.append(
            {
                "start_lon": start[0],
                "start_lat": start[1],
                "start_alt": start[2] * 1000,  # meters
                "end_lon": end[0],
                "end_lat": end[1],
                "end_alt": end[2] * 1000,  # meters
            }
        )

df = pd.DataFrame(ray_segments)

# Define the pydeck layer
arc_layer = pdk.Layer(
    "LineLayer",
    data=df,
    get_source_position="[start_lon, start_lat, start_alt]",
    get_target_position="[end_lon, end_lat, end_alt]",
    get_width=2,
    get_color=[255, 140, 0],
    pickable=False,
    auto_highlight=False,
)

# Set the deck.gl view
view_state = pdk.ViewState(
    longitude=-95, latitude=40, zoom=3, min_zoom=2, max_zoom=15, pitch=60, bearing=0
)

# Render the deck
r = pdk.Deck(
    layers=[arc_layer],
    initial_view_state=view_state,
    map_style="https://basemaps.cartocdn.com/gl/positron-gl-style/style.json",
    tooltip={"text": "HF Ray Segment"},
)

r.to_html("hf_rays_3d.html")
