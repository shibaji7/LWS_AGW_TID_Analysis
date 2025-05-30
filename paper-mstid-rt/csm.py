# cesiumpy_polyline_rays_v4.py

import cesiumpy
# import math # Not needed if Rectangle constructor takes degrees

# Single HF ray data: (longitude, latitude, altitude_km)
polyline_rays_data = [
    [(-95, 39, 0), (-94, 40, 100), (-93, 42, 250), (-92, 45, 0)],
    [(-97, 38, 0), (-96, 39, 80), (-95.5, 40, 180), (-95, 41, 220), (-94, 43, 0)],
    [(-93, 37, 0), (-92.5, 37.5, 50), (-92, 38, 90), (-91, 40, 240), (-90.5, 41, 150), (-90, 42, 0)],
    [(-100, 40, 10), (-98, 40.5, 10), (-96, 40, 10), (-94, 39.5, 10)]
]

def visualize_polyline_rays(rays_data):
    viewer = cesiumpy.CesiumWidget(height="900px", width="100%")

    available_colors = [
        cesiumpy.color.DEEPSKYBLUE, cesiumpy.color.ORANGERED, cesiumpy.color.LIMEGREEN,
        cesiumpy.color.GOLD, cesiumpy.color.DARKVIOLET, cesiumpy.color.CYAN,
        cesiumpy.color.CRIMSON
    ]

    all_lons, all_lats = [], []

    for ray_index, single_ray_coords in enumerate(rays_data):
        cesium_points_for_current_polyline = []
        current_ray_color = available_colors[ray_index % len(available_colors)]
        max_alt_in_ray = max(p[2] for p in single_ray_coords if len(p) > 2) if single_ray_coords else 0

        for point_idx, (lon, lat, alt_km) in enumerate(single_ray_coords):
            alt_m = alt_km * 1000
            cartesian_point = cesiumpy.Cartesian3.fromDegrees(lon, lat, alt_m)
            cesium_points_for_current_polyline.append(cartesian_point)
            all_lons.append(lon)
            all_lats.append(lat)

            create_marker_flag = False
            point_description_html = (
                f"<b>Ray {ray_index + 1}, Point {point_idx + 1}</b><br>"
                f"Longitude: {lon:.2f}°<br>Latitude: {lat:.2f}°<br>Altitude: {alt_km:.1f} km"
            )
            
            marker_pixel_size = 7
            marker_color = current_ray_color.withAlpha(0.8)
            marker_outline_color = cesiumpy.color.BLACK
            marker_outline_width = 1
            marker_name = f"Ray {ray_index+1} Pt {point_idx+1}"

            if point_idx == 0:
                marker_name, marker_color, marker_pixel_size, create_marker_flag = f"Ray {ray_index+1} Start", current_ray_color, 10, True
            elif point_idx == len(single_ray_coords) - 1:
                marker_name, marker_color, marker_pixel_size, create_marker_flag = f"Ray {ray_index+1} End", current_ray_color, 10, True
            elif alt_km == 0:
                marker_name, marker_color, marker_pixel_size, create_marker_flag = f"Ray {ray_index+1} Bounce", cesiumpy.color.WHITE, 8, True
            elif alt_km == max_alt_in_ray and alt_km > 0:
                marker_name, marker_color, marker_pixel_size, create_marker_flag = f"Ray {ray_index+1} Apex", cesiumpy.color.YELLOW, 9, True

            if create_marker_flag:
                point_graphics_options = {
                    "color": marker_color, "pixelSize": marker_pixel_size,
                    "outlineColor": marker_outline_color, "outlineWidth": marker_outline_width,
                }
                # This structure seems to be what your cesiumpy version needs for Entity/Polyline
                try: 
                    marker_entity = cesiumpy.entities.entity.Entity(
                        position=cartesian_point, name=marker_name,
                        description=point_description_html, point=point_graphics_options
                    )
                    viewer.entities.add(marker_entity)
                except (AttributeError, TypeError): # Catch if .Entity is not there or .entity is a module
                    try:
                        marker_entity = cesiumpy.entities.entity.entity( # Try lowercase class name
                            position=cartesian_point, name=marker_name,
                            description=point_description_html, point=point_graphics_options
                        )
                        viewer.entities.add(marker_entity)
                    except Exception as e_innermost:
                         print(f"CRITICAL: Failed to create Entity. Error: {e_innermost}")

        if cesium_points_for_current_polyline:
            polyline_description_html = (
                f"<b>Polyline Ray {ray_index + 1}</b><br>{len(single_ray_coords)} points.<br>"
                f"Max Alt: {max_alt_in_ray:.1f} km"
            )
            try:
                current_polyline = cesiumpy.entities.polyline.Polyline(
                    positions=cesium_points_for_current_polyline, material=current_ray_color.withAlpha(0.75),
                    width=3, name=f"Polyline Ray {ray_index + 1}", description=polyline_description_html,
                )
                viewer.entities.add(current_polyline)
            except (AttributeError, TypeError):
                try:
                    current_polyline = cesiumpy.entities.polyline.polyline( # Try lowercase class name
                        positions=cesium_points_for_current_polyline, material=current_ray_color.withAlpha(0.75),
                        width=3, name=f"Polyline Ray {ray_index + 1}", description=polyline_description_html,
                    )
                    viewer.entities.add(current_polyline)
                except Exception as e_innermost_poly:
                    print(f"CRITICAL: Failed to create Polyline. Error: {e_innermost_poly}")

    if all_lons and all_lats:
        min_lon, max_lon = min(all_lons), max(all_lons)
        min_lat, max_lat = min(all_lats), max(all_lats)
        lon_buffer, lat_buffer = abs(max_lon - min_lon) * 0.15 + 2.0, abs(max_lat - min_lat) * 0.15 + 2.0
        west, south = max(-180.0, min_lon - lon_buffer), max(-90.0, min_lat - lat_buffer)
        east, north = min(180.0, max_lon + lon_buffer), min(90.0, max_lat + lat_buffer)
        if west > east: west, east = east, west
        
        # --- MODIFIED HERE ---
        # Use the Rectangle constructor directly with degree values
        try:
            bounding_rectangle = cesiumpy.Rectangle(west, south, east, north)
            viewer.camera.flyTo(destination=bounding_rectangle, duration=2.0)
        except Exception as e_rect:
            print(f"CRITICAL: Failed to create Rectangle or flyTo. Error: {e_rect}")
            # If direct constructor also fails or flyTo has issues, fallback to home view
            try:
                viewer.fly_home(duration=0)
            except AttributeError:
                print("WARNING: viewer.fly_home() not found. Camera may not reset.")

    else:
        try:
            viewer.fly_home(duration=0)
        except AttributeError:
            print("WARNING: viewer.fly_home() not found. Camera may not reset if no data.")
            pass

    return viewer

if __name__ == "__main__":
    cesium_visualization = visualize_polyline_rays(polyline_rays_data)
    output_filename = "cesium_polyline_rays_v2.html"
    html_string = cesium_visualization.to_html()
    with open(output_filename, "w") as file:
        file.write(html_string)
    print(f"Cesium visualization saved to: {output_filename}")
    print("Open this file in a web browser to view the visualization.")