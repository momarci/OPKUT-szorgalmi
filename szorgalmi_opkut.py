# Könyvtárak importálása: térinformatikai, gráfkezelő és vizualizációs eszközök
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import contextily as ctx
import osmnx as ox
import networkx as nx
from shapely.geometry import Point, LineString

# Bubi adatok és állomások beolvasása Excel fájlokból
bubi = pd.read_excel('C:/Users/molna/Documents/OPKUT/szorg/futar.xlsx')
allomas = pd.read_excel('C:/Users/molna/Documents/OPKUT/szorg/allomas.xlsx')

# Útvonal és állomás adatok GeoDataFrame-mé alakítása
gdf_ut = gpd.GeoDataFrame(bubi, geometry=gpd.points_from_xy(bubi.Szél, bubi.Hossz))
gdf_ut.set_crs(epsg=4326, inplace=True)

gdf_allomas = gpd.GeoDataFrame(allomas, geometry=gpd.points_from_xy(allomas.Szél, allomas.Hossz))
gdf_allomas.set_crs(epsg=4326, inplace=True)

# Hálózat lekérése az útvonal középpontja körüli 10 km-es körzetben
graph = ox.graph_from_point((gdf_ut.geometry.y.mean(), gdf_ut.geometry.x.mean()), dist=10000, network_type='drive')

# Eredeti útvonal távolságának kiszámítása a hálózaton
original_distance = 0
for i in range(len(gdf_ut) - 1):
    origin_node = ox.distance.nearest_nodes(graph, gdf_ut.iloc[i].geometry.x, gdf_ut.iloc[i].geometry.y)
    destination_node = ox.distance.nearest_nodes(graph, gdf_ut.iloc[i + 1].geometry.x, gdf_ut.iloc[i + 1].geometry.y)
    if nx.has_path(graph, origin_node, destination_node):
        original_distance += nx.shortest_path_length(graph, origin_node, destination_node, weight='length')
    else:
        print(f"No path found between original points {i} and {i + 1}.")

print(f"Original route distance: {original_distance:.2f} meters")

# Menetidő és átlagsebesség számítása
total_time_original = bubi['Idő'].sum()
print(f"Total time for the original route: {total_time_original:.2f} minutes")

original_distance_km = original_distance / 1000
total_time_original_hours = total_time_original / 60

if total_time_original_hours > 0:
    average_speed_original = original_distance_km / total_time_original_hours
    print(f"Average speed of the original route: {average_speed_original:.2f} km/h")
else:
    print("Cannot calculate average speed: Total time is zero.")
    average_speed_original = 0

# Útvonal újragenerálása pihenőállomások beiktatásával (30 percenként)
cumulative_time = 0
updated_points = []
last_inserted_index = -1

for i in range(len(gdf_ut) - 1):
    updated_points.append(gdf_ut.iloc[i])  # hozzáadjuk az aktuális pontot
    cumulative_time += gdf_ut.iloc[i]["Idő"]  # menetidő felhalmozása

    # Hálózati csomópontok lekérdezése a gráfból
    origin_node = ox.distance.nearest_nodes(graph, gdf_ut.iloc[i].geometry.x, gdf_ut.iloc[i].geometry.y)
    next_node = ox.distance.nearest_nodes(graph, gdf_ut.iloc[i + 1].geometry.x, gdf_ut.iloc[i + 1].geometry.y)
    
    if nx.has_path(graph, origin_node, next_node):
        distance_to_next = nx.shortest_path_length(graph, origin_node, next_node, weight='length')
    else:
        distance_to_next = 0
        print(f"No path found between points {i} and {i + 1} for time estimation.")

    # Időbecslés a következő ponthoz
    if average_speed_original > 0:
        time_to_next_hours = (distance_to_next / 1000) / average_speed_original
        time_to_next_minutes = time_to_next_hours * 60
    else:
        time_to_next_minutes = 0

    # Ha túlléptük a 30 percet, állomást keresünk beszúrásra
    if cumulative_time + time_to_next_minutes > 30:
        best_station = None
        min_added_distance = float('inf')

        for _, station in gdf_allomas.iterrows():
            station_node = ox.distance.nearest_nodes(graph, station.geometry.x, station.geometry.y)

            if nx.has_path(graph, origin_node, station_node) and nx.has_path(graph, station_node, next_node):
                distance_to_station = nx.shortest_path_length(graph, origin_node, station_node, weight='length')
                distance_from_station = nx.shortest_path_length(graph, station_node, next_node, weight='length')
                total_distance_with_station = distance_to_station + distance_from_station

                original_distance_segment = nx.shortest_path_length(graph, origin_node, next_node, weight='length')
                added_distance = total_distance_with_station - original_distance_segment

                # Legkisebb plusz távolságot eredményező állomás kiválasztása
                if added_distance < min_added_distance:
                    min_added_distance = added_distance
                    best_station = station

        if best_station is not None:
            updated_points.append(best_station)
            last_inserted_index = i
        else:
            print(f"No reachable station found for point {i}.")

        cumulative_time = 0  # idő újraindítása az állomás után

# Utolsó pont hozzáadása
updated_points.append(gdf_ut.iloc[-1])

# Frissített útvonal tárolása GeoDataFrame-ben
gdf_updated = gpd.GeoDataFrame(updated_points, crs=gdf_ut.crs)

# Beszúrt állomások külön kiszűrése
inserted_stations = []
for point in updated_points:
    if point.geometry in gdf_allomas.geometry.values:
        inserted_stations.append(point)

gdf_inserted_stations = gpd.GeoDataFrame(inserted_stations, crs=gdf_allomas.crs)

# Új útvonal szakaszok lekérdezése és vonalak létrehozása
street_distances = []
lines = []
for i in range(len(gdf_updated) - 1):
    origin_node = ox.distance.nearest_nodes(graph, gdf_updated.iloc[i].geometry.x, gdf_updated.iloc[i].geometry.y)
    destination_node = ox.distance.nearest_nodes(graph, gdf_updated.iloc[i + 1].geometry.x, gdf_updated.iloc[i + 1].geometry.y)

    if nx.has_path(graph, origin_node, destination_node):
        distance = nx.shortest_path_length(graph, origin_node, destination_node, weight='length')
        street_distances.append(distance)

        path_nodes = nx.shortest_path(graph, origin_node, destination_node, weight='length')
        path_coords = [(graph.nodes[node]['x'], graph.nodes[node]['y']) for node in path_nodes]

        if len(path_coords) > 1:
            lines.append(LineString(path_coords))
        else:
            print(f"Skipping segment {i} due to insufficient path coordinates.")
    else:
        print(f"No path found between points {i} and {i + 1}. Skipping segment.")

# Vonalas útvonal geoadat létrehozása
gdf_lines = gpd.GeoDataFrame(geometry=lines, crs="EPSG:4326")

# Teljes megtett távolság számítása
total_distance = sum(street_distances)
print(f"Total distance including stops: {total_distance:.2f} meters")

# Beszúrt állomások nevének listázása
inserted_station_names = []
for station in inserted_stations:
    inserted_station_names.append(station['Station'])

print("Inserted Station Names:", inserted_station_names)

# Eredeti útvonal megjelenítése
fig, ax = plt.subplots(figsize=(12, 12))
gdf_lines.plot(ax=ax, color='blue', linewidth=1, label='Original Route')
gdf_ut.plot(ax=ax, color='blue', markersize=10, label='Original Route Points')
gdf_ut.iloc[[0]].plot(ax=ax, color='red', markersize=50, label='Starting Point')
ctx.add_basemap(ax, crs=gdf_ut.crs, source=ctx.providers.CartoDB.Voyager)
ax.set_title('Original Route')
ax.legend()
plt.show()

# Frissített útvonal megjelenítése beszúrt állomásokkal
fig, ax = plt.subplots(figsize=(12, 12))
gdf_lines.plot(ax=ax, color='blue', linewidth=1, label='Updated Route')
gdf_inserted_stations.plot(ax=ax, color='green', markersize=20, label='Inserted Stations')
gdf_ut.plot(ax=ax, color='black', markersize=10, label='Original Route Points')
gdf_ut.iloc[[0]].plot(ax=ax, color='red', markersize=50, label='Starting Point')
ctx.add_basemap(ax, crs=gdf_updated.crs, source=ctx.providers.CartoDB.Voyager)
ax.set_title('Updated Route with Stations')
ax.legend()
plt.show()
