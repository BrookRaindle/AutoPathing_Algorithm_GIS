import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QLineEdit, QPushButton, QColorDialog
from PyQt5.QtGui import QIcon, QColor
from PyQt5.QtCore import Qt
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors as mcolors
import numpy as np
import pandas as pd
import networkx as nx
from shapely.geometry import Point, LineString
import threading


class GeoPackageViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.start_location = None
        self.end_location = None 
        self.mode = "distance"
        self.total_time = 0  # Initialize total time
        self.total_distance = 0  # Initialize total distance
        self.average_fuel_economy = 0  # Initialize average fuel economy
        self.setWindowTitle("AutoPathing Algorithm")
        self.setGeometry(100, 100, 1500, 1000) # Set the window size
        self.ax = None # Initialize the subplot reference
    
        self.setStyleSheet("background-color: white; color: black;")# Set background color
        self.setupUI() # Setup the user interface

    def setupUI(self): # Setup the user interface
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#  Start/End Location inputs

        self.start_input = QLineEdit() # Create a QLineEdit widget for the start location
        self.start_input.setPlaceholderText("Start Location")
        layout.addWidget(self.start_input)

        self.end_input = QLineEdit() # Create a QLineEdit widget for the end location
        self.end_input.setPlaceholderText("End Location")
        layout.addWidget(self.end_input)
        
        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-# Reset button

        self.reset_button = QPushButton("Reset") # Create a QPushButton widget to reset the locations
        self.reset_button.clicked.connect(self.reset_locations)
        layout.addWidget(self.reset_button)
        # Map Visualization
        map_layout = QVBoxLayout() # Create a QVBoxLayout layout for the map visualization

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-# display figure for map

        self.figure = Figure() 
        self.canvas = FigureCanvas(self.figure) # Create a FigureCanvas widget to display the figure
        map_layout.addWidget(self.canvas) # Add the FigureCanvas widget to the layout

        self.addToolBar(Qt.BottomToolBarArea, NavigationToolbar(self.canvas, self))

        map_buttons_layout = QVBoxLayout()

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-# Optimisation buttons

        self.distance_button = QPushButton("Distance") # Create a QPushButton widget to calculate the shortest path distance
        self.distance_button.clicked.connect(self.calculate_shortest_path_distance)
        map_buttons_layout.addWidget(self.distance_button)

        self.time_button = QPushButton("Time") # Create a QPushButton widget to calculate the shortest path time
        self.time_button.clicked.connect(self.calculate_shortest_path_time)
        map_buttons_layout.addWidget(self.time_button)

        self.fuel_efficiency_button = QPushButton("Fuel Efficiency") # Create a QPushButton widget to calculate the shortest path fuel efficiency
        self.fuel_efficiency_button.clicked.connect(self.calculate_shortest_path_fuel)
        map_buttons_layout.addWidget(self.fuel_efficiency_button)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-# Load data and stuff

        map_layout.addLayout(map_buttons_layout)

        layout.addLayout(map_layout)

        self.load_geopackage_data()

        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

    def set_mode(self, mode):
        self.mode = mode

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def calculate_shortest_path_fuel(self):
        # Remove existing route line
        if self.start_location != None and self.end_location != None:
            self.reset_locations()
            self.replace_locations()
            # Calculate shortest path based on distance
            self.set_mode("fuel")
            self.calculate_shortest_path("fuel")
        else:
            self.set_mode("fuel")

    def calculate_shortest_path_distance(self):
        # Remove existing route line
        if self.start_location != None and self.end_location != None:
            self.reset_locations()
            self.replace_locations()
        # Calculate shortest path based on distance
            self.set_mode("distance")
            self.calculate_shortest_path("distance")
        else:
            self.set_mode("distance")

    def calculate_shortest_path_time(self):
        # Remove existing route line
        if self.start_location != None and self.end_location != None:
            self.reset_locations()
            self.replace_locations()
        # Calculate shortest path based on time
            self.set_mode("time")
            self.calculate_shortest_path("time")
        else:
            self.set_mode("time")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

    def load_geopackage_data(self):
        # Load GeoPackage data using GeoPandas ORDER EFFECTS RENDER ORDER
        gpkg_files = ['Clipped_Contour_2m.gpkg', 'clipped_waterlinks.gpkg', 'clipped_roadlinks.gpkg']
        colors = {'Clipped_Contour_2m.gpkg': 'green', 'clipped_roadlinks.gpkg': 'red', 'clipped_waterlinks.gpkg': 'blue'}
        
        merged_geometry = gpd.GeoDataFrame()
        ax = self.figure.add_subplot(111)  # Create a subplot on the figure
        self.ax = ax  # Store the subplot reference in the class variable
        for gpkg_file in gpkg_files:
            data = gpd.read_file(gpkg_file)  # Read the GeoPackage file using GeoPandas
            # Set the geometry column
            data.set_geometry('geometry', inplace=True)
            # Set color based on the file name
            color = colors.get(gpkg_file, 'color')  
            
            data.plot(ax=ax, color=color)  # Plot the data on the subplot

        self.figure.tight_layout()  # Adjust the layout of the figure
        self.figure.subplots_adjust(top=0.9)  # Adjust the spacing between subplots
        self.figure.suptitle("GeoPackage Visualization")  # Set the title of the figure

        # Connect click event to start and end location 
        ax.figure.canvas.mpl_connect('button_press_event', self.on_click)  # Connect the click event to the on_click method

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

    def calculate_shortest_path(self, weight): ### ONLY ROAD
        if self.start_location is None or self.end_location is None:
            print("Start and end locations are required.")
            return
        
        # Load clipped_roadlinks.gpkg
        road_links = gpd.read_file("clipped_roadlinks.gpkg")
        road_links = road_links.set_index('id')  # Set 'id' column as index
        #print(road_links.head())

        if road_links.crs != 'EPSG:27700':
            road_links = road_links.to_crs('EPSG:27700')                            ## EPSG SETTER

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

        # Create a NetworkX graph
        G = nx.Graph()

        node_coordinates = {
        'start': self.start_location,
        'end': self.end_location
        }
        

        # Add road links as edges to the graph
        for idx, road_link in road_links.iterrows():
            start_coord = (road_links.at[idx, 'geometry'].centroid.x, road_links.at[idx, 'geometry'].centroid.y)
            end_coord = (road_links.at[idx, 'geometry'].centroid.x, road_links.at[idx, 'geometry'].centroid.y)
            G.add_edge(road_link['start_node'], road_link['end_node'], id=road_link.name, length=road_link['length'], speed_limit=road_link['speed_limit'])
            node_coordinates[road_link['start_node']] = start_coord
            node_coordinates[road_link['end_node']] = end_coord

        # Add start and end points as nodes to the graph
        start_node = 'start'
        end_node = 'end'

        G.add_node(start_node, geometry=Point(self.start_location))
        G.add_node(end_node, geometry=Point(self.end_location))
        node_coordinates['start'] = (self.start_location[0], self.start_location[1])
        node_coordinates['end'] = (self.end_location[0], self.end_location[1])


        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

        # Find the nearest road link nodes to the start and end locations
        nearest_start_node = self.find_nearest_road_link(road_links, Point(self.start_location))
        nearest_end_node = self.find_nearest_road_link(road_links, Point(self.end_location))


        # Add edges connecting start and end points to the nearest road link nodes
        G.add_edge(start_node, nearest_start_node['start_node'], id=-1, length=nearest_start_node['distance'], speed_limit=0)
        G.add_edge(nearest_end_node['end_node'], end_node, id=-1, length=nearest_end_node['distance'], speed_limit=0)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


        # Calculate fuel efficiency for each edge based on speed limit

        for u, v, d in G.edges(data=True):
            length = d['length']
            speed_limit = d['speed_limit']
            if speed_limit == 0:
                d['speed_limit'] = 30
            else:
                d['fuel'] = self.calculate_fuel_consumption(speed_limit, length)

        # Calculate time for each edge based on speed limit

        for u, v, d in G.edges(data=True):
            length = d['length']
            speed_limit = d['speed_limit']
            if speed_limit == 0:
                d['speed_limit'] = 30
            else:
                d['time'] = length / (speed_limit * 0.44704) 


        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

        #print(G.nodes)
        print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
        #print(G.edges)

        shortest_path_distance = (nx.dijkstra_path(G, source=start_node, target=end_node, weight="length"))
        shortest_path_distance = [node for node in shortest_path_distance if node not in ['start', 'end']]
        print(len(shortest_path_distance))

        shortest_path_time = (nx.dijkstra_path(G, source=start_node, target=end_node, weight="time"))
        shortest_path_time = [node for node in shortest_path_time if node not in ['start', 'end']]
        print(len(shortest_path_time)) 

        shortest_path_fuel = (nx.dijkstra_path(G, source=start_node, target=end_node, weight="fuel"))
        shortest_path_fuel = [node for node in shortest_path_fuel if node not in ['start', 'end']]
        print(len(shortest_path_fuel)) 

        self.shortest_path_list = [shortest_path_distance, shortest_path_time, shortest_path_fuel]

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

        if weight == "distance":
            shortest_path = shortest_path_distance
        if weight == "time": 
            shortest_path = shortest_path_time
        if weight == "fuel":
            shortest_path = shortest_path_fuel


        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#    
        
        total_time = 0
        total_distance = 0
        total_fuel = 0
        for j in range(len(shortest_path) - 1):
            node_u = shortest_path[j]
            node_v = shortest_path[j + 1]
            total_distance += G[node_u][node_v]['length']
            total_time += G[node_u][node_v]['time'] / 60
            total_fuel += G[node_u][node_v]['fuel']
        total_distance_miles = total_distance * 0.000621371
        hours = int(total_time // 60) 
        minutes = int(total_time % 60)
        seconds = int((total_time % 1) * 60)
        print(f"Mode: {self.mode}  \nTime taken: {hours}:, {minutes}:, {seconds}\nDistance: {total_distance_miles} Mile \nFuel Consumed: {total_fuel} Gallons")

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

        for i in range(len(shortest_path) - 1):
            node_u = shortest_path[i]
            node_v = shortest_path[i + 1]
            x = [node_coordinates[node_u][0], node_coordinates[node_v][0]]
            y = [node_coordinates[node_u][1], node_coordinates[node_v][1]]
            self.ax.plot(x, y, color='yellow', linewidth=4)
        self.canvas.draw()

        self.plot_bar_graphs(self.shortest_path_list, G)

        #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def calculate_total_fuel_time_distance(self, list_of_paths, G):
        paths_info = {}  # Dictionary to store information about each path
        path_names = ["Shortest Path Distance", "Shortest Path Time", "Shortest Path Fuel"]
        for index, path in enumerate(list_of_paths):
            total_time = 0
            total_distance = 0
            total_fuel = 0

            for i in range(len(path) - 1):
                node_u = path[i]
                node_v = path[i + 1]
                total_distance += G[node_u][node_v]['length']
                total_time += G[node_u][node_v]['time'] / 60
                total_fuel += G[node_u][node_v]['fuel']

            total_distance_miles = total_distance * 0.000621371
            hours = int(total_time // 60) 
            minutes = int(total_time % 60)
            seconds = int((total_time % 1) * 60)

            paths_info[f"{path_names[index]}"] = {
                'total_time': total_time,
                'total_distance': total_distance_miles,
                'total_fuel': total_fuel
            }

        return paths_info

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

    def plot_bar_graphs(self, list_of_paths, G):
        paths_info = self.calculate_total_fuel_time_distance(list_of_paths, G)

        paths = list(paths_info.keys())  # Labels for paths
        total_fuel = [paths_info[path]['total_fuel'] for path in paths]
        total_time = [paths_info[path]['total_time'] for path in paths]
        total_distance = [paths_info[path]['total_distance'] for path in paths]

        # Plotting graphs
        fig, axs = plt.subplots(3, 1, figsize=(10, 15))  # Create subplots

        # Plot total fuel consumed
        bars = axs[0].bar(paths, total_fuel, color=['#CC4A6D', '#31AFBE', '#B08A2A'])
        axs[0].set_title('Total Fuel Consumed')
        axs[0].set_ylabel('Fuel (Gallons)')
        self.annotate_bars(axs[0], bars)

        # Plot total time taken
        bars = axs[1].bar(paths, total_time, color=['#BF4868', '#2DA6B4', '#A88528'])
        axs[1].set_title('Total Time Taken')
        axs[1].set_ylabel('Time (minutes)')
        self.annotate_bars(axs[1], bars)

        # Plot total distance
        bars = axs[2].bar(paths, total_distance, color=['#AD405E', '#289AA8', '#9C7B25'])
        axs[2].set_title('Total Distance')
        axs[2].set_ylabel('Distance (miles)')
        self.annotate_bars(axs[2], bars)
        plt.show()

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def annotate_bars(self, ax, bars): # adds labels near the top of each bar
        for bar in bars:
            height = bar.get_height()
            ax.annotate('{}'.format(round(height, 2)),
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, -30),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')


    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def plot_gpck(self, dir): # adds shortest path to display
        ax = self.ax
        data = gpd.read_file(dir)
        # Set the geometry column
        data.set_geometry('geometry', inplace=True)
        data.plot(ax=ax, color='yellow')

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
    

    def find_nearest_road_link(self, road_links, point):
        nearest_node = None
        min_distance = float('inf')
        for idx, road_link in road_links.iterrows():
            distance = road_link['geometry'].distance(point)
            if distance < min_distance:
                min_distance = distance
                nearest_node = {'start_node': road_link['start_node'],
                                 'end_node': road_link['end_node'], 'distance': distance}
        return nearest_node

#### Managing UI Methods

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def on_click(self, event):
        x = event.xdata
        y = event.ydata
        if x is not None and y is not None:
            if self.start_location is None:
                self.start_location = (x, y)
                self.start_input.setText(f"({x}, {y})")
            elif self.end_location is None:
                self.end_location = (x, y)
                self.end_input.setText(f"({x}, {y})")
                # Start the path calculation in a separate thread
                self.calculate_shortest_path(self.mode)
                
    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
    
    def remove_route_line(self):
        if self.route_line:
            self.route_line.remove()
            self.canvas.draw()
            self.route_line = None
        print("line removed successfully")

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def reset_locations(self):
        self.prev_start_location = self.start_location
        self.prev_end_location = self.end_location
        self.start_location = None
        self.end_location = None
        self.start_input.clear()
        self.end_input.clear()

        if hasattr(self, 'route_line'): # Remove the yellow line if it exists
            self.route_line.remove()
            del self.route_line
            self.canvas.draw()

        for line in self.ax.lines:
            line.remove()
        self.canvas.draw()

    def replace_locations(self):
        self.start_location = self.prev_start_location
        self.end_location = self.prev_end_location
        self.start_input.setText(f"({self.start_location[0]}, {self.start_location[1]})")
        self.end_input.setText(f"({self.end_location[0]}, {self.end_location[1]})")

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def update_route(self):
        if self.start_location is not None and self.end_location is not None:   
            self.route_line, = self.ax.plot([self.start_location[0], self.end_location[0]],
                                             [self.start_location[1], self.end_location[1]], color='yellow')
            self.canvas.draw()

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


    def calculate_fuel_efficiency(self, speed_mph):
        a = np.where(speed_mph < 45, 24 / 2025, 0.003)
        mpg = a * (speed_mph - 45)**2 + 36
        return mpg
    
    def calculate_fuel_consumption(self, speed, length):
        # speed is in miles per hour
        # fuel consumption rate is in gallons per mile
        
        fuel_consumption_rate_mpg = self.calculate_fuel_efficiency(speed)

        # Convert fuel consumption rate from gallons per mile to gallons per meter
        fuel_consumption_rate_meter_per_gallon = (fuel_consumption_rate_mpg * 1609.34)**-1
        fuel_consumption = length * fuel_consumption_rate_meter_per_gallon # Calculate fuel consumption
        #print(fuel_consumption)
        return fuel_consumption

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


def display_data(dir):
    # Load GeoPackage data using GeoPandas
    gpkg_file = dir
    road_links = gpd.read_file(gpkg_file)

    print(road_links.head())  # Print the first few rows of the GeoDataFrame
    print(road_links.columns)  # Print the column names 
    print(road_links['geometry'])
    # Creates a plot to visualize the road links
    ax = road_links.plot()
    ax.set_title("Road Links")
    plt.show()

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#


if __name__ == "__main__":

    #display_data("clipped_roadlinks.gpkg")
    app = QApplication(sys.argv)
    window = GeoPackageViewer()
    window.show()
    sys.exit(app.exec_())
