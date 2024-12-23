import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import astropy.units as u
from matplotlib import cm
from matplotlib import colors
from PIL import Image
import requests
from io import BytesIO

# Function to calculate altitude and azimuth of a celestial object
def get_altaz(observer_location, ra, dec):
    # Observer's location
    location = EarthLocation(lat=observer_location['latitude'], lon=observer_location['longitude'], height=observer_location['altitude'])
    # Time of observation
    time_now = Time.now()
    # Celestial object (RA/Dec in degrees)
    sky_object = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
    # Transform to AltAz
    altaz_frame = AltAz(obstime=time_now, location=location)
    altaz = sky_object.transform_to(altaz_frame)
    return altaz.alt.degree, altaz.az.degree

# Function to convert AltAz to a direction vector
def calculate_direction_vector(altitude, azimuth):
    alt_rad = np.radians(altitude)
    az_rad = np.radians(azimuth)
    x = np.cos(alt_rad) * np.sin(az_rad)
    y = np.cos(alt_rad) * np.cos(az_rad)
    z = np.sin(alt_rad)
    return np.array([x, y, z])

# Function to draw the Earth with a texture (local map)
def draw_earth_with_texture(ax, texture_image, radius=1):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # Map the texture to the 3D sphere
    ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=texture_image, linewidth=0, antialiased=False)

# Function to plot the observer's location and the line to the celestial object
def plot_location_and_object(ax, observer_lat, observer_lon, direction_vector, radius=1):
    # Convert latitude and longitude to 3D coordinates
    lat_rad = np.radians(observer_lat)
    lon_rad = np.radians(observer_lon)
    observer_x = radius * np.cos(lat_rad) * np.cos(lon_rad)
    observer_y = radius * np.cos(lat_rad) * np.sin(lon_rad)
    observer_z = radius * np.sin(lat_rad)

    ax.scatter(observer_x, observer_y, observer_z, color='red', label='Your Location')

    # Calculate the line endpoint in 3D space
    line_endpoint = np.array([observer_x, observer_y, observer_z]) + 2 * direction_vector
    ax.plot([observer_x, line_endpoint[0]], [observer_y, line_endpoint[1]], [observer_z, line_endpoint[2]], color='orange', label='Line to Celestial Object')

# Function to get a 5-mile map image from OpenStreetMap (or other sources)
def get_local_map_image(center_lat, center_lon, zoom=14, size='640x640', scale=2):
    # OpenStreetMap or similar service URL to get the map
    url = f'https://static-maps.yandex.ru/1.x/?ll={center_lon},{center_lat}&spn=0.1,0.1&size={size}&z={zoom}&l=map&scale={scale}'
    
    # Download the map image using requests
    response = requests.get(url)
    img = Image.open(BytesIO(response.content))
    
    return img

# Observer's location (Astoria, NYC) on Earth's surface
observer_location = {'latitude': 40.768, 'longitude': -73.918, 'altitude': 10}  # Latitude, Longitude, Altitude in meters
earth_radius = 1  # Normalize Earth's radius to 1 for the plot

# Celestial object's coordinates (RA/Dec in degrees) - Sirius as an example
ra = 101.287  # Right Ascension in degrees
dec = -16.716  # Declination in degrees

# Calculate altitude and azimuth
altitude, azimuth = get_altaz(observer_location, ra, dec)
print(f"Altitude: {altitude:.2f}°, Azimuth: {azimuth:.2f}°")

# Convert to direction vector
direction = calculate_direction_vector(altitude, azimuth)

# Get a local map image of the area (5-mile radius around Astoria, NYC)
local_map = get_local_map_image(observer_location['latitude'], observer_location['longitude'])

# Convert the image to an array suitable for plotting
texture_image = np.array(local_map)

# Create the 3D plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Draw Earth with the local map texture and plot the data
draw_earth_with_texture(ax, texture_image, radius=earth_radius)
plot_location_and_object(ax, observer_location['latitude'], observer_location['longitude'], direction, radius=earth_radius)

# Add cardinal directions to the plot
ax.text(earth_radius * 1.2, 0, 0, 'East', color='black', fontsize=12)
ax.text(-earth_radius * 1.2, 0, 0, 'West', color='black', fontsize=12)
ax.text(0, earth_radius * 1.2, 0, 'North', color='black', fontsize=12)
ax.text(0, -earth_radius * 1.2, 0, 'South', color='black', fontsize=12)

# Configure the plot
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

# Set aspect ratio to make the Earth look spherical
ax.set_box_aspect([1, 1, 1])

# Show the plot
plt.show()
