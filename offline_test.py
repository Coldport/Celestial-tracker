import tkinter as tk
from tkinter import ttk
from skyfield.api import load, Topos

def get_object_alt_az_offline(object_name, observer_lat, observer_lon, observer_elevation, observation_time):
    """
    Retrieve the Altitude and Azimuth of a celestial object (planet, star, or satellite) using Skyfield.
    
    :param object_name: Name of the object (e.g., 'Saturn', 'Sirius', 'ISS')
    :param observer_lat: Latitude of the observer in degrees
    :param observer_lon: Longitude of the observer in degrees
    :param observer_elevation: Elevation of the observer in meters
    :param observation_time: Skyfield time object (e.g., ts.now())
    :return: A dictionary with Altitude and Azimuth values in degrees
    """
    # Mapping of planet names to their barycenter IDs
    planet_barycenter_ids = {
        'mercury': 0,
        'venus': 1,
        'earth': 2,
        'mars': 3,
        'jupiter': 4,
        'saturn': 6,
        'uranus': 7,
        'neptune': 8,
        'pluto': 9,
    }
    
    # Load Skyfield ephemeris
    eph = load('de440s.bsp')
    
    # Check if object is a planet
    object_name_lower = object_name.lower()
    if object_name_lower in planet_barycenter_ids:
        # Planet: Use planet barycenter IDs
        planet_barycenter_id = planet_barycenter_ids[object_name_lower]
        planet = eph[planet_barycenter_id]
        
        # Define the observer's location (NYC hardcoded coordinates)
        observer_lat = 40.7725  # Latitude in degrees (NYC)
        observer_lon = -73.9301  # Longitude in degrees (NYC)
        observer_elevation = 10  # Elevation in meters (NYC)
        observer = Topos(latitude_degrees=observer_lat, longitude_degrees=observer_lon, elevation_m=observer_elevation)
        
        # Attach observer to Earth frame
        earth = eph['earth']
        observer_at_time = earth + observer
        
        # Compute the astrometric position of the planet relative to the observer
        astrometric = observer_at_time.at(observation_time).observe(planet)
        alt, az, _ = astrometric.apparent().altaz()
        
        return {"altitude": alt.degrees, "azimuth": az.degrees}
    
    else:
        raise ValueError(f"Object name '{object_name}' is not available in offline mode.")

def convert_to_degrees_from_north(azimuth):
    """Convert azimuth to degrees from North (0° is North, 90° is East, etc.)."""
    return (azimuth + 180) % 360

def convert_to_degrees_from_flat_ground(altitude):
    """Convert altitude to degrees from the flat ground (measuring from the horizon)."""
    return 90 - altitude

def show_result():
    """Fetch and display the celestial object’s altitude and azimuth based on user input."""
    object_name = object_name_combobox.get()
    
    # Load Skyfield timescale and define the observation time
    ts = load.timescale()
    observation_time = ts.now()  # You can change this if you need a specific time
    
    try:
        position = get_object_alt_az_offline(object_name, None, None, None, observation_time)
        alt = position['altitude']
        az = position['azimuth']
        
        azimuth_from_north = convert_to_degrees_from_north(az)
        altitude_from_flat_ground = convert_to_degrees_from_flat_ground(alt)
        
        result_label.config(text=f"At {observation_time.utc_iso()}, {object_name} is located at:\n"
                                f"Altitude: {altitude_from_flat_ground:.2f}° (from flat ground)\n"
                                f"Azimuth: {azimuth_from_north:.2f}° (from North)")
    except Exception as e:
        result_label.config(text=f"Error: {str(e)}")

    # Refresh every second
    window.after(1000, show_result)

# Create the GUI window
window = tk.Tk()
window.title("Celestial Tracker")

# Label and combobox for object selection
tk.Label(window, text="Select Object (e.g., 'Saturn', 'Sirius')").grid(row=0, column=0)
object_name_combobox = ttk.Combobox(window, values=["Saturn", "Mars", "Venus", "Jupiter", "Neptune", "Pluto", "Mercury", "Uranus"])
object_name_combobox.grid(row=0, column=1)
object_name_combobox.set("Saturn")  # Default selection

# Button to calculate result
calculate_button = tk.Button(window, text="Start Tracking", command=show_result)
calculate_button.grid(row=1, column=0, columnspan=2)

# Result label
result_label = tk.Label(window, text="Altitude and Azimuth will appear here.", justify="left")
result_label.grid(row=2, column=0, columnspan=2)

# Start the GUI event loop
window.mainloop()
