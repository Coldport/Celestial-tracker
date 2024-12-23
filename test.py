from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import astropy.units as u

# Define the location for NYC Astoria
nyc_location = EarthLocation(lat=40.7725 * u.deg, lon=-73.9301 * u.deg, height=10 * u.m)

def get_object_alt_az(ra, dec, observer_location, time):
    """
    Calculate the Altitude and Azimuth of a celestial object for a given observer location and time.

    :param ra: Right Ascension of the object in degrees
    :param dec: Declination of the object in degrees
    :param observer_location: An astropy.coordinates.EarthLocation object
    :param time: The time of observation as an astropy.time.Time object
    :return: A dictionary with Altitude and Azimuth values
    """
    # Create SkyCoord object for the celestial object
    sky_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

    # Define the AltAz frame for the observer's location
    altaz_frame = AltAz(location=observer_location, obstime=time)

    # Transform the celestial coordinates to AltAz
    altaz = sky_coord.transform_to(altaz_frame)

    return {
        "altitude": altaz.alt.deg,  # Altitude in degrees
        "azimuth": altaz.az.deg   # Azimuth in degrees
    }

# Example celestial object: Sirius (RA: 101.287, Dec: -16.716)
ra_sirius = 101.287
dec_sirius = -16.716

# Current time
current_time = Time.now()

# Get the Altitude and Azimuth for Sirius
sirius_position = get_object_alt_az(ra_sirius, dec_sirius, nyc_location, current_time)

# Output the results
print(f"At {current_time}, Sirius is located at:")
print(f"Altitude: {sirius_position['altitude']:.2f} degrees")
print(f"Azimuth: {sirius_position['azimuth']:.2f} degrees")
