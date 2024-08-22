import numpy as np

# Calculate Large Scale Clouds
def Large_Cloud_func(Pressure, Specific_Humidity, Relative_Humidity):
    """
    Calculate the large-scale cloud fraction based on pressure, specific humidity, and relative humidity.

    Parameters:
    Pressure (float or array): Atmospheric pressure in hPa.
    Specific_Humidity (float or array): Specific humidity in kg/kg.
    Relative_Humidity (float or array): Relative humidity in percentage (0-100%).

    Returns:
    float or array: Large-scale cloud fraction, adjusted for freeze-dry conditions.
    """
    a = find_a(Pressure)  # Determine the coefficient 'a' based on pressure.
    f = find_freeze_dry(Specific_Humidity, Pressure)  # Apply the freeze-dry adjustment.
    Large_Cloud = np.clip(a * (Relative_Humidity - 1.0) + 1.0, 0.0, 1.0)  # Calculate cloud fraction.
    return Large_Cloud * f  # Apply the freeze-dry factor.

def find_a(Pressure):
    """
    Calculate the coefficient 'a' used in cloud fraction calculation based on pressure.

    Parameters:
    Pressure (float or array): Atmospheric pressure in hPa.

    Returns:
    float or array: Coefficient 'a' for cloud fraction calculation.
    """
    a_t = 13  # Coefficient 'a' at the top of the atmosphere.
    a_s = 36  # Coefficient 'a' at the surface.
    n = 12  # Exponent for the pressure dependency.
    ps = 1000  # Surface pressure in hPa.
    a = a_t + (a_s - a_t) * np.exp(1 - (ps / Pressure)**n)  # Calculate 'a' based on pressure.
    return a

def find_freeze_dry(Specific_Humidity, Pressure):
    """
    Apply a freeze-dry adjustment to the cloud fraction based on specific humidity and pressure.

    Parameters:
    Specific_Humidity (float or array): Specific humidity in kg/kg.
    Pressure (float or array): Atmospheric pressure in hPa.

    Returns:
    float or array: Freeze-dry adjustment factor.
    """
    f = np.clip(Specific_Humidity / qv(Pressure), 0.15, 1.0)  # Calculate and clip the freeze-dry factor.
    return f

def qv(Pressure):
    """
    Calculate the threshold specific humidity for the freeze-dry adjustment.

    Parameters:
    Pressure (float or array): Atmospheric pressure in hPa.

    Returns:
    float or array: Threshold specific humidity in kg/kg.
    """
    q0 = 0.006  # Reference specific humidity at the surface.
    ps = 1000  # Surface pressure in hPa.
    n = 2.5  # Exponent for the pressure dependency.
    qv = q0 * (Pressure / ps)**n  # Calculate the threshold specific humidity.
    return qv

# Calculate Stratocumulus Clouds
def find_freeze_dry_stratocumulus(Specific_Humidity, Pressure):
    """
    Apply a freeze-dry adjustment specific to stratocumulus clouds based on specific humidity and pressure.

    Parameters:
    Specific_Humidity (float or array): Specific humidity in kg/kg.
    Pressure (float or array): Atmospheric pressure in hPa.

    Returns:
    float or array: Stratocumulus-specific freeze-dry adjustment factor.
    """
    qv_stratocumulus = 0.003  # Specific humidity threshold for stratocumulus clouds.
    fsc = np.clip(Specific_Humidity / qv_stratocumulus, 0.15, 1.0)  # Calculate and clip the freeze-dry factor.
    return fsc

def find_ELF(Pressure, Specific_Humidity, Relative_Humidity, Temperature):
    """
    Calculate the Estimated Low-level Cloud Fraction (ELF) for stratocumulus clouds.

    Parameters:
    Pressure (float or array): Atmospheric pressure in hPa.
    Specific_Humidity (float or array): Specific humidity in kg/kg.
    Relative_Humidity (float or array): Relative humidity in percentage (0-100%).
    Temperature (float or array): Atmospheric temperature in Kelvin.

    Returns:
    float or array: Estimated low-level cloud fraction (ELF).
    """
    delta_z = 2750  # Constant scale height for cloud formation.
    Zinv = Find_inversion_height(Pressure, Temperature)  # Find the inversion height.
    Zlcl = find_LCL(Temperature, Relative_Humidity)  # Calculate the lifting condensation level.
    fsc = find_freeze_dry_stratocumulus(Specific_Humidity, Pressure)  # Apply freeze-dry adjustment.
    sqrt_argument = np.maximum(Zinv * Zlcl, 0) # No sqrt of negative numbers.
    ELF = fsc * (1.0 - (np.sqrt(sqrt_argument) / delta_z))  # Calculate ELF based on inversion height and LCL.
    
    return ELF

def Find_inversion_height(Pressure, Temperature):
    """
    Determine the height of the temperature inversion layer based on pressure and temperature profiles.

    Parameters:
    Pressure (float or array): Atmospheric pressure in hPa.
    Temperature (float or array): Atmospheric temperature in Kelvin.

    Returns:
    float: Height of the inversion layer in meters.
    """
    p0 = 1013.25  # Reference pressure in hPa.
    Altitude = 44330 * (1 - (Pressure / p0)**(1 / 5.255))  # Calculate altitude from pressure.
    lapse_rate = find_lapse_rate(Temperature, Pressure)  # Calculate the lapse rate.

    most_negative = np.min(lapse_rate)  # Find the most negative lapse rate.
    index = np.where(lapse_rate == most_negative)[0]  # Find the index of the most negative lapse rate.
    Zinv = Altitude[index[0]] if len(index) > 0 else 0  # Determine the inversion height.

    return Zinv

def find_lapse_rate(Temperature, Pressure):
    """
    Calculate the lapse rate (rate of temperature decrease with height) based on temperature and pressure profiles.

    Parameters:
    Temperature (float or array): Atmospheric temperature in Kelvin.
    Pressure (float or array): Atmospheric pressure in hPa.

    Returns:
    float or array: Lapse rate in K/hPa.
    """
    R = 287.05  # Gas constant for dry air in J/(kg*K).
    cp = 1004  # Specific heat of air at constant pressure in J/(kg*K).
    const = R / cp  # Constant ratio of R to cp.
    theta = Temperature * (1000 / Pressure)**const  # Calculate potential temperature.
    lapse_rate = np.gradient(theta, Pressure)  # Calculate the lapse rate as the gradient of theta with pressure.
    return lapse_rate

def find_LCL(Temperature, Relative_Humidity):
    """
    Calculate the lifting condensation level (LCL) based on temperature and relative humidity.

    Parameters:
    Temperature (float or array): Atmospheric temperature in Kelvin.
    Relative_Humidity (float or array): Relative humidity in percentage (0-100%).

    Returns:
    float or array: Lifting condensation level in meters.
    """
    Td = find_dew_point(Temperature, Relative_Humidity)  # Calculate the dew point temperature.
    LCL = 125 * (Temperature - Td)  # Calculate LCL based on temperature difference.
    return LCL

def find_dew_point(Temperature, Relative_Humidity):
    """
    Calculate the dew point temperature based on ambient temperature and relative humidity.

    Parameters:
    Temperature (float or array): Atmospheric temperature in Kelvin.
    Relative_Humidity (float or array): Relative humidity in percentage (0-100%).

    Returns:
    float or array: Dew point temperature in Kelvin.
    """
    b = 17.625  # Constant used in the dew point calculation.
    c = 243.04  # Constant used in the dew point calculation.
    Temp_Celsius = Temperature - 273.15  # Convert temperature from Kelvin to Celsius.
    gamma = np.log(np.clip(Relative_Humidity / 100, 1e-10, 1.0)) + (b * Temp_Celsius) / (c + Temp_Celsius)  # Calculate gamma.
    Td = (c * gamma) / (b - gamma)  # Calculate the dew point temperature in Celsius.
    return Td + 273.15  # Return the dew point temperature in Kelvin.

# Calculate Stratocumulus Clouds
def Stratocumulus_func(Pressure, Specific_Humidity, Relative_Humidity, Temperature, Vertical_velocity):
    """
    Calculate the stratocumulus cloud fraction based on atmospheric conditions.

    Parameters:
    Pressure (float or array): Atmospheric pressure in hPa.
    Specific_Humidity (float or array): Specific humidity in kg/kg.
    Relative_Humidity (float or array): Relative humidity in percentage (0-100%).
    Temperature (float or array): Atmospheric temperature in Kelvin.
    Vertical_velocity (float or array): Vertical velocity in the atmosphere in Pa/s.

    Returns:
    float or array: Stratocumulus cloud fraction.
    """
    ELF = find_ELF(Pressure, Specific_Humidity, Relative_Humidity, Temperature)  # Calculate ELF.
    b = 1.3  # Coefficient for ELF.
    c = -0.1  # Constant for ELF.
    lapse_rate = find_lapse_rate(Temperature, Pressure)  # Calculate the lapse rate.
    most_negative = np.min(lapse_rate)  # Find the most negative lapse rate.
    index = np.where(lapse_rate == most_negative)[0]  # Find the index of the most negative lapse rate.
    
    if most_negative < -0.08 and Vertical_velocity[index[0]] > 0:
        Csc = np.clip(b * ELF + c, 0.0, 1.0)  # Calculate stratocumulus cloud fraction
    else:
        Csc = 0

    if np.shape(Csc) == ():
        Csc = np.zeros_like(Pressure)
    for i in range(len(Pressure)):
        if Pressure[i] < 750:
            Csc[i] = 0
    # Retain only the first non-zero stratocumulus cloud value for each vertical profile.
    Csc = retain_first_non_zero_1d(Csc)
    return Csc


def retain_first_non_zero_1d(arr):
    """
    Retain only the first non-zero element in a 1D array, setting all other elements to zero.

    Parameters:
    arr (numpy.ndarray): A 1D array from which only the first non-zero element is to be retained.

    Returns:
    numpy.ndarray: A 1D array where all elements are zero except the first non-zero element of the input array.
    """
    # Create an array of zeros with the same shape as the input array.
    result = np.zeros_like(arr)
    
    # Find the indices of all non-zero elements in the input array.
    non_zero_indices = np.nonzero(arr)[0]
    
    # If there is at least one non-zero element, retain the first one.
    if non_zero_indices.size > 0:
        result[non_zero_indices[0]] = arr[non_zero_indices[0]]
    
    return result

def C_total_func(Pressure, Specific_Humidity, Relative_Humidity, Temperature, Vertical_velocity):
    """
    Calculate the total cloud cover by combining large-scale clouds and stratocumulus clouds.

    Parameters:
    Pressure (numpy.ndarray): Atmospheric pressure in hPa.
    Specific_Humidity (numpy.ndarray): Specific humidity in kg/kg.
    Relative_Humidity (numpy.ndarray): Relative humidity in percentage (0-100%).
    Temperature (numpy.ndarray): Atmospheric temperature in Kelvin.
    Vertical_velocity (numpy.ndarray): Vertical velocity in the atmosphere in Pa/s.

    Returns:
    numpy.ndarray: Total cloud cover, which is the maximum of large-scale cloud cover and stratocumulus cloud cover.
    """
    # Calculate the large-scale cloud fraction.
    Large_Cloud = Large_Cloud_func(Pressure, Specific_Humidity, Relative_Humidity)
    
    # Calculate the stratocumulus cloud fraction.
    Stratocumulus = Stratocumulus_func(Pressure, Specific_Humidity, Relative_Humidity, Temperature, Vertical_velocity)
    #print(np.shape(Stratocumulus))
    # if np.shape(Stratocumulus) == ():
    #     Stratocumulus = np.zeros_like(Pressure)
    # for i in range(len(Pressure)):
    #     if Pressure[i] > 750:
    #         Stratocumulus[i] = 0
    # # Retain only the first non-zero stratocumulus cloud value for each vertical profile.
    # Stratocumulus = retain_first_non_zero_1d(Stratocumulus)

    # if Pressure > 750:
    #     Stratocumulus = 0
    # Calculate the total cloud cover by taking the maximum value between large-scale and stratocumulus clouds.
    C_total = np.maximum(Large_Cloud, Stratocumulus)
    
    return C_total
