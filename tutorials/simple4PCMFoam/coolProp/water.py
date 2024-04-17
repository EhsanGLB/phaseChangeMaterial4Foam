#------------------------------- PCM4Foam project -------------------------------#
#Author
    #Ehsan Golab, SUT. All rights reserved.
    #Ehsan1996Golab@gmail.com

#---------------------------------------------------------------------------------#

import CoolProp.CoolProp as CP

temperature = 25  # Temperature in degree Celsius
pressure = 1.01325  # Pressure in bar

# Get thermal properties of water
rho = CP.PropsSI('D', 'T', temperature + 273.15, 'P', pressure * 1e5, 'Water')
kappa = CP.PropsSI('L', 'T', temperature + 273.15, 'P', pressure * 1e5, 'Water')
mu = CP.PropsSI('V', 'T', temperature + 273.15, 'P', pressure * 1e5, 'Water')
Cp = CP.PropsSI('C', 'T', temperature + 273.15, 'P', pressure * 1e5, 'Water')
print(rho, kappa, mu, Cp)
