#------------------------------- PCM4Foam project -------------------------------#
#Author
    #Ehsan Golab, SUT. All rights reserved.
    #Ehsan1996Golab@gmail.com

#---------------------------------------------------------------------------------#

# importing the requests library
import requests
import time
import numpy as np
import matplotlib.pyplot as plt



def geographicData(appid, city):
    # api-endpoint
    URL = "https://api.openweathermap.org/data/2.5/weather"
      
    # defining a params dict for the parameters to be sent to the API
    PARAMS = {'q':city, 'appid':appid}
      
    # sending get request and saving the response as response object
    r = requests.get(url = URL, params = PARAMS)
      
    # extracting latitude and longitude data
    lat = r.json()['coord']['lat']
    lon = r.json()['coord']['lon']
    return lat, lon




def irradianceData(appid, city, date):

    # extracting latitude and longitude data
    lat, lon = geographicData(appid, city)

    # api-endpoint
    URL = "https://demo.openweathermap.org/energy/1.0/solar/data"
      
    # defining a params dict for the parameters to be sent to the API
    PARAMS = {'lat':lat, 'lon':lon, 'date':date, 'appid':appid}
      
    # sending get request and saving the response as response object
    r = requests.get(url = URL, params = PARAMS)
      
    # extracting data in json format
    q = []
    for i in range(0, 24):
        q.append(r.json()['irradiance']['hourly'][i]['clear_sky']['dni'])

    return q




def convCoeff(length, speed):
    rho = 1.184
    kappa = 0.0262
    Cp = 1006
    mu = 1.844e-5
    nu = mu/rho
    alpha = kappa/(rho*Cp)
    Pr = nu/alpha
    Re = (length*speed)/nu

    if Re<5e5:
        Nu = 0.664*(Re**0.5)*(Pr**(1/3))

    else:
        Nu = 0.033*(Re**(4/5))*(Pr**(1/3))

    h = (Nu*kappa)/length
    return h




def makingTable(time, data, name):

    with open(name + ".txt", "w") as file:
        file.write("(" + "\n")

        for i in range(len(time)):
            file.write("    " + "(" + str(time[i]) + " " + str(data[i]) + ")" + "\n")

        file.write(")" + "\n")





appid = '162789149decd6e5e92320af3598515b'
city = 'London'
date='2023-03-30'
#q = irradianceData(appid, city, date)



time = []
s = []
for i in range(0, 24):
    time.append(i*3600)


q = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 445.72, 693.38, 811.18, 877.28, 914.23, 930.77, 930.11, 912.09, 873.2, 804.05, 679.69, 413.08, 0.0, 0.0, 0.0, 0.0]
T = [27.2, 30.5, 31.6, 32.7, 31.6, 31.6, 29.4, 30, 31.1, 33.3, 35.5, 36.6, 37.7, 38.3, 38.3, 39, 38.3, 37.2, 33.8, 32.2, 30, 28.8, 28, 27.5]
s = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3]
T = [(element + 273) for element in T]
h = [(convCoeff(0.2, element)) for element in s]


makingTable(time, q, "qovsTime")
makingTable(time, T, "TovsTime")
makingTable(time, h, "hovsTime")


