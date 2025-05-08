#------------------------------- PCM4Foam project -------------------------------#
#Author
    #Ehsan Golab, SUT. All rights reserved.
    #Ehsan1996Golab@gmail.com

#---------------------------------------------------------------------------------#

# importing the requests library
import requests
import time



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



def currentWeatherData(appid, city):
    # api-endpoint
    URL = "https://api.openweathermap.org/data/2.5/weather"
      
    # defining a params dict for the parameters to be sent to the API
    PARAMS = {'q':city, 'appid':appid}
      
    # sending get request and saving the response as response object
    r = requests.get(url = URL, params = PARAMS)

    # extracting data in json format
    return r.json()



def currentIrradianceData(appid, city):

    # extracting latitude and longitude data
    lat, lon = geographicData(appid, city)

    # api-endpoint
    URL = "https://api.openweathermap.org/data/2.5/solar_radiation"
      
    # defining a params dict for the parameters to be sent to the API
    PARAMS = {'lat':lat, 'lon':lon, 'appid':appid}

    # sending get request and saving the response as response object
    r = requests.get(url = URL, params = PARAMS)

    return r.json()



while True:
    appid = '162789149decd6e5e92320af3598515b'
    city = 'London'
    cwd = currentWeatherData(appid, city)
    print(cwd['main']['temp'])
    print(cwd['main']['humidity'])
    print(cwd['wind']['speed'])
    irData = currentIrradianceData(appid, city)
    time.sleep(1)



