#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''----------------------------------------------------------------------------
 Nazwa narzędzia: Czas Infiltracji 
 Nazwa źródła:    CzasInfiltracji.py
 Wersja:          ArcGIS 10.3.1
 Autor:           Filip Pietruszyński - f.pietruszynski@gmail.com
 Opis:            Narzędzie służące do obliczania czasu infiltracji
                  przez strefę aeracji na dowolnym obszarze za pomocą
                  wzoru Witczaka i Żurek, Bindemana lub Macioszczyka
----------------------------------------------------------------------------'''

''' Importowanie paczki ArcPy, modułu ArcGIS Spatial Analyst extension module
oraz klasy env '''

import arcpy
from arcpy import env
from arcpy.sa import *

''' Włączanie nadpisywania danych w środowisku pracy '''

arcpy.env.overwriteOutput = True

''' Sprawdzanie licencji rozszerzenia ArcGIS Spatial Analyst '''

class BladLicencji(Exception):
    pass

try:
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")

    else:
        raise BladLicencji

except BladLicencji:
    print("Rozszerzenie ArcGIS Spatial Analyst jest niedostępne!")

''' Definiowanie lokalnych zmiennych na podstawie
parametrów wprowadzonych przez użytkownika '''

formulaObl = arcpy.GetParameterAsText(0)
topoRaster = arcpy.GetParameterAsText(1)
ppwRaster = arcpy.GetParameterAsText(2)
wydzieleniaShp = arcpy.GetParameterAsText(3)
wilgObjField = arcpy.GetParameterAsText(4)
wskInfEfField = arcpy.GetParameterAsText(5)
wspolFiltrField = arcpy.GetParameterAsText(6)
porowEfekField = arcpy.GetParameterAsText(7)
opadShp = arcpy.GetParameterAsText(8)
opadField = arcpy.GetParameterAsText(9)
metodaInter = arcpy.GetParameterAsText(10)
numberPoints = arcpy.GetParameterAsText(11)
outputRaster = arcpy.GetParameterAsText(12)

''' Ustalanie wielkości komórki plików rastrowych, generowanych podczas
wykonywania skryptu, na podstawie najmniej rozdzielczego pliku rastrowego '''

cellsizeTopo = arcpy.GetRasterProperties_management(topoRaster, "CELLSIZEX")
cellsizePPW = arcpy.GetRasterProperties_management(ppwRaster, "CELLSIZEX")

if cellsizeTopo >= cellsizePPW:
    cellsize = cellsizeTopo
    
else:
    cellsize = cellsizePPW
    
''' Definiowanie funkcji obliczającej miąższość strefy aeracji oraz
zamieniającej jej ujemne wartości na wartość *0* '''
    
def miazszosc():
    
    miazszoscMinus = Minus(topoRaster, ppwRaster)
    miazszosc0 = Con(miazszoscMinus < 0, 0, miazszoscMinus)
    
    return miazszosc0
    
''' Definiowanie funkcji generującej pliki rastrowe na podstawie parametrów
charakteryzujących poszczególne wydzielenia litologiczne '''

def macierzeRastrowe(x):
    
    if x == "wilgObjRaster":
        
        wilgObjRaster = "wilgObjRaster"
        arcpy.FeatureToRaster_conversion(
                                         wydzieleniaShp,
                                         wilgObjField,
                                         wilgObjRaster,
                                         cellsize)
        
        return wilgObjRaster
    
    elif x == "wskInfEfRaster":
        
        wskInfEfRaster = "wskInfEfRaster"
        arcpy.FeatureToRaster_conversion(
                                         wydzieleniaShp,
                                         wskInfEfField,
                                         wskInfEfRaster,
                                         cellsize)
        
        return wskInfEfRaster
    
    elif x == "wspolFiltrRaster":
        
        wspolFiltrRaster = "wspolFiltrRaster"
        arcpy.FeatureToRaster_conversion(
                                         wydzieleniaShp,
                                         wspolFiltrField,
                                         wspolFiltrRaster,
                                         cellsize)
        
        return wspolFiltrRaster

    else:
        
        porowEfekRaster = "porowEfekRaster"
        arcpy.FeatureToRaster_conversion(
                                         wydzieleniaShp,
                                         porowEfekField,
                                         porowEfekRaster,
                                         cellsize)
        
        return porowEfekRaster

''' Definiowanie funkcji interpolującej dane punktowe przedstawiające roczną
wysokość opadu według metody wybranej przez użytkownika '''

def interpolacja():
    
    if metodaInter == "Spline":
        roczWysOpad = "roczWysOpadRaster"
        roczWysOpad = Spline(
                             opadShp,
                             opadField,
                             cellsize,
                             "REGULARIZED",
                             0.1)

    elif metodaInter == "IDW":
        roczWysOpad = "roczWysOpadRaster"
        roczWysOpad = Idw(
                          opadShp,
                          opadField,
                          cellsize,
                          2,
                          RadiusVariable(numberPoints))

    else:
        roczWysOpad = "roczWysOpadRaster"
        roczWysOpad = NaturalNeighbor(
                                      opadShp,
                                      opadField,
                                      cellsize)
        
    return roczWysOpad

''' Definiowanie funkcji generującej plik rastrowy przedstawiający
czas infiltracji w strefie aeracji za pomocą wzoru Witczaka i Żurek '''

def witczakZurek():

    czasInfiltracji = Divide(
                              Times(
                                    miazszosc(),
                                    macierzeRastrowe(
                                                     "wilgObjRaster")
                                    ),
                              Times(
                                    interpolacja(),
                                    macierzeRastrowe(
                                                     "wskInfEfRaster")
                                    ))
                                    
    return czasInfiltracji

''' Definiowanie funkcji generującej plik rastrowy przedstawiający
czas infiltracji w strefie aeracji za pomocą wzoru Bindemana '''

def bindeman():

    czasInfiltracji = Divide(
                             Times(
                                   miazszosc(),
                                   macierzeRastrowe(
                                                    "porowEfekRaster")
                                   ),
                             Times(
                                   Power(
                                         Times(
                                               interpolacja(),
                                               macierzeRastrowe(
                                                       "wskInfEfRaster")
                                               ),
                                   2),
                             macierzeRastrowe(
                                              "wspolFiltrRaster")
                             )
                      ** (1./3))
                             
    return czasInfiltracji

''' Definiowanie funkcji generującej plik rastrowy przedstawiający
czas infiltracji w strefie aeracji za pomocą wzoru Macioszczyka '''

def macioszczyk():

    czasInfiltracji = Divide(
                             Times(
                                   miazszosc(),
                                   macierzeRastrowe(
                                                    "wilgObjRaster")
                                   ),
                             Times(
                                   Power(
                                         Times(
                                               interpolacja(),
                                               macierzeRastrowe(
                                                       "wskInfEfRaster")
                                               ),
                                   2),
                             macierzeRastrowe(
                                              "wspolFiltrRaster")
                             )
                      ** (1./3))
                             
    return czasInfiltracji

''' Wykonywanie wcześniej zdefiniowanych funkcji w konfiguracji
ustalonej na podstawie parametrów wybranych przez użytkownika '''

if formulaObl == "Wzor Macioszczyka":
    macioszczyk().save(outputRaster)

elif formulaObl == "Wzor Bindemana":
    bindeman().save(outputRaster)

else:
    witczakZurek().save(outputRaster)
    
""" Koniec """