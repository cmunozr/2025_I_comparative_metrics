import cdsapi
import os 

# This script downloads climatic covariates from the Copernicus Climate Data Store (CDS)
# for the Nordic region, specifically for the years 2008 to 2021.
dataset = "insitu-gridded-observations-nordic"
request = {
    "variable": [
        "precipitation",
        "mean_temperature",
        "maximum_temperature",
        "minimum_temperature"
    ],
    "product_type": ["consolidated"],
    "spatial_interpolation_method": ["type_1"],
    "year": [
        "2008", "2009", "2010",
        "2011", "2012", "2013",
        "2014", "2015", "2016",
        "2017", "2018", "2019",
        "2020", "2021"
    ],
    "month": [
        "05", "06", "07",
        "08"
    ],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    "version": ["25_03"]
}


output_directory = "D:/nordic_climate"
output_filename = "may_to_aug_2008-2021.zip"
full_path = os.path.join(output_directory, output_filename)

client = cdsapi.Client()
client.retrieve(dataset, request, full_path)