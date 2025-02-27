# LTN-Detection

![title](https://github.com/user-attachments/assets/6a3af717-7988-401d-a0a2-94a1a077b9a0)

## Introduction

LTN-Detection is an **experimental** project to quantitatively identify plausible Low Traffic Neighbourhoods in the UK using open data.

Low Traffic Neighbourhoods refer to areas in which *filtered pereablity* and *traffic calming* are deployed to *reduce motorised through-traffic* in *residential areas* (definition sourced from [Transport for London](https://madeby.tfl.gov.uk/2020/12/15/low-traffic-neighbourhoods/)). 

To see an example output, try [clicking here](https://froguin99.github.io/newcastleLTNexample/)!

## Methodology

#### Data and tools

Processing is performed in Python, relieing predominalty on [GeoPandas](https://geopandas.org/en/stable/), [OSMnx](https://osmnx.readthedocs.io/en/stable/) and [NetworkX](https://networkx.org/). 

Data used in this project is accessed from [OpenStreetMap](https://www.openstreetmap.org/#map=10/51.6547/-4.0883) using OSMnx, [German Aerospace Centre (DLR)](https://geoservice.dlr.de/web/maps/eoc:guf:4326) using API calls and from the [Ordnance Survey](https://www.ordnancesurvey.co.uk/products/os-open-roads) through manual downloading of data. 

| OpenStreetMap | DLR | Ordnance Survey |
| :-------------------: | :----------: | :----------: |
| Landuse type, Street networks, Public transport routes, Geographic features  | Urban footprint      | Street networks      |

> **Note** you will need to download the OS Open Roads dataset manually and place it in the correct folder (link to folder) to run this project. 

#### Process

At a high level, processing is split into 5  sections:
1. Define neighbourhoods
2. Identify modal filtering density
3. Measure neighbourhood permiablity
4. Identify neighbourhood through traffic
5. Output results as webmaps

Examples of each of these components can be found in the "Examples" folder.

The aim of this processing is to provide a data-driven indication of how plausiablie a given neighbourhood is to be a Low Traffic Neighbourhood. Each neighbourhood is given scores based on modal filtering, neighbourhood permiamblity and through traffic which are combined to provide an overall score. Processing is designed to work at the Local Authority District level, and works best on cities and highly urbanised areas. However it can be applied to much smaller (or larger, depending on how powerful your computer is!) areas. 

## Installation

To run this reposetory, the best option currently is to clone the repository and create a python enviroment via the [ox_151.yaml](https://github.com/Froguin99/LTN-Detection/blob/2e366da18db97ec65d3cec681aeeb974a6e7e4f3/envs/ox_151.yaml) file. 

Code is designed to run in a Jupyter notebook, using vscode. 

As highlighted prior, you will also need to download the [OS Open Roads](https://www.ordnancesurvey.co.uk/products/os-open-roads) dataset and place into the `\data` folder. An alternative using only OpenStreetMap can be found in the `\code` folder which does not require this download, but please be aware that this code is less robust and may not be as up-to-date. 


## Usage

To obtain an output, you must first set the placename(s) in the `\data\placenames.txt` to locations which can be obtained using OSM [nominatim boundary search](https://nominatim.openstreetmap.org/ui/search.html). For example `Newcastle Upon Tyne, United Kingdom` or `London Borough of Lambeth`. Multiple locations should be sepearted with a new line. 

Next, open the file `ltn_scoring_3_12_3_mass_process.ipynb` and replace the variable `file_path = ` with the location of your placenames to be input. Then hit run and let the script do its magic! 

Example scripts for each component of the processing can be found in the `Examples` folder.

## Further information and Contributing

This GitHub repo is very much a living document, hence it can be quite messy at times! If you have any questions or need anything explaining, please do get in touch!

Contributions of any kind are more than welcome, including any tests, bug reports or ideas for improvements. Feel free to open an issue or new discussion on GitHub.

## Acknowledgements and Credits

This project is part of my PhD research at the Geospatial Systems CDT, Newcastle University. Project supervision is from Dr [@Craig-Robson](www.github.com/craig-robson) and Dr Alistair Ford.

This work is made possible by the Engineering and Physical Sciences Research Council (EPSRC)
of Great Britain (grant number EP/S023577/1).



