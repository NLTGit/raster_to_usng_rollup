### Requirements
1. ArcGIS 10.1+ w/ Arcpy  
2. Spatial Analyst License


### Notes
Tool development was performed and tested on a clean install of ArcGIS Desktop 10.6 and 10.7. The only other dependency was that 64 bit geoprocessing was required on win7x64 to create python toolboxes.

To determine if there are issues with python toolboxes in your current environment, open up ArcGIS desktop, open the catalog panel, right click on any folder, click "New" , and select "Python Toolbox". If an error message appears and no toolbox is created, or if a toolbox appears with a red X icon then you may need to install the arcgis 64 bit background geoprocessing. 

