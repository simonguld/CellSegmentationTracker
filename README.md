

# Purpose of implementation

# Guides to using, notable methods, method list


# List of dependencies

* Cellpose (in virtual environmate)
* Fiji -> trackmate-cellpose
* Java 8, Jython
* whatever packages I'm using (in a nice format? Possibly to allow pip install ...)

# Installation guide


# Pretrained models

For hver model:
Billede af segmentering, info om billerne det er trænet på, pixelradius etc

For class:
Gennemgå alle input-argumenter

Gennemgå relevante attributes
Gennemgå alle metoder



Henvis til example_notebook


# CellSegmentationTracker Documentation

## Importing the class
CellSegmentationTracker can be imported as follows:

```
from cellsegmentationtracker import CellSegmentationTracker
```

```
class CellSegmentationTracker.CellSegmentationTracker(self, imagej_filepath, cellpose_python_filepath, image_folder = None, xml_path = None, output_folder = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = False, cellpose_dict = dict(), trackmate_dict = dict()):
```

```
 __init__(self, imagej_filepath, cellpose_python_filepath, image_folder = None, xml_path = None, output_folder = None, use_model = 'CYTO', custom_model_path = None,
show_segmentation = False, cellpose_dict = dict(), trackmate_dict = dict()):
```

**Parameters:**    * **imagej_filepath**(str) - 
                   * **cellpose_python_filepath**(str) - 
                   * **image_folder**(str, default = None) -


**Parameters:**    * **imagej_filepath**(str) - The file path to the ImageJ/Fiji executable. This is required for launching ImageJ macros.
                    * **cellpose_python_filepath**(str) - The file path to the Cellpose Python script or module. This script/module is used to perform cell segmentation.
                    * **image_folder**(str, default=None) - The folder containing input images for processing. If provided, the class will process all images within this folder.
                    * **xml_path**(str, default=None) - The path to an XML file containing additional configuration or metadata. This can be used to customize the processing behavior.
                    * **output_folder**(str, default=None) - The folder where output files will be saved. If not specified, a default output folder will be used.
                    * **use_model**(str, default='CYTO') - Specifies the Cellpose model to use for segmentation. Options include 'CYTO' for cytoplasmic segmentation and others as                             defined by Cellpose.
                    * **custom_model_path**(str, default=None) - If a custom Cellpose model is to be used, provide the path to the model here.
                    * **show_segmentation**(bool, default=False) - Determines whether to display the segmentation results interactively during processing.
                    * **cellpose_dict**(dict, default=dict()) - A dictionary containing additional parameters to pass to the Cellpose segmentation algorithm. Refer to the Cellpose                             documentation for supported parameters.
                    * **trackmate_dict**(dict, default=dict()) - A dictionary containing parameters for configuring TrackMate, an ImageJ plugin for tracking cells over time.


                
                    

**Methods: **

**Attributes:**
