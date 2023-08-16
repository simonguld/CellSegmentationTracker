

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




```
class CellSegmentationTracker:

    def __init__(self, imagej_filepath, cellpose_python_filepath, image_folder = None, xml_path = None, output_folder = None,
                  use_model = 'CYTO', custom_model_path = None, show_segmentation = False, cellpose_dict = dict(), trackmate_dict = dict(),):

        self.img_folder = image_folder
        self.img_path = None
        self.xml_path = xml_path
        self.imagej_filepath = imagej_filepath
        self.fiji_folder_path = os.path.abspath(os.path.join(self.imagej_filepath, os.pardir))
        self.cellpose_python_filepath = cellpose_python_filepath
        self.custom_model_path = custom_model_path
        self.output_folder = output_folder
        self.current_path = os.getcwd()

        self.show_segmentation = show_segmentation
        self.cellpose_dict = cellpose_dict
        self.trackmate_dict = trackmate_dict
        self.use_model = use_model
        self.jython_dict = {}

        self.cellpose_default_values = {
            'TARGET_CHANNEL' : 0,
            'OPTIONAL_CHANNEL_2': 0,
            'CELLPOSE_PYTHON_FILEPATH': self.cellpose_python_filepath,
            'CELLPOSE_MODEL': 'CYTO',
            'CELLPOSE_MODEL_FILEPATH': self.cellpose_python_filepath,
            'CELL_DIAMETER': 30.0,
            'USE_GPU': False,
            'SIMPLIFY_CONTOURS': True
            }
        self.trackmate_default_values = {'LINKING_MAX_DISTANCE': 15.0,
                                         'GAP_CLOSING_MAX_DISTANCE': 15.0,
                                         'MAX_FRAME_GAP': 2,
                                         'ALLOW_TRACK_SPLITTING': True,
                                         'ALLOW_TRACK_MERGING': True,
             }

        self.csv = None
        
        # If no xml file is provided (or present in image folder), cell segmentation, tracking and generation of an xml file is performed
        if self.xml_path is None:
            # If no image folder is provided, raise error
            if self.img_folder is None:
                raise ValueError("No image folder nor xml file provided!")

            if self.img_folder[-4:] == ".tif":
                self.img_path = self.img_folder
                self.img_folder = os.path.abspath(os.path.join(self.img_folder, os.pardir))
            # If more than one .tif file is provided, merge them into one
            else:
                # Make so that such a file has not already been created    
                self.img_path = os.path.join(self.img_folder, "merged.tif")
                if os.path.isfile(self.img_path):
                    os.remove(self.img_path)
                merge_tiff(self.img_folder, img_out_name = os.path.join(self.img_folder, "merged.tif"))
                print("Merged tif files into one file: ", self.img_path)

            # Ensure that xml file has not already been created
            if os.path.isfile(self.img_path.strip(".tif") + ".xml"):
                pass
            else:
                # Generate dictionary for jython script, if not already provided
                self.jython_dict.update({"IMG_PATH": self.img_path, "FIJI_FOLDER_PATH": self.fiji_folder_path, \
                                    "CELLPOSE_PYTHON_FILEPATH": self.cellpose_python_filepath, "CUSTOM_MODEL_PATH": self.custom_model_path,
                                    "OUTPUT_PATH": self.output_folder, "SHOW_SEGMENTATION": self.show_segmentation})
                
                # Specify Cellpose and TrackMate parameters not set by user
                if self.cellpose_dict == {}:
                    self.cellpose_dict = self.cellpose_default_values
                else:
                    for key in self.cellpose_default_values.keys():
                        if key not in self.cellpose_dict.keys():
                            self.cellpose_dict[key] = self.cellpose_default_values[key]
                if self.trackmate_dict == {}:
                    self.trackmate_dict = self.trackmate_default_values
                else:
                    for key in self.trackmate_default_values.keys():
                        if key not in self.trackmate_dict.keys():
                            self.trackmate_dict[key] = self.trackmate_default_values[key]

                self.jython_dict.update({'CELLPOSE_DICT': self.cellpose_dict})
                self.jython_dict.update({'TRACKMATE_DICT': self.trackmate_dict})


                # Save dictionary to json file
                with open(os.path.join(self.current_path,"jython_dict.json"), 'w') as fp:
                    json.dump(self.jython_dict, fp)


                # Run jython script
                os.chdir(self.fiji_folder_path)
                executable = list(os.path.split(self.imagej_filepath))[-1]
                jython_path = os.path.join(self.current_path, "jython_cellpose.py")

                if self.show_segmentation: 
                    command = executable + ' --ij2 --console --run "' + jython_path + '"'
                    os.system(command)
                else:
                    self.xml_path = self.img_path.strip(".tif") + ".xml"
                
                    # call jython script as a subprocess
                    pipe = Popen([executable,'--ij2','--headless', '--run', f"{jython_path}"], stdout=PIPE)
                    # wait for xml file to be created
                    try:
                        while not os.path.isfile(self.xml_path):
                            time.sleep(3)
                            continue
                    except: 
                        KeyboardInterrupt
                        sys.exit(0)
                    # kill subprocess
                    pipe.kill()
                    # print output
                    print(pipe.communicate()[0].decode('ascii'))
  
            # Set xml path to image path
            self.xml_path = self.img_path.strip(".tif") + ".xml"

        # Generate csv file from xml file
        self.xml = self.img_folder.strip(".tif") + ".xml"
  
        self.csv = trackmate_xml_to_csv(self.xml_path, get_tracks = True)

```
