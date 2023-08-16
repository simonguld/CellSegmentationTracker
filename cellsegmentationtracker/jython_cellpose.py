from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import LogRecorder;
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory

from fiji.plugin.trackmate.util import TMUtils
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.action import CaptureOverlayAction
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel

from ij import IJ


from datetime import datetime as dt
import sys, os
import json 

reload(sys)
sys.setdefaultencoding('utf-8')

# ------------------------------------------------------
# 	LOAD SETTINGS FROM FILE.
# ------------------------------------------------------

current_path = os.getcwd()
file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'jython_dict.json')
with open(file_path, 'r') as file:
        param_dict = json.load(file)

# Load paths from settings file.
image_path = param_dict['IMG_PATH']
cellpose_python_filepath = param_dict['CELLPOSE_PYTHON_FILEPATH']
model_directory = param_dict['CUSTOM_MODEL_PATH']
output_path = param_dict['OUTPUT_PATH']

cellpose_dict = param_dict['CELLPOSE_DICT']
trackmate_dict = param_dict['TRACKMATE_DICT']




# Shall we display the results each time?
show_output = param_dict['SHOW_SEGMENTATION']


# ------------------------------------------------------
# 	ACTUAL CODE.
# ------------------------------------------------------

# Open image.
imp  = IJ.openImage( image_path)
dims = imp.getDimensions()

# Make sure dimensions are correct
imp.setDimensions( dims[ 2 ], dims[ 4 ], dims[ 3 ] )
cal = imp.getCalibration()

# Logger -> content will be saved in the XML file.
logger = LogRecorder( Logger.VOID_LOGGER )

logger.log( 'TrackMate-Cellpose analysis script\n' )

dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
logger.log( dt_string + '\n\n' )

#------------------------
# Prepare settings object
#------------------------

settings = Settings(imp)
setup = settings.toStringImageInfo() 

# Configure Cellpose default detector.

settings.detectorFactory = CellposeDetectorFactory()
settings.detectorSettings = cellpose_dict

if cellpose_dict['CELLPOSE_MODEL'] =='CYTO':
     cellpose_dict['CELLPOSE_MODEL'] = PretrainedModel.CYTO
elif cellpose_dict['CELLPOSE_MODEL'] =='CYTO2':
        cellpose_dict['CELLPOSE_MODEL'] = PretrainedModel.CYTO2
elif cellpose_dict['CELLPOSE_MODEL'] == 'NUCLEI':
        cellpose_dict['CELLPOSE_MODEL'] = PretrainedModel.NUCLEI
else:
        cellpose_dict['CELLPOSE_MODEL'] = PretrainedModel.CUSTOM

# Configure tracker
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()

settings.trackerSettings['LINKING_MAX_DISTANCE'] = trackmate_dict['LINKING_MAX_DISTANCE']
settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = trackmate_dict['GAP_CLOSING_MAX_DISTANCE']
settings.trackerSettings['MAX_FRAME_GAP'] = trackmate_dict['MAX_FRAME_GAP']
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = trackmate_dict['ALLOW_TRACK_SPLITTING']
settings.trackerSettings['ALLOW_TRACK_MERGING'] = trackmate_dict['ALLOW_TRACK_MERGING']

# Analyzers 
settings.addAllAnalyzers()


#-------------------
# Instantiate plugin
#-------------------

trackmate = TrackMate( settings )
trackmate.computeSpotFeatures( True )
trackmate.computeTrackFeatures( True )
trackmate.getModel().setLogger( logger )

#--------
# Process
#--------

ok = trackmate.checkInput()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

ok = trackmate.process()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

#----------------
# Save results
#----------------

saveFile = TMUtils.proposeTrackMateSaveFile( settings, logger )
writer = TmXmlWriter( saveFile, logger )
writer.appendLog( logger.toString() )
writer.appendModel( trackmate.getModel() )
writer.appendSettings( trackmate.getSettings() )
writer.writeToFile()
print( "Results saved to: " + saveFile.toString() + '\n' )

#----------------
# Display results
#----------------

if show_output:
    model = trackmate.getModel()
    selectionModel = SelectionModel( model )
    ds = DisplaySettings()
    ds = DisplaySettingsIO.readUserDefault()
    ds.spotDisplayedAsRoi = True
    displayer =  HyperStackDisplayer( model, selectionModel, imp, ds )
    displayer.render()
    displayer.refresh()

print("Tracking and segmentation completed.")

# ------------------------------------------------------
