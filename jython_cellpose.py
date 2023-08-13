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
import sys 

reload(sys)
sys.setdefaultencoding('utf-8')

# ------------------------------------------------------
# 	EDIT FILE PATHS BELOW.
# ------------------------------------------------------

# Shall we display the results each time?
show_output = False

# Channel to process? 
channel_to_process = 1

# Image files to analyse.
file_paths = []
image_file = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\merged.tif"
input_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
model_directory = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0'
cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'
# ------------------------------------------------------
# 	ACTUAL CODE.
# ------------------------------------------------------

# Open image.
imp  = IJ.openImage( image_file)
dims = imp.getDimensions()

# Make sure dimensions are correct
print("DIM", dims)
imp.setDimensions( dims[ 2 ], dims[ 4 ], dims[ 3 ] )
print("DIM", dims)
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

settings.detectorSettings = {
        'TARGET_CHANNEL' : 0,
        'OPTIONAL_CHANNEL_2':0,
        'CELLPOSE_PYTHON_FILEPATH': cellpose_python_filepath,
        'CELLPOSE_MODEL': PretrainedModel.CYTO,
     'CELLPOSE_MODEL_FILEPATH': model_directory,#only use if custom
        'CELL_DIAMETER': 30.0,
        'USE_GPU': True,
        'SIMPLIFY_CONTOURS': True,
}	



# Configure tracker
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
settings.trackerSettings['LINKING_MAX_DISTANCE'] = 15.0
settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 15.0
settings.trackerSettings['MAX_FRAME_GAP'] = 2
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
settings.trackerSettings['ALLOW_TRACK_MERGING'] = True


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

    # capture overlay - RGB file
    image = trackmate.getSettings().imp
    capture = CaptureOverlayAction.capture(image, -1, imp.getNFrames(), logger)
    capture.setTitle("TracksOverlay")
    capture.show()

# ------------------------------------------------------

