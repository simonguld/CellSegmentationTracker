from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import LogRecorder, TMUtils
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory


from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer


import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings, DisplaySettingsIO


from fiji.plugin.trackmate.action import CaptureOverlayAction
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory

from fiji.plugin.trackmate.features.spot import SpotContrastAndSNRAnalyzerFactory, SpotAnalyzerFactory, \
    SpotFitEllipseAnalyzerFactory, SpotIntensityMultiCAnalyzerFactory, SpotMorphologyAnalyzerFactory, SpotShapeAnalyzerFactory, AbstractSpotFeatureAnalyzer
from fiji.plugin.trackmate.features.edges import EdgeSpeedAnalyzer, AbstractEdgeAnalyzer, DirectionalChangeAnalyzer, EdgeTargetAnalyzer, EdgeTimeLocationAnalyzer, EdgeAnalyzer
from fiji.plugin.trackmate.features.track import TrackDurationAnalyzer, TrackAnalyzer, TrackIndexAnalyzer, \
    TrackLocationAnalyzer, TrackSpeedStatisticsAnalyzer, TrackSpotQualityFeatureAnalyzer, AbstractTrackAnalyzer, TrackBranchingAnalyzer, TrackMotilityAnalyzer

from ij import IJ


from datetime import datetime as dt
import sys
import os
import json 
import time

reload(sys)
sys.setdefaultencoding('utf-8')

# ------------------------------------------------------
# 	LOAD SETTINGS FROM FILE.
# ------------------------------------------------------

class_path = os.path.dirname(os.path.realpath(__file__))
file_path = os.path.join(class_path, 'jython_dict.json')
with open(file_path, 'r') as file:
        param_dict = json.load(file)

# Load paths from settings file.
image_path = param_dict['IMG_PATH']
cellpose_python_filepath = param_dict['CELLPOSE_PYTHON_FILEPATH']
output_path = param_dict['OUTPUT_PATH']

base_path = os.path.splitext(image_path)[0]
overlay_path = base_path + '_overlays.avi'

cellpose_dict = param_dict['CELLPOSE_DICT']
trackmate_dict = param_dict['TRACKMATE_DICT']


# Shall we display the results each time?
show_output = param_dict['SHOW_SEGMENTATION']


# ------------------------------------------------------
# 	ACTUAL CODE.
# ------------------------------------------------------

# Open image.
imp  = IJ.openImage(image_path)
dims = imp.getDimensions()

# Make sure dimensions are correct (Format must be: [x,y,c,z,t])
if dims[3] > dims[4]:
        imp.setDimensions(dims[ 2 ], dims[ 4 ], dims[ 3 ] )
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

# Add all analyzers.
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
    time.sleep(2)

    # capture overlay - RGB file
    image = trackmate.getSettings().imp
    capture = CaptureOverlayAction.capture(image, -1, imp.getNFrames(), logger)
    capture.setTitle("TracksOverlay")
    IJ.save(capture, overlay_path)
    capture.show()

    

    

print("Tracking and segmentation completed.")


        
## Generate txt file confirming that the analysis is completed.
with open(os.path.join(class_path, 'cst_analysis_completed.txt'), 'w') as file:
        file.write('Analysis completed')
        file.close()


# ------------------------------------------------------

