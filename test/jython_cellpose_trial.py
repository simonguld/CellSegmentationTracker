
## IMPORTS
import sys, os

from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.util import LogRecorder;
from fiji.plugin.trackmate.detection import DogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate.features.spot import SpotContrastAndSNRAnalyzerFactory, SpotAnalyzerFactory, \
    SpotFitEllipseAnalyzerFactory, SpotIntensityMultiCAnalyzerFactory, SpotMorphologyAnalyzerFactory, SpotShapeAnalyzerFactory
from fiji.plugin.trackmate.features.edges import EdgeSpeedAnalyzer, AbstractEdgeAnalyzer, DirectionalChangeAnalyzer, EdgeTargetAnalyzer, EdgeTimeLocationAnalyzer, EdgeAnalyzer
from fiji.plugin.trackmate.features.track import TrackDurationAnalyzer, TrackAnalyzer, TrackIndexAnalyzer, \
    TrackLocationAnalyzer, TrackSpeedStatisticsAnalyzer, TrackSpotQualityFeatureAnalyzer, AbstractTrackAnalyzer, TrackBranchingAnalyzer, TrackMotilityAnalyzer
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import TMUtils

from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.action import CaptureOverlayAction

from java.io import File
from ij import IJ
from datetime import datetime as dt

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')

### PATHS

file_paths = []
image_path = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648\\merged.tif"
input_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
output_directory = "C:\\Users\\Simon Andersen\\Documents\\Uni\\SummerProject\\valeriias mdck data for simon\\24.08.22_698x648"
model_directory = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0'
cellpose_python_filepath = 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\python.exe'

show_output = False


### FUNCTIONS
def get_imlist(path, format = '.jpg'):

    """
    returns a list of filenames for all png images in a directory
    """
    return [os.path.join(path,f) for f in os.listdir(path) if f.endswith(format)]



# Open image in IJ
imp = IJ.openImage(image_path)



#-------------------------
# Instantiate model object
#-------------------------

model = Model()

# Set logger
# Logger -> content will be saved in the XML file.
logger = LogRecorder( Logger.VOID_LOGGER )

logger.log( 'TrackMate-Cellpose analysis script\n' )

dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
logger.log( dt_string + '\n\n' )

model.setLogger(Logger.IJ_LOGGER)

#------------------------
# Prepare settings object
#------------------------

settings = Settings(imp)

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

# Add the analyzers for some spot features.
# Here we decide brutally to add all of them.
#net.imagej.patcher.LegacyInjector.preinit()


spot_analyzers = [ SpotContrastAndSNRAnalyzerFactory, SpotFitEllipseAnalyzerFactory, SpotIntensityMultiCAnalyzerFactory, SpotShapeAnalyzerFactory]
edge_analyzers = [ EdgeSpeedAnalyzer, DirectionalChangeAnalyzer, EdgeTargetAnalyzer, EdgeTimeLocationAnalyzer] #, AbstractEdgeAnalyzer]
track_analyzers = [ TrackDurationAnalyzer, TrackIndexAnalyzer, TrackLocationAnalyzer, TrackSpeedStatisticsAnalyzer, TrackSpotQualityFeatureAnalyzer, TrackBranchingAnalyzer, TrackMotilityAnalyzer] #AbstractTrackAnalyzer

for spot_analyzer, edge_analyzer, track_analyzer in zip(spot_analyzers, edge_analyzers, track_analyzers):
    settings.addSpotAnalyzerFactory(spot_analyzer())
    settings.addEdgeAnalyzer(edge_analyzer())
    settings.addTrackAnalyzer(track_analyzer())

#settings.addAllAnalyzers()


settings.initialSpotFilterValue = 1

print(str(settings))

#----------------------
# Instantiate trackmate
#----------------------

trackmate = TrackMate(model, settings)
trackmate.computeSpotFeatures( True )

#------------
# Execute all
#------------


ok = trackmate.checkInput()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

ok = trackmate.process()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))


#----------------
# Display results
#----------------

model.getLogger().log('Found ' + str(model.getTrackModel().nTracks(True)) + ' tracks.')

if 0:
    # The feature model, that stores edge and track features.
    fm = model.getFeatureModel()

    # Iterate over all the tracks that are visible.
    for id in model.getTrackModel().trackIDs(True):

        # Fetch the track feature from the feature model.
        v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
        model.getLogger().log('')
        model.getLogger().log('Track ' + str(id) + ': mean velocity = ' + str(v) + ' ' + model.getSpaceUnits() + '/' + model.getTimeUnits())

        # Get all the spots of the current track.
        track = model.getTrackModel().trackSpots(id)
        for spot in track:
            sid = spot.ID()
            # Fetch spot features directly from spot.
            # Note that for spots the feature values are not stored in the FeatureModel
            # object, but in the Spot object directly. This is an exception; for tracks
            # and edges, you have to query the feature model.
            x=spot.getFeature('POSITION_X')
            y=spot.getFeature('POSITION_Y')
            t=spot.getFeature('AREA')
            q=spot.getFeature('QUALITY')
            snr=spot.getFeature('SNR_CH1')
            mean=spot.getFeature('MEAN_INTENSITY_CH1')
            M = model.getLogger().log('\tspot ID = ' + str(sid) + ': x='+str(x)+', y='+str(y)+', t='+str(t)+', q='+str(q) + ', snr='+str(snr) + ', mean = ' + str(mean))

## EXPORT TO XML: 


outFile = File(output_directory, "exportModel.xml")

writer = TmXmlWriter(outFile) 
writer.appendModel(model)
writer.appendSettings(settings)
writer.writeToFile()

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


 
if 0:

    ## EXPORT TO XML: 


    outFile = File(output_directory, "exportModel.xml")

    writer = TmXmlWriter(outFile) 

    writer.appendModel(model)
    writer.appendSettings(settings)
    writer.writeToFile()
        
    