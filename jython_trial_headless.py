import sys, os

from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import Logger
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

from java.io import File
from ij import IJ

### PATHS
output_directory = "C:\\Users\\Simon Andersen\\Projects\\Projects\\CellSegmentationTracking"



# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')


cmd_arg = getArgument()
print(cmd_arg)

# Get currently selected image
# imp = WindowManager.getCurrentImage()
imp = IJ.openImage('https://fiji.sc/samples/FakeTracks.tif')
#imp.show()


#-------------------------
# Instantiate model object
#-------------------------

model = Model()

# Set logger
model.setLogger(Logger.IJ_LOGGER)

#------------------------
# Prepare settings object
#------------------------

settings = Settings(imp)

# Configure detector
settings.detectorFactory = DogDetectorFactory()
settings.detectorSettings = {
    'DO_SUBPIXEL_LOCALIZATION' : True,
    'RADIUS' : 2.5,
    'TARGET_CHANNEL' : 1,
    'THRESHOLD' : 4.,
    'DO_MEDIAN_FILTERING' : False,
}

### CELLPOSE
if 0:
    settings.detectorFactory = CellposeDetectorFactory()
    settings.detectorSettings = {
        'TARGET_CHANNEL' : 0,
        'OPTIONAL_CHANNEL' : 0,
        'CELLPOSE_PYTHON_FILEPATH': 'C:\\Users\\Simon Andersen\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\python.exe',
        'CELLPOSE_MODEL': PretrainedModel.CYTO,
       # 'CELLPOSE_MODEL_FILEPATH': 'C:\\Users\\"Simon Andersen"\\miniconda3\\envs\\cellpose\\lib\\site-packages\\cellpose\\models\\cyto_0',#only use if custom
        'CELL_DIAMETER': 30.0,
        'USE_GPU': False,
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
     
print("you actually made it bro")