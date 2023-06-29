import sys



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

from ij import IJ

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')


# Get currently selected image
# imp = WindowManager.getCurrentImage()
imp = IJ.openImage('https://fiji.sc/samples/FakeTracks.tif')
imp.show()


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
    'THRESHOLD' : 2.,
    'DO_MEDIAN_FILTERING' : False,
}

# Configure tracker
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
settings.trackerSettings['LINKING_MAX_DISTANCE'] = 10.0
settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 10.0
settings.trackerSettings['MAX_FRAME_GAP'] = 3

# Add the analyzers for some spot features.
# Here we decide brutally to add all of them.
#net.imagej.patcher.LegacyInjector.preinit()

if 0:
    from fiji.plugin.trackmate.features.spot import SpotContrastAndSNRAnalyzerFactory, SpotAnalyzerFactory, \
        SpotFitEllipseAnalyzerFactory, SpotIntensityMultiCAnalyzerFactory, SpotMorphologyAnalyzerFactory, SpotShapeAnalyzerFactory
    from fiji.plugin.trackmate.features.edges import EdgeSpeedAnalyzer, AbstractEdgeAnalyzer, DirectionalChangeAnalyzer, EdgeTargetAnalyzer, EdgeTimeLocationAnalyzer, EdgeAnalyzer
    from fiji.plugin.trackmate.features.track import TrackDurationAnalyzer, TrackAnalyzer, TrackIndexAnalyzer, \
        TrackLocationAnalyzer, TrackSpeedStatisticsAnalyzer, TrackSpotQualityFeatureAnalyzer, AbstractTrackAnalyzer, TrackBranchingAnalyzer, TrackMotilityAnalyzer

    spot_analyzers = [ SpotContrastAndSNRAnalyzerFactory, SpotFitEllipseAnalyzerFactory, SpotIntensityMultiCAnalyzerFactory, SpotShapeAnalyzerFactory]
    edge_analyzers = [ EdgeSpeedAnalyzer, DirectionalChangeAnalyzer, EdgeTargetAnalyzer, EdgeTimeLocationAnalyzer] #, AbstractEdgeAnalyzer]
    track_analyzers = [ TrackDurationAnalyzer, TrackIndexAnalyzer, TrackLocationAnalyzer, TrackSpeedStatisticsAnalyzer, TrackSpotQualityFeatureAnalyzer, TrackBranchingAnalyzer, TrackMotilityAnalyzer] #AbstractTrackAnalyzer

    for spot_analyzer, edge_analyzer, track_analyzer in zip(spot_analyzers, edge_analyzers, track_analyzers):
        settings.addSpotAnalyzerFactory(spot_analyzer())
        settings.addEdgeAnalyzer(edge_analyzer())
        settings.addTrackAnalyzer(track_analyzer())

settings.addAllAnalyzers()

print(settings.getSpotAnalyzerFactories())
print(settings.getEdgeAnalyzers())
print(settings.getTrackAnalyzers())
#raise SystemExit("you made it")

#settings.addAllAnalyzers()

# We configure the initial filtering to discard spots 
# with a quality lower than 1.
settings.initialSpotFilterValue = 1.

print(str(settings))

#----------------------
# Instantiate trackmate
#----------------------

trackmate = TrackMate(model, settings)

#------------
# Execute all
#------------


ok = trackmate.checkInput()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

ok = trackmate.process()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

raise SystemExit("you made it")
print("you made it")

#----------------
# Display results
#----------------

model.getLogger().log('Found ' + str(model.getTrackModel().nTracks(True)) + ' tracks.')

# A selection.
sm = SelectionModel( model )

# Read the default display settings.
ds = DisplaySettingsIO.readUserDefault()

# The viewer.
displayer =  HyperStackDisplayer( model, sm, imp, ds ) 
displayer.render()

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
        t=spot.getFeature('FRAME')
        q=spot.getFeature('QUALITY')
        snr=spot.getFeature('SNR_CH1')
        mean=spot.getFeature('MEAN_INTENSITY_CH1')
        model.getLogger().log('\tspot ID = ' + str(sid) + ': x='+str(x)+', y='+str(y)+', t='+str(t)+', q='+str(q) + ', snr='+str(snr) + ', mean = ' + str(mean))
