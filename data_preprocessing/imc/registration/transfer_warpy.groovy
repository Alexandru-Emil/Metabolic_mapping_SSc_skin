/**
 * Native Warpy transfer of Xenium detections to the currently open IMC image.
 * Adapted from BIOP's Transfer_annotations_and_detections_to_current_entry.groovy
 * by Nicolas Chiaruttini and Olivier Burri (Apache-2.0).
 * Source: https://github.com/BIOP/qupath-extension-warpy
 * License: ../../third_party_licenses/Warpy_Apache_2_LICENSE.txt
 */
import qupath.ext.warpy.Warpy
import qupath.lib.io.PathIO

// Set these in the editor, or supply the corresponding environment variables.
def sourceImageId = System.getenv('METABOLIC_REGISTRATION_SOURCE_ID') ?: ''
def outputGeojson = System.getenv('METABOLIC_REGISTRATION_OUTPUT') ?: ''
if (!outputGeojson) throw new IllegalArgumentException('Set outputGeojson to the ROI .geojson output path.')
def targetEntry = getProjectEntry()
if (targetEntry == null) throw new IllegalStateException('Open the target IMC image in its registered QuPath project.')
def sources = Warpy.getCandidateSourceEntries(targetEntry).toList()
def sourceEntry
if (sourceImageId) sourceEntry = sources.find { it.getID() == sourceImageId }
else if (sources.size() == 1) sourceEntry = sources[0]
else throw new IllegalArgumentException('Set sourceImageId explicitly when there are multiple registered source images.')
if (sourceEntry == null) throw new IllegalStateException('No saved Warpy transform matches the requested source and target.')
def transform = Warpy.getRealTransform(sourceEntry, targetEntry)
if (transform == null) throw new IllegalStateException('The source-to-target transform is unavailable.')
def objects = Warpy.getPathObjectsFromEntry(sourceEntry).findAll {
    it.isDetection() && it.getROI() != null && it.getName() != null && !it.getName().isEmpty()
}
if (objects.isEmpty()) throw new IllegalStateException('Import the named Xenium detections into the source image and save it first.')
def transferred = Warpy.transformPathObjects(objects, transform)
if (transferred.size() != objects.size()) throw new IllegalStateException('Warpy did not transfer every source detection.')
def destination = new File(outputGeojson)
destination.getAbsoluteFile().getParentFile().mkdirs()
PathIO.exportObjectsAsGeoJSON(destination, transferred, PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
println 'Exported ' + transferred.size() + ' registered cells to ' + destination
