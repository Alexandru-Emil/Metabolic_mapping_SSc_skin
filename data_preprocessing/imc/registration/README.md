# Xenium–IMC segmentation transfer

1. Run `xenium_to_geojson.py --xenium-output RAW_XENIUM_DIR --output cells.geojson`. The original morphology pixel size is 0.2125 µm; set `--pixel-size-um` if the source image differs.
2. Import the named detections into the Xenium image in QuPath and save the project. Register Xenium and IMC with Warpy, or use the saved project transforms.
3. Open the target IMC image and run `transfer_warpy.groovy`, setting its source ID and output path. Name the output `ROI_NAME.geojson` to match the image CSV. This uses native Warpy affine/TPS transforms, not the experimental Python RBF approximation.
4. Run `polygons_to_masks.py --polygons-dir imc_poly --images-csv images.csv --output-dir imc_mask --mapping-csv cellid_mapping.csv`. The CSV requires `image,width_px,height_px`.

Use these masks for Steinbock's post-Xenium IMC intensity extraction. Dependencies: spatialdata-io, shapely, numpy, affine, rasterio, tifffile; QuPath with Warpy. Registration projects, transforms, raw images and manual gating exports are external inputs.

Adapted from the project notebooks `XeniumSegJson.ipynb` and `JsonToMask.ipynb` (Daria Buszta, Senbai Kang and Yi-Nan Li). The Warpy adapter derives from BIOP's Apache-2.0 example; its license is retained under `third_party_licenses`.
