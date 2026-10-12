"""Rasterize registered QuPath cell polygons and retain the Xenium/IMC ID map.
Adapted from the project's XeniumSegTransfer/JsonToMask.ipynb (Yi-Nan Li).
Uses the original uint16 labels, unit pixel grid, and vertical-flip convention.
"""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polygons-dir", required=True, type=Path)
    parser.add_argument("--images-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mapping-csv", required=True, type=Path)
    args = parser.parse_args()

    import numpy as np
    from affine import Affine
    from rasterio.features import rasterize
    import tifffile

    with args.images_csv.open(newline="", encoding="utf-8-sig") as handle:
        images = list(csv.DictReader(handle))
    if not images or not {"image", "width_px", "height_px"}.issubset(images[0]):
        raise ValueError("Image CSV requires image, width_px, height_px.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    mappings = []
    seen_rois = set()
    for image in images:
        sample_id = Path(image["image"]).stem
        if sample_id in seen_rois:
            raise ValueError(f"Duplicate ROI name: {sample_id}")
        seen_rois.add(sample_id)
        width, height = int(image["width_px"]), int(image["height_px"])
        if width <= 0 or height <= 0:
            raise ValueError("Image dimensions must be positive integers.")
        polygon_path = args.polygons_dir / f"{sample_id}.geojson"
        with polygon_path.open(encoding="utf-8") as handle:
            cells = json.load(handle)["features"]
        if not cells:
            raise ValueError(f"No cell polygons in {polygon_path}")
        geometries, ids, names = [], [], []
        for index, cell in enumerate(cells, start=1):
            props = cell.get("properties", {})
            name = props.get("name")
            cell_id = props.get("cell_id", index)
            if name is None or str(name) == "":
                raise ValueError("Every exported cell requires its original Xenium name.")
            if isinstance(cell_id, bool) or int(cell_id) != float(cell_id):
                raise ValueError("Cell labels must be positive integers.")
            cell_id = int(cell_id)
            if not 0 < cell_id <= np.iinfo(np.uint16).max:
                raise ValueError("Cell label exceeds the original uint16 mask range.")
            if cell["geometry"]["type"] not in {"Polygon", "MultiPolygon"}:
                raise ValueError("Masks require polygon cell geometries.")
            geometries.append(cell["geometry"])
            ids.append(cell_id)
            names.append(str(name))
        if len(set(ids)) != len(ids) or len(set(names)) != len(names):
            raise ValueError(f"Cell IDs/names are not unique within {sample_id}.")
        mask = rasterize(zip(geometries, ids), out_shape=(height, width),
                         transform=Affine(1, 0, 0, 0, -1, height), fill=0, dtype="uint16")
        tifffile.imwrite(args.output_dir / f"{sample_id}.tiff", np.flipud(mask))
        mappings.extend({"sample_id": sample_id, "cell_name": name, "cell_id": label}
                        for name, label in zip(names, ids))
    args.mapping_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.mapping_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "cell_name", "cell_id"])
        writer.writeheader()
        writer.writerows(mappings)
    print(f"Saved {len(images)} masks and {len(mappings)} cell mappings.")


if __name__ == "__main__":
    main()
