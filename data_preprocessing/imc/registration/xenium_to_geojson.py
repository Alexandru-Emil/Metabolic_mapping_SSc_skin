"""Export Xenium boundaries in QuPath pixel coordinates.
Adapted from the project's XeniumSegTransfer/XeniumSegJson.ipynb
(Daria Buszta, Senbai Kang and Yi-Nan Li). No registration is fitted here.
"""
import argparse
import json
from pathlib import Path
from uuid import uuid4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xenium-output", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pixel-size-um", type=float, default=0.2125)
    args = parser.parse_args()
    if not 0 < args.pixel_size_um < float("inf"):
        parser.error("--pixel-size-um must be finite and positive.")

    import spatialdata_io as sdio
    from shapely.affinity import scale
    from shapely.geometry import mapping

    data = sdio.xenium(args.xenium_output, nucleus_boundaries=False,
                       cells_as_circles=False, cells_labels=False)
    polygons = data.shapes["cell_boundaries"]
    if not polygons.index.is_unique:
        raise ValueError("Xenium cell names must be unique.")
    features = []
    for cell_name, row in polygons.iterrows():
        geometry = row.geometry
        if geometry.is_empty:
            continue
        if geometry.geom_type not in {"Polygon", "MultiPolygon"}:
            raise ValueError(f"Unsupported cell geometry: {geometry.geom_type}")
        geometry = scale(geometry, xfact=1 / args.pixel_size_um,
                         yfact=1 / args.pixel_size_um, origin=(0, 0))
        features.append({"type": "Feature", "id": str(uuid4()),
                         "geometry": mapping(geometry),
                         "properties": {"name": str(cell_name), "objectType": "detection"}})
    if not features:
        raise ValueError("No nonempty Xenium cell boundaries were found.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump({"type": "FeatureCollection", "features": features}, handle)
    print(f"Exported {len(features)} cells to {args.output}")


if __name__ == "__main__":
    main()
