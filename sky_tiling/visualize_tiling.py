"""Interactive Plotly visualizations for HEALPix and telescope-tile coverage.

The public plotting functions preserve the behavior of the corresponding
functions in the original module:

* ``plot_coverage_with_boundaries_plotly`` plots every HEALPix pixel in a
  two-dimensional Mollweide projection and optionally overlays tile outlines.
* ``plot_coverage_with_boundaries_plotly_globe`` plots every pixel center,
  highlights covered centers, and adds tile outlines on a rotatable globe.
* ``plot_coverage_healpix_with_boundaries_plotly`` additionally draws every
  HEALPix boundary in the Mollweide projection.
* ``plot_coverage_healpix_with_boundaries_plotly_globe`` optionally draws
  every HEALPix boundary on the three-dimensional globe.

Coverage fractions are summed exactly as in the original code, but are used
only to decide whether a pixel is covered (total coverage greater than zero).
They are not mapped to marker color or opacity.
"""

from __future__ import annotations

from collections.abc import Mapping

import healpy as hp
import numpy as np
import plotly.graph_objects as go
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix
from regions import PolygonSkyRegion, Region

from M4OPT_tiling import footprint


CoverageMap = Mapping[object, Mapping[int, float]]


def _build_coverage_map(
    hpx: HEALPix,
    tile_to_healpix_map: CoverageMap,
) -> np.ndarray:
    """Sum all tile-coverage fractions for every HEALPix pixel."""
    coverage = np.zeros(hp.nside2npix(hpx.nside), dtype=float)

    # Keep the original accumulation order and semantics. This is normally
    # inexpensive compared with constructing and rendering Plotly traces.
    for tile_pixels in tile_to_healpix_map.values():
        for pixel, fraction in tile_pixels.items():
            coverage[int(pixel)] += float(fraction)

    return coverage


def _healpix_lonlat(
    hpx: HEALPix,
    pixel_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert HEALPix indices to the longitude convention used by the plots."""
    if pixel_indices.size == 0:
        empty = np.array([], dtype=float)
        return empty, empty.copy()

    ra, dec = hp.pix2ang(
        hpx.nside,
        pixel_indices,
        nest=(hpx.order == "nested"),
        lonlat=True,
    )

    # This deliberately preserves the original map rotation instead of using
    # the conventional wrap-to-[-180, 180] transformation.
    return ra - 180.0, dec


def _healpix_xyz(
    hpx: HEALPix,
    pixel_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert HEALPix pixel centers to Cartesian unit-sphere coordinates."""
    if pixel_indices.size == 0:
        empty = np.array([], dtype=float)
        return empty, empty.copy(), empty.copy()

    theta, phi = hp.pix2ang(
        hpx.nside,
        pixel_indices,
        nest=(hpx.order == "nested"),
    )

    sin_theta = np.sin(theta)
    return (
        sin_theta * np.cos(phi),
        sin_theta * np.sin(phi),
        np.cos(theta),
    )


def _unit_vector(coord: SkyCoord) -> np.ndarray:
    """Return the normalized Cartesian direction of a scalar SkyCoord."""
    vector = np.asarray(coord.cartesian.xyz.value, dtype=float)
    return vector / np.linalg.norm(vector)


def _vectors_to_plotly_lonlat(
    vectors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert Cartesian directions to the original Plotly map convention."""
    x = vectors[:, 0]
    y = vectors[:, 1]
    z = vectors[:, 2]

    ra = np.mod(np.degrees(np.arctan2(y, x)), 360.0)
    dec = np.degrees(np.arctan2(z, np.hypot(x, y)))
    return ra - 180.0, dec


def _interpolate_rodrigues(
    vector1: np.ndarray,
    vector2: np.ndarray,
    *,
    num_points: int,
    angle_tolerance: float,
    axis_tolerance: float,
) -> np.ndarray:
    """Interpolate a great-circle edge using Rodrigues' rotation formula."""
    dot = np.dot(vector1, vector2)
    angle = np.arccos(np.clip(dot, -1.0, 1.0))

    if angle < angle_tolerance:
        return vector1[None, :]

    axis = np.cross(vector1, vector2)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < axis_tolerance:
        return vector1[None, :]

    axis /= axis_norm
    theta = np.linspace(0.0, angle, num_points)
    cos_theta = np.cos(theta)[:, None]
    sin_theta = np.sin(theta)[:, None]

    # Vectorized form of the rotation matrices used by the original code.
    axis_cross_vector1 = np.cross(axis, vector1)
    axis_dot_vector1 = np.dot(axis, vector1)
    vectors = (
        cos_theta * vector1
        + sin_theta * axis_cross_vector1
        + (1.0 - cos_theta) * axis_dot_vector1 * axis
    )

    vectors /= np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors


def _interpolate_slerp(
    vector1: np.ndarray,
    vector2: np.ndarray,
    *,
    num_points: int,
    angle_tolerance: float,
) -> np.ndarray:
    """Interpolate a great-circle edge with spherical linear interpolation."""
    vector1 = vector1 / np.linalg.norm(vector1)
    vector2 = vector2 / np.linalg.norm(vector2)

    dot = np.clip(np.dot(vector1, vector2), -1.0, 1.0)
    angle = np.arccos(dot)
    if angle < angle_tolerance:
        return vector1[None, :]

    t = np.linspace(0.0, 1.0, num_points)[:, None]
    vectors = (
        np.sin((1.0 - t) * angle) / np.sin(angle) * vector1
        + np.sin(t * angle) / np.sin(angle) * vector2
    )
    vectors /= np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors


def _polygon_boundary_lonlat(
    polygon: PolygonSkyRegion,
    *,
    num_points: int,
    angle_tolerance: float,
    axis_tolerance: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate a closed polygon boundary for a two-dimensional plot."""
    vertices = polygon.vertices
    longitude_parts: list[np.ndarray] = []
    latitude_parts: list[np.ndarray] = []

    for index in range(len(vertices)):
        vector1 = _unit_vector(vertices[index])
        vector2 = _unit_vector(vertices[(index + 1) % len(vertices)])
        vectors = _interpolate_rodrigues(
            vector1,
            vector2,
            num_points=num_points,
            angle_tolerance=angle_tolerance,
            axis_tolerance=axis_tolerance,
        )
        longitude, latitude = _vectors_to_plotly_lonlat(vectors)
        longitude_parts.append(longitude)
        latitude_parts.append(latitude)

    return np.concatenate(longitude_parts), np.concatenate(latitude_parts)


def _polygon_boundary_xyz(
    polygon: PolygonSkyRegion,
    *,
    num_points: int,
    angle_tolerance: float,
) -> np.ndarray:
    """Interpolate a closed polygon boundary on the Cartesian unit sphere."""
    vertices = polygon.vertices
    parts: list[np.ndarray] = []

    for index in range(len(vertices)):
        parts.append(
            _interpolate_slerp(
                _unit_vector(vertices[index]),
                _unit_vector(vertices[(index + 1) % len(vertices)]),
                num_points=num_points,
                angle_tolerance=angle_tolerance,
            )
        )

    return np.vstack(parts)


def _healpix_boundary_lonlat(
    hpx: HEALPix,
    pixel: int,
    *,
    boundary_step: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return sampled points on the actual HEALPix pixel boundary.

    ``healpy.boundaries`` returns ``4 * boundary_step`` points on the
    curvilinear HEALPix boundary. Unlike great-circle interpolation between
    the four corners, these points follow the HEALPix boundary definition.
    """
    boundary = hp.boundaries(
        hpx.nside,
        int(pixel),
        step=boundary_step,
        nest=(hpx.order == "nested"),
    ).T
    closed_boundary = np.vstack((boundary, boundary[0]))
    return _vectors_to_plotly_lonlat(closed_boundary)


def _validate_boundary_step(boundary_step: int) -> int:
    """Validate and normalize the HEALPix samples-per-side setting."""
    if not isinstance(boundary_step, (int, np.integer)) or boundary_step < 1:
        raise ValueError("boundary_step must be a positive integer.")
    return int(boundary_step)


def _limit_pixel_indices(
    pixel_indices: np.ndarray,
    max_boundaries: int | None,
) -> np.ndarray:
    """Apply an optional boundary-count cap without changing pixel order."""
    if max_boundaries is None:
        return pixel_indices
    if not isinstance(max_boundaries, (int, np.integer)) or max_boundaries < 0:
        raise ValueError("max_boundaries must be nonnegative or None.")
    return pixel_indices[: int(max_boundaries)]


def _sphere_grid_trace() -> go.Scatter3d:
    """Create the latitude/longitude wireframe used by both globe plots."""
    longitudes = np.linspace(-np.pi, np.pi, 181)
    latitudes = np.linspace(-np.pi / 2.0, np.pi / 2.0, 91)

    x_values: list[float] = []
    y_values: list[float] = []
    z_values: list[float] = []

    def append_line(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> None:
        x_values.extend(x)
        y_values.extend(y)
        z_values.extend(z)
        x_values.append(np.nan)
        y_values.append(np.nan)
        z_values.append(np.nan)

    for latitude in latitudes:
        append_line(
            np.cos(longitudes) * np.cos(latitude),
            np.sin(longitudes) * np.cos(latitude),
            np.full_like(longitudes, np.sin(latitude)),
        )

    # The original grid selects every fifteenth entry from the longitude
    # sample array, which corresponds to a 30-degree interval.
    for longitude in longitudes[::15]:
        append_line(
            np.cos(longitude) * np.cos(latitudes),
            np.sin(longitude) * np.cos(latitudes),
            np.sin(latitudes),
        )

    return go.Scatter3d(
        x=x_values,
        y=y_values,
        z=z_values,
        mode="lines",
        line=dict(color="gray", width=1),
        hoverinfo="skip",
        showlegend=False,
    )


def _append_xyz_polyline(
    x_values: list[float],
    y_values: list[float],
    z_values: list[float],
    points: np.ndarray,
) -> None:
    """Append one 3D polyline followed by a Plotly line-break marker."""
    x_values.extend(points[:, 0])
    y_values.extend(points[:, 1])
    z_values.extend(points[:, 2])
    x_values.append(np.nan)
    y_values.append(np.nan)
    z_values.append(np.nan)


def _write_show(
    figure: go.Figure,
    output_html: str,
    *,
    message_prefix: str,
) -> None:
    """Write the standalone HTML plot and display it, as in the original."""
    figure.write_html(output_html)
    print(f"{message_prefix}{output_html}")
    figure.show()


def plot_tiling_mollweide(
    hpx: HEALPix,
    tile_coords: list[SkyCoord] | None = None,
    base_region: Region | None = None,
    output_html: str = "tilings_mollweide.html",
):
    """Plot all HEALPix centers and optional tile boundaries on a Mollweide map.

    Parameters
    ----------
    hpx
        HEALPix grid whose pixel centers will be displayed.
    tile_coords
        Optional tile-center coordinates. Tile outlines are drawn only when
        both ``tile_coords`` and ``base_region`` are supplied.
    base_region
        Base field-of-view region translated to each tile center by
        ``footprint``.
    output_html
        Destination for the standalone interactive HTML plot.
    """
    all_pixels = np.arange(hp.nside2npix(hpx.nside), dtype=np.int64)
    longitude, latitude = _healpix_lonlat(hpx, all_pixels)

    figure = go.Figure()
    figure.add_trace(
        go.Scattergeo(
            lon=longitude,
            lat=latitude,
            mode="markers",
            marker=dict(size=3, color="green"),
            name=f"All pixels (NSIDE={hpx.nside})",
        )
    )

    if tile_coords is not None and base_region is not None:
        for tile_id, center in enumerate(tile_coords):
            polygon = footprint(base_region, center)
            if not isinstance(polygon, PolygonSkyRegion):
                continue

            boundary_lon, boundary_lat = _polygon_boundary_lonlat(
                polygon,
                num_points=24,
                angle_tolerance=1e-12,
                axis_tolerance=1e-14,
            )
            figure.add_trace(
                go.Scattergeo(
                    lon=boundary_lon,
                    lat=boundary_lat,
                    mode="lines",
                    line=dict(color="crimson", width=1),
                    name=f"Tile {tile_id} boundary",
                    showlegend=False,
                )
            )

    figure.update_layout(
        title=(
            f"All HEALPix pixels (NSIDE={hpx.nside}) "
            "with optional tile boundaries"
        ),
        geo=dict(
            projection_type="mollweide",
            showcoastlines=False,
            showland=False,
            showocean=False,
            lataxis_showgrid=True,
            lonaxis_showgrid=True,
            lataxis_gridcolor="gray",
            lonaxis_gridcolor="gray",
            bgcolor="rgba(0,0,0,0)",
        ),
        height=800,
        width=1500,
    )
    _write_show(
        figure,
        output_html,
        message_prefix="Interactive plot saved as ",
    )


def plot_tiling_globe_interactive(
    hpx: HEALPix,
    tile_to_healpix_map: dict,
    tile_coords: list[SkyCoord],
    base_region: Region,
    output_html: str = "tilings_globe_interactive.html",
    pixel_size: int = 2,
):
    """Plot all HEALPix centers and tile outlines on a rotatable globe.  """
    coverage = _build_coverage_map(hpx, tile_to_healpix_map)
    all_pixels = np.arange(coverage.size, dtype=np.int64)
    covered_pixels = np.flatnonzero(coverage > 0)
    all_x, all_y, all_z = _healpix_xyz(hpx, all_pixels)
    covered_x, covered_y, covered_z = _healpix_xyz(hpx, covered_pixels)

    figure = go.Figure()
    figure.add_trace(_sphere_grid_trace())
    figure.add_trace(
        go.Scatter3d(
            x=all_x,
            y=all_y,
            z=all_z,
            mode="markers",
            marker=dict(
                size=pixel_size,
                color="gray",
                opacity=0.45,
            ),
            name="All HEALPix pixels",
            hoverinfo="skip",
        )
    )
    figure.add_trace(
        go.Scatter3d(
            x=covered_x,
            y=covered_y,
            z=covered_z,
            mode="markers",
            marker=dict(size=pixel_size, color="green", opacity=0.9),
            name="Covered HEALPix pixels",
            hoverinfo="skip",
        )
    )

    tile_x: list[float] = []
    tile_y: list[float] = []
    tile_z: list[float] = []

    for center in tile_coords:
        polygon = footprint(base_region, center)
        if isinstance(polygon, PolygonSkyRegion):
            points = _polygon_boundary_xyz(
                polygon,
                num_points=24,
                angle_tolerance=1e-8,
            )
            _append_xyz_polyline(tile_x, tile_y, tile_z, points)

    if tile_x:
        figure.add_trace(
            go.Scatter3d(
                x=tile_x,
                y=tile_y,
                z=tile_z,
                mode="lines",
                line=dict(color="red", width=4),
                name="Tile boundaries",
                hoverinfo="skip",
                showlegend=False,
            )
        )

    figure.update_layout(
        title="Coverage with tile boundaries (rotatable globe)",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
            bgcolor="white",
        ),
        height=900,
        width=900,
        showlegend=True,
    )
    _write_show(figure, output_html, message_prefix="Saved: ")


def plot_tiling_with_pixel_boundaries_mollweide(
    hpx: HEALPix,
    tile_to_healpix_map: dict,
    tile_coords: list[SkyCoord],
    base_region: Region,
    output_html: str = "tiling_with_pixel_boundaries.html",
    boundary_step: int = 4,
    max_boundaries: int | None = None,
):
    """Plot all HEALPix cells and tile boundaries on a Mollweide map. """
    coverage = _build_coverage_map(hpx, tile_to_healpix_map)
    all_pixels = np.arange(coverage.size, dtype=np.int64)
    covered_pixels = np.flatnonzero(coverage > 0)
    all_longitude, all_latitude = _healpix_lonlat(hpx, all_pixels)
    covered_longitude, covered_latitude = _healpix_lonlat(
        hpx,
        covered_pixels,
    )
    boundary_step = _validate_boundary_step(boundary_step)
    boundary_pixels = _limit_pixel_indices(all_pixels, max_boundaries)

    figure = go.Figure()
    figure.add_trace(
        go.Scattergeo(
            lon=all_longitude,
            lat=all_latitude,
            mode="markers",
            marker=dict(size=3, color="gray", showscale=False),
            name="All HEALPix pixels",
        )
    )
    figure.add_trace(
        go.Scattergeo(
            lon=covered_longitude,
            lat=covered_latitude,
            mode="markers",
            marker=dict(size=3, color="green", showscale=False),
            name="Covered HEALPix pixels",
        )
    )

    # Combine all pixel outlines into one trace. Text labels retain the pixel
    # identity that previously came from one separately named trace per pixel.
    pixel_boundary_lon: list[float | None] = []
    pixel_boundary_lat: list[float | None] = []
    pixel_boundary_text: list[str | None] = []

    for pixel in boundary_pixels:
        boundary_lon, boundary_lat = _healpix_boundary_lonlat(
            hpx,
            int(pixel),
            boundary_step=boundary_step,
        )
        label = f"Pixel {int(pixel)} Boundary"
        pixel_boundary_lon.extend(boundary_lon)
        pixel_boundary_lat.extend(boundary_lat)
        pixel_boundary_text.extend([label] * boundary_lon.size)
        pixel_boundary_lon.append(None)
        pixel_boundary_lat.append(None)
        pixel_boundary_text.append(None)

    if pixel_boundary_lon:
        figure.add_trace(
            go.Scattergeo(
                lon=pixel_boundary_lon,
                lat=pixel_boundary_lat,
                mode="lines",
                line=dict(color="black", width=0.5),
                name="HEALPix pixel boundaries",
                text=pixel_boundary_text,
                hovertemplate=(
                    "%{text}<br>lon=%{lon:.6f}<br>lat=%{lat:.6f}"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        )

    # Keep one trace per tile because the original function gives every tile
    # its own legend entry and interactive visibility toggle.
    for tile_id, center in enumerate(tile_coords):
        polygon = footprint(base_region, center)
        if not isinstance(polygon, PolygonSkyRegion):
            continue

        boundary_lon, boundary_lat = _polygon_boundary_lonlat(
            polygon,
            num_points=20,
            angle_tolerance=1e-6,
            axis_tolerance=1e-10,
        )
        figure.add_trace(
            go.Scattergeo(
                lon=boundary_lon,
                lat=boundary_lat,
                mode="lines",
                line=dict(color="red", width=1),
                name=f"Tile {tile_id} Boundary",
            )
        )

    figure.update_layout(
        geo=dict(
            projection_type="mollweide",
            showland=False,
            showcountries=False,
            showocean=False,
            bgcolor="rgba(0,0,0,0)",
            lataxis_showgrid=True,
            lonaxis_showgrid=True,
            lataxis_gridcolor="gray",
            lonaxis_gridcolor="gray",
        ),
        title="Tile Coverage Map with Boundaries",
        height=600,
        width=1200,
    )
    _write_show(
        figure,
        output_html,
        message_prefix="Interactive plot saved as ",
    )


def plot_tiling_with_pixel_boundaries_globe_interactive(
    hpx: HEALPix,
    tile_to_healpix_map: (
        dict[int, dict[int, float]] | dict[str, dict[int, float]]
    ),
    tile_coords: list[SkyCoord],
    base_region,
    output_html: str = "coverage_globe.html",
    pixel_point_size: int = 2,
    draw_pixel_boundaries: bool = True,
    max_boundaries: int | None = None,
    boundary_step: int = 4,
):
    """Plot all HEALPix cells and tile outlines on a 3D globe."""
    coverage = _build_coverage_map(hpx, tile_to_healpix_map)
    all_pixels = np.arange(coverage.size, dtype=np.int64)
    covered_pixels = np.flatnonzero(coverage > 0)
    all_x, all_y, all_z = _healpix_xyz(hpx, all_pixels)
    covered_x, covered_y, covered_z = _healpix_xyz(hpx, covered_pixels)

    figure = go.Figure()
    figure.add_trace(_sphere_grid_trace())
    figure.add_trace(
        go.Scatter3d(
            x=all_x,
            y=all_y,
            z=all_z,
            mode="markers",
            marker=dict(
                size=pixel_point_size,
                color="gray",
                opacity=0.45,
            ),
            name="All HEALPix pixels",
            hoverinfo="skip",
        )
    )
    figure.add_trace(
        go.Scatter3d(
            x=covered_x,
            y=covered_y,
            z=covered_z,
            mode="markers",
            marker=dict(size=pixel_point_size, color="green", opacity=0.9),
            name="Covered pixels",
            hoverinfo="skip",
        )
    )

    if draw_pixel_boundaries:
        boundary_step = _validate_boundary_step(boundary_step)
        boundary_pixels = _limit_pixel_indices(all_pixels, max_boundaries)

        boundary_x: list[float] = []
        boundary_y: list[float] = []
        boundary_z: list[float] = []

        for pixel in boundary_pixels:
            boundary = hp.boundaries(
                hpx.nside,
                int(pixel),
                step=boundary_step,
                nest=(hpx.order == "nested"),
            ).T
            closed_boundary = np.vstack((boundary, boundary[0]))
            _append_xyz_polyline(
                boundary_x,
                boundary_y,
                boundary_z,
                closed_boundary,
            )

        if boundary_x:
            figure.add_trace(
                go.Scatter3d(
                    x=boundary_x,
                    y=boundary_y,
                    z=boundary_z,
                    mode="lines",
                    line=dict(color="gray", width=2),
                    name="Pixel boundaries",
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

    tile_x: list[float] = []
    tile_y: list[float] = []
    tile_z: list[float] = []

    for center in tile_coords:
        polygon = footprint(base_region, center)
        if isinstance(polygon, PolygonSkyRegion):
            points = _polygon_boundary_xyz(
                polygon,
                num_points=32,
                angle_tolerance=1e-12,
            )
            _append_xyz_polyline(tile_x, tile_y, tile_z, points)

    if tile_x:
        figure.add_trace(
            go.Scatter3d(
                x=tile_x,
                y=tile_y,
                z=tile_z,
                mode="lines",
                line=dict(color="red", width=6),
                name="Tile boundaries",
                hoverinfo="skip",
                showlegend=False,
            )
        )

    figure.update_layout(
        title="Coverage with tile boundaries (3D globe)",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
            bgcolor="white",
        ),
        height=900,
        width=900,
        showlegend=True,
    )
    _write_show(figure, output_html, message_prefix="Saved: ")