"""Create schematic maps of the energy-balance model (EBM)
applied to the Arctic (manuscript Fig. 2) or reduced-domain
version (with two reference latitudes) for the Southern Ocean
(Supplementary Figure S4).

This uses the cartopy module (version 0.22.0 at time of creating
this script) to generate 'NearsidePerspective' projection view
of the polar region, on which sea ice extent climatology data is
plotted, and over which various patches are drawn using
matplotlib (version 3.8.1) to represent the EBM schematic.

For the Southern Ocean, the figure is created upside down
(specifically, rotated through 180 degrees) because it does not
seem possible to rotate the 'satellite' view in cartopy's
NearsidePerspective projection, only displace the position of
said view above the Earth's surface and adjust the height. It
may in fact be possible, but I could not figure it out and it is
easy enough to just rotate the svg output manually (which I did
using inkscape) before exporting to the final format.

-Jake R. Aylmer

"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import plot_style, plot_tools, script_tools
from utils.observations import \
    get_passive_microwave_siconc_climatology as get_siconc


# Z-order of plot elements:
_zorder = {
    "oht_arrows"     : -1,
    "ax_spines"      : -2,
    "ref_lats_labels": -3,
    "ref_lats_lines" : -4,
    "gridlines"      : -5,
    "land"           : -6,
    "rasterization"  : -7,
    "sea_ice"        : -8
}

def get_zorder(plot_element):
    """Get the z-order of a plot element. Uses the function
    defined in plot_style module, replacing the zorder dict with
    the ad-hoc one above.
    """
    return plot_style.get_zorder(plot_element, _zorder)



def setup_globe(fig, hemi="n", central_longitude=330.0,
                central_latitude=65.0,
                grid_latitudes=np.arange(55.0, 90.0, 5.0),
                grid_lon_interval=30.0,
                grid_line_kw={"color": "lightgrey",
                              "linewidth": 0.25,
                              "linestyle": "-",
                              "zorder": get_zorder("gridlines")},
                ocean_color="#009ABB", land_color="lightgrey",
                globe_edge_color="lightgrey",
                land_kw = {"edgecolor": "none",
                           "zorder": get_zorder("land")}):
    """Create the perspective view of the north or south pole
    as a cartopy NearsidePerspective geoaxes on which data can
    be plotted.
    """
    
    # This the default "satellite_height" in meters hard-coded
    # in the cartopy NearsidePerspective projection,
    # corresponding to true geostationary orbit:
    satellite_height = 35785831.0
    
    # It needs adjusting for the schematic, by a different
    # amount for each hemisphere:
    satellite_height /= (4.0 if hemi == "n" else 3.0)
    
    # Create the NearsidePerspective axes on the figure at an
    # arbitrary position (to be adjusted afterwards):
    ax = fig.add_axes([0,0,0.5,0.5],
        projection=ccrs.NearsidePerspective(
            central_longitude=central_longitude,
            central_latitude=central_latitude,
            satellite_height=satellite_height),
        facecolor=ocean_color)
    
    t = 0.2 if hemi == "n" else -0.6
    b = -0.6 if hemi == "n" else 0.2
    
    plot_tools.distribute_subplots(np.array([ax]), l=-0.5,
                                   r=-0.5, t=t, b=b)
    
    ax.set_global()
    
    for spine in ax.spines:
        ax.spines[spine].set(color=globe_edge_color,
            linewidth=mpl.rcParams["lines.linewidth"],
            zorder=get_zorder("ax_spines"))
    
    land_kw["facecolor"] = land_color
    
    ax.add_feature(cfeature.NaturalEarthFeature(
        "physical", "land", "110m", **land_kw))
    
    # Add gridlines (specify appearance):
    gl = ax.gridlines(**grid_line_kw)
    
    # Set longitudes of gridlines at specified interval:
    gl.xlocator = mpl.ticker.FixedLocator(
        np.arange(-180.0, 180.0, grid_lon_interval))
    
    # Set latitudes of gridlines at specified inputs:
    gl.ylocator = mpl.ticker.FixedLocator(
        (-1 if hemi=="s" else 1)*abs(grid_latitudes))
    
    # Everything below this zorder is rasterized in vector
    # output:
    ax.set_rasterization_zorder(get_zorder("rasterization"))
    
    return ax



def setup_guide_axes(fig):
    """Add invisible axes for annotations, aligned to the figure
    canvas. Coordinates of the bottom-left corner are (0,0), and
    of the top-right corner are (fig_ar, 1) where fig_ar is the
    aspect ratio of the figure (this is so that coordinates on
    this set of axes are equal aspect: a rectangle of width
    dx = 1 and height dy = 1 is thus a square in display
    coordinates.
    """
    
    # Add the axes spanning the whole figure canvas:
    ax = fig.add_axes([0,0,1,1])
    
    # Make everything invisible: no background colour, grid-
    # lines, axes ticks/labels, or spines:
    ax.set_facecolor("none")
    
    ax.grid(False)
    
    ax.tick_params(which="both", axis="both", left=False,
        right=False, top=False, bottom=False, labelleft=False,
        labeltop=False, labelbottom=False)
    
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    
    # Set the limits so that dx = dy is equal aspect on display;
    # determine the figure aspect ratio (width/height):
    figsize = fig.get_size_inches()
    fig_ar = figsize[0]/figsize[1]
    
    ax.set_xlim((0.0, fig_ar))
    ax.set_ylim((0.0, 1.0))
    
    return ax



def add_sea_ice_climatology(ax, hemi="n",
        time_period_1=(1980, 2000), time_period_2=(2001, 2021)):
    """Plot sea ice extent climatology for two time periods.
    Returns the colour that sea ice is plotted in for the later
    time period, needed by one of the annotation labels.
    """
    
    def _fill_artefact(siconc, l1=0, l2=360, p1=-90, p2=90):
        """Set lon-lat box siconc values to 1."""
        return np.where((lon > l1) & (lon < l2) & (lat > p1)
                        & (lat < p2), 1, siconc)
    
    for time_period, cmap, set_val in zip(
            [time_period_1, time_period_2],
            ["Blues_r", "Blues"],
            [0.65, 0.0]):
        
        # Load climatology data:
        lon, lat, siconc = get_siconc(hemi=hemi,
                                      time_period=time_period)
        
        # This hard-coded part fills in specific grid cells with
        # sea ice concentration = 1. This is purely to make the
        # figure, which just needs to illustrate the typical sea
        # ice edge, look better and mainly fills in spurious
        # gaps around the Arctic/Antarctic coastlines. It does
        # not alter data along the sea ice edge in the open
        # ocean areas.
        lon = lon % 360  # lon in range 0-360
        
        # Arctic regions:
        siconc = _fill_artefact(siconc, 0, 360, 88, 90)
        siconc = _fill_artefact(siconc, 240, 290, 60, 90)
        siconc = _fill_artefact(siconc, 230, 335, 67, 90)
        siconc = _fill_artefact(siconc, 330, 340, 69, 90)
        siconc = _fill_artefact(siconc, 280, 299, 55, 90)
        siconc = _fill_artefact(siconc, 275, 285, 50, 56)
        siconc = _fill_artefact(siconc, 60, 210, 62, 82)
        siconc = _fill_artefact(siconc, 150, 167, 60, 65)
        siconc = _fill_artefact(siconc, 134, 144, 50, 65)
        siconc = _fill_artefact(siconc, 144, 156, 58.75, 60)
        siconc = _fill_artefact(siconc, 161, 163.7, 57.3, 60)
        siconc = _fill_artefact(siconc, 164.7, 166.7, 59.6, 60)
        siconc = _fill_artefact(siconc, 168, 170.3, 60.1, 60.6)
        siconc = _fill_artefact(siconc, 168, 169.4, 60.3, 60.6)
        siconc = _fill_artefact(siconc, 171, 172, 60.4, 62)
        siconc = _fill_artefact(siconc, 171.8, 173.2, 60.6, 62)
        siconc = _fill_artefact(siconc, 173.3, 174.3, 61.5, 62)
        siconc = _fill_artefact(siconc, 193, 203, 58.5, 62.1)
        siconc = _fill_artefact(siconc, 192.6, 193.3, 59.6, 60.1)
        siconc = _fill_artefact(siconc, 330, 350, 80, 83)
        siconc = _fill_artefact(siconc, 325, 345, 74, 78)
        siconc = _fill_artefact(siconc, 44, 62, 79, 82)
        siconc = _fill_artefact(siconc, 16.3, 33.5, 77.4, 80.5)
        siconc = _fill_artefact(siconc, 10.9, 20, 79.4, 80.4)
        siconc = _fill_artefact(siconc, 15.8, 17, 76.7, 87.6)
        siconc = _fill_artefact(siconc, 298.6, 300.7, 54.95, 56.3)
        siconc = _fill_artefact(siconc, 317, 322.8, 65.2, 66.1)
        siconc = _fill_artefact(siconc, 322, 322.9, 65.7, 66.1)
        siconc = _fill_artefact(siconc, 315.9, 317.8, 61.5, 65)
        siconc = _fill_artefact(siconc, 317.6, 319.2, 62.8, 65)
        siconc = _fill_artefact(siconc, 316.3, 317, 60.6, 61)
        siconc = _fill_artefact(siconc, 318, 319.3, 63.7, 65)
        siconc = _fill_artefact(siconc, 56.2, 61.4, 75, 76.1)
        siconc = _fill_artefact(siconc, 52.8, 61, 69.5, 76.1)
        siconc = _fill_artefact(siconc, 47.9, 60.9, 67.4, 69.5)
        siconc = _fill_artefact(siconc, 32, 48.4, 63.7, 67.86)
        siconc = _fill_artefact(siconc, 43.9, 46.9, 67.9, 68.7)
        siconc = _fill_artefact(siconc, 19.4, 26, 63.6, 66.3)
        siconc = _fill_artefact(siconc, 19.3, 23, 63.4, 65.3)
        siconc = _fill_artefact(siconc, 21.05, 23, 62.78, 63.7)
        siconc = _fill_artefact(siconc, 21.17, 22, 62.7, 62.9)
        
        # Southern Ocean regions:
        siconc = _fill_artefact(siconc, 95, 100, -65, -60)
        siconc = _fill_artefact(siconc, 0, 360, -90, -68)
        siconc = _fill_artefact(siconc, 75, 150, -90, -65)
        # ---------------------------------------------------- #
        
        # Just plot area of sea ice extent; create mask where
        # sea ice concentration exceeds 15% threshold:
        mask = np.where((siconc > 0.15), 0, 1)
        siconc = np.ma.masked_array(
            set_val*np.ones(np.shape(siconc)), mask=mask)
        
        cmap_use = mpl.colormaps.get_cmap(cmap)
        norm = mpl.colors.BoundaryNorm(np.arange(0., 1.01, .05),
                                       ncolors=cmap_use.N,
                                       clip=True)
        
        ax.pcolormesh(lon, lat, siconc[1:,1:], cmap=cmap_use,
                      norm=norm, zorder=get_zorder("sea_ice"),
                      transform=ccrs.PlateCarree())
    
    # Return colour of sea ice plotted (second time period):
    return cmap_use(1)



def translate(xy, dr=(0.0, 0.0)):
    """Cartesian translation transformation of points
    xy = [(x1, y1), (x2, y2), ...] each by the vector
    dr = (dx, dy).
    """
    return xy + np.reshape(dr, (1,2))



def rotate(xy, theta=180.0):
    """Rotation transformation of points xy = [(x1, y1), ...]
    about the center of mass of all points in xy by angle theta
    anticlockwise.
    """
    theta_r = np.pi*theta/180.0
    
    rot_matrix = np.array([[np.cos(theta_r), np.sin(theta_r)],
                           [-np.sin(theta_r), np.cos(theta_r)]])
    
    xy_cm = np.array([[np.mean(xy[:,0]), np.mean(xy[:,1])]])
    
    return np.dot(rot_matrix, translate(xy, -xy_cm).T).T + xy_cm



def vertical_arrow(xy, size=1.0, scale=0.02, height=1.0,
                   rotation=None, slant_factor=0.5,
                   color="tab:blue"):
    """Draw a filled arrow centered at location xy = (x, y).
    
    
    Optional parameters
    -------------------
    size : float, default = 1.0
        Overall size of the arrow relative to default (set by
        scale).
    
    scale : float, default = 0.02
        Overall size of the default arrow (depends on the
        coordinate system).
    
    height : float, default = 1.0
        Height of the arrow tail relative to the default.
    
    rotation : float or None (default)
        Angle in degrees by which to rotate the arrow. Default
        is None, equivalently zero, and the arrow points
        vertically upwards.
    
    slant_factor : float, default = 0.5
        Degree to which arrow vertices are distorted to give it
        a 'slanted'/perspective/3D/angled appearance. With
        slant_factor = 0, arrows are symmetrical/2D/flat in
        appearance.
    
    color : matplotlib colour identifier (string, float, etc.)
        Colour of the arrow; default = "tab:blue".
    
    
    Returns
    -------
    matplotlib.patches.Polygon instance
    
    """
    
    # Define some metrics which define the shape of the arrow
    # (units are arbitrary here):
    tailW = 1.5         # arrow tail width
    tailH = 7.5*height  # arrow tail height
    headH = 2.5         # arrow head height
    headW = 3.5         # arrow head width
    
    totalH = tailH + headH
    
    # Coordinates of the arrow shape nodes (unclosed path), 
    # pointing upward by default:
    verts = np.array([
        [-tailW/2, -totalH/2 + slant_factor],
        [tailW/2, -totalH/2 - 0.75*slant_factor],
        [tailW/2, totalH/2 - headH - slant_factor],
        [headW/2, totalH/2 - headH - 2*slant_factor],
        [-0.1*slant_factor, totalH/2],
        [-headW/2, totalH/2 - headH + slant_factor],
        [-tailW/2, totalH/2 - headH]
    ])
    
    # Center the arrow vertices (before scaling):
    verts = translate(verts, dr=-np.mean(verts, axis=0))
    
    # Apply stretch transformation:
    verts *= scale*size
    
    # Move the arrow to the specified location:
    verts = translate(verts, dr=np.array(xy))
    
    # Rotate as required:
    if rotation is not None:
        verts = rotate(verts, theta=rotation)
    
    # Return the Polygon patch:
    return mpl.patches.Polygon(verts, closed=True,
        edgecolor="k", facecolor=color,
        linewidth=1.5*mpl.rcParams["lines.linewidth"])



def add_arctic_annotations(ax, fig_ar=1.2,
                           sea_ice_color="white"):
    """Add Arctic-specific annotations (arrows and text labels)
    to the figure. ax should be the invisible/cartesian axes
    spanning the figure canvas with size aspect ratio fig_ar.
    Also requires sea_ice_color (matplotlib color) because one
    of the text labels needs a background patch and is
    positioned over sea ice.
    """
    
    ax.add_patch(vertical_arrow((0.56*fig_ar, 0.81),
                 color="tab:orange"))
    ax.annotate("Outgoing longwave\nradiation ("
                + r"$F_\mathrm{OLR}$)", (0.51*fig_ar, 0.95),
                ha="right", va="top")
    
    ax.add_patch(vertical_arrow((0.67*fig_ar, 0.5),
                 color="tab:olive", rotation=180.0, height=3.0))
    ax.annotate("Net shortwave\nradiation ("
                + r"$F_\mathrm{sw}$)", (0.7*fig_ar, 0.88),
                ha="left", va="top")
    
    ax.add_patch(vertical_arrow((0.4*fig_ar, 0.62),
                 color="tab:green"))
    ax.add_patch(vertical_arrow((0.34*fig_ar, 0.65),
                 color="tab:green", rotation=180.0, height=0.8))
    ax.annotate("Surface radiation/\nturbulent heat\n"
                + "fluxes (" + r"$F_\mathrm{up}$, "
                + r"$F_\mathrm{down}$)", (0.43*fig_ar, 0.57),
                ha="left", va="top",
                fontsize=plot_style.fs_0-4.5,
                bbox={"facecolor": sea_ice_color,
                      "edgecolor": "none", "pad": 2})
    
    ax.annotate("Atmospheric\nheat transport\n(AHT)",
        (0.22*fig_ar, 0.2), ha="right", va="top",
        fontsize=plot_style.fs_0-3,
        bbox={"facecolor": "white", "alpha": 0.75,
              "edgecolor": "none", "pad": 2})
    
    ax.annotate("Ocean heat\ntransport (OHT)",
        (0.71*fig_ar, 0.14), ha="left", va="top",
        fontsize=plot_style.fs_0-3,
        bbox={"facecolor": "white", "alpha": 0.75,
              "edgecolor": "none", "pad": 2})



def add_southern_ocean_annotations(ax, fig_ar=1.2,
        central_longitude=150.0, ref_lat_label_lon=135):
    """As in add_arctic_annotations but for the southern
    ocean figure. Here, need the central_longitude of the
    projection and longitude of reference latitude label to add
    one additional text label ("zonal-mean Antarctic coast")
    underneath at the correct rotation.
    """
    
    ax.add_patch(vertical_arrow((0.36*fig_ar, 0.175),
                 color="tab:orange", rotation=180.0))
    ax.annotate("Outgoing longwave\nradiation ("
                + r"$F_\mathrm{OLR}$)", (0.78*fig_ar, 0.14),
                ha="right", va="top", rotation=180.0)
    
    ax.add_patch(vertical_arrow((0.2*fig_ar, 0.5),
                 color="tab:olive", height=3.0))
    ax.annotate("Net shortwave\nradiation ("
                + r"$F_\mathrm{sw}$)", (0.16*fig_ar, 0.1),
                ha="center", va="center", rotation=180.0)
    
    ax.add_patch(vertical_arrow((0.73*fig_ar, 0.25),
                 color="tab:green", rotation=180.0))
    ax.add_patch(vertical_arrow((0.79*fig_ar, 0.35),
                 color="tab:green", height=0.8))
    ax.annotate("Surface radiation/\nturbulent heat\n"
                + "fluxes (" + r"$F_\mathrm{up}$, "
                + r"$F_\mathrm{down}$)", (0.78*fig_ar, 0.18),
                ha="left", va="center", rotation=180.0,
                fontsize=plot_style.fs_0-4.5)
    
    ax.annotate("OHT and AHT away from\nSouthern Ocean and\n"
                + "towards the south pole", (0.5*fig_ar, 0.46),
                rotation=180.0, ha="center", va="center",
                fontsize=plot_style.fs_0-4.5,
                bbox={"facecolor": "white", "edgecolor": "none",
                      "alpha": 0.75, "pad": 2})
    
    ax.annotate("(zonal-median\nAntarctic coast)",
                (0.41*fig_ar, 0.82),
                rotation=central_longitude - ref_lat_label_lon
                         + 180.0, ha="center", va="center",
                fontsize=plot_style.fs_0-4.5)



def add_heat_transports(ax, ref_lats=[65.0],
                        aht_lons=[np.arange(40, 360, 60)],
                        oht_lons=[np.array([191, 302, 330, 358])],
                        offset_lat=0.75, aht_color="tab:blue",
                        oht_color="tab:red"):
    """Add the heat transport arrows at each references latitude
    ref_lats (list of float) at specified longitudes.
    
    
    Optional parameters
    -------------------
    aht_lons : list of arrays (or lists) of float
        Longitudes at which to draw AHT arrows for each
        reference latitude.
    
    oht_lons : list of arrays (or lists) of float
        As in aht_lons but for OHT arrows.
    
    offset_lat : float
        Draw arrows slightly off-center in the meridional
        direction by offset_lat (degrees_north).
    
    aht_color : str or matplotlib color
        Face color of AHT arrows.
    
    oht_color : str or matplotlib color
        Face color of OHT arrows.
    """
    
    shemi = ref_lats[0] < 0
    
    if shemi:
        offset_lat *= -1.0
    
    for j in range(len(ref_lats)):
        for lons, col in zip([aht_lons[j], oht_lons[j]],
                             [aht_color, oht_color]):
            ax.quiver(np.array(lons),
                (offset_lat+ref_lats[j])*np.ones(len(lons)),
                np.zeros(len(lons)),
                (-1 if shemi else 1)*np.ones(len(lons)),
                color=col, transform=ccrs.PlateCarree(),
                pivot="mid", edgecolor="k", units="inches",
                scale_units="inches", scale=4.0,
                width=0.0325, headlength=3.5,
                headaxislength=3.5, headwidth=3,
                zorder=get_zorder("oht_arrows"),
                linewidth=1.5*mpl.rcParams["lines.linewidth"],
                joinstyle="miter", capstyle="projecting")



def add_reference_latitudes(ax, ref_lats=[65],
        ref_lats_labels_lons=[315], central_longitude=330,
        hemi="n", lon_labels_facecolors=["lightgrey"]):
    """Add reference latitude dashed lines and labels at
    specified longitudes (ref_lats_labels_lons).
    
    
    Optional parameters
    -------------------
    central_longitude : int or float
        Central longitude of plot (that which is also passed to
        cartopy.ccrs.NearsidePerspective).
    
    hemi : str, "n" or "s"
        Which hemisphere is being drawn (needed to determine
        rotation of labels).
    
    lon_labels_facecolors : list of matplotlib colors
        Facecolors (backgrounds) of bbox patches for text
        labels.
    
    """
    
    line_kw = {"linewidth": 1.5*mpl.rcParams["lines.linewidth"],
               "color": "k", "linestyle": (0, (5, 5)),
               "zorder": get_zorder("ref_lats_lines"),
               "transform": ccrs.PlateCarree()}
    
    annotate_kw = {"fontsize": plot_style.fs_0 - 2.5,
                   "color": "k", "ha": "center", "va": "center",
                   "bbox": {"alpha": 1, "edgecolor": "none",
                            "pad": 0},
                   "zorder": get_zorder("ref_lats_lines"),
                   "transform": ccrs.PlateCarree()}
    
    nlabs = len(ref_lats_labels_lons)
    
    for j in range(len(ref_lats)):
        
        # Reference latitudes are lines of constant latitude;
        # need to specify longitudes of a fairly large number
        # for ax.plot() to 'join up' along:
        X_circ = np.arange(0, 361, 2)
        Y_circ = ref_lats[j]*np.ones(np.shape(X_circ))
        ax.plot(X_circ, Y_circ, **line_kw)
        
        if hemi == "n":
            rotation = (ref_lats_labels_lons[j%nlabs]
                        - central_longitude)
        else:
            # Need to account for rotation of figure itself and
            # southern ocean inverted perspective...
            rotation = (central_longitude
                        - ref_lats_labels_lons[j%nlabs]
                        - 180.0)
        
        annotate_kw["bbox"]["facecolor"] = \
            lon_labels_facecolors[j%len(lon_labels_facecolors)]
        
        ax.annotate(
            f"{abs(ref_lats[j]):.0f}" + u"\u00b0" + hemi.upper(),
            (ref_lats_labels_lons[j%nlabs], ref_lats[j]),
            rotation=rotation, **annotate_kw)



def main():
    
    # Decrease standard font size slightly:
    plot_style.fs_0 -= 1
    plot_style.set_mpl_rcParams()
    
    land_color = "lightgray"
    ocean_color = "#009ABB"
    
    prsr = script_tools.argument_parser("Create EBM schematics")
    script_tools.add_cmip_selection_cmd_args(prsr)
    script_tools.add_plotting_cmd_args(prsr, n_panels=0,
        n_cbars=0, text_labels=[])
    
    cmd, _ = script_tools.get_args(prsr,
        suppress_options=["n_models", "yravg1", "yravg2",
                          "exp"])
    
    ref_lats = list(cmd.reflats)
    hemi = "s" if (ref_lats[0] < 0 and ref_lats[1] < 0) else "n"
    
    # Dotted lines are drawn at references latitudes; only in
    # Arctic case is second needed / if second reference
    # latitude is not the pole:
    if abs(ref_lats[1]) > 85.0:
        del ref_lats[1]
    if abs(ref_lats[0]) > 85.0:
        del ref_lats[0]
    
    # Various specifications differ per hemisphere:
    if hemi == "s":
        
        # Arguments passed to cartopy.ccrs.NearsidePerspective
        # and also required to determine rotation of some text
        # labels (hence these not hardcoded in the setup_globe
        # function):
        central_longitude = 150
        central_latitude = -65
        
        # Longitudes at which reference latitude labels are
        # written, list of values for each reference latitude:
        ref_lats_labels_lons = [181, 137]
        
        # Background colours of those text labels (depends
        # where they are drawn!):
        lon_labels_facecolors = [ocean_color, land_color]
        
        # Longitudes at which to draw AHT and OHT arrows, a list
        # of arrays or lists for each reference latitude:
        aht_lons = [np.arange(-135, 180, 60)]*2
        oht_lons = [np.arange(-165, 180, 60), [-165, -45]]
        
    else:
        # As above but for the Arctic (hemi == "n"):
        central_longitude = 330
        central_latitude = 65
        ref_lats_labels_lons = [313.7]
        lon_labels_facecolors = [land_color]
        
        aht_lons = [np.arange(40, 360, 60)]
        oht_lons = [[191, 302, 330, 358]]
    
    # ------------------------------------------------------- #
    
    fig_ar = 1.2  # aspect ratio (width/height)
    
    fig = plt.figure(figsize=(plot_style.fig_width_single_in,
        plot_style.fig_width_single_in/fig_ar))
    fig.set_facecolor("none")
    
    # Setup geoAxes for the globe, to which sea ice, reference
    # latitude, and OHT/AHT arrows are drawn:
    ax_globe = setup_globe(fig, hemi=hemi,
        central_longitude=central_longitude,
        central_latitude=central_latitude,
        ocean_color=ocean_color, land_color=land_color)
    
    # Setup a 'guide' axes for other annotations:
    ax_guide = setup_guide_axes(fig)
    
    # Add sea ice climatology data:
    sea_ice_color = add_sea_ice_climatology(ax_globe, hemi=hemi)
    
    add_reference_latitudes(ax_globe, ref_lats=ref_lats,
        ref_lats_labels_lons=ref_lats_labels_lons, hemi=hemi,
        central_longitude=central_longitude,
        lon_labels_facecolors=lon_labels_facecolors)
    
    if hemi == "n":
        add_arctic_annotations(ax_guide, fig_ar=fig_ar,
            sea_ice_color=sea_ice_color)
    else:
        add_southern_ocean_annotations(ax_guide, fig_ar=fig_ar,
            central_longitude=central_longitude,
            ref_lat_label_lon=ref_lats_labels_lons[1])
    
    add_heat_transports(ax_globe, ref_lats=ref_lats,
                        oht_lons=oht_lons, aht_lons=aht_lons)
    
    # Only save png copy for Arctic (Fig. 2); for Southern Ocean
    # (Fig. S4), needs to be rotated first. These figures must
    # also have embedded raster graphics for the sea ice
    # specifically, which is set by being below a rasterization
    # zorder (see _zorder dictionary):
    plot_tools.finish_fig(fig, savefig=cmd.savefig,
        file_fmts=[".svg"] + ([".png"] if hemi=="n" else []),
        set_raster_level=True, file_name=cmd.savefigname,
        fig_metadata={"Title": cmd.savefigtitle})



if __name__ == "__main__":
    main()
