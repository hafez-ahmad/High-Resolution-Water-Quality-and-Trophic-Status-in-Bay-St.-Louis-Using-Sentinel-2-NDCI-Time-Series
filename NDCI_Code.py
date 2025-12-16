# %%
import pandas as pd
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
import geemap,ee

from datetime import datetime

# %%
# # authenticate to the Earth Engine API
ee.Authenticate()
ee.Initialize()

# %%
import geopandas as gpd

# %%
# read from geodatabase C:\Users\hafez\MSU\GISProjects\BayStLouise\BayStLouise\BayStLouise.gdb
gdb_path = "projects/ee-wfahydrologymsu100/assets/BayStLouise"
# feature_collection = ee.FeatureCollection(gdb_path)
aoi = ee.FeatureCollection(gdb_path)

# Cloud masking using Scene Classification Layer (SCL)
def mask_s2_clouds(image):
    scl = image.select('SCL')
    mask = scl.neq(3).And(scl.neq(7)).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))  # Remove cloud, shadows, saturated
    return image.updateMask(mask).copyProperties(image, ["system:time_start"])

# Function to calculate NDCI
def compute_ndci(image):
    red = image.select('B4')      # 665 nm
    red_edge = image.select('B5') # 705 nm
    ndci = red_edge.subtract(red).divide(red_edge.add(red)).rename('NDCI')
    return image.addBands(ndci).clip(aoi)

# Define time range
start = ee.Date('2017-03-01')
end = ee.Date('2025-07-01')  # July exclusive

# Build monthly composites
def make_monthly_composites(start_date, end_date):
    composites = []
    date = start_date
    while date.millis().lt(end_date.millis()).getInfo():
        next_month = date.advance(1, 'month')
        filtered = ee.ImageCollection('COPERNICUS/S2_SR') \
            .filterBounds(aoi) \
            .filterDate(date, next_month) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) \
            .map(mask_s2_clouds) \
            .map(compute_ndci)

        composite = filtered.median().set({
            'system:time_start': date.millis(),
            'date': date.format('YYYY-MM')
        })
        composites.append(composite)
        date = next_month
    return composites

# Create list of monthly composites
monthly_list = make_monthly_composites(start, end)

# Export each image to Google Drive
for img in monthly_list:
    date_str = img.get('date').getInfo()
    task = ee.batch.Export.image.toDrive(
        image=img.select('NDCI'),
        description=f'NDCI_{date_str}',
        folder='NDCI_Exports',
        fileNamePrefix=f'NDCI_{date_str}',
        region=aoi.geometry().bounds().getInfo()['coordinates'],
        scale=10,
        fileFormat='GeoTIFF',
        maxPixels=1e13
    )
    task.start()
    print(f"Exporting NDCI for {date_str}...")

# %%
# NDCI with classification
def compute_ndci(image):
    red = image.select('B4')
    red_edge = image.select('B5')
    ndci = red_edge.subtract(red).divide(red_edge.add(red)).rename('NDCI')
    return image.addBands(ndci).clip(aoi)

# Load and filter Sentinel-2 collection
s2 = ee.ImageCollection('COPERNICUS/S2_SR') \
    .filterBounds(aoi) \
    .filterDate('2017-03-01', '2025-07-01') \
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) \
    .map(mask_s2_clouds) \
    .map(compute_ndci)

# Generate average image per month (all years)
def monthly_average(month):
    month = ee.Number(month).int()
    filtered = s2.filter(ee.Filter.calendarRange(month, month, 'month'))
    mean_img = filtered.select('NDCI').mean() \
        .set({'month': month.format('%02d')})
    return mean_img.clip(aoi)

# Create ImageCollection of 12 monthly averages
months = ee.List.sequence(1, 12)
monthly_avg_ic = ee.ImageCollection(months.map(monthly_average))

# Export each monthly average image to Drive
for m in range(1, 13):
    month_str = f"{m:02d}"
    img = monthly_avg_ic.filter(ee.Filter.eq('month', month_str)).first()
    
    task = ee.batch.Export.image.toDrive(
        image=img,
        description=f'Avg_NDCI_Month_{month_str}',
        folder='NDCI_Monthly_Averages',
        fileNamePrefix=f'Avg_NDCI_{month_str}',
        region=aoi.geometry().bounds().getInfo()['coordinates'],
        scale=10,
        fileFormat='GeoTIFF',
        maxPixels=1e13
    )
    task.start()
    print(f"Export started for month {month_str}")

# %%
file=r"C:\Users\hafez\MSU\GISProjects\BayStLouise\BayStLouise\BayStLouise.gdb"
name="SamplingPoint"
data= gpd.read_file(file, layer=name)
# remove the geometry column
data = data.drop(columns='geometry')
# Convert wide format to tidy/long format
df_long = data.melt(
    id_vars=['Site'],                # keep the 'Site' column fixed
    var_name='Date',                 # move wide columns to 'Date'
    value_name='NDCI'                # measurement value column
)

# OPTIONAL: Convert 'Date' column from e.g., 'NDCI_2019_01' to datetime
df_long['Date'] = df_long['Date'].str.replace('NDCI_', '', regex=False)
df_long['Date'] = pd.to_datetime(df_long['Date'], format='%Y_%m')
# make tidy
df_long.head( )

# %%
# save the tidy data to a CSV file
output_csv = r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\NDCI_Tidy.csv"
df_long.to_csv(output_csv, index=False)

# %%
df_long.describe()

# %%

# Ensure proper datetime format
df_long['Date'] = pd.to_datetime(df_long['Date'])

# Get first 5 unique sites
sites = df_long['Site'].unique()[:5]

# Create subplots (2 columns layout)
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 5))
axes = axes.flatten()

# Plot histograms for each site
for i, site in enumerate(sites):
    ax = axes[i]
    site_data = df_long[df_long['Site'] == site]['NDCI']
    
    ax.hist(site_data, bins=50, color='steelblue', edgecolor='black')
    ax.set_title(f"({chr(97+i)}) Transect {site}", loc='left', fontsize=12, weight='bold')
    ax.set_xlabel("NDCI", fontsize=12, weight='bold')
    ax.set_ylabel("Frequency", fontsize=12, weight='bold')

    stats = site_data.describe()
    annotation = (
        f"Min: {stats['min']:.2f}\n"
        f"Max: {stats['max']:.2f}\n"
        f"Mean: {stats['mean']:.2f}\n"
        f"Std: {stats['std']:.2f}"
    )
    ax.text(0.98, 0.95, annotation, transform=ax.transAxes,
            ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Remove unused subplot if odd number
if len(sites) < len(axes):
    fig.delaxes(axes[-1])

#fig.suptitle("NDCI Distribution for First Five Sites", fontsize=14, weight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save figure
output_path = "C:\\Users\\hafez\\MSU\\Research\\algal\\baystlouise\\paper\\figure\\NDCI_Histograms_First5Sites_500dpi.jpg"
plt.savefig(output_path, dpi=500)
plt.show()

output_path


# %%
from scipy.interpolate import make_interp_spline
# Group by Site and Date, then average NDCI (to remove duplicate values per site-date)
df_avg = df_long.groupby(['Site', 'Date'])['NDCI'].mean().reset_index()

# Get the first 5 unique sites
sites_avg = df_avg['Site'].unique()[:5]

# Set up the plot
plt.figure(figsize=(12, 6))
sns.set(style="whitegrid")

# Define line styles and colors
line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
colors = sns.color_palette("tab10", n_colors=len(sites_avg))

# Loop and plot smoothed splines
for i, site in enumerate(sites_avg):
    site_data = df_avg[df_avg['Site'] == site].sort_values('Date')
    x = site_data['Date']
    y = site_data['NDCI']

    x_numeric = x.map(datetime.toordinal)

    if len(x_numeric) > 3:
        x_smooth = np.linspace(x_numeric.min(), x_numeric.max(), 300)
        spline = make_interp_spline(x_numeric, y, k=3)
        y_smooth = spline(x_smooth)
        x_smooth_dates = [datetime.fromordinal(int(val)) for val in x_smooth]

        plt.plot(
            x_smooth_dates,
            y_smooth,
            label=f"{i+1}",
            linewidth=2,
            linestyle=line_styles[i % len(line_styles)],
            color=colors[i % len(colors)]
        )

# Customize plot
#plt.title("Smoothed NDCI Time Series (Site-Averaged)\nDistinct Colors and Line Styles", fontsize=14, weight='bold')
plt.xlabel("Date", fontsize=12, weight='bold')
plt.ylabel("NDCI", fontsize=12, weight='bold')
plt.legend(title="Transect", loc='upper right')
plt.tight_layout()

# Save the figure
styled_legend_path = "C:\\Users\\hafez\\MSU\\Research\\algal\\baystlouise\\paper\\figure\\NDCI_TimeSeries_SiteAveraged_Spline_Styled_Legend_First5Sites_500dpi.jpg"
plt.savefig(styled_legend_path, dpi=500)
plt.show()

styled_legend_path


# %%
output_csv = r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\NDCI_Tidy.csv"
df= pd.read_csv(output_csv)
# find max and min and mean of NDCI from 2019 july and august
df['Date'] = pd.to_datetime(df['Date'])
# Filter for July and August 2019
july_august_2019 = df[(df['Date'].dt.year == 2019) & (df['Date'].dt.month.isin([7, 8]))]
# Calculate min, max, and mean NDCI  first july and then august
july_stats = july_august_2019[july_august_2019['Date'].dt.month == 7]['NDCI'].agg(['min', 'max', 'mean'])
august_stats = july_august_2019[july_august_2019['Date'].dt.month == 8]['NDCI'].agg(['min', 'max', 'mean'])
print("July 2019 NDCI Stats:")
print(july_stats)
print("\nAugust 2019 NDCI Stats:")
print(august_stats)

# %%
df_avg.head(10)

# %%
# apply normality test for each site
from scipy.stats import shapiro
def test_normality(df, site):
    site_data = df[df['Site'] == site]['NDCI']
    stat, p_value = shapiro(site_data)
    return stat, p_value
# Test normality for each of the first 5 sites
normality_results = {}
for site in sites:
    stat, p_value = test_normality(df_avg, site)
    normality_results[site] = {
        'Shapiro-Wilk Statistic': stat,
        'p-value': p_value
    }
# Convert results to DataFrame for better readability
normality_df = pd.DataFrame(normality_results).T
normality_df

# %% [markdown]
# All five sites reject the null hypothesis of normality (p < 0.05). This confirms the earlier choice of the non-parametric Kruskal–Wallis test for comparing NDCI distributions across sites is appropriate.

# %%
import pandas as pd
import numpy as np
from scipy.stats import kruskal, f_oneway
from itertools import combinations

# Reload averaged data
df_avg['Month'] = df_avg['Date'].dt.month
df_avg['Year'] = df_avg['Date'].dt.year

# 1. Site with the most variability
site_variability = df_avg.groupby('Site')['NDCI'].std().sort_values(ascending=False)
most_variable_site = site_variability.idxmax()

# 2. Date(s) with most variability across sites
date_variability = df_avg.groupby('Date')['NDCI'].std().sort_values(ascending=False)
most_variable_date = date_variability.idxmax()
most_variable_month = pd.to_datetime(most_variable_date).strftime('%B %Y')

# 3. Sites with highest and lowest NDCI (on average)
mean_ndci_by_site = df_avg.groupby('Site')['NDCI'].mean().sort_values()
lowest_mean_site = mean_ndci_by_site.idxmin()
highest_mean_site = mean_ndci_by_site.idxmax()

# 4. Kruskal-Wallis test for similarity (non-parametric due to unknown distribution)
groups = [group['NDCI'].values for name, group in df_avg.groupby('Site')]
kw_stat, kw_p = kruskal(*groups)

# 5. Classification function
def classify_ndci(ndci):
    if ndci < -0.1:
        return 'Clean / Oligotrophic'
    elif ndci < 0.0:
        return 'Low algae concentration'
    elif ndci < 0.1:
        return 'Moderate eutrophic'
    elif ndci < 0.2:
        return 'High eutrophic'
    elif ndci < 0.4:
        return 'Very high, bloom likely'
    elif ndci < 0.5:
        return 'Algal bloom'
    else:
        return 'Severe bloom (Hypereutrophic)'

df_avg['NDCI_Class'] = df_avg['NDCI'].apply(classify_ndci)

# Count classification frequency by site
classification_counts = df_avg.groupby(['Site', 'NDCI_Class']).size().unstack(fill_value=0)



# Prepare summary
summary = {
    "Most Variable Site (STD)": most_variable_site,
    "Date with Most Variability Across Sites": most_variable_month,
    "Highest Mean NDCI Site": highest_mean_site,
    "Lowest Mean NDCI Site": lowest_mean_site,
    "Kruskal-Wallis H-statistic": kw_stat,
    "Kruskal-Wallis p-value": kw_p
}

summary


# %%
classification_counts 

# %%
# save classification counts to CSV
classification_counts.to_csv(r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\NDCI_Classification_Counts.csv")

# %%
import xarray as xr
import rioxarray

# %%

import rasterio
from rasterio.plot import show

from matplotlib.colors import ListedColormap, BoundaryNorm

# Directory containing classified monthly rasters
folder = r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\ndci\Classify"  # Adjusted to match upload; update path if local

# Ordered list of months for consistent labeling
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Corresponding filenames
filenames = [f"CopyRaster_OutRasterdataset_Avg_NDCI_{i:02d}_Classified.tif" for i in range(1, 13)]

# Define classification colormap and labels
class_colors = ['#003366', '#3366cc', '#66cc66', '#ffff66', '#ff9900', '#ff3300', '#990099']
class_labels = [
    '1: Oligotrophic',
    '2: Low eutrophic',
    '3: Moderate eutrophic',
    '4: High eutrophic',
    '5: Very high (bloom likely)',
    '6: Algal bloom',
    '7: Severe bloom'
]
# Re-plot using imshow without norm, since this is classified integer raster data
fig, axes = plt.subplots(3, 4, figsize=(16, 10))
axes = axes.flatten()

last_img = None

for idx, (ax, fname, month) in enumerate(zip(axes, filenames, months)):
    filepath = os.path.join(folder, fname)
    if os.path.exists(filepath):
        with rasterio.open(filepath) as src:
            data = src.read(1)
            # masked out nan
            data = np.where(data == 0, np.nan, data)
            img = ax.imshow(data, cmap=cmap, vmin=1, vmax=7)  # direct range for classified categories
            last_img = img
            ax.set_title(month, fontsize=10)
            ax.axis('off')
    else:
        ax.set_visible(False)

# Colorbar using discrete values
if last_img is not None:
    cbar = fig.colorbar(last_img, ax=axes, orientation='vertical', fraction=0.015, pad=0.02)
    cbar.set_ticks(np.arange(1, 8))  # 1 to 7
    cbar.set_ticklabels([label.split(":")[1] for label in class_labels])
    cbar.set_label("NDCI Class")


# %%
import glob
import os
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Set font style
plt.rcParams.update({'font.size': 12, 'font.weight': 'bold'})

# Directory containing monthly average NDCI rasters (March 2017–June 2025)
ndci_folder = r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\ndci\monthl"
titles = ['Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'Jun.',
          'Jul.', 'Aug.', 'Sep.', 'Oct.', 'Nov.', 'Dec.']

# List all monthly GeoTIFFs
ndci_tifs = sorted(glob.glob(os.path.join(ndci_folder, "*.tif")))

# Create 4×3 subplot grid
fig, axs = plt.subplots(4, 3, figsize=(15, 8))

# Normalization range for NDCI
vmin, vmax = -0.2, 0.25
norm = Normalize(vmin=vmin, vmax=vmax)

# Colormap
cmap = plt.colormaps['jet']

# Store image handles for shared colorbar
images = []

# Plot each month's raster
for i, raster_path in enumerate(ndci_tifs):
    with rasterio.open(raster_path) as src:
        array = src.read(1)
        array = np.ma.masked_where(np.isnan(array) | (array == 0), array)  # Mask invalid
        ax = axs[i // 3, i % 3]
        ax.set_facecolor('white')
        im = ax.imshow(array, cmap=cmap, norm=norm, interpolation='none')
        images.append(im)
        ax.set_title(titles[i], fontweight='bold')
        ax.axis('off')

# Add colorbar (continuous)
cbar = fig.colorbar(images[0], ax=axs.ravel().tolist(), orientation='vertical',
                    label='Normalized Difference Chlorophyll Index (NDCI)')
cbar.ax.set_ylabel('NDCI', fontweight='bold')

# Save figure
output_path = r'C:\Users\hafez\MSU\Research\Figures\Chapter2_NDCI_monthly_map.jpg'
plt.savefig(output_path, dpi=500)
plt.show()
# Caption: Average Monthly NDCI from March 2017 to June 2025.Remarks: The NDCI values are normalized between -0.2 and 0.25, with a color gradient from blue (low NDCI) to red (high NDCI). The maps show seasonal variations in chlorophyll concentration, indicating potential algal blooms during warmer months. 

# %%
import glob
import os
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Set font style
plt.rcParams.update({'font.size': 12, 'font.weight': 'bold'})

# Directory and years
ndci_folder = r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\ndci\NDCI_Exports-20250614T044459Z-1-001"
titles = [2019, 2020, 2021, 2022, 2023, 2024]

# List and filter GeoTIFFs
ndci_tifs = sorted(glob.glob(os.path.join(ndci_folder, "*.tif")))
ndci_tifs = [f for f in ndci_tifs if any(str(year) in f for year in titles)]

# Setup plot
fig, axs = plt.subplots(3, 2, figsize=(18, 10))
axs = axs.flatten()

vmin, vmax = -0.2, 0.25
norm = Normalize(vmin=vmin, vmax=vmax)
cmap = plt.colormaps['jet']
images = []

for idx, year in enumerate(titles):
    year_tifs = [f for f in ndci_tifs if str(year) in f]
    if not year_tifs:
        continue

    arrays = []
    for tif in year_tifs:
        with rasterio.open(tif) as src:
            data = src.read(1).astype('float32')
            data[data == 0] = np.nan
            arrays.append(data)

    # Calculate yearly average
    avg_array = np.nanmean(arrays, axis=0)
    avg_array = np.ma.masked_invalid(avg_array)
    # print yearly max and min
    print(f"Year: {year}, Max NDCI: {np.nanmax(avg_array):.4f}, Min NDCI: {np.nanmin(avg_array):.4f}")
    ax = axs[idx]
    ax.set_facecolor('white')
    im = ax.imshow(avg_array, cmap=cmap, norm=norm)
    images.append(im)
    ax.set_title(f"{year}", fontweight='bold')
    ax.axis('off')

# Hide any unused subplots
for j in range(len(titles), len(axs)):
    axs[j].axis('off')

# Add shared colorbar
cbar = fig.colorbar(images[0], ax=axs.ravel().tolist(), orientation='vertical')
cbar.ax.set_ylabel('NDCI', fontweight='bold')
plt.suptitle('Yearly Average NDCI (2019–2024)', fontsize=16, fontweight='bold')
# Save figure
output_path = r'C:\Users\hafez\MSU\Research\Figures\Chapter2_NDCI_yearly_map.jpg'
plt.savefig(output_path, dpi=500)
plt.show()
# Caption: Yearly Average NDCI from 2019 to 2024. Remarks: The NDCI values are normalized between -0.2 and 0.25, with a color gradient from blue (low NDCI) to red (high NDCI). The maps show annual variations in chlorophyll concentration, indicating potential algal blooms during warmer months.

# %%
import os
import glob
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
# Add colorbar to the right of all subplots, not overlapping
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.tight_layout(rect=[0, 0, 0.93, 1])  # leave space for colorbar on the right



# Set font style
plt.rcParams.update({'font.size': 12, 'font.weight': 'bold'})

# Directory containing NDCI files
ndci_folder = r"C:\Users\hafez\MSU\Research\algal\baystlouise\paper\data\ndci\NDCI_Exports-20250614T044459Z-1-001"

# Select July and August 2019 files
target_months = ['NDCI_2019_07', 'NDCI_2019_08']
ndci_files = [f for f in glob.glob(os.path.join(ndci_folder, "*.tif")) if any(month in f for month in target_months)]

# Plot
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
vmin, vmax = -0.2, 0.25
norm = Normalize(vmin=vmin, vmax=vmax)
cmap = plt.colormaps['jet']
images = []
annotation =['(a)', '(b)']

for i, tif in enumerate(ndci_files):
    with rasterio.open(tif) as src:
        array = src.read(1).astype('float32')
        array = np.ma.masked_where((array == 0) | np.isnan(array), array)
        print(f"File: {os.path.basename(tif)}, Max NDCI: {np.nanmax(array):.4f}, Min NDCI: {np.nanmin(array):.4f}")
        ax = axs[i]
        ax.set_facecolor('white')
        im = ax.imshow(array, cmap=cmap, norm=norm)
        images.append(im)
        ax.set_title(os.path.basename(tif).split('.')[0], fontweight='bold')
        ax.text(0.01, 0.99, annotation[i], transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        ax.axis('off')

# Add colorbar
# Create a new axis for the colorbar
divider = make_axes_locatable(axs[-1])
cax = divider.append_axes("right", size="5%", pad=0.15)

cbar = fig.colorbar(images[0], cax=cax, orientation='vertical')
cbar.ax.set_ylabel('NDCI', fontweight='bold')
# save the figure
output_path = r'C:\Users\hafez\MSU\Research\Figures\Chapter2_NDCI_July_August_2019.jpg'
plt.savefig(output_path, dpi=500, bbox_inches='tight')
plt.tight_layout()
plt.show()
# Caption: NDCI Maps for July and August 2019. Remarks: The NDCI values are normalized between -0.2 and 0.25, with a color gradient from blue (low NDCI) to red (high NDCI). The maps show seasonal variations in chlorophyll concentration, indicating potential algal blooms during warmer months.


