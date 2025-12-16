# NDCI-BSL: Sentinel-2 NDCI Time Series for Bay St. Louis

NDCI-BSL is a set of Google Earth Engine (GEE) and Python tools designed to compute, analyze, and visualize the **Normalized Difference Chlorophyll Index (NDCI)** from Sentinel-2 MSI imagery over **Bay St. Louis, Mississippi (USA)**.

The workflow supports:
- Generation of **monthly NDCI composites** (March 2017 – June 2025),
- **Trophic status classification** across multiple estuarine transects,
- Extraction of **time series and statistics** for hydrologically distinct zones.

This repository accompanies the publication:

> Ahmad, H. (2025). *High‐Resolution Spatiotemporal Monitoring of Water Quality and Trophic Status in Bay St. Louis Using Sentinel‐2 NDCI Time Series on Google Earth Engine.* Transactions in GIS, 29(8), e70166.

---

## Requirements and Usage

The core workflow is implemented in **Google Earth Engine** with optional **Python** scripts for local analysis and plotting.

### Google Earth Engine

You need a GEE account to run the scripts.

1. Log in to the GEE Code Editor and set up authenticationtication.  
2. Load Python script in Text Editor.  
3. Set your:
   - Study region (Bay St. Louis AOI),
   - Time range,
   - Output folder on Google Drive or Earth Engine assets.
4. Run the script to generate monthly NDCI composites and zonal statistics.
5.  Download and run post analysis
