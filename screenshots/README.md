# Screenshots

This folder contains screenshots of the AIQSAR application for use in the project README and documentation.

## Naming Convention

| Filename | Description |
|---|---|
| `01_data_explorer.png` | Data Explorer tab — file upload panel and missing-value heatmap |
| `02_correlation_heatmap.png` | Interactive correlation heatmap of molecular descriptors |
| `03_qsar_pipeline.png` | Full QSAR Pipeline results — model comparison table and validation metrics |
| `04_williams_plot.png` | Applicability Domain — Williams Plot |
| `05_model_grader.png` | AI-powered Model Grader output with letter grade and recommendation |
| `06_analysis_dashboard.png` | Consolidated Analysis Dashboard showing all key diagnostic plots |

## How to Contribute Screenshots

1. Run AIQSAR locally:
   ```r
   shiny::runApp("AIQSAR_1.R")
   ```
2. Navigate to each tab listed above.
3. Capture a screenshot at **1280 × 800** resolution (or similar 16:10 / 16:9 ratio) in **PNG** format.
4. Save the file using the exact filename from the table above.
5. Open a pull request adding your screenshot(s) to this folder.

## Guidelines

- Use a clean browser window with no personal data loaded.
- Prefer a light theme or the application's default theme.
- Crop out browser chrome (address bar, tabs) — capture only the app viewport.
- Keep file sizes reasonable (compress with tools like [Squoosh](https://squoosh.app/) if needed).
