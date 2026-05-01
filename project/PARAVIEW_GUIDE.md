# ParaView Visual Guide (Phase 1)

This guide explains how to turn exported CSV files into report figures.

## Files You Already Have

- Path bundle data: `output/path_bundle.csv`
- Barrier helper data: `output/barrier_plane.csv`

## A) Monte Carlo Path Bundle Figure

1. Open ParaView.
2. `File -> Open` and select `output/path_bundle.csv`.
3. Click `Apply`.
4. In `Properties`:
   - set `X Column = time`
   - set `Y Column = stock_price`
   - set `Z Column = path_id`
5. Convert points to lines:
   - `Filters -> Alphabetical -> Table To Points` (if needed)
   - then `Filters -> Alphabetical -> Plot Data Over Time` is not needed here
   - use `Filters -> Alphabetical -> Tube` after generating line representation
6. For connected path lines:
   - apply `Filters -> Temporal/Pathline` only if data is temporal.
   - for CSV static table, use `Filters -> Alphabetical -> Programmable Filter` or `Plot Over Line` alternatives if needed.
   - easiest workflow: use `Table To Points`, then `Glyph` or `Tube` to visually show ensemble structure.
7. Color by `path_id` (or solid color), adjust tube radius.
8. Save screenshot:
   - `File -> Save Screenshot`
   - save to `visuals/paraview_path_bundle.png`

## B) Barrier Geometry Figure

1. Load `output/path_bundle.csv` as above.
2. Load `output/barrier_plane.csv` in a second source.
3. For barrier data:
   - map `X=time`, `Y=barrier`, `Z=path_id`
   - display as lines (or tube with thin radius) so the barrier level appears as a sheet/guide.
4. Overlay both sources in one view.
5. Set barrier color to a strong contrasting color (e.g., red).
6. Save screenshot:
   - `visuals/paraview_barrier_overlay.png`

## C) Optional Animation (Path Evolution)

1. After loading path bundle, ensure animation controls are visible.
2. If using time-indexed data, use `Play` to preview.
3. `File -> Save Animation`
4. Recommended:
   - format: `.mp4`
   - frame rate: 20-30
   - save to `visuals/paraview_paths.mp4`

## D) Suggested Figure Captions

- Path bundle: "Simulated GBM path ensemble in \((t, S_t, i)\) space."
- Barrier overlay: "Barrier option geometry showing path trajectories relative to barrier level."

## Practical Notes

- Keep camera angle fixed across screenshots for consistency.
- Use a white background for print readability.
- If rendering is dense, downsample paths (e.g., 100-200) for cleaner visuals.
