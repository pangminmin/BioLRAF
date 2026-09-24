#!/usr/bin/env python3

"""BioLRAF GSE226365 3D ternary visualization notebook.

Input: Example/GSE226365/processed/GSE226365_gficf_scores.csv
Expected score columns: Os, Ch, Ad, MSC, and celltype.
The notebook calculates normalized lineage activities, MSC_ratio,
Diff_score, and lineage_dominant, then exports plots colored by celltype
and lineage_dominant.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


FIG_WIDTH = 300
FIG_HEIGHT = 300
EXPORT_SCALE = 2
DEFAULT_CENTER_DELTA = 0.10

LINEAGE_ORDER = [
    "Os-dominant",
    "Ch-dominant",
    "Ad-dominant",
    "Non-dominant",
]

LINEAGE_COLORS = {
    "Os-dominant": "#597EB7",
    "Ch-dominant": "#F29696",
    "Ad-dominant": "#F0CC63",
    "Non-dominant": "#9CD379",
}

# Four fixed views used by the BioLRAF result browser.
CAMERA_LIST = [
    dict(
        eye=dict(x=-0.092, y=-2.097, z=0.539),
        center=dict(x=-0.0106, y=-0.023, z=-0.078),
        up=dict(x=0.0, y=0.0, z=1.0),
        projection=dict(type="perspective"),
    ),
    dict(
        eye=dict(x=1.681, y=1.032, z=0.571),
        center=dict(x=0.0, y=0.0, z=0.0),
        up=dict(x=0.0, y=0.0, z=1.0),
        projection=dict(type="perspective"),
    ),
    dict(
        eye=dict(x=-1.658, y=1.053, z=0.598),
        center=dict(x=0.0, y=0.0, z=0.0),
        up=dict(x=0.0, y=0.0, z=1.0),
        projection=dict(type="perspective"),
    ),
    dict(
        eye=dict(x=0.001, y=-0.312, z=2.030),
        center=dict(x=0.0, y=0.0, z=0.0),
        up=dict(x=0.0, y=0.0, z=1.0),
        projection=dict(type="perspective"),
    ),
]


def read_score_table(path: Path) -> pd.DataFrame:
    """Read a BioLRAF score table and preserve a cell identifier index."""
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")

    df = pd.read_csv(path)

    if "cell_id" in df.columns:
        df["cell_id"] = df["cell_id"].astype(str)
        df = df.set_index("cell_id", drop=False)
    elif len(df.columns) > 0 and str(df.columns[0]).startswith("Unnamed:"):
        first_col = df.columns[0]
        df[first_col] = df[first_col].astype(str)
        df = df.set_index(first_col, drop=True)
        df.index.name = "cell_id"
        df["cell_id"] = df.index.astype(str)
    else:
        df.index = df.index.astype(str)
        df.index.name = "cell_id"
        df["cell_id"] = df.index

    return df


def calculate_biolraf_metrics(
    df: pd.DataFrame,
    center_delta: float = DEFAULT_CENTER_DELTA,
) -> pd.DataFrame:
    """Calculate BioLRAF normalized metrics and dominant lineage states."""
    required = ["Os", "Ch", "Ad", "MSC", "celltype"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    result = df.copy()
    for column in ["Os", "Ch", "Ad", "MSC"]:
        result[column] = pd.to_numeric(result[column], errors="coerce")

    result["total_score"] = result[["Os", "Ch", "Ad", "MSC"]].sum(axis=1)
    result["MSC_ratio"] = np.where(
        result["total_score"] > 0,
        result["MSC"] / result["total_score"],
        np.nan,
    )
    result["Diff_score"] = 1.0 - result["MSC_ratio"]

    result["sum_tri"] = result[["Os", "Ch", "Ad"]].sum(axis=1)
    for lineage in ["Os", "Ch", "Ad"]:
        result[f"{lineage}_n"] = np.where(
            result["sum_tri"] > 0,
            result[lineage] / result["sum_tri"],
            np.nan,
        )

    result["is_center"] = (
        (np.abs(result["Os_n"] - 1.0 / 3.0) <= center_delta)
        & (np.abs(result["Ch_n"] - 1.0 / 3.0) <= center_delta)
        & (np.abs(result["Ad_n"] - 1.0 / 3.0) <= center_delta)
    )

    conditions = [
        result["is_center"],
        (result["Os_n"] > result["Ch_n"])
        & (result["Os_n"] > result["Ad_n"]),
        (result["Ch_n"] > result["Os_n"])
        & (result["Ch_n"] > result["Ad_n"]),
        (result["Ad_n"] > result["Os_n"])
        & (result["Ad_n"] > result["Ch_n"]),
    ]
    choices = [
        "Non-dominant",
        "Os-dominant",
        "Ch-dominant",
        "Ad-dominant",
    ]
    result["lineage_dominant"] = np.select(
        conditions,
        choices,
        default="Non-dominant",
    )
    result["lineage_dominant"] = pd.Categorical(
        result["lineage_dominant"],
        categories=LINEAGE_ORDER,
        ordered=True,
    )

    required_complete = [
        "Os",
        "Ch",
        "Ad",
        "MSC",
        "Diff_score",
        "Os_n",
        "Ch_n",
        "Ad_n",
        "lineage_dominant",
        "celltype",
    ]
    return result.dropna(subset=required_complete).copy()


def barycentric_to_xyz(
    a,
    b,
    c,
    d,
    s0: float = 1.0,
    area_linear: bool = False,
    clip_d: bool = True,
):
    """Convert Os, Ch, Ad, and Diff_score values into 3D coordinates."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    d = np.asarray(d, dtype=float)

    total = a + b + c
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(total == 0, 1.0 / 3.0, a / total)
        q = np.where(total == 0, 1.0 / 3.0, b / total)
        r = np.where(total == 0, 1.0 / 3.0, c / total)

    d_plot = np.clip(d, 0, 1) if clip_d else d
    s = s0 * np.sqrt(d_plot) if area_linear else s0 * d_plot
    height = np.sqrt(3.0) / 2.0 * s

    vertex_ch = np.stack([np.zeros_like(s), 2.0 / 3.0 * height], axis=-1)
    vertex_os = np.stack([-s / 2.0, -1.0 / 3.0 * height], axis=-1)
    vertex_ad = np.stack([s / 2.0, -1.0 / 3.0 * height], axis=-1)

    weights = np.stack([q, p, r], axis=-1)[..., None]
    vertices = np.stack([vertex_ch, vertex_os, vertex_ad], axis=-2)
    xy = (weights * vertices).sum(axis=-2)

    return xy[..., 0], xy[..., 1], d_plot


def plot_3d_ternary(
    df: pd.DataFrame,
    color_by: str = "Diff_score",
    color_type: str = "auto",
    color_map: dict[str, str] | None = None,
    colorscale=None,
    s0: float = 1.2,
    area_linear: bool = False,
    draw_wireframes: bool = True,
    point_size: float = 4,
    legend_marker_size: float = 10,
    opacity: float = 0.75,
    cmin=None,
    cmax=None,
    na_color: str = "gray",
    colorbar_title: str | None = None,
    colorbar_len: float = 0.45,
    colorbar_thickness: float = 14,
    title: str | None = None,
    show: bool = False,
) -> go.Figure:
    """Create a BioLRAF 3D ternary projection using Plotly."""
    required = ["Os", "Ch", "Ad", "Diff_score"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if color_by not in df.columns:
        raise ValueError(f"color_by='{color_by}' is not present in the input data.")
    if color_type not in {"auto", "continuous", "discrete"}:
        raise ValueError("color_type must be auto, continuous, or discrete.")

    plot_df = df.copy()
    plot_df["x"], plot_df["y"], plot_df["z"] = barycentric_to_xyz(
        a=plot_df["Os"],
        b=plot_df["Ch"],
        c=plot_df["Ad"],
        d=plot_df["Diff_score"],
        s0=s0,
        area_linear=area_linear,
        clip_d=True,
    )
    if "cell_id" not in plot_df.columns:
        plot_df["cell_id"] = plot_df.index.astype(str)

    if color_type == "auto":
        color_type_use = (
            "continuous"
            if pd.api.types.is_numeric_dtype(plot_df[color_by])
            else "discrete"
        )
    else:
        color_type_use = color_type

    if colorscale is None:
        colorscale = [
            [0.00, "#000004"],
            [0.25, "#3B0F70"],
            [0.50, "#B63679"],
            [0.75, "#FB8861"],
            [1.00, "#FCFDBF"],
        ]

    hover_columns = []
    for column in [
        "Os",
        "Ch",
        "Ad",
        "MSC",
        "Diff_score",
        "MSC_ratio",
        "dataset",
        "celltype",
        color_by,
    ]:
        if column in plot_df.columns and column not in hover_columns:
            hover_columns.append(column)

    fig = go.Figure()

    if color_type_use == "continuous":
        valid = plot_df[plot_df[color_by].notna()].copy()
        missing_color = plot_df[plot_df[color_by].isna()].copy()
        customdata = valid[hover_columns].to_numpy() if hover_columns else None
        hovertemplate = _build_hover_template(plot_df, hover_columns)

        colorbar = dict(
            title=dict(
                text=colorbar_title or str(color_by),
                side="right",
                font=dict(size=12),
            ),
            len=colorbar_len,
            thickness=colorbar_thickness,
            x=0.95,
            y=0.50,
            tickfont=dict(size=11),
        )
        if cmin == 0 and cmax == 1:
            colorbar.update(
                tickmode="array",
                tickvals=[0, 0.25, 0.50, 0.75, 1.00],
                ticktext=["0", "0.25", "0.50", "0.75", "1.00"],
            )

        fig.add_trace(
            go.Scatter3d(
                x=valid["x"],
                y=valid["y"],
                z=valid["z"],
                mode="markers",
                text=valid["cell_id"],
                customdata=customdata,
                hovertemplate=hovertemplate,
                showlegend=False,
                marker=dict(
                    size=point_size,
                    color=valid[color_by],
                    colorscale=colorscale,
                    cmin=cmin,
                    cmax=cmax,
                    opacity=opacity,
                    colorbar=colorbar,
                ),
            )
        )

        if not missing_color.empty:
            fig.add_trace(
                go.Scatter3d(
                    x=missing_color["x"],
                    y=missing_color["y"],
                    z=missing_color["z"].fillna(0),
                    mode="markers",
                    text=missing_color["cell_id"],
                    customdata=missing_color[hover_columns].to_numpy(),
                    hovertemplate=hovertemplate,
                    name="NA",
                    marker=dict(size=point_size, color=na_color, opacity=opacity),
                )
            )
    else:
        values = plot_df[color_by].astype(str)
        observed = list(pd.unique(values))
        if color_by == "lineage_dominant":
            groups = [group for group in LINEAGE_ORDER if group in observed]
            groups += [group for group in observed if group not in groups]
        else:
            groups = observed

        palette = px.colors.qualitative.Set1 + px.colors.qualitative.Dark24
        colors = {group: palette[i % len(palette)] for i, group in enumerate(groups)}
        if color_map:
            colors.update(color_map)

        for group in groups:
            subset = plot_df[values == group].copy()
            hovertemplate = _build_hover_template(
                plot_df,
                hover_columns,
                group_label=f"{color_by}: {group}<br>",
            )
            fig.add_trace(
                go.Scatter3d(
                    x=subset["x"],
                    y=subset["y"],
                    z=subset["z"],
                    mode="markers",
                    text=subset["cell_id"],
                    customdata=subset[hover_columns].to_numpy(),
                    hovertemplate=hovertemplate,
                    name=group,
                    showlegend=False,
                    marker=dict(
                        size=point_size,
                        color=colors.get(group, "gray"),
                        opacity=opacity,
                    ),
                )
            )

        for group in groups:
            fig.add_trace(
                go.Scatter3d(
                    x=[None],
                    y=[None],
                    z=[None],
                    mode="markers",
                    name=group,
                    showlegend=True,
                    hoverinfo="skip",
                    marker=dict(
                        size=legend_marker_size,
                        color=colors.get(group, "gray"),
                        opacity=1,
                    ),
                )
            )

    if draw_wireframes:
        _add_wireframes(fig, s0=s0, area_linear=area_linear)
    _add_vertex_labels(fig, s0=s0)

    fig.update_layout(
        title=title or f"3D ternary projection colored by {color_by}",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Diff_score",
            zaxis=dict(range=[0, 1]),
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=0.8),
        ),
        legend=dict(title=str(color_by), x=1.02, y=1, itemsizing="constant"),
        margin=dict(l=0, r=0, b=0, t=50),
    )

    if show:
        fig.show()
    return fig


def _build_hover_template(
    plot_df: pd.DataFrame,
    hover_columns: list[str],
    group_label: str = "",
) -> str:
    template = (
        "cell: %{text}<br>"
        + group_label
        + "x: %{x:.4f}<br>y: %{y:.4f}<br>z: %{z:.4f}<br>"
    )
    for index, column in enumerate(hover_columns):
        if pd.api.types.is_numeric_dtype(plot_df[column]):
            template += f"{column}: %{{customdata[{index}]:.4f}}<br>"
        else:
            template += f"{column}: %{{customdata[{index}]}}<br>"
    return template + "<extra></extra>"


def _add_wireframes(fig: go.Figure, s0: float, area_linear: bool) -> None:
    for level in np.linspace(0.0, 1.0, 6):
        scale = s0 * np.sqrt(level) if area_linear else s0 * level
        height = np.sqrt(3.0) / 2.0 * scale
        vertex_ch = np.array([0.0, 2.0 / 3.0 * height, level])
        vertex_os = np.array([-scale / 2.0, -1.0 / 3.0 * height, level])
        vertex_ad = np.array([scale / 2.0, -1.0 / 3.0 * height, level])
        vertices = [vertex_ch, vertex_os, vertex_ad, vertex_ch]
        fig.add_trace(
            go.Scatter3d(
                x=[vertex[0] for vertex in vertices],
                y=[vertex[1] for vertex in vertices],
                z=[vertex[2] for vertex in vertices],
                mode="lines",
                line=dict(color="gray", width=3),
                showlegend=False,
                hoverinfo="skip",
            )
        )


def _add_vertex_labels(fig: go.Figure, s0: float) -> None:
    fig.add_trace(
        go.Scatter3d(
            x=[0],
            y=[0],
            z=[0],
            mode="markers",
            marker=dict(size=4, color="black"),
            showlegend=False,
            hoverinfo="skip",
        )
    )

    height = np.sqrt(3.0) / 2.0 * s0
    vertex_ch = np.array([0.0, 2.0 / 3.0 * height, 1.0])
    vertex_os = np.array([-s0 / 2.0, -1.0 / 3.0 * height, 1.0])
    vertex_ad = np.array([s0 / 2.0, -1.0 / 3.0 * height, 1.0])
    fig.add_trace(
        go.Scatter3d(
            x=[vertex_ch[0], vertex_os[0], vertex_ad[0]],
            y=[vertex_ch[1] + 0.05, vertex_os[1] - 0.05, vertex_ad[1] - 0.05],
            z=[1.0, 1.0, 1.0],
            mode="text",
            text=["Ch", "Os", "Ad"],
            textfont=dict(size=14, color="black"),
            showlegend=False,
            hoverinfo="skip",
        )
    )


def configure_static_layout(fig: go.Figure) -> go.Figure:
    """Return a main-plot-only figure for fixed PNG export."""
    static_fig = go.Figure(fig)
    hidden_axis = dict(
        showbackground=False,
        showgrid=False,
        zeroline=False,
        showline=False,
        showticklabels=False,
        ticks="",
        title="",
    )
    static_fig.update_layout(
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        autosize=False,
        margin=dict(l=0, r=0, t=0, b=0),
        title=None,
        showlegend=False,
        scene=dict(
            domain=dict(x=[0.0, 1.0], y=[0.0, 1.0]),
            bgcolor="rgba(0,0,0,0)",
            aspectmode="manual",
            aspectratio=dict(x=1.0, y=1.0, z=0.75),
            xaxis=hidden_axis,
            yaxis=hidden_axis,
            zaxis=hidden_axis,
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    return static_fig


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "group"


def export_total_dataset_plots(
    df: pd.DataFrame,
    output_dir: Path,
    dataset_label: str,
) -> None:
    """Export total-dataset plots colored by celltype and lineage_dominant."""
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_dataset = safe_filename(dataset_label)

    print("Cell counts by cell type:")
    print(df["celltype"].value_counts())
    print("\nDominant lineage states by cell type:")
    print(pd.crosstab(df["celltype"], df["lineage_dominant"]))

    plot_settings = [
        ("celltype", None),
        ("lineage_dominant", LINEAGE_COLORS),
    ]

    for color_by, color_map in plot_settings:
        print(f"\nProcessing total dataset colored by: {color_by}")
        plot_df = df.copy()
        if color_by == "lineage_dominant":
            plot_df = plot_df.sort_values("lineage_dominant")

        interactive_fig = plot_3d_ternary(
            df=plot_df,
            color_by=color_by,
            color_type="discrete",
            color_map=color_map,
            s0=1.2,
            area_linear=True,
            point_size=3,
            opacity=0.5,
            legend_marker_size=10,
            title=f"{dataset_label}: colored by {color_by}",
            show=False,
        )

        html_path = output_dir / (
            f"{safe_dataset}_3D_ternary_{safe_filename(color_by)}.html"
        )
        interactive_fig.write_html(
            html_path,
            include_plotlyjs=True,
            full_html=True,
        )
        print(f"Saved: {html_path}")

        static_fig = configure_static_layout(interactive_fig)
        for view_number, camera in enumerate(CAMERA_LIST, start=1):
            static_fig.update_layout(scene_camera=camera)
            png_path = output_dir / (
                f"{safe_dataset}_3D_ternary_{safe_filename(color_by)}_"
                f"fixed_{view_number}.png"
            )
            try:
                static_fig.write_image(
                    png_path,
                    width=FIG_WIDTH,
                    height=FIG_HEIGHT,
                    scale=EXPORT_SCALE,
                )
            except Exception as exc:
                raise RuntimeError(
                    "PNG export failed. Install a Plotly-compatible Kaleido "
                    "version in the current Python environment."
                ) from exc
            print(f"Saved: {png_path}")


def export_celltype_plots(
    df: pd.DataFrame,
    output_dir: Path,
    dataset_label: str,
    celltypes: list[str] | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    available = [str(value) for value in pd.unique(df["celltype"])]
    selected = available if not celltypes else celltypes

    print("Cell counts by cell type:")
    print(df["celltype"].value_counts())
    print("\nDominant lineage states by cell type:")
    print(pd.crosstab(df["celltype"], df["lineage_dominant"]))

    for celltype_name in selected:
        subset = df[df["celltype"].astype(str) == str(celltype_name)].copy()
        if subset.empty:
            print(f"Skipped {celltype_name}: no cells.")
            continue

        subset = subset.sort_values("lineage_dominant")
        safe_celltype = safe_filename(celltype_name)
        safe_dataset = safe_filename(dataset_label)
        print(f"\nProcessing cell type: {celltype_name} ({len(subset)} cells)")
        print(subset["lineage_dominant"].value_counts(sort=False))

        interactive_fig = plot_3d_ternary(
            df=subset,
            color_by="lineage_dominant",
            color_type="discrete",
            color_map=LINEAGE_COLORS,
            s0=1.2,
            area_linear=True,
            point_size=3,
            opacity=0.5,
            legend_marker_size=10,
            title=f"{dataset_label}: {celltype_name}",
            show=False,
        )

        html_path = output_dir / (
            f"{safe_dataset}_{safe_celltype}_3D_ternary_lineage_dominant.html"
        )
        interactive_fig.write_html(
            html_path,
            include_plotlyjs=True,
            full_html=True,
        )
        print(f"Saved: {html_path}")

        static_fig = configure_static_layout(interactive_fig)
        for view_number, camera in enumerate(CAMERA_LIST, start=1):
            static_fig.update_layout(scene_camera=camera)
            png_path = output_dir / (
                f"{safe_dataset}_{safe_celltype}_3D_ternary_"
                f"lineage_dominant_fixed_{view_number}.png"
            )
            try:
                static_fig.write_image(
                    png_path,
                    width=FIG_WIDTH,
                    height=FIG_HEIGHT,
                    scale=EXPORT_SCALE,
                )
            except Exception as exc:
                raise RuntimeError(
                    "PNG export failed. Install a Plotly-compatible Kaleido "
                    "version in the current Python environment."
                ) from exc
            print(f"Saved: {png_path}")


# ============================================================
# Notebook execution cell: load GSE226365 scores, calculate BioLRAF
# metrics, and export total-dataset ternary plots
# ============================================================

from pathlib import Path

# Find the repository root whether Jupyter was started from the repo root
# or from the Methods directory.
REPO_ROOT = Path.cwd()
if not (REPO_ROOT / "Example").exists():
    if (REPO_ROOT.parent / "Example").exists():
        REPO_ROOT = REPO_ROOT.parent

INPUT_CSV = (
    REPO_ROOT
    / "Example"
    / "GSE226365"
    / "processed"
    / "GSE226365_gficf_scores.csv"
)

OUTPUT_DIR = (
    REPO_ROOT
    / "Example"
    / "GSE226365"
    / "results"
    / "3D_ternary"
)

DATASET_LABEL = "GSE226365"
CENTER_DELTA = 0.10

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input CSV not found: {INPUT_CSV}\n"
        "Check the file path and filename."
    )

print(f"Input CSV: {INPUT_CSV}")
print(f"Output directory: {OUTPUT_DIR}")

# Read the cell-level score table generated by BioLRAF_analysis.R.
scores_raw = read_score_table(INPUT_CSV)

# The plotting functions require these four lineage/MSC score columns.
required_columns = ["Os", "Ch", "Ad", "MSC", "celltype"]
missing_columns = [
    column for column in required_columns
    if column not in scores_raw.columns
]
if missing_columns:
    raise ValueError(
        f"Missing required columns: {missing_columns}\n"
        f"Columns found in CSV: {scores_raw.columns.tolist()}\n"
        "Check that the gene-set score columns exported by the R script "
        "are named Os, Ch, Ad, and MSC, and that celltype is present."
    )

# Calculate normalized lineage activities, differentiation score,
# and dominant lineage state.
scores = calculate_biolraf_metrics(
    scores_raw,
    center_delta=CENTER_DELTA,
)

if scores.empty:
    raise ValueError(
        "No cells remain after metric calculation and missing-value filtering. "
        "Check the score columns and metadata."
    )

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Save a processed table containing the original scores and derived metrics.
processed_csv = (
    OUTPUT_DIR
    / f"{DATASET_LABEL}_BioLRAF_scores_with_states.csv"
)
scores.to_csv(processed_csv, index=False)
print(f"Saved processed score table: {processed_csv}")
print(f"Cells retained for plotting: {len(scores)}")

# Export two full-dataset views: colored by cell type and by dominant lineage.
export_total_dataset_plots(
    df=scores,
    output_dir=OUTPUT_DIR,
    dataset_label=DATASET_LABEL,
)

print("GSE226365 visualization completed.")
