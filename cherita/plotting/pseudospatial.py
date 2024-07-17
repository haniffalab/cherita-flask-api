from __future__ import annotations
from typing import Literal
import base64
import json
import zarr
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from cherita.utils.adata_utils import (
    parse_data,
    get_index_in_array,
    get_group_index,
    get_row_from_zarr_df,
)
from cherita.resources.errors import BadRequest


def pseudospatial_gene(
    adata_group: zarr.Group,
    marker_id: str = None,
    marker_name: str = None,
    mask: str = "spatial",
    var_names_col: str = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):

    if not marker_id and not marker_name:
        raise BadRequest("Either 'varId' or 'varName' must be provided")

    if "masks" not in adata_group.uns or mask not in adata_group.uns["masks"]:
        raise KeyError(f"Mask {mask} not found in adata")

    if "polygons" not in adata_group.uns["masks"][mask] or not len(
        adata_group.uns["masks"][mask]["polygons"]
    ):
        raise ValueError(f"No polygons found in mask {mask}")

    if plot_format not in ["png", "svg", "html"]:
        raise BadRequest("Invalid format. Must be one of 'png', 'svg', 'html'")

    if marker_name:
        if var_names_col:
            marker_idx = get_index_in_array(adata_group.var[var_names_col], marker_name)
            marker_id = get_group_index(adata_group.var)[marker_idx]
        else:
            marker_id = marker_name

    if "varm" in adata_group.uns["masks"][mask].keys():  # precomputed mean expression
        cols = parse_data(
            adata_group.obs[adata_group.uns["masks"][mask]["varm"][()]]["categories"]
        )
        values_dict = get_row_from_zarr_df(
            adata_group.varm[adata_group.uns["masks"][mask]["varm"][()]],
            marker_id,
            cols,
        )
    else:
        raise ValueError(
            "No precomputed mean expression found in dataset for requested masks"
        )

    return plot_polygons(
        adata_group, values_dict, mask, plot_format=plot_format, **kwargs
    )


def plot_polygons(
    adata_group: zarr.Group,
    values_dict: pd.DataFrame,
    mask: str = "spatial",
    colormap: str = "viridis",
    show_colorbar: bool = True,
    min_value: float = None,
    max_value: float = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    width: int = 500,
    height: int = 500,
    full_html: bool = False,
):
    min_value = float(min_value) if min_value else min(values_dict.values())
    max_value = float(max_value) if max_value else max(values_dict.values())

    normalized_values = {
        k: (v - min_value) / (max_value - min_value) for k, v in values_dict.items()
    }
    color_values = {
        k: sample_colorscale(colormap, [v], colortype="rgb")[0]
        for k, v in normalized_values.items()
    }
    fig = go.Figure()

    for polygon in adata_group.uns["masks"][mask]["polygons"].keys():
        fig.add_trace(
            go.Scatter(
                x=list(
                    *adata_group.uns["masks"][mask]["polygons"][polygon][:, :, 0, 0]
                ),
                y=list(
                    *adata_group.uns["masks"][mask]["polygons"][polygon][:, :, 0, 1]
                ),
                line=dict(
                    color=(color_values[polygon]),
                    width=1,
                ),
                mode="lines",
                showlegend=False,
                fill="toself",
                fillcolor=(color_values[polygon]),
                text="<br>".join(
                    [
                        f"<b>{polygon}</b>",
                        (
                            "Mean expression: 0"
                            if values_dict[polygon] == 0
                            else (
                                f"Mean expression: {values_dict[polygon]:.3e}"
                                if values_dict[polygon] < 0.001
                                else f"Mean expression: {values_dict[polygon]:,.3f}"
                            )
                        ),
                    ]
                ),
                hoverinfo="text",
            )
        )

    if isinstance(show_colorbar, str):
        show_colorbar = show_colorbar.lower() in ["true", "1"]
    if show_colorbar:
        colorbar_trace = go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            showlegend=False,
            marker=dict(
                colorscale=colormap,
                showscale=True,
                colorbar=dict(thickness=10, tickmode="auto"),
                cmin=min_value,
                cmax=max_value,
            ),
        )
        fig.add_trace(colorbar_trace)

    fig.update_xaxes(visible=False, fixedrange=True)
    fig.update_yaxes(
        visible=False, fixedrange=True, autorange="reversed", scaleanchor="x"
    )
    fig.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        width=int(width),
        height=int(height),
        modebar_remove=["select2d", "lasso2d"],
        margin=dict(pad=0, l=10, r=10, t=10, b=10),
    )

    if plot_format in ["png", "svg"]:
        img_bytes = fig.to_image(format=plot_format)
        img_str = base64.b64encode(img_bytes).decode()
        return img_str
    elif plot_format == "html":
        if isinstance(full_html, str):
            full_html = full_html.lower() in ["true", "1"]
        return fig.to_html(
            config=dict(displaylogo=False), include_plotlyjs="cdn", full_html=full_html
        )
    else:  # json
        json.loads(fig.to_json())
