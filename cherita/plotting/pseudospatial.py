from __future__ import annotations
from typing import Literal
import base64
import json
import zarr
import logging
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from cherita.utils.adata_utils import (
    parse_data,
    get_index_in_array,
    get_group_index,
    get_row_from_zarr_df,
)
from cherita.resources.errors import BadRequest, NotInData, InvalidKey


def validate_pseudospatial(adata_group: zarr.Group, mask: str):
    if "masks" not in adata_group.uns or mask not in adata_group.uns["masks"]:
        raise NotInData(f"Mask '{mask}' not found in adata")

    if "polygons" not in adata_group.uns["masks"][mask] or not len(
        adata_group.uns["masks"][mask]["polygons"]
    ):
        raise NotInData(f"No polygons found in mask {mask}")


def validate_format(format: str):
    if format not in ["png", "svg", "html"]:
        raise BadRequest(
            f"Invalid format '{format}'. Must be one of 'png', 'svg', 'html'"
        )


def pseudospatial_gene(
    adata_group: zarr.Group,
    marker_id: str = None,
    marker_name: str = None,
    mask: str = "spatial",
    var_names_col: str = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask)
    validate_format(plot_format)
    if not marker_id and not marker_name:
        raise BadRequest("Either 'varId' or 'varName' must be provided")

    try:
        if marker_name:
            if var_names_col:
                marker_idx = get_index_in_array(
                    adata_group.var[var_names_col], marker_name
                )
                marker_id = get_group_index(adata_group.var)[marker_idx]
            else:
                marker_id = marker_name
                marker_idx = get_index_in_array(
                    get_group_index(adata_group.var), marker_id
                )
    except (KeyError, InvalidKey):
        raise NotInData(f"Marker '{marker_name}' not found in dataset")

    if "varm" in adata_group.uns["masks"][mask].keys():
        # precomputed mean expression
        logging.info(f"Using precomputed mean expression of {marker_id}")
        obs_colname = adata_group.uns["masks"][mask]["obs"][()]
        cats = parse_data(adata_group.obs[obs_colname]).categories
        values_dict = get_row_from_zarr_df(
            adata_group.varm[adata_group.uns["masks"][mask]["varm"][()]],
            marker_id,
            cats,
        )
    else:
        # compute mean expression
        logging.info(f"Computing mean expression for {marker_id}")
        obs_colname = adata_group.uns["masks"][mask]["obs"][()]
        obs_col = parse_data(adata_group.obs[obs_colname])
        cats = obs_col.categories
        values_dict = {
            cat: adata_group.X[np.flatnonzero(obs_col.isin([cat])), marker_idx].mean(0)
            for cat in cats
        }

    return plot_polygons(
        adata_group, values_dict, mask, plot_format=plot_format, **kwargs
    )


def pseudospatial_categorical(
    adata_group: zarr.Group,
    obs_colname: str,
    obs_values: list[str] = None,
    mode: Literal["counts", "across", "within"] = "counts",
    mask: str = "spatial",
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask)
    validate_format(plot_format)
    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in adata")
    if mode not in ["counts", "across", "within"]:
        raise BadRequest(
            f"Invalid mode '{mode}'. Must be one of 'counts', 'across', 'within'"
        )

    mask_obs_colname = adata_group.uns["masks"][mask]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])

    obs_col = parse_data(adata_group.obs[obs_colname])

    if not len(obs_values):
        values_dict = {cat: None for cat in mask_obs_col.categories}
    else:
        crosstab = pd.crosstab(mask_obs_col, obs_col)

        if mode == "counts":
            values_dict = {
                cat: crosstab[obs_values].loc[cat].sum()
                for cat in mask_obs_col.categories
            }
        elif mode == "across":
            values_dict = {
                cat: crosstab[obs_values].loc[cat].sum() / crosstab.loc[cat].sum()
                for cat in mask_obs_col.categories
            }
        elif mode == "within":
            values_dict = {
                cat: crosstab[obs_values].loc[cat].sum()
                / crosstab[obs_values].sum().sum()
                for cat in mask_obs_col.categories
            }

    return plot_polygons(
        adata_group, values_dict, mask, plot_format=plot_format, **kwargs
    )


def pseudospatial_continuous(
    adata_group: zarr.Group,
    obs_colname: str,
    mask: str = "spatial",
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask)
    validate_format(plot_format)
    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in adata")

    mask_obs_colname = adata_group.uns["masks"][mask]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])

    obs_col = parse_data(adata_group.obs[obs_colname])

    df = pd.DataFrame({mask_obs_colname: mask_obs_col, obs_colname: obs_col})
    mean_table = df.pivot_table(
        index=mask_obs_colname, values=obs_colname, aggfunc="mean"
    )

    values_dict = mean_table.to_dict()[obs_colname]

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

    min_value = (
        float(min_value)
        if min_value is not None
        else min([v for v in values_dict.values() if v is not None] or [0])
    )
    max_value = (
        float(max_value)
        if max_value is not None
        else max([v for v in values_dict.values() if v is not None] or [1])
    )

    normalized_values = {
        k: (((v - min_value) / (max_value - min_value)) if v is not None else None)
        for k, v in values_dict.items()
    }
    color_values = {
        k: (
            (
                sample_colorscale(
                    colormap,
                    [min(max(v, 0), 1)],
                    colortype="rgb",
                )[0]
            )
            if v is not None
            else "rgba(0,0,0,0)"
        )
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
                    color=(
                        color_values[polygon]
                        if values_dict[polygon] is not None
                        else "rgb(0,0,0)"
                    ),
                    width=1,
                ),
                mode="lines",
                showlegend=False,
                fill="toself",
                fillcolor=(color_values[polygon]),
                hoverinfo="text",
                text="<br>".join(
                    [
                        f"<b>{polygon}</b>",
                        (
                            "Mean expression: 0"
                            if values_dict[polygon] == 0
                            else (
                                f"Mean expression: {values_dict[polygon]:.3e}"
                                if values_dict[polygon] is not None
                                and values_dict[polygon] < 0.001
                                else (
                                    f"Mean expression: {values_dict[polygon]:,.3f}"
                                    if values_dict[polygon] is not None
                                    else "Mean expression: N/A"
                                )
                            )
                        ),
                    ]
                ),
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
            config=dict(displaylogo=False),
            include_plotlyjs="cdn",
            full_html=full_html,
        )
    else:  # json
        json.loads(fig.to_json())
