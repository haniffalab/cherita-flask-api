from __future__ import annotations
from typing import Union, Literal
import base64
import json
import zarr
import logging
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from cherita.utils.adata_utils import parse_data
from cherita.utils.models import Marker
from cherita.resources.errors import BadRequest, NotInData


def validate_pseudospatial(adata_group: zarr.Group, mask_set: str):
    if "masks" not in adata_group.uns or mask_set not in adata_group.uns["masks"]:
        raise NotInData(f"Mask '{mask_set}' not found in adata")

    if "polygons" not in adata_group.uns["masks"][mask_set] or not len(
        adata_group.uns["masks"][mask_set]["polygons"]
    ):
        raise NotInData(f"No polygons found in mask {mask_set}")


def validate_format(format: str):
    if format not in ["png", "svg", "html", "json"]:
        raise BadRequest(
            f"Invalid format '{format}'. Must be one of 'png', 'svg', 'html', 'json'"
        )


def pseudospatial_gene(
    adata_group: zarr.Group,
    var_key: Union[int, str, dict],
    mask_set: str = "spatial",
    mask_values: list[str] = None,
    var_names_col: str = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)

    marker = Marker.from_any(adata_group, var_key, var_names_col)

    # compute mean expression
    logging.info(f"Computing mean expression for {marker.name}")
    mask_obs_colname = adata_group.uns["masks"][mask_set]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])
    masks = mask_obs_col.categories
    values_dict = {
        m: (
            None
            if mask_values and m not in mask_values
            else marker.X[np.flatnonzero(mask_obs_col.isin([m]))].mean(0)
        )
        for m in masks
    }

    return plot_polygons(
        adata_group,
        values_dict,
        mask_set,
        text="Mean expression",
        plot_format=plot_format,
        **kwargs,
    )


def pseudospatial_categorical(
    adata_group: zarr.Group,
    obs_colname: str,
    obs_values: list[str] = None,
    mode: Literal["counts", "across", "within"] = "counts",
    mask_set: str = "spatial",
    mask_values: list[str] = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)
    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in adata")
    if mode not in ["counts", "across", "within"]:
        raise BadRequest(
            f"Invalid mode '{mode}'. Must be one of 'counts', 'across', 'within'"
        )

    mask_obs_colname = adata_group.uns["masks"][mask_set]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])
    masks = mask_obs_col.categories

    obs_col = parse_data(adata_group.obs[obs_colname])

    if not len(obs_values):
        values_dict = {m: None for m in masks}
    else:
        crosstab = pd.crosstab(mask_obs_col, obs_col)
        if mask_values:
            crosstab = crosstab.loc[mask_values]

        if mode == "counts":
            values_dict = {
                m: (
                    None
                    if mask_values and m not in mask_values
                    else crosstab[obs_values].loc[m].sum()
                )
                for m in masks
            }
            text = "Counts"
        elif mode == "across":
            values_dict = {
                m: (
                    None
                    if mask_values and m not in mask_values
                    else (
                        crosstab[obs_values].loc[m].sum()
                        / crosstab[obs_values].sum().sum()
                    )
                    * 100
                )
                for m in masks
            }
            text = "% across regions"
        elif mode == "within":
            values_dict = {
                m: (
                    None
                    if mask_values and m not in mask_values
                    else (crosstab[obs_values].loc[m].sum() / crosstab.loc[m].sum())
                    * 100
                )
                for m in masks
            }
            text = "% within region"

    return plot_polygons(
        adata_group, values_dict, mask_set, text=text, plot_format=plot_format, **kwargs
    )


def pseudospatial_continuous(
    adata_group: zarr.Group,
    obs_colname: str,
    mask_set: str = "spatial",
    mask_values: list[str] = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)
    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in adata")

    mask_obs_colname = adata_group.uns["masks"][mask_set]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])

    obs_col = parse_data(adata_group.obs[obs_colname])

    df = pd.DataFrame({mask_obs_colname: mask_obs_col, obs_colname: obs_col})
    if mask_values:
        df = df[df[mask_obs_colname].isin(mask_values)]
    mean_table = df.pivot_table(
        index=mask_obs_colname, values=obs_colname, aggfunc="mean"
    )

    values_dict = mean_table.to_dict()[obs_colname]

    return plot_polygons(
        adata_group,
        values_dict,
        mask_set,
        text="Mean value",
        plot_format=plot_format,
        **kwargs,
    )


def pseudospatial_masks(
    adata_group: zarr.Group,
    mask_set: str = "spatial",
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    width: int = 500,
    height: int = 500,
    full_html: bool = False,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)

    values_dict = {
        k: None for k in adata_group.uns["masks"][mask_set]["polygons"].keys()
    }

    return plot_polygons(
        adata_group,
        values_dict,
        mask_set,
        show_colorbar=False,
        plot_format=plot_format,
        width=width,
        height=height,
        full_html=full_html,
    )


def plot_polygons(
    adata_group: zarr.Group,
    values_dict: pd.DataFrame,
    mask_set: str = "spatial",
    colormap: str = "viridis",
    show_colorbar: bool = True,
    min_value: float = None,
    max_value: float = None,
    text: str = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    width: int = 500,
    height: int = 500,
    full_html: bool = False,
):
    min_value = (
        float(min_value)
        if min_value is not None
        else min(
            [v for v in values_dict.values() if v is not None and not np.isnan(v)]
            or [0]
        )
    )
    max_value = (
        float(max_value)
        if max_value is not None
        else max(
            [v for v in values_dict.values() if v is not None and not np.isnan(v)]
            or [1]
        )
    )

    normalized_values = {
        k: ((v - min_value) / (max_value - min_value))
        for k, v in values_dict.items()
        if v is not None and not np.isnan(v)
    }
    color_values = {
        k: (
            sample_colorscale(
                colormap,
                [min(max(v, 0), 1)],
                colortype="rgb",
            )[0]
        )
        for k, v in normalized_values.items()
        if v is not None
    }
    text = (text + ": ") if text is not None else ""

    fig = go.Figure()

    for polygon in adata_group.uns["masks"][mask_set]["polygons"].keys():
        line_color = (
            values_dict.get(polygon) is not None
            and color_values.get(polygon)
            or "rgb(0,0,0)"
        )
        fill_color = color_values.get(polygon, "rgba(0,0,0,0)")
        value = values_dict.get(polygon)

        fig.add_trace(
            go.Scatter(
                x=list(
                    *adata_group.uns["masks"][mask_set]["polygons"][polygon][:, :, 0, 0]
                ),
                y=list(
                    *adata_group.uns["masks"][mask_set]["polygons"][polygon][:, :, 0, 1]
                ),
                line=dict(
                    color=line_color,
                    width=1,
                ),
                mode="lines",
                showlegend=False,
                fill="toself",
                fillcolor=fill_color,
                hoverinfo="text",
                text="<br>".join(
                    [
                        f"<b>{polygon}</b>",
                        (
                            f"{text}0"
                            if value == 0
                            else (
                                f"{text}{value:.3e}"
                                if value is not None and value < 0.001
                                else (
                                    f"{text}{value:,.3f}"
                                    if value is not None
                                    else f"{text}N/A"
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
        return json.loads(fig.to_json())
