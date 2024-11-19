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
from cherita.utils.adata_utils import parse_data, to_categorical
from cherita.utils.models import Marker
from cherita.resources.errors import BadRequest, InvalidObs, InvalidVar, NotInData


def validate_pseudospatial(adata_group: zarr.Group, mask_set: str):
    if "masks" not in adata_group.uns or mask_set not in adata_group.uns["masks"]:
        raise NotInData(f"Mask '{mask_set}' not found in AnnData")

    if "polygons" not in adata_group.uns["masks"][mask_set] or not len(
        adata_group.uns["masks"][mask_set]["polygons"]
    ):
        raise NotInData(f"No polygons found in mask {mask_set}")


def validate_format(format: str):
    if format not in ["png", "svg", "html", "json"]:
        raise BadRequest(
            f"Invalid format '{format}'. Must be one of 'png', 'svg', 'html', 'json'"
        )


# @TODO: support obs_indices


def pseudospatial_gene(
    adata_group: zarr.Group,
    var_key: Union[int, str, dict],
    obs_col: dict = None,
    obs_values: list[str] = None,
    mask_set: str = "spatial",
    mask_values: list[str] = None,
    var_names_col: str = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)

    try:
        marker = Marker.from_any(adata_group, var_key, var_names_col)
    except InvalidVar:
        raise NotInData(f"{var_key} not found in data")

    # compute mean expression
    logging.info(f"Computing mean expression for {marker.name}")
    mask_obs_colname = adata_group.uns["masks"][mask_set]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])
    masks = mask_obs_col.categories

    if mask_values is None:
        mask_values = masks

    obs_indices = None
    if obs_col and obs_values is not None:
        obs_colname = obs_col["name"]
        try:
            obs = parse_data(adata_group.obs[obs_colname])
        except KeyError as e:
            raise InvalidObs(f"Invalid observation {e}")

        categorical_obs, _ = to_categorical(obs, **obs_col)

        if obs_values is not None:
            obs_indices = np.flatnonzero(categorical_obs.isin(obs_values))

    mask_obs_col = (
        mask_obs_col[obs_indices] if obs_indices is not None else mask_obs_col
    )

    values_dict = {}
    add_text = {}
    X = marker.get_X_at(obs_indices)
    for m in masks:
        if m not in mask_values:
            values_dict[m] = None
            add_text[m] = "0 cells"
        else:
            mask_indx = np.flatnonzero(mask_obs_col.isin([m]))
            vals = X[mask_indx] if not marker.isSet else X[:, mask_indx]
            values_dict[m] = np.mean(vals) if len(vals) else None
            add_text[m] = f"{len(mask_indx):,} cells"

    return plot_polygons(
        adata_group,
        values_dict,
        mask_set,
        text="Mean expression",
        add_text=add_text,
        plot_format=plot_format,
        **kwargs,
    )


def pseudospatial_categorical(
    adata_group: zarr.Group,
    obs_col: dict,
    obs_values: list[str] = None,
    mode: Literal["across", "within"] = "across",
    mask_set: str = "spatial",
    mask_values: list[str] = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)

    obs_colname = obs_col["name"]

    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in AnnData")
    if mode not in ["across", "within"]:
        raise BadRequest(f"Invalid mode '{mode}'. Must be 'across' or 'within'")

    mask_obs_colname = adata_group.uns["masks"][mask_set]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])
    masks = mask_obs_col.categories

    if mask_values is None:
        mask_values = masks

    try:
        obs = parse_data(adata_group.obs[obs_colname])
    except KeyError as e:
        raise InvalidObs(f"Invalid observation {e}")

    # could be removed as obs is supposed to be categorical?
    cat_obs, _ = to_categorical(obs, **obs_col)

    if obs_values is None:
        obs_values = cat_obs.categories

    if not len(obs_values):
        values_dict = {m: None for m in masks}
        text = None
        add_text = {m: None for m in masks}
    else:
        crosstab = pd.crosstab(mask_obs_col, cat_obs)
        crosstab = crosstab.loc[mask_values]

        values_dict = {}
        add_text = {}
        for m in masks:
            if m not in mask_values:
                values_dict[m] = None
                add_text[m] = "0 cells"
            else:
                if mode == "across":
                    s = crosstab[obs_values].loc[m].sum()
                    total = crosstab[obs_values].sum().sum()
                elif mode == "within":
                    s = crosstab[obs_values].loc[m].sum()
                    total = crosstab.loc[m].sum()
                values_dict[m] = s / total * 100
                add_text[m] = f"{s:,} out of {total:,}"

        text = "% across masks" if mode == "across" else "% within mask"
    return plot_polygons(
        adata_group,
        values_dict,
        mask_set,
        text=text,
        add_text=add_text,
        plot_format=plot_format,
        **kwargs,
    )


def pseudospatial_continuous(
    adata_group: zarr.Group,
    obs_col: dict,
    obs_values: list[str] = None,
    mask_set: str = "spatial",
    mask_values: list[str] = None,
    plot_format: Literal["png", "svg", "html", "json"] = "png",
    **kwargs,
):
    validate_pseudospatial(adata_group, mask_set)
    validate_format(plot_format)

    obs_colname = obs_col["name"]

    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in AnnData")

    mask_obs_colname = adata_group.uns["masks"][mask_set]["obs"][()]
    mask_obs_col = parse_data(adata_group.obs[mask_obs_colname])

    try:
        obs = parse_data(adata_group.obs[obs_colname])
    except KeyError as e:
        raise InvalidObs(f"Invalid observation {e}")

    categorical_obs, _ = to_categorical(obs, **obs_col)
    df = pd.DataFrame(
        {mask_obs_colname: mask_obs_col, obs_colname: obs, "_cat": categorical_obs}
    )

    if mask_values is None:
        mask_values = mask_obs_col.categories

    if obs_values is not None:
        df = df[df["_cat"].isin(obs_values)]

    df = df[df[mask_obs_colname].isin(mask_values)]

    mean_table = df.pivot_table(
        index=mask_obs_colname, values=obs_colname, aggfunc="mean"
    )
    count_table = df.pivot_table(
        index=mask_obs_colname, values=obs_colname, aggfunc="count"
    )

    values_dict = mean_table.to_dict()[obs_colname]
    add_text = {
        k: f"{v:,} cells" for k, v in count_table.to_dict()[obs_colname].items()
    }

    return plot_polygons(
        adata_group,
        values_dict,
        mask_set,
        text="Mean value",
        add_text=add_text,
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
    add_text: dict = {},
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

    color_values = {
        k: (
            sample_colorscale(
                colormap,
                [
                    min(
                        max(
                            (v - min_value)
                            / (
                                (max_value - min_value)
                                + (np.spacing(1) if max_value == min_value else 0)
                            ),
                            0,
                        ),
                        1,
                    )
                ],  # use np.spacing to avoid division by zero
                colortype="rgb",
            )[0]
        )
        for k, v in values_dict.items()
        if v is not None and not np.isnan(v)
    }
    text = (text + ": ") if text is not None else ""

    fig = go.Figure()

    for polygon in adata_group.uns["masks"][mask_set]["polygons"].keys():
        line_color = (
            values_dict.get(polygon) is not None
            and color_values.get(polygon)
            or "rgb(0,0,0,0.5)"
        )
        fill_color = color_values.get(polygon)
        value = values_dict.get(polygon)

        x_coords = list(
            *adata_group.uns["masks"][mask_set]["polygons"][polygon][:, :, 0, 0]
        )
        y_coords = list(
            *adata_group.uns["masks"][mask_set]["polygons"][polygon][:, :, 0, 1]
        )

        # Add first coord to end to close the polygon
        x_coords.append(x_coords[0])
        y_coords.append(y_coords[0])

        fig.add_trace(
            go.Scatter(
                x=x_coords,
                y=y_coords,
                line=dict(
                    color=line_color,
                    width=1,
                ),
                mode="lines",
                showlegend=False,
                fill="toself",
                fillcolor=fill_color or "rgba(0,0,0,0)",
                fillpattern=(
                    dict(
                        shape="/",
                        bgcolor="rgba(0,0,0,0)",
                        fgcolor="rgba(0,0,0,0.5)",
                        size=5,
                        solidity=0.2,
                    )
                    if not fill_color
                    else None
                ),
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
                        add_text.get(polygon, "") or "",
                    ]
                ),
                meta=dict(
                    value=value,
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
                coloraxis="coloraxis",
                showscale=True,
                colorbar=dict(thickness=10, tickmode="auto"),
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
        coloraxis=dict(
            colorscale=colormap,
            cmin=min_value,
            cmax=max_value,
        ),
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
        fig.update_layout(width=None, height=None)
        response_obj = json.loads(fig.to_json())
        response_obj["range"] = [float(min_value), float(max_value)]
        return response_obj
