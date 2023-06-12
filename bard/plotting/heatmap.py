import json
import zarr
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.validator_cache import ValidatorCache

from bard.utils.adata_utils import get_group_index, get_index_in_array


def heatmap(adata_group: zarr.Group, marker: str = None):
    marker_idx = get_index_in_array(get_group_index(adata_group.var), marker)

    df = pd.DataFrame(adata_group.X[:, marker_idx], columns=[marker]).transpose()
    fig = go.Figure(go.Heatmap(z=df, y=df.index.tolist()))

    return json.loads(fig.to_json())
